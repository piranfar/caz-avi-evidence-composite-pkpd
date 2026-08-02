"""Analyses added in response to peer review.

Three questions the reviewed draft could not answer:

1. The avibactam target. The draft used fCss >= 4 mg/L, the target adopted by the
   continuous-infusion TDM literature, and reported that avibactam limits joint
   attainment. But 4 mg/L is the fixed avibactam concentration used in
   susceptibility testing, whereas the registrational PK/PD analyses used a
   threshold of 1 mg/L. Since the one-at-a-time analysis already ranked this
   threshold above every pharmacokinetic parameter, the conclusion has to be
   shown across the range of thresholds actually used in the literature.

2. The lowest renal class. The scenario excludes renal replacement therapy, yet
   renal function in the 0-30 class is drawn uniformly down to zero, so a
   substantial fraction of that class would in practice be dialysed. Because
   those subjects have the lowest clearance they dominate the exposure ceiling,
   so the permissibility result for that class needs to be shown with and
   without them.

3. Monte Carlo precision. Attainment was reported to two decimals; the binomial
   standard error at 20,000 subjects per class is around 0.3 percentage points,
   so the reported precision needs quantifying.

Usage:
    python reviewer_response_analyses.py
"""

from __future__ import annotations

import csv
import os

import numpy as np

from cazavi_analyses import (
    REGIMENS, SELECTED_REGIMENS, Scenario, clearances, compute_cfr,
    draw_population, evaluate, load_mic_distributions, write_csv,
)
from reproduce_primary_run import (
    AVI_FRACTION, CAZ_FRACTION, EKFC_CLASSES, MIC_GRID, N_PER_CLASS,
    OMEGA_AVI, OMEGA_CAZ, PRIMARY_SEED, RHO, TOX_THRESHOLD,
)

from cazavi_analyses import DEFAULT_OUT as OUT

HERE = os.path.dirname(os.path.abspath(__file__))

# Thresholds spanning the published range: 1 mg/L is the registrational PK/PD
# threshold used for dose selection; 4 mg/L is both the fixed concentration in
# susceptibility testing and the steady-state target used in the CI/TDM work.
AVI_THRESHOLDS = (1.0, 2.0, 4.0, 6.0, 8.0)

# Lowest renal function retained in the sensitivity analysis. Below roughly
# 15 mL/min/1.73 m2 a critically ill patient would ordinarily be receiving renal
# replacement therapy, which the primary scenario excludes.
RRT_FLOOR = 15.0


def avibactam_threshold_sweep(dists):
    """Joint attainment and CFR across the published range of avibactam targets."""
    population = draw_population(N_PER_CLASS, PRIMARY_SEED)
    rows, cfr_rows = [], []
    for ct in AVI_THRESHOLDS:
        sc = Scenario(name=f"AVI_CT_{ct:g}", avi_ct=ct)
        evaluated = evaluate(population, sc)
        for r in evaluated:
            if r["regimen"] in SELECTED_REGIMENS and r["mic_mg_l"] in (2, 4, 8):
                rows.append({"avi_ct_mg_l": ct, **{k: r[k] for k in (
                    "regimen", "ekfc_class", "mic_mg_l", "caz_pta_pct",
                    "avi_attainment_pct", "joint_pta_pct")}})
        for c in compute_cfr([r for r in evaluated
                              if r["regimen"] in SELECTED_REGIMENS], dists):
            if c["distribution_id"] == "LEE2022_KPC_KP":
                cfr_rows.append({"avi_ct_mg_l": ct, "regimen": c["regimen"],
                                 "ekfc_class": c["ekfc_class"],
                                 "caz_cfr_pct": c["caz_cfr_pct"],
                                 "avi_weighted_pct": c["avi_weighted_pct"],
                                 "joint_cfr_pct": c["joint_cfr_pct"]})
    return rows, cfr_rows


def limiting_component_by_threshold(rows):
    """How often each component binds, as the avibactam threshold moves."""
    out = []
    for ct in AVI_THRESHOLDS:
        for mic in (2, 4, 8):
            sel = [r for r in rows if r["avi_ct_mg_l"] == ct and r["mic_mg_l"] == mic]
            avi = sum(1 for r in sel if r["avi_attainment_pct"] < r["caz_pta_pct"] - 0.5)
            caz = sum(1 for r in sel if r["caz_pta_pct"] < r["avi_attainment_pct"] - 0.5)
            joint = [r["joint_pta_pct"] for r in sel]
            out.append({"avi_ct_mg_l": ct, "mic_mg_l": mic, "n_regimens": len(sel),
                        "avibactam_limiting": avi, "ceftazidime_limiting": caz,
                        "joint_pta_min_pct": min(joint), "joint_pta_max_pct": max(joint)})
    return out


def rrt_floor_sensitivity():
    """Exposure and the toxicity screen with and without sub-dialysis renal function."""
    rng = np.random.default_rng(PRIMARY_SEED)
    cov = np.array([[OMEGA_CAZ**2, RHO * OMEGA_CAZ * OMEGA_AVI],
                    [RHO * OMEGA_CAZ * OMEGA_AVI, OMEGA_AVI**2]])
    rows = []
    for cls, (lo, hi) in EKFC_CLASSES.items():
        ekfc = rng.uniform(lo, hi, N_PER_CLASS)
        eta = rng.multivariate_normal(np.zeros(2), cov, size=N_PER_CLASS)
        keep = ekfc >= RRT_FLOOR
        for regimen, (rc, dose_g, interval_h) in REGIMENS.items():
            if rc != cls or regimen not in SELECTED_REGIMENS:
                continue
            from reproduce_primary_run import CL0_CAZ, EXP_CAZ, EKFC_REF, FU_CAZ
            cl = CL0_CAZ * (ekfc / EKFC_REF) ** EXP_CAZ * np.exp(eta[:, 0])
            css = dose_g * 1000.0 * CAZ_FRACTION / interval_h / cl
            tox = css > TOX_THRESHOLD
            rows.append({
                "ekfc_class": cls, "regimen": regimen,
                "pct_below_floor": 100.0 * (~keep).mean(),
                "toxicity_all_pct": 100.0 * tox.mean(),
                "toxicity_above_floor_pct": 100.0 * tox[keep].mean() if keep.any() else float("nan"),
                "share_of_exceedances_below_floor_pct":
                    100.0 * tox[~keep].sum() / tox.sum() if tox.sum() else 0.0,
                "median_css_all": float(np.median(css)),
                "median_css_above_floor": float(np.median(css[keep])) if keep.any() else float("nan"),
            })
    return rows


def monte_carlo_precision(n_per_class=N_PER_CLASS):
    """Binomial standard error and half-width of the 95% interval for a proportion."""
    rows = []
    for p in (0.10, 0.30, 0.50, 0.70, 0.90):
        se = np.sqrt(p * (1 - p) / n_per_class)
        rows.append({"attainment_pct": 100 * p,
                     "standard_error_pp": round(100 * se, 3),
                     "half_width_95ci_pp": round(100 * 1.96 * se, 3)})
    return rows


def full_regimen_grid(dists):
    """Every evaluated regimen, not only the selected one per class."""
    population = draw_population(N_PER_CLASS, PRIMARY_SEED)
    pta = evaluate(population)
    cfr = {(c["regimen"], c["distribution_id"]): c for c in compute_cfr(pta, dists)}
    rows = []
    for regimen, (cls, dose_g, interval_h) in REGIMENS.items():
        at = {r["mic_mg_l"]: r for r in pta if r["regimen"] == regimen}
        primary = cfr[(regimen, "LEE2022_KPC_KP")]
        rows.append({
            "regimen": regimen, "ekfc_class": cls,
            "product_dose_g": dose_g, "interval_h": interval_h,
            "daily_dose_g": round(dose_g * 24.0 / interval_h, 3),
            "joint_pta_mic2_pct": round(at[2]["joint_pta_pct"], 1),
            "joint_pta_mic4_pct": round(at[4]["joint_pta_pct"], 1),
            "joint_pta_mic8_pct": round(at[8]["joint_pta_pct"], 1),
            "caz_cfr_pct": round(primary["caz_cfr_pct"], 1),
            "joint_cfr_pct": round(primary["joint_cfr_pct"], 1),
            "toxicity_pct": round(at[4]["toxicity_pct"], 1),
            "permissible": "yes" if at[4]["toxicity_pct"] <= 15 else "no",
            "selected": "yes" if regimen in SELECTED_REGIMENS else "no",
        })
    return rows


def main():
    dists = load_mic_distributions()

    sweep, sweep_cfr = avibactam_threshold_sweep(dists)
    write_csv(sweep, os.path.join(OUT, "avi_threshold_sweep_pta.csv"))
    write_csv(sweep_cfr, os.path.join(OUT, "avi_threshold_sweep_cfr.csv"))
    limiting = limiting_component_by_threshold(sweep)
    write_csv(limiting, os.path.join(OUT, "avi_threshold_limiting_component.csv"))

    floor = rrt_floor_sensitivity()
    write_csv(floor, os.path.join(OUT, "rrt_floor_sensitivity.csv"))
    write_csv(monte_carlo_precision(), os.path.join(OUT, "monte_carlo_precision.csv"))
    write_csv(full_regimen_grid(dists), os.path.join(OUT, "full_regimen_grid.csv"))

    print("\navibactam threshold: limiting component and joint PTA")
    print(f"  {'CT':>4} {'MIC':>4} {'AVI lim':>8} {'CAZ lim':>8} {'joint PTA range':>20}")
    for r in limiting:
        print(f"  {r['avi_ct_mg_l']:>4g} {r['mic_mg_l']:>4g} {r['avibactam_limiting']:>8}"
              f" {r['ceftazidime_limiting']:>8}"
              f"   {r['joint_pta_min_pct']:>6.1f}–{r['joint_pta_max_pct']:.1f}")

    print("\njoint CFR by avibactam threshold (KPC-KP, selected regimens)")
    for ct in AVI_THRESHOLDS:
        vals = [r["joint_cfr_pct"] for r in sweep_cfr if r["avi_ct_mg_l"] == ct]
        avi = [r["avi_weighted_pct"] for r in sweep_cfr if r["avi_ct_mg_l"] == ct]
        print(f"  CT {ct:g} mg/L: joint CFR {min(vals):5.1f}–{max(vals):5.1f}%"
              f"   avibactam attainment {min(avi):5.1f}–{max(avi):5.1f}%")

    print("\nsub-dialysis renal function in the lowest class")
    for r in floor:
        if r["ekfc_class"] != "0–30":
            continue
        print(f"  {r['regimen']}: {r['pct_below_floor']:.1f}% below {RRT_FLOOR:g} mL/min/1.73 m2;"
              f" they contribute {r['share_of_exceedances_below_floor_pct']:.1f}% of exceedances")
        print(f"      exposure-screen exceedance  all {r['toxicity_all_pct']:.2f}%"
              f"  ->  restricted {r['toxicity_above_floor_pct']:.2f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
