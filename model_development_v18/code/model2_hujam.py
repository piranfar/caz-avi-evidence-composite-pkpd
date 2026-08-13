"""MODEL 2 — Hierarchical Uncertainty-integrated Joint Attainment Model (HU-JAM).

Decision-analytic layer over the verified engine in `model2_engine.py`.

WHAT IS NEW HERE
    The primary model reports a point estimate of attainment and explores
    sensitivity one parameter at a time. This module propagates the published
    uncertainty jointly and converts it into decisions:

      * attainment as a distribution, not a point estimate
      * P(avibactam is the limiting component | MIC) — a probability, not a label
      * P(each regimen is the optimal one)
      * P(a single fixed avibactam target selects the WRONG regimen)
      * expected regret of the fixed-threshold decision rule
      * expected value of perfect information, overall and per parameter

WHAT IT IS NOT
    Not a pharmacokinetic model. The target distributions are NOT established
    exposure-response relationships. The rho scenarios are NOT a pooled posterior
    across populations. No output is a bedside dosing recommendation.

Run:  python model2_hujam.py            (writes CSVs to ../outputs)
      python model2_hujam.py --quick    (smaller run for a smoke test)
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import replace

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import model2_engine as E   # noqa: E402
import reproduce_primary_run as P   # noqa: E402

OUT = os.path.join(os.path.dirname(E.HERE), "outputs")
FIG = os.path.join(os.path.dirname(E.HERE), "figures")

MASTER_SEED = 20260811
TOX_CEILING = 15.0          # exposure-screen permissibility ceiling, analyst-specified
PRIMARY_DIST = "LEE2022_KPC_KP"


# ---------------------------------------------------------------- scenarios --
# LAYER 4 — the avibactam target. Every distribution is labelled by KIND.
# None of these is an established clinical exposure-response relationship.

TARGET_SCENARIOS = {
    "T1_point_1": dict(
        kind="point mass", label="1 mg/L (regulatory target)",
        note="Animal-derived: Berkhout et al. 2015/2016, murine neutropenic thigh/lung "
             "infection, dose-fractionation (Antimicrob Agents Chemother 60(1):368-75, "
             "doi:10.1128/AAC.01269-15, PMC4704241). Expressed in the source as %fT>C_T "
             "(percentage of the dosing interval free avibactam exceeds 1 mg/L), a "
             "time-above-threshold index, NOT a steady-state concentration -- and, per "
             "R3 (NOVELTY_ROUTES.md), the source study never tested continuous infusion, "
             "so there is no published data point at which this index was actually "
             "evaluated under the exposure this scenario borrows it for. Used here as a "
             "steady-state target only for comparability with T2-T7, not because the "
             "source validates that use."),
    "T2_point_4": dict(
        kind="point mass", label="4 mg/L (EUCAST testing concentration)",
        note="The manuscript's current primary scenario. A susceptibility-testing "
             "convention repurposed as a clinical target."),
    "T3_discrete": dict(
        kind="discrete, equal weights", label="{0.5, 1, 4} mg/L, 1/3 each",
        note="The three distinct constructs in the evidence: hollow-fibre regrowth "
             "threshold, regulatory target, EUCAST convention. Equal weights are an "
             "ANALYST CHOICE expressing no preference, not an evidence synthesis."),
    "T4_uniform": dict(
        kind="uniform", label="U(0.5, 4) mg/L",
        note="Agnostic across the defensible range. ANALYST-SPECIFIED bounds."),
    "T5_triangular": dict(
        kind="triangular", label="Tri(0.25, 4; mode 1) mg/L",
        note="Mass concentrated at the regulatory value. ANALYST-SPECIFIED."),
    "T6_lognormal": dict(
        kind="lognormal", label="median 1, 95th centile 4 mg/L",
        note="Smooth and right-skewed. ANALYST-SPECIFIED."),
    "T7_coleman_evidence": dict(
        kind="discrete, equal weights, EVIDENCE-DERIVED",
        label="{0.15, 0.22, 0.28} mg/L, 1/3 each",
        note="NOT analyst-specified. Coleman et al. 2014 (Antimicrob Agents Chemother "
             "58(6):3366-72, doi:10.1128/AAC.00080-14, PMC4068505), Table 2: the avibactam "
             "concentration below which regrowth commenced, in the continuous-infusion-"
             "ceftazidime / single-bolus-avibactam hollow-fibre experiment. THREE strains "
             "(not eight -- eight strains were used elsewhere in the same paper, for a "
             "separate single-dose killing experiment that did not estimate this threshold): "
             "E. cloacae 293HT96 (CT <=0.15 mg/L at 12 h), K. pneumoniae 283CF5 (CT <=0.22 "
             "mg/L), K. pneumoniae Tunisie K4 (CT <=0.28 mg/L). A fourth Table 2 entry, a "
             "re-estimate of the SAME E. cloacae strain at 18-20 h (CT ~0.2 mg/L), is excluded "
             "here to avoid double-counting one strain; it is a rough within-strain check, not "
             "a fourth independent data point (0.15 vs ~0.2 at the two timepoints, same "
             "direction, same order of magnitude). "
             "This is deliberately NOT presented as a random-effects synthesis with an "
             "estimated between-strain variance: three points (one of them not fully "
             "independent) cannot support one. It is the plain empirical distribution over the "
             "three measured strains -- a new strain from the same source population is modelled "
             "as equally likely to resemble any of the three actually measured, no smoother "
             "and no more precise than that. All three source values are upper bounds "
             "('<=', extrapolated from exponential-decline curves at the last pre-regrowth "
             "timepoint sampled), so using them as point values is itself a conservative "
             "(upward-biased) simplification of an already-small dataset. "
             "This is a distribution over an IN-VITRO REGROWTH THRESHOLD in a hollow-fibre "
             "model, not over the clinical target -- the gap between the two is exactly the "
             "problem T1-T6 already document, and this scenario does not close it."),
}
# 2.5 mg/L is EXCLUDED throughout: it is the aztreonam-avibactam target and does
# not belong to a ceftazidime model.

# LAYER 3 — the clearance correlation. Two populations, never pooled.
RHO_SCENARIOS = {
    "C1_cojutti": dict(
        label="0.94, published RSE 23.8%",
        note="Cojutti 2024, non-RRT ICU on continuous infusion. THE PRIMARY "
             "ANALYSIS — the manuscript's own assumption, with its own uncertainty."),
    "C2_model1": dict(
        label="0.703, profile-likelihood 0.380-0.874",
        note="Model 1, CRRT cohort on intermittent infusion. A DIFFERENT "
             "POPULATION. Sensitivity analysis only; not a replacement value."),
    "C3_agnostic": dict(
        label="U(0.38, 0.98)",
        note="Spans the union of both intervals. ANALYST-SPECIFIED."),
}


def sample_target(rng, scenario):
    if scenario == "T1_point_1":
        return 1.0
    if scenario == "T2_point_4":
        return 4.0
    if scenario == "T3_discrete":
        return float(rng.choice([0.5, 1.0, 4.0]))
    if scenario == "T4_uniform":
        return float(rng.uniform(0.5, 4.0))
    if scenario == "T5_triangular":
        return float(rng.triangular(0.25, 1.0, 4.0))
    if scenario == "T6_lognormal":
        return float(np.exp(rng.normal(0.0, np.log(4.0) / 1.6448536)))
    if scenario == "T7_coleman_evidence":
        return float(rng.choice([0.15, 0.22, 0.28]))
    raise KeyError(scenario)


def sample_rho(rng, scenario):
    if scenario == "C1_cojutti":
        # RSE 23.8% on a correlation of 0.94, applied on the Fisher-z scale
        # rather than on rho directly. The source does not state the scale on
        # which the standard error was computed, so applying it is a scenario
        # assumption either way -- but sampling Normal(0.94, 0.238*0.94) on the
        # rho scale directly and clipping to a valid correlation was tested
        # (test_model2.py) and found to clip 39.5% of draws to the 0.999
        # boundary, pulling the sampled MEAN down to 0.877 -- a full 0.06 below
        # the value this scenario is meant to represent. Sampling on the Fisher-z
        # scale, as C2 already does below, removes the boundary clipping
        # entirely and recovers a median of 0.940 (mean 0.919, pulled down only
        # by Jensen's inequality, not by truncation).
        z0 = np.arctanh(0.94)
        sd = 0.238 * z0
        return float(np.tanh(rng.normal(z0, sd)))
    if scenario == "C2_model1":
        # Fisher z, calibrated so the 95% interval matches the profile likelihood
        z0 = np.arctanh(0.703)
        sd = (np.arctanh(0.874) - np.arctanh(0.380)) / (2 * 1.959964)
        return float(np.tanh(rng.normal(z0, sd)))
    if scenario == "C3_agnostic":
        return float(rng.uniform(0.38, 0.98))
    raise KeyError(scenario)


# ------------------------------------------------------------------- utility --

LAMBDA = 1.0     # exchange rate, see below. ANALYST-SPECIFIED.


def utility(res, weights, lam=LAMBDA):
    """Net benefit: joint CFR minus an exposure penalty.

        U = joint CFR  -  lambda * exceedance

    lambda is the exchange rate between a percentage point of joint attainment and
    a percentage point of subjects above the 104 mg/L exposure screen. lambda = 1
    treats them as equally weighted. It is ANALYST-SPECIFIED and varied in
    sensitivity analysis.

    An earlier version used the manuscript's 15% ceiling as a hard feasibility
    constraint, giving zero utility above it. That was discarded: the source
    describes the 15% cut-off as arbitrary, and a hard constraint turns an
    acknowledged arbitrary choice into a cliff where 14.9% exceedance scores 85 and
    15.1% scores 0. Expected regret then measures the cliff rather than the decision,
    which is what the first run showed — regret of 20 to 33 percentage points with a
    95th centile near 95, driven entirely by draws that crossed the threshold.
    A linear penalty keeps the trade-off explicit and the regret interpretable.
    """
    return float(res["joint_pta"] @ weights) - lam * res["exceedance"]


def utility_constrained(res, weights, ceiling=TOX_CEILING):
    """The hard-constraint variant, retained and reported as a secondary analysis."""
    return 0.0 if res["exceedance"] > ceiling else float(res["joint_pta"] @ weights)


# ------------------------------------------------------------------ main loop --

def run(n_draws, n_per_class, rho_scenario, target_scenario, seed=MASTER_SEED,
        dist_id=PRIMARY_DIST, tau_between=0.0):
    """One outer loop over parameter, correlation and target uncertainty.

    The virtual population is drawn ONCE and reused across draws (common random
    numbers), so a difference between draws reflects the parameter change and not
    Monte Carlo noise. This is the same device the v16 code uses for its
    one-at-a-time analyses, applied here to the full joint uncertainty.
    """
    pop = E.draw_population(n_per_class, seed)
    dists = E.load_mic_distributions()
    weights = dists[dist_id]["w"]
    mics = np.asarray(P.MIC_GRID, float)
    regimens = list(P.REGIMENS)
    rng = np.random.default_rng(seed + 7919)

    n_r = len(regimens)
    U = np.zeros((n_draws, n_r))
    cfr_all = np.zeros((n_draws, n_r))
    exceed = np.zeros((n_draws, n_r))
    jpta = np.zeros((n_draws, n_r, len(mics)))
    avi_lim = np.zeros((n_draws, n_r, len(mics)), dtype=bool)
    drawn = {"rho": np.zeros(n_draws), "target": np.zeros(n_draws),
             "cl0_avi": np.zeros(n_draws), "omega_avi": np.zeros(n_draws),
             "cl0_caz": np.zeros(n_draws), "fu_avi": np.zeros(n_draws)}

    for m in range(n_draws):
        pr = E.draw_parameters(rng, tau_between=tau_between)
        pr = replace(pr, rho=sample_rho(rng, rho_scenario),
                     avi_target=sample_target(rng, target_scenario))
        for key in drawn:
            drawn[key][m] = pr.avi_target if key == "target" else getattr(pr, key)

        res = E.evaluate(pop, pr, regimens)
        for j, reg in enumerate(regimens):
            d = res[reg]
            U[m, j] = utility(d, weights)
            cfr_all[m, j] = float(d["joint_pta"] @ weights)
            exceed[m, j] = d["exceedance"]
            jpta[m, j] = d["joint_pta"]
            avi_lim[m, j] = d["avi_pta"] < d["caz_pta"]

    return dict(regimens=regimens, mics=mics, U=U, cfr=cfr_all, exceed=exceed,
                jpta=jpta, avi_lim=avi_lim, drawn=drawn, weights=weights,
                rho_scenario=rho_scenario, target_scenario=target_scenario,
                n_draws=n_draws, n_per_class=n_per_class, tau_between=tau_between)


# -------------------------------------------------------------------- outputs --

def by_class(r):
    """Group regimen indices by renal-function class."""
    out = {}
    for j, reg in enumerate(r["regimens"]):
        out.setdefault(P.REGIMENS[reg][0], []).append(j)
    return out


def integrated_attainment(r):
    rows = []
    for j, reg in enumerate(r["regimens"]):
        c = r["cfr"][:, j]
        rows.append({
            "rho_scenario": r["rho_scenario"], "target_scenario": r["target_scenario"],
            "regimen": reg, "ekfc_class": P.REGIMENS[reg][0],
            "joint_cfr_median": round(float(np.median(c)), 2),
            "joint_cfr_p2.5": round(float(np.percentile(c, 2.5)), 2),
            "joint_cfr_p97.5": round(float(np.percentile(c, 97.5)), 2),
            "prediction_interval_width_pp": round(
                float(np.percentile(c, 97.5) - np.percentile(c, 2.5)), 2),
            "p_permissible_pct": round(100.0 * float(np.mean(
                r["exceed"][:, j] <= TOX_CEILING)), 1),
            "p_joint_cfr_ge_80_pct": round(100.0 * float(np.mean(c >= 80.0)), 1),
        })
    return rows


def limiting_probability(r):
    rows = []
    for j, reg in enumerate(r["regimens"]):
        for k, mic in enumerate(r["mics"]):
            p = 100.0 * float(r["avi_lim"][:, j, k].mean())
            rows.append({
                "rho_scenario": r["rho_scenario"],
                "target_scenario": r["target_scenario"],
                "regimen": reg, "ekfc_class": P.REGIMENS[reg][0], "mic_mg_l": mic,
                "p_avibactam_limiting_pct": round(p, 1),
                "p_ceftazidime_limiting_pct": round(100.0 - p, 1),
                "classification": ("avibactam" if p >= 90 else
                                   "ceftazidime" if p <= 10 else "uncertain"),
            })
    return rows


def optimality_and_regret(r):
    """P(regimen optimal), P(misselection under a fixed target), expected regret, EVPI."""
    groups = by_class(r)
    U = r["U"]

    # the fixed-threshold decision: point estimates, avibactam target 4 mg/L
    pop = E.draw_population(r["n_per_class"], MASTER_SEED)
    fixed = E.evaluate(pop, replace(E.BASE, avi_target=4.0), r["regimens"])
    u_fixed = np.array([utility(fixed[reg], r["weights"]) for reg in r["regimens"]])

    rows, summary = [], []
    for cls, idx in groups.items():
        idx = np.array(idx)
        Uc = U[:, idx]
        best = idx[np.argmax(Uc, axis=1)]
        r_fixed = int(idx[int(np.argmax(u_fixed[idx]))])

        for j in idx:
            rows.append({
                "rho_scenario": r["rho_scenario"],
                "target_scenario": r["target_scenario"],
                "ekfc_class": cls, "regimen": r["regimens"][j],
                "p_optimal_pct": round(100.0 * float(np.mean(best == j)), 1),
                "is_fixed_threshold_choice": "yes" if j == r_fixed else "no",
                "mean_utility": round(float(U[:, j].mean()), 2),
            })

        u_best = Uc.max(axis=1)
        regret = u_best - U[:, r_fixed]
        evpi = float(u_best.mean() - Uc.mean(axis=0).max())
        summary.append({
            "rho_scenario": r["rho_scenario"],
            "target_scenario": r["target_scenario"], "ekfc_class": cls,
            "fixed_threshold_regimen": r["regimens"][r_fixed],
            "p_misselection_pct": round(100.0 * float(np.mean(best != r_fixed)), 1),
            "expected_regret_pp": round(float(regret.mean()), 3),
            "regret_p95_pp": round(float(np.percentile(regret, 95)), 3),
            "max_regret_pp": round(float(regret.max()), 3),
            "evpi_pp": round(evpi, 3),
            "modal_optimal_regimen": r["regimens"][int(np.bincount(best).argmax())],
        })
    return rows, summary


def evppi(r, parameter, degree=4):
    """Expected value of partial perfect information (Strong & Oakley, 2014).

    Regress each regimen's utility on the parameter of interest across the existing
    single-loop sample, then take the expectation of the maximum of the fitted
    conditional expectations. No nested simulation is required.
    """
    x = r["drawn"][parameter]
    if np.std(x) < 1e-12:
        return {"parameter": parameter, "evppi_pp": 0.0,
                "note": "parameter held fixed in this scenario"}
    xs = (x - x.mean()) / x.std()
    groups = by_class(r)
    out = []
    for cls, idx in groups.items():
        idx = np.array(idx)
        Uc = r["U"][:, idx]
        fitted = np.column_stack([
            np.polyval(np.polyfit(xs, Uc[:, k], degree), xs)
            for k in range(Uc.shape[1])])
        val = float(fitted.max(axis=1).mean() - Uc.mean(axis=0).max())
        out.append({"rho_scenario": r["rho_scenario"],
                    "target_scenario": r["target_scenario"],
                    "ekfc_class": cls, "parameter": parameter,
                    "evppi_pp": round(max(val, 0.0), 3)})
    return out


def write(rows, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {os.path.relpath(path, os.path.dirname(E.HERE))} ({len(rows)} rows)")


def convergence_check(n_per_class, seeds=(20260811, 20260812, 20260813)):
    """Stability of the slowest-converging outputs across independent seeds."""
    rows = []
    for m in (250, 500, 1000, 2000):
        vals = []
        for s in seeds:
            r = run(m, n_per_class, "C1_cojutti", "T4_uniform", seed=s)
            _, summ = optimality_and_regret(r)
            vals.append((float(np.mean([x["evpi_pp"] for x in summ])),
                         float(np.mean([x["p_misselection_pct"] for x in summ]))))
        evpi = np.array([v[0] for v in vals])
        mis = np.array([v[1] for v in vals])
        rows.append({"n_draws": m, "n_seeds": len(seeds),
                     "evpi_mean_pp": round(float(evpi.mean()), 3),
                     "evpi_range_pp": round(float(evpi.max() - evpi.min()), 3),
                     "misselection_mean_pct": round(float(mis.mean()), 2),
                     "misselection_range_pct": round(float(mis.max() - mis.min()), 2)})
        print(f"    M={m:5}  EVPI {evpi.mean():6.3f} (range {evpi.max()-evpi.min():.3f})"
              f"   misselection {mis.mean():5.2f}% (range {mis.max()-mis.min():.2f})")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    n_draws = 300 if args.quick else 2000
    n_per_class = 2000 if args.quick else 4000

    print("=" * 76)
    print("MODEL 2 — HU-JAM")
    print("=" * 76)
    ok, _ = E.verify_against_frozen(verbose=True)
    if not ok:
        print("  ABORT: the engine no longer reproduces v16.")
        return 1
    print(f"\n  outer draws {n_draws}, virtual subjects {n_per_class*5:,}, "
          f"seed {MASTER_SEED}")
    print(f"  primary MIC distribution {PRIMARY_DIST}\n")

    attain, limiting, optimal, summary, voi = [], [], [], [], []
    combos = [(c, t) for c in RHO_SCENARIOS for t in TARGET_SCENARIOS]
    for i, (c, t) in enumerate(combos, 1):
        r = run(n_draws, n_per_class, c, t)
        attain += integrated_attainment(r)
        limiting += limiting_probability(r)
        o, s = optimality_and_regret(r)
        optimal += o
        summary += s
        for par in ("rho", "target", "cl0_avi", "omega_avi", "cl0_caz", "fu_avi"):
            res = evppi(r, par)
            if isinstance(res, list):
                voi += res
        print(f"  [{i:2}/{len(combos)}] {c:12} x {t:14}  "
              f"misselection {np.mean([x['p_misselection_pct'] for x in s]):5.1f}%   "
              f"EVPI {np.mean([x['evpi_pp'] for x in s]):6.3f} pp")

    write(attain, "model2_integrated_attainment.csv")
    write(limiting, "model2_limiting_probability.csv")
    write(optimal, "model2_regimen_optimality.csv")
    write(summary, "model2_misselection_regret.csv")
    write(voi, "model2_evppi.csv")

    # --- exchange-rate sensitivity -------------------------------------------
    # The regimen chosen depends on how a percentage point of exposure-screen
    # exceedance is traded against a percentage point of joint attainment. The
    # manuscript makes that trade implicitly through a 15% ceiling its own source
    # calls arbitrary. Making it explicit shows how much rides on it.
    print("\n  Exchange-rate (lambda) sensitivity")
    lam_rows = []
    for lam in (0.0, 0.5, 1.0, 2.0, 4.0):
        r = run(n_draws, n_per_class, "C1_cojutti", "T2_point_4")
        pop = E.draw_population(n_per_class, MASTER_SEED)
        base = E.evaluate(pop, replace(E.BASE, avi_target=4.0), r["regimens"])
        for cls, idx in by_class(r).items():
            u = np.array([utility(base[r["regimens"][j]], r["weights"], lam)
                          for j in idx])
            uc = np.array([utility_constrained(base[r["regimens"][j]], r["weights"])
                           for j in idx])
            lam_rows.append({
                "lambda": lam, "ekfc_class": cls,
                "chosen_regimen_net_benefit": r["regimens"][idx[int(np.argmax(u))]],
                "chosen_regimen_hard_15pct_ceiling": r["regimens"][idx[int(np.argmax(uc))]],
                "agree": "yes" if int(np.argmax(u)) == int(np.argmax(uc)) else "NO",
            })
            print(f"    lambda {lam:4.1f}  {cls:9}  net benefit picks "
                  f"{lam_rows[-1]['chosen_regimen_net_benefit']:4}   "
                  f"15% ceiling picks {lam_rows[-1]['chosen_regimen_hard_15pct_ceiling']:4}"
                  f"   {'' if lam_rows[-1]['agree']=='yes' else '<-- differ'}")
    write(lam_rows, "model2_lambda_sensitivity.csv")

    print("\n  Convergence check (three independent seeds)")
    write(convergence_check(n_per_class), "model2_convergence.csv")

    meta = [{"scenario_id": k, "layer": "target", "kind": v["kind"],
             "definition": v["label"], "note": v["note"]}
            for k, v in TARGET_SCENARIOS.items()]
    meta += [{"scenario_id": k, "layer": "clearance correlation", "kind": "scenario",
              "definition": v["label"], "note": v["note"]}
             for k, v in RHO_SCENARIOS.items()]
    meta += [{"scenario_id": k, "layer": "parameter uncertainty", "kind": "distribution",
              "definition": "", "note": v} for k, v in E.UNCERTAINTY_NOTES.items()]
    write(meta, "model2_scenario_register.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
