"""What dose would have been needed, and what it would cost.

The paper currently reports that no licensed regimen reaches target at the
clinical breakpoint. That is a qualitative negative. These analyses turn it into
a number by asking three questions the existing model can answer.

1. Dose escalation. Sweep total daily dose well past the licensed maximum and
   find two crossing points at MIC 8 mg/L: the dose at which joint PTA reaches
   90%, and the dose at which exposure-screen exceedance passes 15%. Whichever
   comes first decides whether the breakpoint is reachable at any dose, and by
   how much the two constraints miss each other.

2. Resistance suppression. The targets used throughout are killing targets, not
   mutant-suppression targets. Re-running the grid against a stricter
   ceftazidime target approximates what suppression would demand.

3. Time to steady state. Exposure was computed algebraically, so the first day
   or two of therapy sits outside the model. The accumulation half-life follows
   from the clearance already in the model, so the size of that blind spot can
   be stated rather than merely acknowledged.

Usage:
    python dose_escalation_analyses.py
"""

from __future__ import annotations

import os

import numpy as np

from cazavi_analyses import (
    DEFAULT_OUT as OUT, Scenario, clearances, compute_cfr, draw_population,
    evaluate, load_mic_distributions, write_csv,
)
from reproduce_primary_run import (
    AVI_CT, AVI_FRACTION, CAZ_FRACTION, CAZ_TARGET, EKFC_CLASSES, FU_AVI, FU_CAZ,
    N_PER_CLASS, PRIMARY_SEED, TOX_THRESHOLD,
)

# Licensed maximum is 2.5 g of product every 6 h, i.e. 10 g/day. The sweep runs
# to three times that so the crossing points, if any, are bracketed.
DAILY_DOSES_G = [1.25, 1.875, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0,
                 17.5, 20.0, 22.5, 25.0, 27.5, 30.0]
LICENSED_MAX_G = 10.0
BREAKPOINT_MIC = 8.0
PTA_TARGET = 90.0
TOX_CEILING = 15.0

# Volume of distribution from the source model, used only for the accumulation
# half-life; it does not enter the steady-state concentration.
VSS_L = 18.0 + 18.1


def dose_sweep(population, mic=BREAKPOINT_MIC, avi_ct=AVI_CT, caz_target=CAZ_TARGET):
    """Joint PTA and exposure-screen exceedance against total daily dose."""
    rows = []
    for cls, (ekfc, z) in population.items():
        cl_caz, cl_avi = clearances(ekfc, z, Scenario())
        for daily in DAILY_DOSES_G:
            rate_caz = daily * 1000.0 * CAZ_FRACTION / 24.0
            rate_avi = daily * 1000.0 * AVI_FRACTION / 24.0
            css_caz, css_avi = rate_caz / cl_caz, rate_avi / cl_avi
            caz_ok = (css_caz * FU_CAZ) / mic >= caz_target
            avi_ok = css_avi * FU_AVI >= avi_ct
            rows.append({
                "ekfc_class": cls, "daily_dose_g": daily,
                "licensed": "yes" if daily <= LICENSED_MAX_G else "no",
                "mic_mg_l": mic,
                "caz_pta_pct": 100.0 * caz_ok.mean(),
                "avi_attainment_pct": 100.0 * avi_ok.mean(),
                "joint_pta_pct": 100.0 * (caz_ok & avi_ok).mean(),
                "exceedance_pct": 100.0 * np.mean(css_caz > TOX_THRESHOLD),
                "median_css_caz": float(np.median(css_caz)),
            })
    return rows


def crossings(rows):
    """Where each constraint is first met, by linear interpolation in dose."""
    out = []
    for cls in EKFC_CLASSES:
        sel = sorted((r for r in rows if r["ekfc_class"] == cls),
                     key=lambda r: r["daily_dose_g"])
        dose = [r["daily_dose_g"] for r in sel]
        pta = [r["joint_pta_pct"] for r in sel]
        tox = [r["exceedance_pct"] for r in sel]

        def cross(y, level, rising=True):
            for i in range(1, len(y)):
                a, b = y[i - 1], y[i]
                if (rising and a < level <= b) or (not rising and a > level >= b):
                    f = (level - a) / (b - a)
                    return dose[i - 1] + f * (dose[i] - dose[i - 1])
            return None

        d_pta = cross(pta, PTA_TARGET)
        d_tox = cross(tox, TOX_CEILING)
        out.append({
            "ekfc_class": cls,
            "dose_for_90pct_joint_pta_g": round(d_pta, 2) if d_pta else "not reached by 30 g",
            "dose_at_15pct_exceedance_g": round(d_tox, 2) if d_tox else "not reached by 30 g",
            "safety_limit_reached_first": ("yes" if d_tox and (not d_pta or d_tox < d_pta)
                                           else "no"),
            "exceedance_at_90pct_pta": (round(np.interp(d_pta, dose, tox), 1)
                                        if d_pta else None),
            "joint_pta_at_licensed_max": round(np.interp(LICENSED_MAX_G, dose, pta), 1),
            "exceedance_at_licensed_max": round(np.interp(LICENSED_MAX_G, dose, tox), 1),
        })
    return out


def resistance_suppression(dists, population):
    """A stricter ceftazidime target, as a proxy for mutant suppression."""
    rows = []
    for target, label in ((CAZ_TARGET, "killing (fCss/MIC >= 4)"),
                          (6.0, "intermediate (fCss/MIC >= 6)"),
                          (8.0, "suppression proxy (fCss/MIC >= 8)")):
        sc = Scenario(name=f"CAZ_TARGET_{target:g}", caz_target=target)
        evaluated = evaluate(population, sc)
        for c in compute_cfr(evaluated, dists):
            if c["distribution_id"] != "LEE2022_KPC_KP":
                continue
            at4 = next(r for r in evaluated
                       if r["regimen"] == c["regimen"] and r["mic_mg_l"] == 4)
            rows.append({"caz_target": target, "target_label": label,
                         "regimen": c["regimen"], "ekfc_class": c["ekfc_class"],
                         "joint_pta_mic4_pct": round(at4["joint_pta_pct"], 1),
                         "joint_cfr_pct": round(c["joint_cfr_pct"], 1),
                         "toxicity_pct": round(c["toxicity_pct"], 1)})
    return rows


def time_to_steady_state(population):
    """Accumulation half-life and time to 95% of steady state, by renal class."""
    rows = []
    for cls, (ekfc, z) in population.items():
        cl_caz, _ = clearances(ekfc, z, Scenario())
        for label, cl in (("median", float(np.median(cl_caz))),
                          ("5th percentile (slowest clearers)",
                           float(np.quantile(cl_caz, 0.05)))):
            half_life = np.log(2) * VSS_L / cl
            rows.append({"ekfc_class": cls, "clearance_stratum": label,
                         "clearance_l_h": round(cl, 3),
                         "half_life_h": round(half_life, 1),
                         "time_to_95pct_steady_state_h": round(4.32 * half_life, 1)})
    return rows


def main():
    dists = load_mic_distributions()
    population = draw_population(N_PER_CLASS, PRIMARY_SEED)

    sweep = dose_sweep(population)
    write_csv(sweep, os.path.join(OUT, "dose_escalation_mic8.csv"))
    cross = crossings(sweep)
    write_csv(cross, os.path.join(OUT, "dose_escalation_crossings.csv"))

    sweep4 = dose_sweep(population, mic=4.0)
    write_csv(sweep4, os.path.join(OUT, "dose_escalation_mic4.csv"))

    suppression = resistance_suppression(dists, population)
    write_csv(suppression, os.path.join(OUT, "resistance_suppression_targets.csv"))

    tss = time_to_steady_state(population)
    write_csv(tss, os.path.join(OUT, "time_to_steady_state.csv"))

    print("\ndose needed at the breakpoint (MIC 8 mg/L)")
    print(f"  {'class':10}{'90% joint PTA':>16}{'15% exceedance':>17}"
          f"{'safety first?':>15}{'exceed. at 90% PTA':>20}")
    for r in cross:
        print(f"  {r['ekfc_class']:10}{str(r['dose_for_90pct_joint_pta_g']):>16}"
              f"{str(r['dose_at_15pct_exceedance_g']):>17}"
              f"{r['safety_limit_reached_first']:>15}"
              f"{str(r['exceedance_at_90pct_pta']):>20}")

    print("\nat the licensed maximum (10 g/day)")
    for r in cross:
        print(f"  {r['ekfc_class']:10} joint PTA {r['joint_pta_at_licensed_max']:5.1f}%"
              f"   exceedance {r['exceedance_at_licensed_max']:5.1f}%")

    print("\nstricter ceftazidime target (KPC-KP, selected regimens)")
    for target in (4.0, 6.0, 8.0):
        v = [r["joint_cfr_pct"] for r in suppression if r["caz_target"] == target
             and r["regimen"] in ("R1", "R8", "R10", "R12", "R13")]
        print(f"  fCss/MIC >= {target:g}: joint CFR {min(v):5.1f}–{max(v):5.1f}%")

    print("\ntime to 95% of steady state")
    for r in tss:
        if r["clearance_stratum"].startswith("5th"):
            print(f"  {r['ekfc_class']:10} slowest clearers: t1/2 {r['half_life_h']:6.1f} h"
                  f"  ->  {r['time_to_95pct_steady_state_h']:6.1f} h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
