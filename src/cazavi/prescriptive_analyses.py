"""Three prescriptive outputs: what dose, up to what MIC, and what to measure.

Everything the paper has established so far is a constraint. These three turn
the same model around and ask what it permits.

Safety-optimal dose. For each renal class, the daily dose that maximises joint
attainment subject to the exposure ceiling, rather than the licensed dose
applied uniformly. In the two highest classes this is higher than the licensed
maximum and gains attainment; in the three lowest it is far lower, because the
licensed maximum is not safe there.

Monotherapy ceiling. The highest MIC at which that dose still reaches target,
per renal class — the point at which the combination should give way to
something else.

Dose targets. For an individual whose concentration is measured, the efficacy
target and the exposure screen bound the concentration to aim at, and the dose
that reaches it follows from that subject's clearance. What is reported from
this is the dose itself and whether it lies within the licensed range, both of
which vary with clearance.

WITHDRAWN: earlier drafts also reported the proportion of subjects "placeable"
inside those bounds. That statistic is an identity, not a simulation result:
the concentration a subject needs and the concentration at which they breach the
screen are both proportional to their clearance, so clearance cancels and the
test reduces to MIC <= TOX_THRESHOLD * FU_CAZ / CAZ_TARGET. Every subject
returns the same answer and the proportion is 100% or 0% for the whole cohort at
once. The `in_window_any_dose_pct` columns below are retained only so the
withdrawal can be checked; see critique_response.py, test A.

Usage:
    python prescriptive_analyses.py
"""

from __future__ import annotations

import os

import numpy as np

from cazavi_analyses import (
    DEFAULT_OUT as OUT, SELECTED_REGIMENS, _cholesky, compute_cfr,
    draw_population, evaluate, load_mic_distributions, write_csv,
)
from reproduce_primary_run import (
    AVI_CT, AVI_FRACTION, CAZ_FRACTION, CAZ_TARGET, CL0_AVI, CL0_CAZ, EKFC_REF,
    EXP_AVI, EXP_CAZ, FU_AVI, FU_CAZ, MIC_GRID, N_PER_CLASS, OMEGA_AVI,
    OMEGA_CAZ, PRIMARY_SEED, RHO, TOX_THRESHOLD,
)

TOX_CEILING = 15.0
LICENSED_MAX_DAILY_G = 10.0          # product, i.e. 8 g ceftazidime + 2 g avibactam
DOSE_GRID_G = np.arange(0.625, 30.01, 0.125)


def clearances(ekfc, z):
    l = _cholesky(OMEGA_CAZ, OMEGA_AVI, RHO)
    eta = z @ l.T
    return (CL0_CAZ * (ekfc / EKFC_REF) ** EXP_CAZ * np.exp(eta[:, 0]),
            CL0_AVI * (ekfc / EKFC_REF) ** EXP_AVI * np.exp(eta[:, 1]))


def safety_optimal_dose(population):
    """The dose that maximises attainment without breaching the exposure ceiling."""
    rows = []
    for cls, (ekfc, z) in population.items():
        cl_caz, cl_avi = clearances(ekfc, z)
        best = {}
        for daily in DOSE_GRID_G:
            css_c = daily * 1000.0 * CAZ_FRACTION / 24.0 / cl_caz
            tox = 100.0 * np.mean(css_c > TOX_THRESHOLD)
            if tox > TOX_CEILING:
                continue
            css_a = daily * 1000.0 * AVI_FRACTION / 24.0 / cl_avi
            avi_ok = css_a * FU_AVI >= AVI_CT
            for mic in (4.0, 8.0):
                caz_ok = (css_c * FU_CAZ) / mic >= CAZ_TARGET
                j = 100.0 * np.mean(caz_ok & avi_ok)
                if mic not in best or j > best[mic][0]:
                    best[mic] = (j, daily, tox)
        for mic, (j, daily, tox) in sorted(best.items()):
            rows.append({"ekfc_class": cls, "mic_mg_l": mic,
                         "optimal_daily_dose_g": round(float(daily), 3),
                         "vs_licensed_max": ("above" if daily > LICENSED_MAX_DAILY_G
                                             else "at" if abs(daily - LICENSED_MAX_DAILY_G) < 0.2
                                             else "below"),
                         "joint_pta_pct": round(j, 1),
                         "exceedance_pct": round(tox, 1)})
    return rows


def monotherapy_ceiling(population, thresholds=(90.0, 80.0)):
    """Highest MIC at which the safety-optimal dose still reaches target."""
    rows = []
    for cls, (ekfc, z) in population.items():
        cl_caz, cl_avi = clearances(ekfc, z)
        for level in thresholds:
            ceiling = None
            for mic in MIC_GRID:
                reached = False
                for daily in DOSE_GRID_G:
                    css_c = daily * 1000.0 * CAZ_FRACTION / 24.0 / cl_caz
                    if 100.0 * np.mean(css_c > TOX_THRESHOLD) > TOX_CEILING:
                        continue
                    css_a = daily * 1000.0 * AVI_FRACTION / 24.0 / cl_avi
                    j = 100.0 * np.mean(((css_c * FU_CAZ) / mic >= CAZ_TARGET)
                                        & (css_a * FU_AVI >= AVI_CT))
                    if j >= level:
                        reached = True
                        break
                if reached:
                    ceiling = mic
            rows.append({"ekfc_class": cls, "attainment_level_pct": level,
                         "highest_mic_reached_mg_l": ceiling if ceiling else "none"})
    return rows


def therapeutic_window():
    """The concentration window an individual must sit in, and where it closes."""
    rows = []
    for mic in MIC_GRID:
        floor_total = CAZ_TARGET * mic / FU_CAZ          # total Css needed for efficacy
        rows.append({"mic_mg_l": mic,
                     "caz_css_floor_mg_l": round(floor_total, 1),
                     "caz_css_ceiling_mg_l": TOX_THRESHOLD,
                     "window_open": "yes" if floor_total <= TOX_THRESHOLD else "no",
                     "window_width_mg_l": round(max(TOX_THRESHOLD - floor_total, 0), 1),
                     "avi_css_floor_mg_l": round(AVI_CT / FU_AVI, 2)})
    closes_at = CAZ_TARGET / FU_CAZ
    return rows, TOX_THRESHOLD * FU_CAZ / CAZ_TARGET


def individualised_attainment(population):
    """Dose needed per subject, and whether it is licensed.

    `in_window_any_dose_pct` is an identity in the MIC (see module docstring)
    and is retained only for the withdrawal check; the informative columns are
    the required dose and the proportion within the licensed range.
    """
    rows = []
    for cls, (ekfc, z) in population.items():
        cl_caz, cl_avi = clearances(ekfc, z)
        for mic in (4.0, 8.0, 16.0):
            # rate each subject needs for efficacy, and the rate that would
            # breach the ceiling for that same subject
            need_caz = CAZ_TARGET * mic / FU_CAZ * cl_caz            # mg/h
            cap_caz = TOX_THRESHOLD * cl_caz
            need_avi = AVI_CT / FU_AVI * cl_avi
            feasible = need_caz <= cap_caz
            daily_caz = need_caz * 24.0 / 1000.0
            daily_avi = need_avi * 24.0 / 1000.0
            daily_product = np.maximum(daily_caz / CAZ_FRACTION, daily_avi / AVI_FRACTION)
            within_licensed = feasible & (daily_product <= LICENSED_MAX_DAILY_G)
            rows.append({
                "ekfc_class": cls, "mic_mg_l": mic,
                "in_window_any_dose_pct": round(100.0 * feasible.mean(), 1),
                "in_window_within_licensed_pct": round(100.0 * within_licensed.mean(), 1),
                "median_daily_product_needed_g": round(float(np.median(daily_product)), 2),
                "p90_daily_product_needed_g": round(float(np.quantile(daily_product, 0.9)), 2)})
    return rows


# Epithelial lining fluid penetration, as applied elsewhere in the analysis:
# the effective free concentration at the site is the total plasma steady-state
# concentration times the unbound fraction times the penetration ratio.
ELF_SCENARIOS = {"central estimate": (0.52, 0.42), "conservative": (0.30, 0.30)}


def lung_window():
    """The concentration bounds restated at the site of infection.

    The efficacy floor rises by the reciprocal of the penetration ratio; the
    exposure screen does not move, because it is a plasma quantity. As in
    plasma, where these bounds cross is arithmetic in the targets and the
    unbound fraction and contains no simulated quantity.
    """
    rows = []
    for label, (pc, pa) in [("plasma", (1.0, 1.0))] + list(ELF_SCENARIOS.items()):
        closes_at = TOX_THRESHOLD * FU_CAZ * pc / CAZ_TARGET
        rows.append({
            "compartment": label,
            "caz_penetration": pc, "avi_penetration": pa,
            "caz_css_floor_multiple_of_mic": round(CAZ_TARGET / (FU_CAZ * pc), 2),
            "caz_css_ceiling_mg_l": TOX_THRESHOLD,
            "avi_css_floor_mg_l": round(AVI_CT / (FU_AVI * pa), 2),
            "window_closes_at_mic_mg_l": round(closes_at, 1),
            "multiple_of_clinical_breakpoint": round(closes_at / 8.0, 2),
            "open_at_breakpoint": "yes" if closes_at >= 8.0 else "no",
        })
    return rows


def lung_individualised(population, mics=(4.0, 8.0)):
    """Proportion placeable inside the window at the site of infection."""
    rows = []
    for label, (pc, pa) in [("plasma", (1.0, 1.0))] + list(ELF_SCENARIOS.items()):
        for cls, (ekfc, z) in population.items():
            cl_caz, cl_avi = clearances(ekfc, z)
            for mic in mics:
                need_caz = CAZ_TARGET * mic / (FU_CAZ * pc) * cl_caz
                cap_caz = TOX_THRESHOLD * cl_caz
                need_avi = AVI_CT / (FU_AVI * pa) * cl_avi
                feasible = need_caz <= cap_caz
                daily = np.maximum(need_caz / CAZ_FRACTION, need_avi / AVI_FRACTION) * 24.0 / 1000.0
                rows.append({
                    "compartment": label, "ekfc_class": cls, "mic_mg_l": mic,
                    "in_window_any_dose_pct": round(100.0 * float(feasible.mean()), 1),
                    "in_window_within_licensed_pct":
                        round(100.0 * float((feasible & (daily <= LICENSED_MAX_DAILY_G)).mean()), 1),
                    "median_daily_product_needed_g": round(float(np.median(daily)), 2)})
    return rows


def prescriptive_grid(population):
    """Renal class against MIC: what the licensed dose does, and what to do.

    The earlier version of this grid reported only what the licensed regimen
    attains in each cell. That leaves the reader with a verdict and no action.
    This version adds the quantity a clinician actually needs — the daily dose
    that places a patient inside the window, and whether it lies within the
    licensed range — so each cell names a step rather than a shortfall.
    """
    licensed = {r["ekfc_class"]: r for r in evaluate(population)
                if r["regimen"] in SELECTED_REGIMENS}
    rows = []
    for cls, (ekfc, z) in population.items():
        cl_caz, cl_avi = clearances(ekfc, z)
        per_mic = {r["mic_mg_l"]: r for r in evaluate(population)
                   if r["ekfc_class"] == cls and r["regimen"] in SELECTED_REGIMENS}
        for mic in MIC_GRID:
            r = per_mic[mic]
            caz, avi = r["caz_pta_pct"], r["avi_attainment_pct"]
            limiting = ("avibactam" if avi < caz - 0.5
                        else "ceftazidime" if caz < avi - 0.5 else "comparable")
            need_caz = CAZ_TARGET * mic / FU_CAZ * cl_caz
            cap_caz = TOX_THRESHOLD * cl_caz
            need_avi = AVI_CT / FU_AVI * cl_avi
            feasible = need_caz <= cap_caz
            daily = np.maximum(need_caz / CAZ_FRACTION, need_avi / AVI_FRACTION) * 24.0 / 1000.0
            within = feasible & (daily <= LICENSED_MAX_DAILY_G)
            pct_within = 100.0 * float(within.mean())
            median_daily = float(np.median(daily))

            if r["joint_pta_pct"] >= 90.0 and r["toxicity_pct"] <= TOX_CEILING:
                action = "licensed regimen attains target"
            elif not feasible.all():
                action = "outside the window at any dose; use a different agent"
            elif pct_within >= 90.0:
                action = "measure and titrate within the licensed range"
            else:
                action = ("measure and titrate; the placing dose exceeds the licensed "
                          "maximum for most patients")
            rows.append({
                "ekfc_class": cls, "regimen": r["regimen"], "mic_mg_l": mic,
                "licensed_joint_pta_pct": round(r["joint_pta_pct"], 1),
                "licensed_exceedance_pct": round(r["toxicity_pct"], 1),
                "limiting_component": limiting,
                "in_window_any_dose_pct": round(100.0 * float(feasible.mean()), 1),
                "median_placing_dose_g_day": round(median_daily, 2),
                "p90_placing_dose_g_day": round(float(np.quantile(daily, 0.9)), 2),
                "placing_dose_within_licensed_pct": round(pct_within, 1),
                "action": action})
    return rows


def main():
    dists = load_mic_distributions()
    population = draw_population(N_PER_CLASS, PRIMARY_SEED)

    dose = safety_optimal_dose(population)
    write_csv(dose, os.path.join(OUT, "safety_optimal_dose.csv"))
    ceiling = monotherapy_ceiling(population)
    write_csv(ceiling, os.path.join(OUT, "monotherapy_ceiling.csv"))
    window, closes_at = therapeutic_window()
    write_csv(window, os.path.join(OUT, "therapeutic_window.csv"))
    indiv = individualised_attainment(population)
    write_csv(indiv, os.path.join(OUT, "individualised_attainment.csv"))
    lw = lung_window()
    write_csv(lw, os.path.join(OUT, "lung_therapeutic_window.csv"))
    li = lung_individualised(population)
    write_csv(li, os.path.join(OUT, "lung_individualised_attainment.csv"))
    grid = prescriptive_grid(population)
    write_csv(grid, os.path.join(OUT, "prescriptive_decision_grid.csv"))

    print("\nsafety-optimal daily dose (product, g/day)")
    print(f"  {'class':10}{'MIC':>5}{'dose':>8}{'vs licensed':>13}{'joint PTA':>11}{'exceed':>8}")
    for r in dose:
        print(f"  {r['ekfc_class']:10}{r['mic_mg_l']:5.0f}{r['optimal_daily_dose_g']:8.2f}"
              f"{r['vs_licensed_max']:>13}{r['joint_pta_pct']:10.1f}%{r['exceedance_pct']:7.1f}%")

    print("\nhighest MIC reached at the safety-optimal dose")
    for level in (90.0, 80.0):
        vals = [f"{r['ekfc_class']}: {r['highest_mic_reached_mg_l']}"
                for r in ceiling if r["attainment_level_pct"] == level]
        print(f"  at >={level:.0f}%   " + " | ".join(vals))

    print(f"\ntherapeutic window: total ceftazidime Css between "
          f"{CAZ_TARGET / FU_CAZ:.1f} x MIC and {TOX_THRESHOLD:.0f} mg/L")
    print(f"  avibactam floor {AVI_CT / FU_AVI:.2f} mg/L")
    print(f"  window closes at MIC {closes_at:.1f} mg/L "
          f"({closes_at / 8:.1f} times the clinical breakpoint)")

    print("\nproportion that can be placed in the window by adjusting dose")
    print(f"  {'class':10}{'MIC':>5}{'any dose':>11}{'within licensed':>17}{'median g/day':>14}")
    for r in indiv:
        print(f"  {r['ekfc_class']:10}{r['mic_mg_l']:5.0f}{r['in_window_any_dose_pct']:10.1f}%"
              f"{r['in_window_within_licensed_pct']:16.1f}%"
              f"{r['median_daily_product_needed_g']:14.2f}")

    print("\nthe same window at the site of infection")
    print(f"  {'compartment':20}{'floor':>10}{'closes at MIC':>15}{'vs breakpoint':>15}{'open there':>12}")
    for r in lw:
        print(f"  {r['compartment']:20}{r['caz_css_floor_multiple_of_mic']:8.2f}x"
              f"{r['window_closes_at_mic_mg_l']:14.1f}"
              f"{r['multiple_of_clinical_breakpoint']:14.2f}x"
              f"{r['open_at_breakpoint']:>12}")

    print("\nproportion placeable in the window at the site of infection")
    print(f"  {'compartment':20}{'MIC':>5}{'any dose':>16}{'within licensed':>18}")
    for comp in ["plasma"] + list(ELF_SCENARIOS):
        for mic in (4.0, 8.0):
            sel = [r for r in li if r["compartment"] == comp and r["mic_mg_l"] == mic]
            any_lo = min(r["in_window_any_dose_pct"] for r in sel)
            any_hi = max(r["in_window_any_dose_pct"] for r in sel)
            lic_lo = min(r["in_window_within_licensed_pct"] for r in sel)
            lic_hi = max(r["in_window_within_licensed_pct"] for r in sel)
            print(f"  {comp:20}{mic:5.0f}{f'{any_lo:.1f}-{any_hi:.1f}%':>16}"
                  f"{f'{lic_lo:.1f}-{lic_hi:.1f}%':>18}")

    print("\nprescriptive grid: what each cell tells a clinician to do")
    from collections import Counter
    for action, n in Counter(r["action"] for r in grid).most_common():
        print(f"  {n:3} of {len(grid)} cells   {action}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
