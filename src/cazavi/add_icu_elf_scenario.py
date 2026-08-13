"""Apply the 2026 ICU pneumonia ELF penetration ratios to the primary simulation.

The primary manuscript originally applied ELF/plasma ratios derived from healthy
volunteers (0.52/0.42) and an earlier conservative scenario (0.30/0.30). This
script adds a fixed, patient-derived scenario using the median ratios reported
by Benitez-Cano et al. (Critical Care 2026;30:305): 0.41 for ceftazidime and
0.44 for avibactam. The reported ratios are used as fixed scenario values, not
as a distribution of between-subject penetration.
"""

from __future__ import annotations

import csv
import os

import numpy as np

from reproduce_primary_run import (
    AVI_CT,
    AVI_FRACTION,
    CAZ_FRACTION,
    CAZ_TARGET,
    FU_AVI,
    FU_CAZ,
    EKFC_CLASSES,
    MIC_GRID,
    N_PER_CLASS,
    PRIMARY_SEED,
    REGIMENS,
    run_primary,
    sample_clearances,
)

HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(HERE))
# In the repository the analysis outputs live in data/processed; run standalone
# from a working folder they live in ./outputs beside the script.
OUT = (os.path.join(_REPO, "data", "processed")
       if os.path.isdir(os.path.join(_REPO, "data", "processed"))
       else os.path.join(HERE, "outputs"))
os.makedirs(OUT, exist_ok=True)
SELECTED_REGIMENS = ("R1", "R8", "R10", "R12", "R13")

# Full-precision class weights specified in Supplementary Table S3c.
ICU_WEIGHTS = dict(zip(
    EKFC_CLASSES,
    (0.150000, 0.159524, 0.178571, 0.309202, 0.202703),
))

# ATLAS 2020 KPC-producing K. pneumoniae distribution (379 isolates), as
# tabulated in Supplementary Table S3b.  The bins are aligned to MIC_GRID.
KPC_COUNTS = {
    0.0625: 5, 0.125: 2, 0.25: 17, 0.5: 66, 1: 149, 2: 102,
    4: 32, 8: 4, 16: 1, 32: 0, 64: 1,
}
KPC_WEIGHTS = {mic: count / 379 for mic, count in KPC_COUNTS.items()}


SCENARIOS = (
    ("plasma", "Plasma", 1.00, 1.00),
    ("icu_trial", "ELF, ICU pneumonia trial", 0.41, 0.44),
    ("healthy_volunteer", "ELF, healthy-volunteer estimate", 0.52, 0.42),
    ("conservative", "ELF, conservative scenario", 0.30, 0.30),
)


def _draw_clearance_population():
    """Recreate the frozen primary draw, in the original class order."""
    rng = np.random.default_rng(PRIMARY_SEED)
    return {
        cls: sample_clearances(rng, N_PER_CLASS, lo, hi)
        for cls, (lo, hi) in EKFC_CLASSES.items()
    }


def _evaluate_scenario(population, caz_ratio, avi_ratio):
    """Evaluate the selected regimens after applying fixed ELF ratios."""
    rows = []
    for regimen in SELECTED_REGIMENS:
        cls, dose_g, interval_h = REGIMENS[regimen]
        cl_caz, cl_avi = population[cls]
        css_caz = dose_g * 1000.0 * CAZ_FRACTION / interval_h / cl_caz
        css_avi = dose_g * 1000.0 * AVI_FRACTION / interval_h / cl_avi
        avi_ok = css_avi * FU_AVI * avi_ratio >= AVI_CT
        for mic in MIC_GRID:
            caz_ok = (css_caz * FU_CAZ * caz_ratio) / mic >= CAZ_TARGET
            rows.append({
                "regimen": regimen,
                "ekfc_class": cls,
                "mic_mg_l": mic,
                "joint_pta_pct": 100.0 * float((caz_ok & avi_ok).mean()),
            })
    return rows


def _cfr(rows, regimen):
    pta = {row["mic_mg_l"]: row["joint_pta_pct"]
           for row in rows if row["regimen"] == regimen}
    return sum(KPC_WEIGHTS[mic] * pta[mic] for mic in KPC_WEIGHTS)


def _write_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _summarize_variable_penetration(population):
    """Stress-test the fixed-ratio assumption over the updated evidence range.

    The ICU trial extends the avibactam upper bound from 0.42 to 0.44.  As in
    the original analysis, these are scenario draws spanning point estimates,
    not an asserted empirical distribution of penetration.
    """
    summaries = []
    for label, dependence in (
        ("comonotonic (used)", "same"),
        ("independent", "independent"),
        ("countermonotonic", "opposite"),
    ):
        rng = np.random.default_rng(PRIMARY_SEED + 101)
        rows = []
        for regimen in SELECTED_REGIMENS:
            cls, dose_g, interval_h = REGIMENS[regimen]
            cl_caz, cl_avi = population[cls]
            css_caz = dose_g * 1000.0 * CAZ_FRACTION / interval_h / cl_caz
            css_avi = dose_g * 1000.0 * AVI_FRACTION / interval_h / cl_avi
            u = rng.uniform(size=css_caz.size)
            v = (u if dependence == "same" else 1.0 - u if dependence == "opposite"
                 else rng.uniform(size=css_caz.size))
            caz_ratio = 0.30 + u * (0.52 - 0.30)
            avi_ratio = 0.30 + v * (0.44 - 0.30)
            avi_ok = css_avi * FU_AVI * avi_ratio >= AVI_CT
            for mic in MIC_GRID:
                caz_ok = (css_caz * FU_CAZ * caz_ratio) / mic >= CAZ_TARGET
                rows.append({
                    "regimen": regimen,
                    "ekfc_class": cls,
                    "mic_mg_l": mic,
                    "joint_pta_pct": 100.0 * float((caz_ok & avi_ok).mean()),
                })
        cfr_by_regimen = {regimen: _cfr(rows, regimen) for regimen in SELECTED_REGIMENS}
        pta8 = [row["joint_pta_pct"] for row in rows if row["mic_mg_l"] == 8]
        summaries.append({
            "penetration_dependence": label,
            "caz_ratio_range": "0.30-0.52",
            "avi_ratio_range": "0.30-0.44",
            "joint_pta_mic8_low_pct": round(min(pta8), 1),
            "joint_pta_mic8_high_pct": round(max(pta8), 1),
            "joint_cfr_low_pct": round(min(cfr_by_regimen.values()), 1),
            "joint_cfr_high_pct": round(max(cfr_by_regimen.values()), 1),
            "population_weighted_joint_cfr_pct": round(sum(
                ICU_WEIGHTS[next(row["ekfc_class"] for row in rows if row["regimen"] == regimen)]
                * cfr_by_regimen[regimen]
                for regimen in SELECTED_REGIMENS
            ), 1),
        })
    return summaries


def _individualized_icu_elf_dose(population):
    """Subject-specific product dose required in plasma and ICU-derived ELF."""
    rows = []
    for scenario, caz_ratio, avi_ratio in (
        ("plasma", 1.00, 1.00),
        ("icu_trial", 0.41, 0.44),
    ):
        for regimen in SELECTED_REGIMENS:
            cls, _, _ = REGIMENS[regimen]
            cl_caz, cl_avi = population[cls]
            for mic in (4.0, 8.0, 16.0):
                # Required infusion rate of product (mg/h) for each component.
                caz_rate = CAZ_TARGET * mic * cl_caz / (FU_CAZ * caz_ratio * CAZ_FRACTION)
                avi_rate = AVI_CT * cl_avi / (FU_AVI * avi_ratio * AVI_FRACTION)
                daily_product = np.maximum(caz_rate, avi_rate) * 24.0 / 1000.0
                rows.append({
                    "scenario": scenario,
                    "ekfc_class": cls,
                    "mic_mg_l": mic,
                    "median_daily_product_needed_g": round(float(np.median(daily_product)), 2),
                    "within_10_g_day_modelled_cap_pct": round(
                        100.0 * float(np.mean(daily_product <= 10.0)), 1),
                })
    return rows


def run():
    population = _draw_clearance_population()
    detail, summary = [], []

    for key, label, caz_ratio, avi_ratio in SCENARIOS:
        evaluated = _evaluate_scenario(population, caz_ratio, avi_ratio)

        weighted_cfr = 0.0
        for regimen in SELECTED_REGIMENS:
            regimen_rows = [row for row in evaluated if row["regimen"] == regimen]
            cfr = _cfr(evaluated, regimen)
            first = regimen_rows[0]
            attainment = {
                row["mic_mg_l"]: row
                for row in regimen_rows
            }
            detail.append({
                "scenario": key,
                "scenario_label": label,
                "caz_elf_to_plasma_ratio": caz_ratio,
                "avi_elf_to_plasma_ratio": avi_ratio,
                "regimen": regimen,
                "ekfc_class": first["ekfc_class"],
                "joint_pta_mic4_pct": round(attainment[4]["joint_pta_pct"], 1),
                "joint_pta_mic8_pct": round(attainment[8]["joint_pta_pct"], 1),
                "joint_cfr_pct": round(cfr, 1),
            })
            weighted_cfr += ICU_WEIGHTS[first["ekfc_class"]] * cfr

        selected = [row for row in detail if row["scenario"] == key]
        summary.append({
            "scenario": key,
            "scenario_label": label,
            "caz_elf_to_plasma_ratio": caz_ratio,
            "avi_elf_to_plasma_ratio": avi_ratio,
            "joint_pta_mic8_low_pct": round(min(row["joint_pta_mic8_pct"] for row in selected), 1),
            "joint_pta_mic8_high_pct": round(max(row["joint_pta_mic8_pct"] for row in selected), 1),
            "joint_cfr_low_pct": round(min(row["joint_cfr_pct"] for row in selected), 1),
            "joint_cfr_high_pct": round(max(row["joint_cfr_pct"] for row in selected), 1),
            "population_weighted_joint_cfr_pct": round(weighted_cfr, 1),
        })

    _write_csv(detail, os.path.join(OUT, "lung_penetration_icu_trial_detail.csv"))
    _write_csv(summary, os.path.join(OUT, "lung_penetration_icu_trial_summary.csv"))
    variable = _summarize_variable_penetration(population)
    _write_csv(variable, os.path.join(OUT, "lung_penetration_icu_trial_variable.csv"))
    dose = _individualized_icu_elf_dose(population)
    _write_csv(dose, os.path.join(OUT, "lung_penetration_icu_trial_dose.csv"))
    for row in summary:
        print(
            f"{row['scenario']:17} CAZ/AVI {row['caz_elf_to_plasma_ratio']:.2f}/"
            f"{row['avi_elf_to_plasma_ratio']:.2f}; MIC 8 PTA "
            f"{row['joint_pta_mic8_low_pct']:.1f}-{row['joint_pta_mic8_high_pct']:.1f}%; "
            f"CFR {row['joint_cfr_low_pct']:.1f}-{row['joint_cfr_high_pct']:.1f}%; "
            f"weighted CFR {row['population_weighted_joint_cfr_pct']:.1f}%"
        )
    for row in variable:
        print(
            f"{row['penetration_dependence']:17} variable CAZ/AVI "
            f"{row['caz_ratio_range']}/{row['avi_ratio_range']}; MIC 8 PTA "
            f"{row['joint_pta_mic8_low_pct']:.1f}-{row['joint_pta_mic8_high_pct']:.1f}%; "
            f"weighted CFR {row['population_weighted_joint_cfr_pct']:.1f}%"
        )
    print("ICU ELF individualized dose at MIC 8")
    for row in dose:
        if row["scenario"] == "icu_trial" and row["mic_mg_l"] == 8.0:
            print(f"  {row['ekfc_class']}: {row['median_daily_product_needed_g']:.2f} g/day; "
                  f"{row['within_10_g_day_modelled_cap_pct']:.1f}% at or below 10 g/day")


if __name__ == "__main__":
    run()
