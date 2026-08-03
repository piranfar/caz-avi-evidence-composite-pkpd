"""Five extensions that make the result usable rather than only defensible.

Population mix. Attainment has been reported per renal class with 20,000
subjects in each, which is a design for comparing classes, not a description of
any ICU. Reweighting the same outputs to the renal-function distribution of the
source cohort gives one number for a population rather than five for strata.

Protein binding as a distribution. The unbound fractions have been fixed at 0.85
and 0.92. In critically ill patients albumin varies widely and protein binding
with it, and free concentration is what both targets are expressed in. Sampling
the fraction per subject rather than fixing it moves variability into the place
it belongs.

Augmented renal clearance. The highest class runs to 150 mL/min/1.73 m², so
augmented clearance is present but diluted by subjects who do not have it.
Defining it explicitly at >=130 isolates the subgroup where ceftazidime exposure
is lowest.

Lung penetration. The principal indication is nosocomial pneumonia and the whole
analysis is in plasma. Applying published epithelial lining fluid penetration
ratios asks whether the window that is narrow in plasma is open at all at the
site of infection.

Decision grid. Renal class against MIC, with the binding component and the
permissibility status in each cell.

Usage:
    python scope_extension_analyses.py
"""

from __future__ import annotations

import os

import numpy as np

from cazavi_analyses import (
    DEFAULT_OUT as OUT, Scenario, compute_cfr, draw_population, evaluate,
    load_mic_distributions, write_csv, SELECTED_REGIMENS, _cholesky,
)
from reproduce_primary_run import (
    AVI_CT, AVI_FRACTION, CAZ_FRACTION, CAZ_TARGET, CL0_AVI, CL0_CAZ,
    EKFC_CLASSES, EKFC_REF, EXP_AVI, EXP_CAZ, FU_AVI, FU_CAZ, MIC_GRID,
    N_PER_CLASS, OMEGA_AVI, OMEGA_CAZ, PRIMARY_SEED, REGIMENS, RHO,
    TOX_THRESHOLD,
)

# Renal-function distribution of the source ICU cohort, which reported eGFR
# 92 mL/min/1.73 m2 with an interquartile range of 50-113. Class weights are
# obtained by interpolating that cumulative distribution; they describe one
# published cohort and are not a general ICU distribution.
ICU_WEIGHTS = {"0–30": 0.150, "31–60": 0.160, "61–90": 0.178,
               "91–120": 0.309, "121–150": 0.203}

# Between-subject variability in protein binding, expressed as the range over
# which the unbound fraction is sampled. The bounds are those the sensitivity
# analysis already used; here they describe variation between patients rather
# than a scenario applied to everyone.
FU_CAZ_RANGE = (0.80, 0.90)
FU_AVI_RANGE = (0.87, 0.97)

# Augmented renal clearance, conventionally defined at or above
# 130 mL/min/1.73 m2.
ARC_FLOOR = 130.0

# Epithelial lining fluid penetration. The compartmental analysis of plasma and
# ELF in healthy volunteers reports 52% for ceftazidime and 42% for avibactam at
# efficacy-relevant concentrations; an earlier estimate averaged across the
# whole concentration range gives about 30% for both, retained as the
# conservative case.
ELF_SCENARIOS = {"central estimate": (0.52, 0.42),
                 "conservative": (0.30, 0.30)}


def population_weighted_cfr(dists):
    """One CFR for a population, from the same per-class outputs."""
    rows = evaluate(draw_population(N_PER_CLASS, PRIMARY_SEED))
    cfr = compute_cfr([r for r in rows if r["regimen"] in SELECTED_REGIMENS], dists)
    out = []
    for did in sorted({c["distribution_id"] for c in cfr}):
        joint = tox = caz = 0.0
        for c in cfr:
            if c["distribution_id"] != did:
                continue
            w = ICU_WEIGHTS[c["ekfc_class"]]
            joint += w * c["joint_cfr_pct"]
            caz += w * c["caz_cfr_pct"]
            tox += w * c["toxicity_pct"]
        out.append({"distribution_id": did, "weighting": "source cohort",
                    "population_joint_cfr_pct": round(joint, 1),
                    "population_caz_cfr_pct": round(caz, 1),
                    "population_exceedance_pct": round(tox, 1)})
    return out


def variable_protein_binding(seed=PRIMARY_SEED):
    """Attainment with the unbound fraction sampled per subject."""
    rng = np.random.default_rng(seed + 7)
    base = draw_population(N_PER_CLASS, seed)
    rows = []
    for regimen, (cls, dose_g, interval_h) in REGIMENS.items():
        if regimen not in SELECTED_REGIMENS:
            continue
        ekfc, z = base[cls]
        l = _cholesky(OMEGA_CAZ, OMEGA_AVI, RHO)
        eta = z @ l.T
        cl_caz = CL0_CAZ * (ekfc / EKFC_REF) ** EXP_CAZ * np.exp(eta[:, 0])
        cl_avi = CL0_AVI * (ekfc / EKFC_REF) ** EXP_AVI * np.exp(eta[:, 1])
        css_caz = dose_g * 1000.0 * CAZ_FRACTION / interval_h / cl_caz
        css_avi = dose_g * 1000.0 * AVI_FRACTION / interval_h / cl_avi
        for label, fu_c, fu_a in (
                ("fixed", np.full(len(ekfc), FU_CAZ), np.full(len(ekfc), FU_AVI)),
                ("sampled", rng.uniform(*FU_CAZ_RANGE, len(ekfc)),
                 rng.uniform(*FU_AVI_RANGE, len(ekfc)))):
            avi_ok = css_avi * fu_a >= AVI_CT
            for mic in (2, 4, 8):
                caz_ok = (css_caz * fu_c) / mic >= CAZ_TARGET
                rows.append({"protein_binding": label, "regimen": regimen,
                             "ekfc_class": cls, "mic_mg_l": mic,
                             "caz_pta_pct": round(100 * caz_ok.mean(), 1),
                             "avi_attainment_pct": round(100 * avi_ok.mean(), 1),
                             "joint_pta_pct": round(100 * (caz_ok & avi_ok).mean(), 1)})
    return rows


def arc_subgroup(seed=PRIMARY_SEED):
    """The augmented-clearance subgroup, isolated from the highest class."""
    rng = np.random.default_rng(seed)
    cov = np.array([[OMEGA_CAZ ** 2, RHO * OMEGA_CAZ * OMEGA_AVI],
                    [RHO * OMEGA_CAZ * OMEGA_AVI, OMEGA_AVI ** 2]])
    rows = []
    for cls, (lo, hi) in EKFC_CLASSES.items():
        ekfc = rng.uniform(lo, hi, N_PER_CLASS)
        eta = rng.multivariate_normal(np.zeros(2), cov, size=N_PER_CLASS)
        if cls != "121–150":
            continue
        cl_caz = CL0_CAZ * (ekfc / EKFC_REF) ** EXP_CAZ * np.exp(eta[:, 0])
        cl_avi = CL0_AVI * (ekfc / EKFC_REF) ** EXP_AVI * np.exp(eta[:, 1])
        dose_g, interval_h = 2.5, 6
        css_caz = dose_g * 1000.0 * CAZ_FRACTION / interval_h / cl_caz
        css_avi = dose_g * 1000.0 * AVI_FRACTION / interval_h / cl_avi
        for label, mask in (("whole class (121–150)", np.ones(len(ekfc), bool)),
                            (f"augmented clearance (>={ARC_FLOOR:g})", ekfc >= ARC_FLOOR)):
            avi_ok = (css_avi * FU_AVI)[mask] >= AVI_CT
            for mic in (2, 4, 8):
                caz_ok = ((css_caz * FU_CAZ)[mask]) / mic >= CAZ_TARGET
                rows.append({"subgroup": label, "n": int(mask.sum()), "mic_mg_l": mic,
                             "caz_pta_pct": round(100 * caz_ok.mean(), 1),
                             "avi_attainment_pct": round(100 * avi_ok.mean(), 1),
                             "joint_pta_pct": round(100 * (caz_ok & avi_ok).mean(), 1),
                             "median_css_caz": round(float(np.median(css_caz[mask])), 1)})
    return rows


def lung_penetration(dists):
    """Attainment at the site of infection, applying published ELF ratios."""
    population = draw_population(N_PER_CLASS, PRIMARY_SEED)
    rows = []
    for label, (pc, pa) in [("plasma", (1.0, 1.0))] + list(ELF_SCENARIOS.items()):
        sc = Scenario(name=label, fu_caz=FU_CAZ * pc, fu_avi=FU_AVI * pa)
        evaluated = [r for r in evaluate(population, sc)
                     if r["regimen"] in SELECTED_REGIMENS]
        for c in compute_cfr(evaluated, dists):
            if c["distribution_id"] != "LEE2022_KPC_KP":
                continue
            at = {r["mic_mg_l"]: r for r in evaluated if r["regimen"] == c["regimen"]}
            rows.append({"compartment": label, "regimen": c["regimen"],
                         "ekfc_class": c["ekfc_class"],
                         "joint_pta_mic4_pct": round(at[4]["joint_pta_pct"], 1),
                         "joint_pta_mic8_pct": round(at[8]["joint_pta_pct"], 1),
                         "joint_cfr_pct": round(c["joint_cfr_pct"], 1)})
    return rows


def decision_grid():
    """Renal class against MIC: what happens, and which component decides."""
    rows = evaluate(draw_population(N_PER_CLASS, PRIMARY_SEED))
    out = []
    for regimen, (cls, _, _) in REGIMENS.items():
        if regimen not in SELECTED_REGIMENS:
            continue
        sel = {r["mic_mg_l"]: r for r in rows if r["regimen"] == regimen}
        for mic in MIC_GRID:
            r = sel[mic]
            j, tox = r["joint_pta_pct"], r["toxicity_pct"]
            caz, avi = r["caz_pta_pct"], r["avi_attainment_pct"]
            limiting = ("avibactam" if avi < caz - 0.5
                        else "ceftazidime" if caz < avi - 0.5 else "comparable")
            status = ("exceeds exposure ceiling" if tox > 15
                      else "attains target" if j >= 90
                      else "marginal" if j >= 80 else "fails target")
            out.append({"ekfc_class": cls, "regimen": regimen, "mic_mg_l": mic,
                        "joint_pta_pct": round(j, 1), "exceedance_pct": round(tox, 1),
                        "limiting_component": limiting, "status": status})
    return out


def main():
    dists = load_mic_distributions()

    pop = population_weighted_cfr(dists)
    write_csv(pop, os.path.join(OUT, "population_weighted_cfr.csv"))
    fu = variable_protein_binding()
    write_csv(fu, os.path.join(OUT, "variable_protein_binding.csv"))
    arc = arc_subgroup()
    write_csv(arc, os.path.join(OUT, "arc_subgroup.csv"))
    elf = lung_penetration(dists)
    write_csv(elf, os.path.join(OUT, "lung_penetration.csv"))
    grid = decision_grid()
    write_csv(grid, os.path.join(OUT, "decision_grid.csv"))

    print("\npopulation-weighted CFR (source-cohort renal mix)")
    for r in pop:
        print(f"  {r['distribution_id']:28} joint {r['population_joint_cfr_pct']:5.1f}%"
              f"   exceedance {r['population_exceedance_pct']:5.1f}%")

    print("\nprotein binding fixed vs sampled (joint PTA, MIC 4)")
    for reg in SELECTED_REGIMENS:
        v = {r["protein_binding"]: r["joint_pta_pct"] for r in fu
             if r["regimen"] == reg and r["mic_mg_l"] == 4}
        print(f"  {reg:5} fixed {v['fixed']:5.1f}%   sampled {v['sampled']:5.1f}%"
              f"   ({v['sampled'] - v['fixed']:+.1f})")

    print("\naugmented renal clearance")
    for r in arc:
        if r["mic_mg_l"] in (4, 8):
            print(f"  {r['subgroup']:34} MIC {r['mic_mg_l']}: joint PTA "
                  f"{r['joint_pta_pct']:5.1f}%   median Css {r['median_css_caz']:5.1f}")

    print("\nlung penetration (joint CFR, KPC-KP)")
    for comp in ["plasma"] + list(ELF_SCENARIOS):
        v = [r["joint_cfr_pct"] for r in elf if r["compartment"] == comp]
        w = [r["joint_pta_mic8_pct"] for r in elf if r["compartment"] == comp]
        print(f"  {comp:18} CFR {min(v):5.1f}–{max(v):5.1f}%   "
              f"joint PTA at MIC 8 {min(w):5.1f}–{max(w):5.1f}%")

    from collections import Counter
    print("\ndecision grid")
    print("  " + str(dict(Counter(r["status"] for r in grid))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
