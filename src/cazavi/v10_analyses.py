"""Two additions for v10, and the figure the second-assay section now needs.

Penetration variability. Applying a single penetration ratio to every subject
treats a quantity measured with error as a constant. There is no published
between-subject distribution for it, so the ratio is instead drawn per subject
across the span of the published point estimates — 30% to 52% for ceftazidime
and 30% to 42% for avibactam — which is a scenario bounding the published
disagreement, not a measured distribution, and is labelled as such.

Second-assay figure. The previous figure reported misclassification alone.
Accuracy is the wrong summary here because attainment is common: the informative
quantities are specificity, which is what governs false reassurance, and how
both degrade with assay error.

Usage:
    python v10_analyses.py
"""

from __future__ import annotations

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cazavi_analyses import (
    DEFAULT_OUT as OUT, SELECTED_REGIMENS, _cholesky, compute_cfr,
    draw_population, load_mic_distributions, write_csv,
)
from reproduce_primary_run import (
    AVI_CT, AVI_FRACTION, CAZ_FRACTION, CAZ_TARGET, CL0_AVI, CL0_CAZ, EKFC_REF,
    EXP_AVI, EXP_CAZ, FU_AVI, FU_CAZ, MIC_GRID, N_PER_CLASS, OMEGA_AVI,
    OMEGA_CAZ, PRIMARY_SEED, REGIMENS, RHO,
)
from scope_extension_analyses import ICU_WEIGHTS

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
# in the repository, figures sit at the top level; standalone they sit beside the script
FIGDIR = (os.path.join(_REPO, "figures")
          if os.path.isdir(os.path.join(_REPO, "data", "processed"))
          else os.path.join(_HERE, "figures"))
os.makedirs(FIGDIR, exist_ok=True)
CAZ_PEN = (0.30, 0.52)
AVI_PEN = (0.30, 0.42)


def clearances(renal, z):
    eta = z @ _cholesky(OMEGA_CAZ, OMEGA_AVI, RHO).T
    return (CL0_CAZ * (renal / EKFC_REF) ** EXP_CAZ * np.exp(eta[:, 0]),
            CL0_AVI * (renal / EKFC_REF) ** EXP_AVI * np.exp(eta[:, 1]))


def penetration_variability(seed=PRIMARY_SEED):
    """Fixed penetration against a per-subject draw across the published span."""
    dists = load_mic_distributions()
    pop = draw_population(N_PER_CLASS, seed)
    rng = np.random.default_rng(seed + 101)
    scenarios = {
        "plasma": None,
        "fixed central (0.52 / 0.42)": (0.52, 0.42),
        "fixed conservative (0.30 / 0.30)": (0.30, 0.30),
        "drawn per subject (0.30-0.52 / 0.30-0.42)": "draw",
    }
    rows = []
    for label, spec in scenarios.items():
        per_regimen = []
        for regimen in SELECTED_REGIMENS:
            cls, dose_g, interval = REGIMENS[regimen]
            renal, z = pop[cls]
            cl_caz, cl_avi = clearances(renal, z)
            c_caz = dose_g * 1000.0 * CAZ_FRACTION / interval / cl_caz
            c_avi = dose_g * 1000.0 * AVI_FRACTION / interval / cl_avi
            if spec is None:
                pc = pa = 1.0
            elif spec == "draw":
                u = rng.uniform(size=c_caz.size)      # shared draw: a subject with
                pc = CAZ_PEN[0] + u * (CAZ_PEN[1] - CAZ_PEN[0])   # good ceftazidime
                pa = AVI_PEN[0] + u * (AVI_PEN[1] - AVI_PEN[0])   # penetration also
            else:                                                  # has good avibactam
                pc, pa = spec
            avi_ok = c_avi * FU_AVI * pa >= AVI_CT
            grid = {}
            for mic in MIC_GRID:
                caz_ok = (c_caz * FU_CAZ * pc) / mic >= CAZ_TARGET
                grid[mic] = {"regimen": regimen, "ekfc_class": cls,
                             "dose_g": dose_g, "interval_h": interval,
                             "mic_mg_l": mic,
                             "caz_pta_pct": 100.0 * float(caz_ok.mean()),
                             "avi_attainment_pct": 100.0 * float(avi_ok.mean()),
                             "joint_pta_pct": 100.0 * float((caz_ok & avi_ok).mean()),
                             "toxicity_pct": 0.0}
            per_regimen += list(grid.values())
        cfr = compute_cfr(per_regimen, dists)
        kpc = [c for c in cfr if c["distribution_id"] == "LEE2022_KPC_KP"]
        pop_cfr = sum(ICU_WEIGHTS[c["ekfc_class"]] * c["joint_cfr_pct"] for c in kpc)
        pta8 = [r["joint_pta_pct"] for r in per_regimen if r["mic_mg_l"] == 8]
        rows.append({"scenario": label,
                     "joint_pta_mic8_low_pct": round(min(pta8), 1),
                     "joint_pta_mic8_high_pct": round(max(pta8), 1),
                     "joint_cfr_low_pct": round(min(c["joint_cfr_pct"] for c in kpc), 1),
                     "joint_cfr_high_pct": round(max(c["joint_cfr_pct"] for c in kpc), 1),
                     "population_joint_cfr_pct": round(pop_cfr, 1)})
    print("\npenetration treated as fixed against drawn per subject")
    print(f"  {'scenario':44}{'joint PTA at MIC 8':>22}{'population CFR':>17}")
    for r in rows:
        span = f"{r['joint_pta_mic8_low_pct']:.1f}–{r['joint_pta_mic8_high_pct']:.1f}%"
        print(f"  {r['scenario']:44}{span:>22}"
              f"{r['population_joint_cfr_pct']:16.1f}%")
    write_csv(rows, os.path.join(OUT, "penetration_variability.csv"))
    return rows


def penetration_dependence(seed=PRIMARY_SEED):
    """Does the dependence imposed between the two penetration ratios matter?

    The per-subject draw above ties the two components together, so a subject
    with good ceftazidime penetration also has good avibactam penetration. That
    is an assumption, not an observation. Drawing them independently and in
    opposition bounds it.
    """
    dists = load_mic_distributions()
    pop = draw_population(N_PER_CLASS, seed)
    rows = []
    for label, mode in (("comonotonic (used)", "same"), ("independent", "indep"),
                        ("countermonotonic", "opposite")):
        rng = np.random.default_rng(seed + 101)
        rec = []
        for regimen in SELECTED_REGIMENS:
            cls, dose_g, interval = REGIMENS[regimen]
            renal, z = pop[cls]
            cl_caz, cl_avi = clearances(renal, z)
            c_caz = dose_g * 1000.0 * CAZ_FRACTION / interval / cl_caz
            c_avi = dose_g * 1000.0 * AVI_FRACTION / interval / cl_avi
            u = rng.uniform(size=c_caz.size)
            v = (u if mode == "same"
                 else 1.0 - u if mode == "opposite"
                 else rng.uniform(size=c_caz.size))
            pc = CAZ_PEN[0] + u * (CAZ_PEN[1] - CAZ_PEN[0])
            pa = AVI_PEN[0] + v * (AVI_PEN[1] - AVI_PEN[0])
            avi_ok = c_avi * FU_AVI * pa >= AVI_CT
            for mic in MIC_GRID:
                caz_ok = (c_caz * FU_CAZ * pc) / mic >= CAZ_TARGET
                rec.append({"regimen": regimen, "ekfc_class": cls, "dose_g": dose_g,
                            "interval_h": interval, "mic_mg_l": mic,
                            "caz_pta_pct": 100.0 * float(caz_ok.mean()),
                            "avi_attainment_pct": 100.0 * float(avi_ok.mean()),
                            "joint_pta_pct": 100.0 * float((caz_ok & avi_ok).mean()),
                            "toxicity_pct": 0.0})
        kpc = [c for c in compute_cfr(rec, dists)
               if c["distribution_id"] == "LEE2022_KPC_KP"]
        pw = sum(ICU_WEIGHTS[c["ekfc_class"]] * c["joint_cfr_pct"] for c in kpc)
        p8 = [r["joint_pta_pct"] for r in rec if r["mic_mg_l"] == 8]
        rows.append({"penetration_dependence": label,
                     "population_joint_cfr_pct": round(pw, 1),
                     "joint_pta_mic8_low_pct": round(min(p8), 1),
                     "joint_pta_mic8_high_pct": round(max(p8), 1)})
    print("\ndependence imposed between the two penetration ratios")
    for r in rows:
        print(f"  {r['penetration_dependence']:22} population CFR "
              f"{r['population_joint_cfr_pct']:5.1f}%   joint PTA at MIC 8 "
              f"{r['joint_pta_mic8_low_pct']:5.1f}-{r['joint_pta_mic8_high_pct']:.1f}%")
    write_csv(rows, os.path.join(OUT, "penetration_dependence.csv"))
    return rows


def second_assay_figure():
    with open(os.path.join(OUT, "critique2_second_assay_operating.csv")) as fh:
        rows = list(csv.DictReader(fh))
    rhos = sorted({float(r["rho"]) for r in rows})
    cvs = sorted({float(r["assay_cv_pct"]) for r in rows})

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    ax = axes[0]
    styles = {0.0: ("-", "o"), 10.0: ("--", "s"), 20.0: (":", "^")}
    for metric, colour, label in (("sensitivity_pct", "#1b3a6b", "Sensitivity"),
                                  ("specificity_pct", "#c1663c", "Specificity"),
                                  ("accuracy_pct", "#7a7a7a", "Accuracy")):
        for cv in cvs:
            y = [float(next(r for r in rows if float(r["rho"]) == p
                            and float(r["assay_cv_pct"]) == cv)[metric]) for p in rhos]
            ls, mk = styles[cv]
            ax.plot(rhos, y, ls, marker=mk, ms=4, lw=1.9, color=colour,
                    label=f"{label}, assay CV {cv:.0f}%" if cv == 0 else None)
    ax.axvline(0.94, color="#b03030", ls=":", lw=1.2)
    ax.text(0.925, 8, "ρ reported\nin the source model", fontsize=7.5,
            color="#b03030", ha="right")
    ax.set_xlabel("Assumed correlation between component clearances (ρ)")
    ax.set_ylabel("Operating characteristic (%)")
    ax.set_ylim(0, 104)
    ax.set_title("(a) Predicting avibactam attainment from one assay",
                 loc="left", fontsize=10.5)
    ax.legend(frameon=False, fontsize=8, loc="center left")
    ax.grid(alpha=0.25, lw=0.6)

    ax = axes[1]
    width = 0.25
    idx = np.arange(len(rhos))
    for k, cv in enumerate(cvs):
        y = [float(next(r for r in rows if float(r["rho"]) == p
                        and float(r["assay_cv_pct"]) == cv)["false_reassurance_pct"])
             for p in rhos]
        ax.bar(idx + (k - 1) * width, y, width,
               color=["#1b3a6b", "#3e8e7e", "#c1663c"][k], label=f"assay CV {cv:.0f}%")
    ax.set_xticks(idx)
    ax.set_xticklabels([f"ρ = {p:.2f}" for p in rhos])
    ax.set_ylabel("Patients wrongly reported as attaining (%)")
    ax.set_title("(b) False reassurance from a single assay", loc="left", fontsize=10.5)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25, lw=0.6, axis="y")

    fig.tight_layout()
    path = os.path.join(FIGDIR, "fig_second_assay.png")
    fig.savefig(path, dpi=300)
    print(f"\n  wrote {os.path.basename(path)}")


if __name__ == "__main__":
    penetration_variability()
    penetration_dependence()
    second_assay_figure()
