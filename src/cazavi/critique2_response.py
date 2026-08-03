"""Test the second round of objections, where computation can settle them.

1  Free versus total avibactam. The 4 mg/L value is a total concentration in a
   broth well. The model applies it to a free concentration. How much does that
   choice move the answer?

2  What does "limiting component" mean. Comparing marginal attainment does not
   establish which component causes joint failure in an individual, especially
   with correlated clearance. Three per-subject definitions are computed and
   compared.

3  What a second assay is worth, done properly. The point-prediction figure
   ignores assay error and reports accuracy alone. Here the avibactam
   concentration is predicted from its conditional distribution given the
   measured ceftazidime concentration, assay error is added, and sensitivity,
   specificity and predictive values are reported.

4  Whether the population CFR rests on synthetic weights. The class weights come
   from the source cohort's reported eGFR distribution, but renal function is
   uniform within each class. The alternative is a continuous draw from the
   reported distribution.

Usage:
    python critique2_response.py
"""

from __future__ import annotations

import os

import numpy as np

from cazavi_analyses import (
    DEFAULT_OUT as OUT, SELECTED_REGIMENS, _cholesky, compute_cfr,
    draw_population, evaluate, load_mic_distributions, write_csv,
)
from reproduce_primary_run import (
    AVI_CT, AVI_FRACTION, CAZ_FRACTION, CAZ_TARGET, CL0_AVI, CL0_CAZ,
    EKFC_CLASSES, EKFC_REF, EXP_AVI, EXP_CAZ, FU_AVI, FU_CAZ, MIC_GRID,
    N_PER_CLASS, OMEGA_AVI, OMEGA_CAZ, PRIMARY_SEED, REGIMENS, RHO,
)
from scope_extension_analyses import ICU_WEIGHTS

LICENSED_MAX = 10.0


def clearances(renal, z, rho=RHO):
    eta = z @ _cholesky(OMEGA_CAZ, OMEGA_AVI, rho).T
    return (CL0_CAZ * (renal / EKFC_REF) ** EXP_CAZ * np.exp(eta[:, 0]),
            CL0_AVI * (renal / EKFC_REF) ** EXP_AVI * np.exp(eta[:, 1]))


def css(regimen, cl_caz, cl_avi):
    _, dose_g, interval = REGIMENS[regimen]
    return (dose_g * 1000.0 * CAZ_FRACTION / interval / cl_caz,
            dose_g * 1000.0 * AVI_FRACTION / interval / cl_avi)


# ------------------------------------------------------------------ 1
def free_vs_total(pop):
    """The avibactam target read as a free concentration, and as a total one."""
    print("\n1. Avibactam target: free versus total concentration")
    rows = []
    for regimen in SELECTED_REGIMENS:
        cls = REGIMENS[regimen][0]
        renal, z = pop[cls]
        cl_caz, cl_avi = clearances(renal, z)
        c_caz, c_avi = css(regimen, cl_caz, cl_avi)
        avi_free = 100.0 * float(np.mean(c_avi * FU_AVI >= AVI_CT))
        avi_total = 100.0 * float(np.mean(c_avi >= AVI_CT))
        row = {"regimen": regimen, "ekfc_class": cls,
               "avi_attainment_free_pct": round(avi_free, 1),
               "avi_attainment_total_pct": round(avi_total, 1),
               "difference_pp": round(avi_total - avi_free, 1)}
        for mic in (4.0, 8.0):
            caz_ok = (c_caz * FU_CAZ) / mic >= CAZ_TARGET
            row[f"joint_free_mic{mic:.0f}"] = round(
                100.0 * float(np.mean(caz_ok & (c_avi * FU_AVI >= AVI_CT))), 1)
            row[f"joint_total_mic{mic:.0f}"] = round(
                100.0 * float(np.mean(caz_ok & (c_avi >= AVI_CT))), 1)
        rows.append(row)
    print(f"   {'reg':5}{'class':10}{'AVI free':>10}{'AVI total':>11}{'Δ pp':>8}"
          f"{'joint free 8':>14}{'joint total 8':>15}")
    for r in rows:
        print(f"   {r['regimen']:5}{r['ekfc_class']:10}{r['avi_attainment_free_pct']:9.1f}%"
              f"{r['avi_attainment_total_pct']:10.1f}%{r['difference_pp']:8.1f}"
              f"{r['joint_free_mic8']:13.1f}%{r['joint_total_mic8']:14.1f}%")
    d = [r["difference_pp"] for r in rows]
    print(f"   Reading 4 mg/L as a total rather than a free concentration raises "
          f"avibactam attainment by {min(d):.1f}–{max(d):.1f} pp.")
    write_csv(rows, os.path.join(OUT, "critique2_free_vs_total.csv"))
    return rows


# ------------------------------------------------------------------ 2
def limiting_definitions(pop):
    """Three per-subject definitions of which component limits attainment."""
    print("\n2. What 'limiting component' means, under three definitions")
    rows = []
    for regimen in SELECTED_REGIMENS:
        cls = REGIMENS[regimen][0]
        renal, z = pop[cls]
        cl_caz, cl_avi = clearances(renal, z)
        c_caz, c_avi = css(regimen, cl_caz, cl_avi)
        avi_ok = c_avi * FU_AVI >= AVI_CT
        for mic in MIC_GRID:
            caz_ok = (c_caz * FU_CAZ) / mic >= CAZ_TARGET
            joint = caz_ok & avi_ok
            fail = ~joint
            n_fail = int(fail.sum())

            # (a) marginal: which component has the lower attainment
            marg = ("avibactam" if avi_ok.mean() < caz_ok.mean() - 0.005
                    else "ceftazidime" if caz_ok.mean() < avi_ok.mean() - 0.005
                    else "comparable")

            # (b) attribution: among subjects who fail, which target they miss
            only_caz = int((fail & ~caz_ok & avi_ok).sum())
            only_avi = int((fail & caz_ok & ~avi_ok).sum())
            both = int((fail & ~caz_ok & ~avi_ok).sum())
            attrib = ("ceftazidime" if only_caz > only_avi else
                      "avibactam" if only_avi > only_caz else "comparable") \
                if n_fail else "none fail"

            # (c) dose: for each subject, which component needs the larger dose
            need_caz = CAZ_TARGET * mic / FU_CAZ * cl_caz / CAZ_FRACTION
            need_avi = AVI_CT / FU_AVI * cl_avi / AVI_FRACTION
            caz_binds = float(np.mean(need_caz > need_avi))
            dose_def = ("ceftazidime" if caz_binds > 0.5 else "avibactam")

            rows.append({
                "regimen": regimen, "ekfc_class": cls, "mic_mg_l": mic,
                "joint_pta_pct": round(100.0 * joint.mean(), 1),
                "def_a_marginal": marg,
                "def_b_attribution": attrib,
                "fail_ceftazidime_only_pct": round(100.0 * only_caz / max(n_fail, 1), 1),
                "fail_avibactam_only_pct": round(100.0 * only_avi / max(n_fail, 1), 1),
                "fail_both_pct": round(100.0 * both / max(n_fail, 1), 1),
                "def_c_dose": dose_def,
                "ceftazidime_needs_more_dose_pct": round(100.0 * caz_binds, 1),
                "definitions_agree": len({marg, attrib, dose_def} - {"comparable",
                                                                    "none fail"}) == 1,
            })
    print(f"   {'MIC':>6}  {'(a) marginal':>13}  {'(b) attribution':>16}  {'(c) dose':>11}"
          f"  {'agree':>6}")
    for mic in MIC_GRID:
        sel = [r for r in rows if r["mic_mg_l"] == mic]
        a = {r["def_a_marginal"] for r in sel}
        b = {r["def_b_attribution"] for r in sel}
        c = {r["def_c_dose"] for r in sel}
        agree = sum(1 for r in sel if r["definitions_agree"])
        print(f"   {mic:6g}  {'/'.join(sorted(a)):>13}  {'/'.join(sorted(b)):>16}"
              f"  {'/'.join(sorted(c)):>11}  {agree}/{len(sel)}")
    n_agree = sum(1 for r in rows if r["definitions_agree"])
    print(f"   The three definitions agree in {n_agree} of {len(rows)} regimen-MIC cells.")
    write_csv(rows, os.path.join(OUT, "critique2_limiting_definitions.csv"))
    return rows


# ------------------------------------------------------------------ 3
def second_assay_proper(seed=PRIMARY_SEED, assay_cv=(0.0, 0.10, 0.20)):
    """Conditional prediction with assay error, reported as a classifier."""
    print("\n3. What a second assay is worth, with assay error and proper operating "
          "characteristics")
    rng = np.random.default_rng(seed + 31)
    rows = []
    for rho in (0.94, 0.75, 0.50, 0.0):
        pop = draw_population(N_PER_CLASS, seed)
        for cv in assay_cv:
            tp = fp = tn = fn = 0
            for regimen in SELECTED_REGIMENS:
                cls = REGIMENS[regimen][0]
                renal, z = pop[cls]
                cl_caz, cl_avi = clearances(renal, z, rho)
                c_caz, c_avi = css(regimen, cl_caz, cl_avi)
                truth = c_avi * FU_AVI >= AVI_CT

                # what the laboratory reports for ceftazidime
                obs = c_caz * np.exp(rng.normal(0.0, np.sqrt(np.log(1 + cv ** 2)),
                                                c_caz.size)) if cv > 0 else c_caz

                # eta_caz implied by the observed concentration, then the
                # conditional distribution of eta_avi given it
                typ_caz = CL0_CAZ * (renal / EKFC_REF) ** EXP_CAZ
                _, dose_g, interval = REGIMENS[regimen]
                cl_caz_hat = dose_g * 1000.0 * CAZ_FRACTION / interval / obs
                eta_caz_hat = np.log(cl_caz_hat / typ_caz)
                mu = rho * (OMEGA_AVI / OMEGA_CAZ) * eta_caz_hat
                sd = OMEGA_AVI * np.sqrt(max(1.0 - rho ** 2, 1e-12))

                typ_avi = CL0_AVI * (renal / EKFC_REF) ** EXP_AVI
                rate_avi = dose_g * 1000.0 * AVI_FRACTION / interval
                # attain  <=>  cl_avi <= rate*fu/AVI_CT  <=>  eta_avi <= log(...)
                thresh = np.log((rate_avi * FU_AVI / AVI_CT) / typ_avi)
                from math import erf, sqrt
                zscore = (thresh - mu) / sd
                p_attain = 0.5 * (1.0 + np.vectorize(erf)(zscore / sqrt(2.0)))
                pred = p_attain >= 0.5

                tp += int(np.sum(pred & truth)); fp += int(np.sum(pred & ~truth))
                tn += int(np.sum(~pred & ~truth)); fn += int(np.sum(~pred & truth))
            n = tp + fp + tn + fn
            rows.append({
                "rho": rho, "assay_cv_pct": round(100 * cv, 0),
                "accuracy_pct": round(100.0 * (tp + tn) / n, 1),
                "sensitivity_pct": round(100.0 * tp / max(tp + fn, 1), 1),
                "specificity_pct": round(100.0 * tn / max(tn + fp, 1), 1),
                "ppv_pct": round(100.0 * tp / max(tp + fp, 1), 1),
                "npv_pct": round(100.0 * tn / max(tn + fn, 1), 1),
                "false_reassurance_pct": round(100.0 * fp / n, 1),
                "prevalence_attaining_pct": round(100.0 * (tp + fn) / n, 1)})
    print(f"   {'rho':>5}{'assay CV':>10}{'accuracy':>10}{'sens':>8}{'spec':>8}"
          f"{'PPV':>8}{'NPV':>8}{'false reassurance':>19}")
    for r in rows:
        print(f"   {r['rho']:5.2f}{r['assay_cv_pct']:9.0f}%{r['accuracy_pct']:9.1f}%"
              f"{r['sensitivity_pct']:7.1f}%{r['specificity_pct']:7.1f}%"
              f"{r['ppv_pct']:7.1f}%{r['npv_pct']:7.1f}%"
              f"{r['false_reassurance_pct']:18.1f}%")
    write_csv(rows, os.path.join(OUT, "critique2_second_assay_operating.csv"))
    return rows


# ------------------------------------------------------------------ 4
def population_weighting(dists, seed=PRIMARY_SEED):
    """Uniform-within-class against a continuous draw from the reported cohort."""
    print("\n4. Does the population CFR depend on the synthetic within-class draw?")
    base = evaluate(draw_population(N_PER_CLASS, seed))
    cfr = compute_cfr([r for r in base if r["regimen"] in SELECTED_REGIMENS], dists)

    # continuous alternative: lognormal eGFR matched to the reported median 92
    # and interquartile range 50-113, truncated to 5-150
    rng = np.random.default_rng(seed + 77)
    mu, sigma = np.log(92.0), (np.log(113.0) - np.log(50.0)) / (2 * 0.6745)
    draw = np.clip(rng.lognormal(mu, sigma, 500_000), 5.0, 150.0)
    edges = [0, 30, 60, 90, 120, 150]
    names = list(EKFC_CLASSES)
    cont = {}
    for i, name in enumerate(names):
        cont[name] = float(np.mean((draw > edges[i]) & (draw <= edges[i + 1])))
    total = sum(cont.values())
    cont = {k: v / total for k, v in cont.items()}

    rows = []
    for did in sorted({c["distribution_id"] for c in cfr}):
        out = {"distribution_id": did}
        for label, weights in (("source cohort (used)", ICU_WEIGHTS),
                               ("continuous lognormal", cont)):
            j = sum(weights[c["ekfc_class"]] * c["joint_cfr_pct"]
                    for c in cfr if c["distribution_id"] == did)
            out[label] = round(j, 1)
        out["difference_pp"] = round(out["continuous lognormal"]
                                     - out["source cohort (used)"], 1)
        rows.append(out)
    print("   class weights: " + "  ".join(
        f"{k} {ICU_WEIGHTS[k]:.3f}/{cont[k]:.3f}" for k in names))
    for r in rows:
        print(f"   {r['distribution_id']:28} source-cohort weights "
              f"{r['source cohort (used)']:5.1f}%   continuous draw "
              f"{r['continuous lognormal']:5.1f}%   Δ {r['difference_pp']:+.1f} pp")
    write_csv(rows, os.path.join(OUT, "critique2_population_weighting.csv"))
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    dists = load_mic_distributions()
    pop = draw_population(N_PER_CLASS, PRIMARY_SEED)
    free_vs_total(pop)
    limiting_definitions(pop)
    second_assay_proper()
    population_weighting(dists)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
