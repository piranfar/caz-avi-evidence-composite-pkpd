"""Test the eight objections raised against the v8 manuscript.

Each is either confirmed, partly confirmed, or refuted by computation rather
than by argument. Nothing here is written into the manuscript; this establishes
what is true first.

A  Is the therapeutic window a simulation result or an identity?
B  Does restricting the cohort to EKFC >= 15 overturn "the ceiling binds first"?
C  Is "the ceiling binds first" a property of the drug or of the variance?
D  Over what range of the exposure ceiling and the exceedance cut-off does the
   ordering survive?
E  Under rho = 0.94, how much does the second assay add?
F  What does avibactam do in renal impairment, given it carries no ceiling?

Usage:
    python critique_response.py
"""

from __future__ import annotations

import os

import numpy as np

from cazavi_analyses import DEFAULT_OUT as OUT, SELECTED_REGIMENS, _cholesky, write_csv
from reproduce_primary_run import (
    AVI_CT, AVI_FRACTION, CAZ_FRACTION, CAZ_TARGET, CL0_AVI, CL0_CAZ, EKFC_CLASSES,
    EKFC_REF, EXP_AVI, EXP_CAZ, FU_AVI, FU_CAZ, N_PER_CLASS, OMEGA_AVI, OMEGA_CAZ,
    PRIMARY_SEED, REGIMENS, RHO, TOX_THRESHOLD,
)
from scope_extension_analyses import ICU_WEIGHTS

TOX_CEILING = 15.0
DOSE_GRID = np.arange(0.625, 40.01, 0.125)
PTA_TARGET = 90.0
RRT_FLOOR = 15.0


def draw(classes, n, seed, omega_caz=OMEGA_CAZ, omega_avi=OMEGA_AVI, rho=RHO):
    rng = np.random.default_rng(seed)
    l = _cholesky(omega_caz, omega_avi, rho)
    pop = {}
    for cls, (lo, hi) in classes.items():
        renal = rng.uniform(lo, hi, n)
        z = rng.standard_normal((n, 2))
        eta = z @ l.T
        pop[cls] = (renal,
                    CL0_CAZ * (renal / EKFC_REF) ** EXP_CAZ * np.exp(eta[:, 0]),
                    CL0_AVI * (renal / EKFC_REF) ** EXP_AVI * np.exp(eta[:, 1]))
    return pop


def crossing(cl_caz, cl_avi, mic=8.0, tox=TOX_THRESHOLD, ceiling=TOX_CEILING,
             target=PTA_TARGET):
    """Daily dose reaching the attainment target, and the dose reaching the ceiling."""
    d_pta = d_ceil = None
    exceed_at_pta = None
    for daily in DOSE_GRID:
        css_c = daily * 1000.0 * CAZ_FRACTION / 24.0 / cl_caz
        css_a = daily * 1000.0 * AVI_FRACTION / 24.0 / cl_avi
        ex = 100.0 * float(np.mean(css_c > tox))
        j = 100.0 * float(np.mean(((css_c * FU_CAZ) / mic >= CAZ_TARGET)
                                  & (css_a * FU_AVI >= AVI_CT)))
        if d_ceil is None and ex > ceiling:
            d_ceil = float(daily)
        if d_pta is None and j >= target:
            d_pta, exceed_at_pta = float(daily), ex
        if d_pta is not None and d_ceil is not None:
            break
    return d_pta, d_ceil, exceed_at_pta


# --------------------------------------------------------------------------- A
def test_a():
    print("\nA. Is the therapeutic window a simulation result?")
    pop = draw(EKFC_CLASSES, N_PER_CLASS, PRIMARY_SEED)
    rows = []
    for cls, (renal, cl_caz, cl_avi) in pop.items():
        for mic in (4.0, 8.0, 16.0, 22.0, 23.0, 32.0):
            need = CAZ_TARGET * mic / FU_CAZ * cl_caz          # mg/h for efficacy
            cap = TOX_THRESHOLD * cl_caz                        # mg/h at the ceiling
            feasible = need <= cap
            rows.append({"ekfc_class": cls, "mic_mg_l": mic,
                         "placeable_pct": round(100.0 * float(feasible.mean()), 4),
                         "distinct_values": int(len(np.unique(feasible)))})
    closes = TOX_THRESHOLD * FU_CAZ / CAZ_TARGET
    degenerate = all(r["distinct_values"] == 1 for r in rows)
    print(f"   need_caz = {CAZ_TARGET}*MIC/{FU_CAZ}*CL ; cap_caz = {TOX_THRESHOLD}*CL")
    print("   clearance appears on both sides and cancels, leaving")
    print(f"     placeable  <=>  MIC <= {TOX_THRESHOLD}*{FU_CAZ}/{CAZ_TARGET} = {closes:.1f} mg/L")
    print(f"   every subject in every class therefore gives the same answer: "
          f"{'CONFIRMED' if degenerate else 'not confirmed'}")
    for r in rows:
        if r["ekfc_class"] == "0–30":
            print(f"     MIC {r['mic_mg_l']:5.1f}  placeable {r['placeable_pct']:6.1f}%  "
                  f"({r['distinct_values']} distinct subject-level value)")
    print("   VERDICT: the 100% figure is an identity, not a population statistic.")
    print("            What is NOT an identity: whether the placing dose is licensed.")
    lic = []
    for cls, (renal, cl_caz, cl_avi) in pop.items():
        need = CAZ_TARGET * 8.0 / FU_CAZ * cl_caz
        need_a = AVI_CT / FU_AVI * cl_avi
        daily = np.maximum(need / CAZ_FRACTION, need_a / AVI_FRACTION) * 24.0 / 1000.0
        lic.append((cls, 100.0 * float(np.mean(daily <= 10.0)),
                    float(np.median(daily)), int(len(np.unique(daily <= 10.0)))))
    for cls, pct, med, nd in lic:
        print(f"     {cls:9} within licensed at MIC 8: {pct:5.1f}%  "
              f"median {med:5.2f} g/day  ({nd} distinct values -> genuine)")
    write_csv(rows, os.path.join(OUT, "critique_a_window_identity.csv"))
    return degenerate


# --------------------------------------------------------------------------- B
def test_b():
    print("\nB. Does excluding EKFC < 15 (who would be on RRT) overturn the ordering?")
    full = draw(EKFC_CLASSES, N_PER_CLASS, PRIMARY_SEED)
    classes_rrt = dict(EKFC_CLASSES)
    classes_rrt["0–30"] = (RRT_FLOOR, 30.0)
    restricted = draw(classes_rrt, N_PER_CLASS, PRIMARY_SEED)

    rows = []
    for label, pop in (("all EKFC", full), ("EKFC >= 15", restricted)):
        for cls, (renal, cl_caz, cl_avi) in pop.items():
            d_pta, d_ceil, ex = crossing(cl_caz, cl_avi)
            # exceedance at the licensed regimen for this class
            reg = next(r for r in SELECTED_REGIMENS if REGIMENS[r][0] == cls)
            _, dose_g, interval = REGIMENS[reg]
            css_c = dose_g * 1000.0 * CAZ_FRACTION / interval / cl_caz
            rows.append({"cohort": label, "ekfc_class": cls,
                         "licensed_exceedance_pct": round(100.0 * float(np.mean(css_c > TOX_THRESHOLD)), 2),
                         "dose_for_90pct_pta_g": d_pta if d_pta else "not reached",
                         "dose_at_ceiling_g": d_ceil if d_ceil else "not reached",
                         "exceedance_at_pta_dose_pct": round(ex, 1) if ex else "",
                         "ceiling_first": "yes" if (d_ceil and (not d_pta or d_ceil < d_pta)) else "no"})
    print(f"   {'cohort':12}{'class':9}{'licensed exceed':>17}{'dose 90%':>11}"
          f"{'dose ceiling':>14}{'exceed there':>14}{'ceiling first':>15}")
    for r in rows:
        print(f"   {r['cohort']:12}{r['ekfc_class']:9}{r['licensed_exceedance_pct']:16.2f}%"
              f"{str(r['dose_for_90pct_pta_g']):>11}{str(r['dose_at_ceiling_g']):>14}"
              f"{str(r['exceedance_at_pta_dose_pct']):>13}%{r['ceiling_first']:>15}")
    for label in ("all EKFC", "EKFC >= 15"):
        w = sum(ICU_WEIGHTS[r["ekfc_class"]] * r["licensed_exceedance_pct"]
                for r in rows if r["cohort"] == label)
        n_first = sum(1 for r in rows if r["cohort"] == label and r["ceiling_first"] == "yes")
        print(f"   {label:12} population exceedance {w:5.2f}%   ceiling binds first in "
              f"{n_first}/5 classes")
    write_csv(rows, os.path.join(OUT, "critique_b_rrt_exclusion.csv"))
    return rows


# --------------------------------------------------------------------------- C
def test_c():
    print("\nC. Is the ordering a drug property or a variance property?")
    rows = []
    for cv in (0.10, 0.20, 0.30, 0.45, 0.6792):
        om_c = float(np.sqrt(np.log(1 + cv ** 2)))
        om_a = float(np.sqrt(np.log(1 + (cv * 0.7691 / 0.6792) ** 2)))
        pop = draw(EKFC_CLASSES, N_PER_CLASS, PRIMARY_SEED, om_c, om_a)
        first = 0
        exs = []
        for cls, (renal, cl_caz, cl_avi) in pop.items():
            d_pta, d_ceil, ex = crossing(cl_caz, cl_avi)
            if d_ceil and (not d_pta or d_ceil < d_pta):
                first += 1
            if ex is not None:
                exs.append(ex)
        rows.append({"cv_caz": cv, "ceiling_first_classes": first,
                     "exceedance_at_pta_dose_min": round(min(exs), 1) if exs else "",
                     "exceedance_at_pta_dose_max": round(max(exs), 1) if exs else ""})
        print(f"   CV {cv*100:5.1f}%  ceiling binds first in {first}/5 classes   "
              f"exceedance at the attaining dose "
              f"{min(exs):5.1f}–{max(exs):5.1f}%" if exs else "")
    print("   VERDICT: if the ordering disappears at low CV it is a variance result,")
    print("            not a statement about ceftazidime–avibactam.")
    write_csv(rows, os.path.join(OUT, "critique_c_variance.csv"))
    return rows


# --------------------------------------------------------------------------- D
def test_d():
    print("\nD. How far can the ceiling and the cut-off move before the ordering fails?")
    pop = draw(EKFC_CLASSES, N_PER_CLASS, PRIMARY_SEED)
    rows = []
    for tox in (80.0, 104.0, 150.0, 200.0, 300.0, 500.0):
        for cut in (5.0, 15.0, 30.0):
            first = 0
            for cls, (renal, cl_caz, cl_avi) in pop.items():
                d_pta, d_ceil, _ = crossing(cl_caz, cl_avi, tox=tox, ceiling=cut)
                if d_ceil and (not d_pta or d_ceil < d_pta):
                    first += 1
            rows.append({"exposure_ceiling_mg_l": tox, "exceedance_cutoff_pct": cut,
                         "ceiling_first_classes": first})
    print(f"   {'ceiling':>9} | " + " | ".join(f"cut {c:.0f}%" for c in (5, 15, 30)))
    for tox in (80.0, 104.0, 150.0, 200.0, 300.0, 500.0):
        vals = [next(r["ceiling_first_classes"] for r in rows
                     if r["exposure_ceiling_mg_l"] == tox and r["exceedance_cutoff_pct"] == c)
                for c in (5.0, 15.0, 30.0)]
        print(f"   {tox:8.0f}  | " + " | ".join(f"  {v}/5  " for v in vals))
    write_csv(rows, os.path.join(OUT, "critique_d_ceiling_sensitivity.csv"))
    return rows


# --------------------------------------------------------------------------- E
def test_e():
    print("\nE. Under rho = 0.94, what does the second assay add?")
    rows = []
    for rho in (0.94, 0.75, 0.50, 0.0):
        pop = draw(EKFC_CLASSES, N_PER_CLASS, PRIMARY_SEED, rho=rho)
        agree = miss = tot = 0
        ratios = []
        for cls, (renal, cl_caz, cl_avi) in pop.items():
            reg = next(r for r in SELECTED_REGIMENS if REGIMENS[r][0] == cls)
            _, dose_g, interval = REGIMENS[reg]
            css_c = dose_g * 1000.0 * CAZ_FRACTION / interval / cl_caz
            css_a = dose_g * 1000.0 * AVI_FRACTION / interval / cl_avi
            truth = css_a * FU_AVI >= AVI_CT
            # predict avibactam from the measured ceftazidime, using typical clearances
            typ_ratio = (CL0_CAZ * (renal / EKFC_REF) ** EXP_CAZ) / \
                        (CL0_AVI * (renal / EKFC_REF) ** EXP_AVI)
            pred_a = css_c * (AVI_FRACTION / CAZ_FRACTION) * typ_ratio
            pred = pred_a * FU_AVI >= AVI_CT
            agree += int(np.sum(pred == truth)); miss += int(np.sum(pred != truth))
            tot += truth.size
            ratios.append(css_a / css_c)
        r = np.concatenate(ratios)
        rows.append({"rho": rho,
                     "avibactam_status_misclassified_pct": round(100.0 * miss / tot, 2),
                     "avi_caz_ratio_median": round(float(np.median(r)), 4),
                     "avi_caz_ratio_cv_pct": round(100.0 * float(np.std(r) / np.mean(r)), 1)})
        print(f"   rho {rho:4.2f}  predicting avibactam status from measured ceftazidime "
              f"misclassifies {100.0*miss/tot:5.2f}% of subjects   "
              f"(AVI:CAZ Css ratio CV {100.0*float(np.std(r)/np.mean(r)):4.1f}%)")
    print("   VERDICT: a low misclassification rate at rho = 0.94 means the paper's own")
    print("            correlation makes the second assay largely predictable.")
    write_csv(rows, os.path.join(OUT, "critique_e_second_assay.csv"))
    return rows


# --------------------------------------------------------------------------- F
def test_f():
    print("\nF. Avibactam in renal impairment, given it carries no ceiling")
    pop = draw(EKFC_CLASSES, N_PER_CLASS, PRIMARY_SEED)
    rows = []
    for cls, (renal, cl_caz, cl_avi) in pop.items():
        reg = next(r for r in SELECTED_REGIMENS if REGIMENS[r][0] == cls)
        _, dose_g, interval = REGIMENS[reg]
        css_c = dose_g * 1000.0 * CAZ_FRACTION / interval / cl_caz
        css_a = dose_g * 1000.0 * AVI_FRACTION / interval / cl_avi
        rows.append({"ekfc_class": cls, "regimen": reg,
                     "median_css_caz": round(float(np.median(css_c)), 1),
                     "median_css_avi": round(float(np.median(css_a)), 1),
                     "p95_css_avi": round(float(np.quantile(css_a, 0.95)), 1),
                     "avi_over_caz_median": round(float(np.median(css_a / css_c)), 3)})
    print(f"   {'class':9}{'reg':6}{'median CAZ':>12}{'median AVI':>12}{'p95 AVI':>10}{'AVI/CAZ':>10}")
    for r in rows:
        print(f"   {r['ekfc_class']:9}{r['regimen']:6}{r['median_css_caz']:12.1f}"
              f"{r['median_css_avi']:12.1f}{r['p95_css_avi']:10.1f}{r['avi_over_caz_median']:10.3f}")
    print("   The licensed regimens are renally adjusted, which is why avibactam does not")
    print("   run away in the low classes; the asymmetry is that no published avibactam")
    print("   exposure limit exists to test against.")
    write_csv(rows, os.path.join(OUT, "critique_f_avibactam_exposure.csv"))
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    test_a(); test_b(); test_c(); test_d(); test_e(); test_f()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
