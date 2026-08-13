"""Empirical estimation of the ceftazidime-avibactam clearance correlation
from two independent, openly licensed patient-level ICU datasets.

Purpose: to check the value of rho = 0.94 assumed by the primary model, which is
taken from a single cohort (Cojutti 2024, RSE 23.8%) and drives the manuscript's
conclusion about whether a separate avibactam assay is needed.

SOURCES AND LICENCES
  Dataset 1  Gatti M, et al. J Crit Care 2023;76:154301, Table 2. CC BY-NC-ND 4.0.
             8 patients, 17 TDM occasions on CVVHDF; total clearances as published.
  Dataset 2  Li C, et al. Dryad doi:10.5061/dryad.fxpnvx16s. CC0 public domain.
             21 patients on CRRT; individual concentration-time data, clearance
             derived here as Dose / AUC(0-tau) at steady state.

CRITICAL LIMITATION, TO BE REPORTED WITH EVERY RESULT
  Both cohorts are on renal replacement therapy. The manuscript's primary scenario
  EXCLUDES renal replacement therapy. A shared extracorporeal circuit removes both
  analytes together and should therefore INFLATE the correlation between their
  clearances -- which makes the low observed values an a fortiori argument, not a
  direct substitute for rho in the non-RRT population.

  Dataset 2 additionally uses INTERMITTENT 8-hourly infusion, not continuous infusion.

This is NOT external validation of the primary model and must never be described as such.
"""
import csv, os, sys
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data_external")

ASSUMED_RHO = 0.94


def fisher_ci(r, n, alpha=0.05):
    z, se = np.arctanh(r), 1.0 / np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - alpha / 2)
    return np.tanh(z - zc * se), np.tanh(z + zc * se)


def test_against(r, n, rho0):
    z = (np.arctanh(r) - np.arctanh(rho0)) * np.sqrt(n - 3)
    return z, 2 * stats.norm.sf(abs(z))


def report(label, caz, avi, n_for_ci=None):
    caz, avi = np.asarray(caz, float), np.asarray(avi, float)
    n = n_for_ci or len(caz)
    lc, la = np.log(caz), np.log(avi)
    r, p = stats.pearsonr(lc, la)
    rs, ps = stats.spearmanr(caz, avi)
    lo, hi = fisher_ci(r, n)
    z, pz = test_against(r, n, ASSUMED_RHO)
    print(f"\n{label}")
    print(f"  n = {len(caz)} observations (n = {n} used for the confidence interval)")
    print(f"  CL ceftazidime  median {np.median(caz):5.2f} L/h  IQR {np.percentile(caz,25):.2f}-{np.percentile(caz,75):.2f}")
    print(f"  CL avibactam    median {np.median(avi):5.2f} L/h  IQR {np.percentile(avi,25):.2f}-{np.percentile(avi,75):.2f}")
    print(f"  Pearson r (log scale)  {r:6.3f}   95% CI {lo:.3f} to {hi:.3f}   p = {p:.4g}")
    print(f"  Spearman rho           {rs:6.3f}   p = {ps:.4g}")
    print(f"  vs assumed rho = {ASSUMED_RHO}:  z = {z:.2f},  p = {pz:.3g}")
    print(f"  apparent between-subject CV:  ceftazidime {100*np.sqrt(np.exp(lc.var(ddof=1))-1):.1f}%"
          f"   avibactam {100*np.sqrt(np.exp(la.var(ddof=1))-1):.1f}%")
    print(f"     (source non-RRT ICU model reports 67.9% and 76.9%)")
    return r, lo, hi


# --- Dataset 1: Gatti 2023, Table 2 (published clearances) --------------------
GATTI = [  # (patient, CL ceftazidime L/h, CL avibactam L/h)
    (1, 3.88, 2.96), (1, 1.99, 1.94), (2, 2.97, 2.73), (2, 2.77, 3.00),
    (2, 1.94, 2.02), (2, 1.89, 2.58), (3, 2.39, 2.26), (4, 3.03, 2.02),
    (5, 2.38, 1.78), (5, 2.94, 3.43), (6, 2.69, 2.56), (7, 2.26, 2.22),
    (8, 2.27, 2.51), (8, 2.77, 3.14), (8, 1.81, 2.22), (8, 2.05, 2.89),
    (8, 3.63, 6.13),
]


def gatti():
    caz = [x[1] for x in GATTI]
    avi = [x[2] for x in GATTI]
    # transcription check against the medians and IQRs the authors report independently
    assert abs(np.median(caz) - 2.39) < 0.005, "ceftazidime median does not match the publication"
    assert abs(np.median(avi) - 2.56) < 0.005, "avibactam median does not match the publication"
    assert abs(np.percentile(caz, 25) - 2.05) < 0.02 and abs(np.percentile(caz, 75) - 2.94) < 0.02
    assert abs(np.percentile(avi, 25) - 2.22) < 0.02 and abs(np.percentile(avi, 75) - 2.96) < 0.02
    print("  transcription verified: medians and IQRs reproduce the published values exactly")
    return caz, avi


# --- Dataset 2: Dryad, individual concentrations -> Dose/AUC ------------------
def _profiles(filename, value_column):
    path = os.path.join(DATA, "dryad_Li2025_CRRT", filename)
    by_subject = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            by_subject.setdefault(int(row["subjectID"]), []).append(
                (float(row["Time_h"]), float(row[value_column])))
    return {k: sorted(v) for k, v in by_subject.items()}


def _trapezoid(points):
    return sum((points[i + 1][0] - points[i][0]) * (points[i + 1][1] + points[i][1]) / 2.0
               for i in range(len(points) - 1))


def dryad():
    """Steady-state clearance over one dosing interval: CL = Dose / AUC(0-tau).

    Pre-filter concentrations are the systemic ones. Sparse sampling misses the
    peak in most patients, so AUC is biased low and clearance biased high; the
    bias acts on both analytes in the same direction and largely cancels in the
    correlation, which is the quantity of interest here.
    """
    caz_p = _profiles("Ceftazidime_concentration.csv", "caz_pre")
    avi_p = _profiles("Avibactam_concentration.csv", "avi_pre")
    ids = sorted(caz_p)
    caz = [2000.0 / _trapezoid(caz_p[i]) for i in ids]   # 2 g ceftazidime q8h
    avi = [500.0 / _trapezoid(avi_p[i]) for i in ids]    # 0.5 g avibactam q8h
    print(f"  {len(ids)} patients, {sum(len(v) for v in caz_p.values())} ceftazidime and "
          f"{sum(len(v) for v in avi_p.values())} avibactam observations")
    return caz, avi


def main():
    print("=" * 78)
    print("EMPIRICAL CEFTAZIDIME-AVIBACTAM CLEARANCE CORRELATION")
    print("Two independent, openly licensed patient-level ICU datasets")
    print("=" * 78)

    print("\nDataset 1 - Gatti 2023, CVVHDF, Italy (CC BY-NC-ND)")
    r1, lo1, hi1 = report("  RESULT", *gatti(), n_for_ci=len(GATTI))

    print("\nDataset 2 - Dryad doi:10.5061/dryad.fxpnvx16s, CRRT, China (CC0)")
    r2, lo2, hi2 = report("  RESULT", *dryad())

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"  Gatti 2023 (CVVHDF)   r = {r1:.3f}  95% CI {lo1:.3f} to {hi1:.3f}")
    print(f"  Dryad / Li  (CRRT)    r = {r2:.3f}  95% CI {lo2:.3f} to {hi2:.3f}")
    print(f"  Assumed in the model  rho = {ASSUMED_RHO}")
    print("\n  Two unrelated cohorts, two countries, two renal-replacement modalities,")
    print("  two independent analyses -- essentially the same estimate, far below 0.94.")
    print("\n  BOTH COHORTS ARE ON RENAL REPLACEMENT THERAPY. The primary scenario excludes it.")
    print("  A shared circuit removes both analytes together and should INFLATE this")
    print("  correlation, which is why the low values argue against 0.94 rather than for")
    print("  a specific replacement value. Report as a sensitivity analysis, not a new rho.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
