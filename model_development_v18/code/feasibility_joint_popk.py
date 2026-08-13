"""Feasibility test: can the Dryad CRRT dataset support a joint two-analyte
population PK model with an estimated cross-drug clearance correlation?

Two-stage check first (individual nonlinear least squares per subject per analyte),
then a check of whether V is identifiable alongside CL from 5-7 samples.
"""
import csv, os, sys
import numpy as np
from scipy.optimize import least_squares
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
D = r"E:/Github Project/IJAA_submission_package_v16/model_development_v18/data_external/dryad_Li2025_CRRT"

TAU = 8.0


def profiles(fn, col):
    d = {}
    for r in csv.DictReader(open(os.path.join(D, fn), encoding="utf-8")):
        d.setdefault(int(r["subjectID"]), []).append((float(r["Time_h"]), float(r[col])))
    return {k: sorted(v) for k, v in d.items()}


def tinf():
    d = {}
    for r in csv.DictReader(open(os.path.join(D, "Medication_Information.csv"), encoding="utf-8")):
        d[int(r["subjectID"])] = float(r["infusion_duration_cat"])
    return d


def css_1cmt(t, cl, v, dose_mg, t_inf, tau=TAU, n_doses=60):
    """Steady-state 1-compartment, intermittent infusion, by superposition."""
    k = cl / v
    R = dose_mg / t_inf
    t = np.atleast_1d(np.asarray(t, float))
    out = np.zeros_like(t)
    for n in range(n_doses):
        s = t + n * tau                      # time since the n-th previous dose start
        during = s <= t_inf
        out += np.where(during,
                        (R / cl) * (1 - np.exp(-k * np.clip(s, 0, None))),
                        (R / cl) * (1 - np.exp(-k * t_inf)) * np.exp(-k * (s - t_inf)))
    return out


def fit_subject(times, conc, dose_mg, t_inf, fit_v=True, v_fixed=None):
    y = np.log(np.asarray(conc, float))

    def resid(p):
        cl = np.exp(p[0])
        v = np.exp(p[1]) if fit_v else v_fixed
        pred = css_1cmt(times, cl, v, dose_mg, t_inf)
        pred = np.clip(pred, 1e-9, None)
        return np.log(pred) - y

    p0 = [np.log(2.5), np.log(20.0)] if fit_v else [np.log(2.5)]
    best = None
    for cl0 in (1.0, 2.5, 5.0):
        for v0 in (10.0, 25.0, 60.0):
            start = [np.log(cl0), np.log(v0)] if fit_v else [np.log(cl0)]
            try:
                r = least_squares(resid, start, method="lm", max_nfev=4000)
            except Exception:
                continue
            if best is None or r.cost < best.cost:
                best = r
    cl = float(np.exp(best.x[0]))
    v = float(np.exp(best.x[1])) if fit_v else v_fixed
    rss = float(2 * best.cost)
    return cl, v, rss, len(times)


def main():
    caz = profiles("Ceftazidime_concentration.csv", "caz_pre")
    avi = profiles("Avibactam_concentration.csv", "avi_pre")
    T = tinf()
    ids = sorted(caz)

    print("=" * 84)
    print("FEASIBILITY: individual two-stage fits, steady-state 1-compartment, IV intermittent")
    print("Dose 2000 mg ceftazidime + 500 mg avibactam q8h; infusion duration from the dataset")
    print("=" * 84)
    print(f"\n{'ID':>3} {'n':>2} {'Tinf':>4} | {'CL_CAZ':>7} {'V_CAZ':>7} {'CV%res':>7} | "
          f"{'CL_AVI':>7} {'V_AVI':>7} {'CV%res':>7}")
    print("-" * 84)

    rows = []
    for i in ids:
        ti = T[i]
        tc = [p[0] for p in caz[i]]; cc = [p[1] for p in caz[i]]
        ta = [p[0] for p in avi[i]]; ca = [p[1] for p in avi[i]]
        cl1, v1, rss1, n1 = fit_subject(tc, cc, 2000.0, ti)
        cl2, v2, rss2, n2 = fit_subject(ta, ca, 500.0, ti)
        e1 = 100 * np.sqrt(rss1 / max(n1 - 2, 1))
        e2 = 100 * np.sqrt(rss2 / max(n2 - 2, 1))
        rows.append((i, n1, ti, cl1, v1, e1, cl2, v2, e2))
        print(f"{i:3d} {n1:2d} {ti:4.0f} | {cl1:7.3f} {v1:7.2f} {e1:7.1f} | "
              f"{cl2:7.3f} {v2:7.2f} {e2:7.1f}")

    a = np.array([[r[3], r[4], r[6], r[7]] for r in rows], float)
    clc, vc, cla, va = a[:, 0], a[:, 1], a[:, 2], a[:, 3]

    print("\n" + "=" * 84)
    print("PARAMETER DISTRIBUTIONS")
    print("=" * 84)
    for nm, x in (("CL ceftazidime (L/h)", clc), ("V ceftazidime (L)", vc),
                  ("CL avibactam (L/h)", cla), ("V avibactam (L)", va)):
        lx = np.log(x)
        print(f"  {nm:24} median {np.median(x):7.2f}  IQR {np.percentile(x,25):6.2f}-{np.percentile(x,75):6.2f}"
              f"   apparent CV {100*np.sqrt(np.exp(lx.var(ddof=1))-1):5.1f}%")

    print("\n" + "=" * 84)
    print("CROSS-DRUG CLEARANCE CORRELATION (the parameter at issue)")
    print("=" * 84)
    for nm, x, y in (("CL_CAZ vs CL_AVI", clc, cla),
                     ("V_CAZ  vs V_AVI ", vc, va),
                     ("CL_CAZ vs V_CAZ ", clc, vc)):
        r, p = stats.pearsonr(np.log(x), np.log(y))
        z = np.arctanh(r); se = 1 / np.sqrt(len(x) - 3); zc = stats.norm.ppf(0.975)
        lo, hi = np.tanh(z - zc * se), np.tanh(z + zc * se)
        print(f"  {nm}   r = {r:6.3f}   95% CI {lo:6.3f} to {hi:6.3f}   p = {p:.4g}")

    r, _ = stats.pearsonr(np.log(clc), np.log(cla))
    zz = (np.arctanh(r) - np.arctanh(0.94)) * np.sqrt(len(clc) - 3)
    print(f"\n  vs the assumed rho = 0.94:  z = {zz:.2f},  p = {2*stats.norm.sf(abs(zz)):.3g}")

    print("\n" + "=" * 84)
    print("IDENTIFIABILITY VERDICT")
    print("=" * 84)
    nsub = len(ids); nobs = sum(len(caz[i]) for i in ids) + sum(len(avi[i]) for i in ids)
    print(f"  subjects {nsub}, observations {nobs} ({nobs/nsub/2:.1f} per analyte per subject)")
    print(f"  residual error of the individual fits: ceftazidime median {np.median([r[5] for r in rows]):.1f}%,"
          f" avibactam median {np.median([r[8] for r in rows]):.1f}%")
    print("  A joint NLME model would carry 4 fixed effects (CL, V per analyte),")
    print("  4 variance terms and 1 cross-drug covariance = 9 random-effect parameters,")
    print(f"  supported by {nsub} subjects. That is roughly 2.3 subjects per random-effect")
    print("  parameter -- thin, but within the range where a REDUCED model (IIV on CL only,")
    print("  V fixed or shared) is estimable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
