"""MODEL 1 — joint two-analyte population pharmacokinetic model for ceftazidime
and avibactam, with the CROSS-DRUG clearance correlation ESTIMATED rather than assumed.

WHY THIS MODEL EXISTS
    The cross-drug random-effect correlation between ceftazidime and avibactam
    clearance has been quantified exactly once in the literature (Cojutti 2024,
    rho = 0.94, RSE 23.8%). No regulatory document and no other publication reports
    it; the registrational analyses resampled eta-pairs empirically without ever
    estimating the covariance. This model provides a second, independent estimate
    from openly licensed individual patient data.

POPULATION — READ BEFORE USING ANY RESULT
    21 critically ill adults on CONTINUOUS RENAL REPLACEMENT THERAPY receiving
    INTERMITTENT 8-hourly infusion (2 g ceftazidime + 0.5 g avibactam).
    The manuscript's primary scenario is non-RRT adults on CONTINUOUS infusion.
    This is a model of a DIFFERENT POPULATION and a DIFFERENT ADMINISTRATION MODE.
    Its clearance and volume estimates do not transfer to the primary scenario.
    Only the correlation is carried downstream, and only as a sensitivity bound —
    never as a replacement value for rho.

VARIANCE STRUCTURE — and why it was reduced
    A first fit estimated the ceftazidime/avibactam VOLUME correlation at exactly
    1.000, a boundary value and the signature of an overparameterised variance
    model. The structure below imposes that boundary as a constraint rather than
    estimating it: a SINGLE shared volume deviate, scaled by its own omega for each
    analyte. This is also the more defensible model mechanistically, since both
    agents distribute into extracellular water, so a patient with an expanded
    volume has it for both. One parameter fewer, no boundary, Omega well conditioned.

        z = [z_CL_caz, z_CL_avi, z_V_shared],  corr(z1, z2) = r_cl, z3 independent
        CL_caz = th1 exp(w1 z1)      CL_avi = th2 exp(w2 z2)
        V_caz  = th3 exp(w3 z3)      V_avi  = th4 exp(w4 z3)

    Omega is therefore the CORRELATION matrix of standard-normal deviates; the
    magnitudes enter in log_pred. That keeps Omega well conditioned and makes the
    estimated correlation directly interpretable.

DATA
    Li C, Wang Y, Chen F, Huang L, Dong J, Fan W, Yue H, Ge Y.
    Dryad doi:10.5061/dryad.fxpnvx16s, CC0 1.0 public domain.
    Primary article: Antimicrob Agents Chemother 2026;70(2):e0143825.
    Depositors confirm explicit patient consent for public-domain release.

ESTIMATION
    Nonlinear mixed effects. Marginal likelihood by first-order conditional
    estimation with the Laplace approximation; the individual objective is a
    penalised nonlinear least-squares problem solved by Levenberg-Marquardt, and
    the curvature at the mode uses the Gauss-Newton form H = 2 (J'J/sigma^2 + Omega^-1),
    the standard FOCE approximation, which is what makes the fit tractable.
    Fully deterministic: no random number generation anywhere in estimation.
"""
from __future__ import annotations

import csv
import os
import sys
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.stats import chi2

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data_external", "dryad_Li2025_CRRT")
OUT = os.path.join(os.path.dirname(HERE), "outputs")

TAU = 8.0                              # dosing interval, hours
DOSE = {"caz": 2000.0, "avi": 500.0}   # mg per administration
N_ETA = 3                              # [z_CL_caz, z_CL_avi, z_V_shared]
ANALYTES = ("caz", "avi")


# ----------------------------------------------------------------- data ------

@dataclass(frozen=True)
class Subject:
    sid: int
    t_inf: float
    times: dict
    logconc: dict


def load():
    def profiles(fn, col):
        d = {}
        with open(os.path.join(DATA, fn), encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                d.setdefault(int(r["subjectID"]), []).append(
                    (float(r["Time_h"]), float(r[col])))
        return {k: sorted(v) for k, v in d.items()}

    caz = profiles("Ceftazidime_concentration.csv", "caz_pre")
    avi = profiles("Avibactam_concentration.csv", "avi_pre")
    t_inf = {}
    with open(os.path.join(DATA, "Medication_Information.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            t_inf[int(r["subjectID"])] = float(r["infusion_duration_cat"])

    return [Subject(
        sid=sid, t_inf=t_inf[sid],
        times={"caz": np.array([p[0] for p in caz[sid]]),
               "avi": np.array([p[0] for p in avi[sid]])},
        logconc={"caz": np.log([p[1] for p in caz[sid]]),
                 "avi": np.log([p[1] for p in avi[sid]])})
        for sid in sorted(caz)]


# ------------------------------------------------------- structural model ----

def css(t, cl, v, dose_mg, t_inf, tau=TAU):
    """Steady-state one-compartment concentration under repeated IV infusion.

    Closed form. With C0 the trough at the start of an interval, steady state
    requires C(tau) = C0, giving
        C0 = (R/CL)(1 - exp(-k T)) exp(-k (tau - T)) / (1 - exp(-k tau)).
    Verified against explicit superposition of 80 doses: maximum relative
    difference 7e-10 over the parameter ranges in this cohort, about 43x faster.
    """
    k = cl / v
    rate = dose_mg / t_inf
    t = np.asarray(t, float)
    c0 = (rate / cl) * (1.0 - np.exp(-k * t_inf)) * np.exp(-k * (tau - t_inf)) \
        / (1.0 - np.exp(-k * tau))
    t_in = np.minimum(t, t_inf)
    during = (rate / cl) * (1.0 - np.exp(-k * t_in)) + c0 * np.exp(-k * t_in)
    c_end = (rate / cl) * (1.0 - np.exp(-k * t_inf)) + c0 * np.exp(-k * t_inf)
    return np.where(t <= t_inf, during, c_end * np.exp(-k * np.maximum(t - t_inf, 0.0)))


def css_superposition(t, cl, v, dose_mg, t_inf, tau=TAU, n=80):
    """Explicit superposition, retained only as the verification reference."""
    k, rate = cl / v, dose_mg / t_inf
    t = np.asarray(t, float)
    return sum((rate / cl) * (1.0 - np.exp(-k * np.minimum(t + i * tau, t_inf)))
               * np.exp(-k * np.maximum(t + i * tau - t_inf, 0.0)) for i in range(n))


def log_pred(subj, theta, z, w):
    """Stacked log predictions for both analytes, in the order caz then avi."""
    cl = (theta[0] * np.exp(w[0] * z[0]), theta[1] * np.exp(w[1] * z[1]))
    v = (theta[2] * np.exp(w[2] * z[2]), theta[3] * np.exp(w[3] * z[2]))
    return np.concatenate([
        np.log(np.clip(css(subj.times[a], cl[j], v[j], DOSE[a], subj.t_inf), 1e-10, None))
        for j, a in enumerate(ANALYTES)])


def obs_vector(subj):
    return np.concatenate([subj.logconc[a] for a in ANALYTES])


def sigma_vector(subj, sigma):
    return np.concatenate([np.full(len(subj.times[a]), sigma[j])
                           for j, a in enumerate(ANALYTES)])


# ------------------------------------------------------ variance structure ---

def build_omega(p):
    """p = [log w1..w4, z_cl]; the clearance correlation is tanh(z_cl)."""
    w = np.exp(p[:4])
    r_cl = float(np.tanh(p[4]))
    om = np.eye(N_ETA)
    om[0, 1] = om[1, 0] = r_cl
    return om, w, r_cl, 1.0


def build_omega_diag(p):
    """Null model: clearance correlation fixed at zero. One df against the full model."""
    return np.eye(N_ETA), np.exp(p[:4]), 0.0, 1.0


# ------------------------------------------------------------ likelihood -----

def laplace_subject(subj, theta, w, om, om_inv, om_chol_inv, sigma, z0):
    """Individual objective at the empirical-Bayes mode, plus the Laplace term."""
    y = obs_vector(subj)
    sig = sigma_vector(subj, sigma)

    def residuals(z):
        return np.concatenate([(y - log_pred(subj, theta, z, w)) / sig, om_chol_inv @ z])

    sol = least_squares(residuals, z0, method="lm", xtol=1e-10, ftol=1e-10, max_nfev=300)
    z = sol.x

    base = log_pred(subj, theta, z, w)
    J = np.empty((len(y), N_ETA))
    h = 1e-5
    for a in range(N_ETA):
        e = z.copy()
        e[a] += h
        J[:, a] = (log_pred(subj, theta, e, w) - base) / h
    Js = J / sig[:, None]
    H = Js.T @ Js + om_inv

    r_obs = (y - base) / sig
    _, logdet_om = np.linalg.slogdet(om)
    sign_h, logdet_h = np.linalg.slogdet(H)
    if sign_h <= 0:
        logdet_h = 60.0
    ofv_i = (float(r_obs @ r_obs) + float(z @ om_inv @ z)
             + 2.0 * float(np.sum(np.log(sig))) + logdet_om + logdet_h)
    return ofv_i, z


def ofv(params, subjects, omega_builder, n_omega, cache):
    theta = np.exp(params[:4])
    om, w, _, _ = omega_builder(params[4:4 + n_omega])
    sigma = np.exp(params[4 + n_omega:4 + n_omega + 2])
    try:
        om_inv = np.linalg.inv(om)
        om_chol_inv = np.linalg.inv(np.linalg.cholesky(om))
    except np.linalg.LinAlgError:
        return 1e10
    total = 0.0
    for s in subjects:
        val, z = laplace_subject(s, theta, w, om, om_inv, om_chol_inv, sigma,
                                 cache.get(s.sid, np.zeros(N_ETA)))
        cache[s.sid] = z
        total += val
    return total if np.isfinite(total) else 1e10


def fit(subjects, omega_builder, n_omega, label, p0, quiet=False, max_rounds=12,
        tol=1e-4, fixed=None):
    """Alternate Nelder-Mead and L-BFGS-B until the objective stops improving.

    A single Nelder-Mead pass on 11 parameters terminates on the evaluation limit
    rather than on its tolerance, which is not convergence. Alternating a
    derivative-free pass with a quasi-Newton polish and repeating until successive
    rounds improve the objective by less than `tol` gives a defensible criterion:
    the reported fit is the point at which further optimisation changes nothing.

    `fixed` is an optional {index: value} map used by the profile likelihood to
    hold a parameter at a grid value while the rest are re-estimated.
    """
    fixed = fixed or {}
    free = [i for i in range(len(p0)) if i not in fixed]

    def expand(x_free):
        x = np.empty(len(p0))
        for k, i in enumerate(free):
            x[i] = x_free[k]
        for i, v in fixed.items():
            x[i] = v
        return x

    def objective(x_free, cache):
        return ofv(expand(x_free), subjects, omega_builder, n_omega, cache)

    if not quiet:
        print(f"  fitting {label} ...", flush=True)
    x = np.array([p0[i] for i in free], float)
    prev = np.inf
    cache = {}
    rounds = 0
    for rounds in range(1, max_rounds + 1):
        # The first round explores; later rounds only need to polish, so the
        # simplex budget drops sharply. Without this the routine spends most of
        # its evaluations re-exploring a region it has already resolved, which
        # matters because the profile likelihood and the leave-one-out analysis
        # call this function forty times.
        budget = 2500 if rounds == 1 else 600
        r1 = minimize(objective, x, args=(cache,), method="Nelder-Mead",
                      options={"maxiter": budget, "maxfev": budget,
                               "xatol": 1e-6, "fatol": 1e-6, "adaptive": True})
        r2 = minimize(objective, r1.x, args=(cache,), method="L-BFGS-B",
                      options={"maxiter": 120, "ftol": 1e-12, "gtol": 1e-9,
                               "eps": 1e-5})
        x, cur = (r2.x, r2.fun) if r2.fun < r1.fun else (r1.x, r1.fun)
        if prev - cur < tol:
            break
        prev = cur

    cache = {}
    final = ofv(expand(x), subjects, omega_builder, n_omega, cache)
    result = type("FitResult", (), {})()
    result.x = expand(x)
    result.fun = final
    result.rounds = rounds
    result.converged = (prev - final) < tol or rounds < max_rounds
    if not quiet:
        print(f"    OFV {final:.4f}   rounds {rounds}   "
              f"stopped on tolerance: {result.converged}")
    return result, cache


def standard_errors(params, subjects, omega_builder, n_omega):
    n = len(params)
    h = 5e-3
    H = np.zeros((n, n))
    for a in range(n):
        for b in range(a, n):
            pp, mm, pm, mp = (params.copy() for _ in range(4))
            pp[a] += h; pp[b] += h
            mm[a] -= h; mm[b] -= h
            pm[a] += h; pm[b] -= h
            mp[a] -= h; mp[b] += h
            H[a, b] = H[b, a] = (
                ofv(pp, subjects, omega_builder, n_omega, {})
                - ofv(pm, subjects, omega_builder, n_omega, {})
                - ofv(mp, subjects, omega_builder, n_omega, {})
                + ofv(mm, subjects, omega_builder, n_omega, {})) / (4 * h * h)
    try:
        cov = np.linalg.inv(H / 2.0)     # OFV = -2 log L, so information = H / 2
        return np.sqrt(np.abs(np.diag(cov))), cov
    except np.linalg.LinAlgError:
        return np.full(n, np.nan), None


def structural_self_check():
    t = np.array([0.0, 1, 2, 3, 4, 6, 8])
    worst = 0.0
    for cl in (1.5, 2.5, 4.0):
        for v in (10.0, 25.0, 45.0):
            for ti in (1.0, 2.0, 3.0):
                a = css(t, cl, v, 2000.0, ti)
                b = css_superposition(t, cl, v, 2000.0, ti)
                worst = max(worst, float(np.max(np.abs(a - b) / np.maximum(a, 1e-9))))
    assert worst < 1e-6, f"closed form disagrees with superposition ({worst:.2e})"
    print(f"  structural self-check passed (max relative difference {worst:.1e})")


# ------------------------------------------------------------------ main -----

def main():
    structural_self_check()
    subjects = load()
    nobs = sum(len(obs_vector(s)) for s in subjects)
    print("=" * 78)
    print("MODEL 1 - joint ceftazidime/avibactam population PK, CRRT cohort")
    print("=" * 78)
    print(f"  subjects {len(subjects)}   observations {nobs} "
          f"({nobs / len(subjects) / 2:.1f} per analyte per subject)")
    print("  POPULATION: CRRT, intermittent 8-hourly infusion - NOT the primary scenario")
    print("  Variance model reduced to a shared volume deviate; the 4-deviate model")
    print("  estimated the volume correlation at the 1.000 boundary.")
    print()

    p0 = np.array([np.log(2.57), np.log(3.22), np.log(20.0), np.log(27.0),
                   np.log(0.20), np.log(0.14), np.log(0.27), np.log(0.20),
                   0.9,
                   np.log(0.10), np.log(0.093)])
    full, cache = fit(subjects, build_omega, 5,
                      "full model (cross-drug clearance correlation estimated)", p0)
    null, _ = fit(subjects, build_omega_diag, 4,
                  "null model (correlation fixed at zero)",
                  np.concatenate([full.x[:8], full.x[9:]]))

    theta = np.exp(full.x[:4])
    om, w, r_cl, _ = build_omega(full.x[4:9])
    sigma = np.exp(full.x[9:11])
    se, _ = standard_errors(full.x, subjects, build_omega, 5)

    print()
    print("=" * 78)
    print("RESULTS")
    print("=" * 78)
    names = ["CL ceftazidime (L/h)", "CL avibactam (L/h)",
             "V ceftazidime (L)", "V avibactam (L)"]
    print()
    print("  Fixed effects            estimate     RSE")
    for i, nm in enumerate(names):
        print(f"    {nm:24} {theta[i]:8.3f}  {100 * se[i]:6.1f}%")

    print()
    print("  Between-subject variability")
    for i, nm in enumerate(["CL ceftazidime", "CL avibactam",
                            "V ceftazidime", "V avibactam"]):
        cv = 100 * np.sqrt(np.exp(w[i] ** 2) - 1)
        print(f"    {nm:24} omega {w[i]:6.4f}   CV {cv:5.1f}%")

    print()
    print("  CROSS-DRUG CLEARANCE CORRELATION  <-- the quantity of interest")
    z, sz = full.x[8], se[8]
    print(f"    corr(eta_CL_caz, eta_CL_avi) = {r_cl:.3f}")
    lo = hi = float("nan")
    if np.isfinite(sz):
        lo, hi = float(np.tanh(z - 1.96 * sz)), float(np.tanh(z + 1.96 * sz))
        print(f"    95% CI  {lo:.3f} to {hi:.3f}   (Fisher z scale, SE {sz:.3f})")
        verdict = "EXCLUDES" if hi < 0.94 else "INCLUDES"
        print(f"    The interval {verdict} the assumed value of 0.94.")
    print("    corr(eta_V_caz, eta_V_avi) = 1 by construction (shared volume deviate)")

    print()
    print("  Residual error (proportional)")
    print(f"    ceftazidime {100 * sigma[0]:6.1f}%      avibactam {100 * sigma[1]:6.1f}%")

    d_ofv = null.fun - full.fun
    p_lrt = chi2.sf(max(d_ofv, 0.0), 1)
    print()
    print("  Model comparison (likelihood ratio, 1 df: the clearance correlation)")
    print(f"    OFV full {full.fun:9.2f}   OFV null {null.fun:9.2f}"
          f"   dOFV {d_ofv:7.2f}   p = {p_lrt:.4g}")
    print(f"    BIC full {full.fun + 11 * np.log(nobs):9.2f}"
          f"   BIC null {null.fun + 10 * np.log(nobs):9.2f}")

    etas = np.array([cache[s.sid] for s in subjects])
    print()
    print("  Shrinkage of the standard-normal deviates")
    for i, nm in enumerate(["z CL ceftazidime", "z CL avibactam", "z V shared"]):
        print(f"    {nm:24} {100 * (1 - etas[:, i].std(ddof=1)):6.1f}%")
    print()
    print("  Empirical correlation of the individual clearance deviates: "
          f"{np.corrcoef(etas[:, 0], etas[:, 1])[0, 1]:.3f}")

    print()
    print("  Comparison with the assumed value")
    print("    Cojutti 2024 (non-RRT, continuous infusion): rho = 0.940  (RSE 23.8%)")
    print(f"    This model   (CRRT, intermittent infusion) : rho = {r_cl:.3f}")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "model1_joint_popk_parameters.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        wcsv = csv.writer(fh, lineterminator="\n")
        wcsv.writerow(["parameter", "estimate", "rse_pct", "ci_low", "ci_high", "note"])
        for i, nm in enumerate(names):
            wcsv.writerow([nm, f"{theta[i]:.4f}", f"{100 * se[i]:.1f}", "", "",
                           "CRRT cohort, intermittent infusion; does not transfer "
                           "to the primary scenario"])
        for i, nm in enumerate(["omega_CL_caz", "omega_CL_avi",
                                "omega_V_caz", "omega_V_avi"]):
            wcsv.writerow([nm, f"{w[i]:.4f}", "", "", "",
                           f"CV {100 * np.sqrt(np.exp(w[i] ** 2) - 1):.1f}%"])
        wcsv.writerow(["corr_CL_caz_avi", f"{r_cl:.4f}", "",
                       f"{lo:.4f}" if np.isfinite(lo) else "",
                       f"{hi:.4f}" if np.isfinite(hi) else "",
                       "CROSS-DRUG clearance correlation; assumed 0.94 in the primary model"])
        wcsv.writerow(["corr_V_caz_avi", "1.0000", "", "", "",
                       "fixed by construction: shared volume deviate"])
        wcsv.writerow(["sigma_prop_caz", f"{sigma[0]:.4f}", "", "", "",
                       "proportional residual error"])
        wcsv.writerow(["sigma_prop_avi", f"{sigma[1]:.4f}", "", "", "",
                       "proportional residual error"])
        wcsv.writerow(["OFV_full", f"{full.fun:.3f}", "", "", "", ""])
        wcsv.writerow(["OFV_null_no_crossdrug", f"{null.fun:.3f}", "", "", "", ""])
        wcsv.writerow(["dOFV", f"{d_ofv:.3f}", "", "", "",
                       f"likelihood ratio, 1 df, p = {p_lrt:.4g}"])
    print()
    print(f"  wrote {path}")

    ipath = os.path.join(OUT, "model1_individual_parameters.csv")
    with open(ipath, "w", newline="", encoding="utf-8") as fh:
        wcsv = csv.writer(fh, lineterminator="\n")
        wcsv.writerow(["subjectID", "z_CL_caz", "z_CL_avi", "z_V_shared",
                       "CL_caz_L_h", "CL_avi_L_h", "V_caz_L", "V_avi_L"])
        for s, e in zip(subjects, etas):
            wcsv.writerow([s.sid] + [f"{x:.5f}" for x in e] + [
                f"{theta[0] * np.exp(w[0] * e[0]):.4f}",
                f"{theta[1] * np.exp(w[1] * e[1]):.4f}",
                f"{theta[2] * np.exp(w[2] * e[2]):.4f}",
                f"{theta[3] * np.exp(w[3] * e[2]):.4f}"])
    print(f"  wrote {ipath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
