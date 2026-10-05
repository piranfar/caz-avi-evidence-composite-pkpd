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

VARIANCE STRUCTURE
    Both drugs are measured in the same pre-filter sample at the same times, so the
    model carries a correlation c between the 2 drugs' residual errors in each sample.
    Without it, a sample-level error shared by both drugs (sampling time, handling,
    dilution) can only be absorbed by the random effects, which biases the clearance
    correlation. With c in the model:

        z = [z_CL_caz, z_CL_avi, z_V_caz, z_V_avi]
        corr(z1, z2) = r_cl (the quantity of interest), corr(z3, z4) = r_v
        CL_caz = th1 exp(w1 z1)      CL_avi = th2 exp(w2 z2)
        V_caz  = th3 exp(w3 z3)      V_avi  = th4 exp(w4 z4)
        log C_obs = log C_pred + eps,  sd(eps) = sigma_drug,  corr(eps_caz, eps_avi) = c

    Reduction rule for the volume deviates: share them (one deviate scaled for each
    drug) only when their estimated correlation sits at the 1.0 boundary. Without c it
    did (r_v = 1.000), which is why an earlier version shared them; with c it does not
    (r_v about 0.70), so they are estimated separately. build_omega (3 deviates, shared
    volume) is kept for that comparison and for the tests.

    Omega is the CORRELATION matrix of standard-normal deviates; the magnitudes enter
    in log_pred. That keeps Omega well conditioned and makes the estimated correlation
    directly interpretable.

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
    the standard FOCE approximation, which is what makes the fit tractable. The paired
    residuals are whitened before they enter the objective (see whiten()), which is the
    exact bivariate normal likelihood for each sample.
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
    """Stacked log predictions for both analytes, in the order caz then avi.

    With 3 deviates the 2 volumes share z[2]; with 4, z[3] is the avibactam volume deviate.
    """
    cl = (theta[0] * np.exp(w[0] * z[0]), theta[1] * np.exp(w[1] * z[1]))
    v = (theta[2] * np.exp(w[2] * z[2]), theta[3] * np.exp(w[3] * z[3 if len(z) > 3 else 2]))
    return np.concatenate([
        np.log(np.clip(css(subj.times[a], cl[j], v[j], DOSE[a], subj.t_inf), 1e-10, None))
        for j, a in enumerate(ANALYTES)])


def obs_vector(subj):
    return np.concatenate([subj.logconc[a] for a in ANALYTES])


def sigma_vector(subj, sigma):
    return np.concatenate([np.full(len(subj.times[a]), sigma[j])
                           for j, a in enumerate(ANALYTES)])


def whiten(v, n_caz, c):
    """Remove the correlation c between the 2 drugs' standardized residuals in the same sample.

    Both drugs are measured in each sample, so the rows are paired: row i of ceftazidime and
    row n_caz + i of avibactam come from the same tube. The avibactam rows are replaced by
    their residual on the ceftazidime rows, (b - c a) / sqrt(1 - c^2), which turns the paired
    bivariate normal into 2 independent standard normals. Works on vectors and on Jacobians.
    """
    if c == 0.0:
        return v
    a, b = v[:n_caz], v[n_caz:]
    return np.concatenate([a, (b - c * a) / np.sqrt(1.0 - c * c)])


def residual_correlation(params, n_omega):
    """Residual correlation between the 2 drugs (tanh scale), or 0 when it is not in the model."""
    k = 4 + n_omega + 2
    return float(np.tanh(params[k])) if len(params) > k else 0.0


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


def build_omega4(p):
    """Separate volume deviates: p = [log w1..w4, z_cl, z_v]; volume correlation tanh(z_v)."""
    w = np.exp(p[:4])
    r_cl, r_v = float(np.tanh(p[4])), float(np.tanh(p[5]))
    om = np.eye(4)
    om[0, 1] = om[1, 0] = r_cl
    om[2, 3] = om[3, 2] = r_v
    return om, w, r_cl, r_v


# ------------------------------------------------------------ likelihood -----

def laplace_subject(subj, theta, w, om, om_inv, om_chol_inv, sigma, z0, c=0.0):
    """Individual objective at the empirical-Bayes mode, plus the Laplace term.

    c is the correlation between the 2 drugs' residual errors in the same sample.
    """
    y = obs_vector(subj)
    sig = sigma_vector(subj, sigma)
    n_caz, n_eta = len(subj.times["caz"]), len(z0)
    if c != 0.0:
        assert np.array_equal(subj.times["caz"], subj.times["avi"]), subj.sid

    def residuals(z):
        return np.concatenate([whiten((y - log_pred(subj, theta, z, w)) / sig, n_caz, c),
                               om_chol_inv @ z])

    sol = least_squares(residuals, z0, method="lm", xtol=1e-10, ftol=1e-10, max_nfev=300)
    z = sol.x

    base = log_pred(subj, theta, z, w)
    J = np.empty((len(y), n_eta))
    h = 1e-5
    for a in range(n_eta):
        e = z.copy()
        e[a] += h
        J[:, a] = (log_pred(subj, theta, e, w) - base) / h
    Js = whiten(J / sig[:, None], n_caz, c)
    H = Js.T @ Js + om_inv

    r_obs = whiten((y - base) / sig, n_caz, c)
    _, logdet_om = np.linalg.slogdet(om)
    sign_h, logdet_h = np.linalg.slogdet(H)
    if sign_h <= 0:
        logdet_h = 60.0
    ofv_i = (float(r_obs @ r_obs) + float(z @ om_inv @ z)
             + 2.0 * float(np.sum(np.log(sig))) + n_caz * np.log(1.0 - c * c)
             + logdet_om + logdet_h)
    return ofv_i, z


def ofv(params, subjects, omega_builder, n_omega, cache):
    theta = np.exp(params[:4])
    om, w, _, _ = omega_builder(params[4:4 + n_omega])
    sigma = np.exp(params[4 + n_omega:4 + n_omega + 2])
    c = residual_correlation(params, n_omega)
    if abs(c) >= 0.999:
        return 1e10
    try:
        om_inv = np.linalg.inv(om)
        om_chol_inv = np.linalg.inv(np.linalg.cholesky(om))
    except np.linalg.LinAlgError:
        return 1e10
    total = 0.0
    for s in subjects:
        val, z = laplace_subject(s, theta, w, om, om_inv, om_chol_inv, sigma,
                                 cache.get(s.sid, np.zeros(om.shape[0])), c)
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

# Primary model: 4 deviates plus the residual correlation (13 parameters).
# Layout: [log th1..th4, log w1..w4, atanh r_cl, atanh r_v, log sigma_caz, log sigma_avi, atanh c]
P0 = np.array([np.log(2.576), np.log(3.228), np.log(20.06), np.log(27.36),
               np.log(0.205), np.log(0.139), np.log(0.267), np.log(0.194),
               np.arctanh(0.59), np.arctanh(0.70),
               np.log(0.0997), np.log(0.0932), np.arctanh(0.63)])
N_OMEGA = 6
IDX_RHO, IDX_RV, IDX_C = 8, 9, 12


def unpack(x):
    """theta, Omega, w, r_cl, r_v, sigma, c for a 13-parameter vector of the primary model."""
    theta = np.exp(x[:4])
    om, w, r_cl, r_v = build_omega4(x[4:4 + N_OMEGA])
    sigma = np.exp(x[4 + N_OMEGA:6 + N_OMEGA])
    return theta, om, w, r_cl, r_v, sigma, residual_correlation(x, N_OMEGA)


def main():
    structural_self_check()
    subjects = load()
    nobs = sum(len(obs_vector(s)) for s in subjects)
    print("=" * 78)
    print("MODEL 1 - joint ceftazidime/avibactam population PK, CRRT cohort")
    print("=" * 78)
    print(f"  subjects {len(subjects)}   observations {nobs} "
          f"({nobs / len(subjects) / 2:.1f} per analyte per subject)")
    print("  4 deviates (separate volume deviates) and a residual correlation between the")
    print("  2 drugs measured in the same sample")
    print()

    full, cache = fit(subjects, build_omega4, N_OMEGA,
                      "full model (cross-drug clearance correlation estimated)", P0)
    null, _ = fit(subjects, build_omega4, N_OMEGA, "null model (clearance correlation fixed at zero)",
                  full.x, fixed={IDX_RHO: 0.0})
    theta, om, w, r_cl, r_v, sigma, c = unpack(full.x)
    se, cov = standard_errors(full.x, subjects, build_omega4, N_OMEGA)

    print()
    print("=" * 78)
    print("RESULTS")
    print("=" * 78)
    names = ["CL ceftazidime (L/h)", "CL avibactam (L/h)", "V ceftazidime (L)", "V avibactam (L)"]
    print()
    print("  Fixed effects            estimate     RSE")
    for i, nm in enumerate(names):
        print(f"    {nm:24} {theta[i]:8.3f}  {100 * se[i]:6.1f}%")
    print()
    print("  Between-subject variability")
    for i, nm in enumerate(["CL ceftazidime", "CL avibactam", "V ceftazidime", "V avibactam"]):
        print(f"    {nm:24} omega {w[i]:6.4f}   CV {100 * np.sqrt(np.exp(w[i] ** 2) - 1):5.1f}%")
    print()
    print("  CROSS-DRUG CLEARANCE CORRELATION  <-- the quantity of interest")
    z, sz = full.x[IDX_RHO], se[IDX_RHO]
    print(f"    corr(eta_CL_caz, eta_CL_avi) = {r_cl:.4f}")
    lo = hi = float("nan")
    if np.isfinite(sz):
        lo, hi = float(np.tanh(z - 1.96 * sz)), float(np.tanh(z + 1.96 * sz))
        print(f"    Wald 95% CI on the Fisher z scale {lo:.3f} to {hi:.3f} (SE {sz:.3f}); "
              f"the profile-likelihood interval is in model1_finalise.py")
    print(f"    corr(eta_V_caz, eta_V_avi) = {r_v:.4f}")
    print(f"    residual correlation in the same sample = {c:.4f}")
    print()
    print("  Residual error (proportional)")
    print(f"    ceftazidime {100 * sigma[0]:6.2f}%      avibactam {100 * sigma[1]:6.2f}%")

    d_ofv = null.fun - full.fun
    p_lrt = chi2.sf(max(d_ofv, 0.0), 1)
    print()
    print("  Model comparison (likelihood ratio, 1 df: the clearance correlation)")
    print(f"    OFV full {full.fun:9.4f}   OFV null {null.fun:9.4f}   dOFV {d_ofv:7.3f}   p = {p_lrt:.4g}")

    etas = np.array([cache[s.sid] for s in subjects])
    print()
    print("  Shrinkage of the standard-normal deviates")
    for i, nm in enumerate(["z CL ceftazidime", "z CL avibactam", "z V ceftazidime", "z V avibactam"]):
        print(f"    {nm:24} {100 * (1 - etas[:, i].std(ddof=1)):6.1f}%")
    print()
    print("  Empirical correlation of the individual clearance deviates: "
          f"{np.corrcoef(etas[:, 0], etas[:, 1])[0, 1]:.3f}")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "model1_joint_popk_parameters.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        wcsv = csv.writer(fh, lineterminator="\n")
        wcsv.writerow(["parameter", "estimate", "rse_pct", "ci_low", "ci_high", "note"])
        for i, nm in enumerate(names):
            wcsv.writerow([nm, f"{theta[i]:.6f}", f"{100 * se[i]:.1f}", "", "", ""])
        for i, nm in enumerate(["omega_CL_caz", "omega_CL_avi", "omega_V_caz", "omega_V_avi"]):
            wcsv.writerow([nm, f"{w[i]:.6f}", "", "", "", f"CV {100 * np.sqrt(np.exp(w[i] ** 2) - 1):.1f}%"])
        wcsv.writerow(["corr_CL_caz_avi", f"{r_cl:.6f}", "", f"{lo:.6f}" if np.isfinite(lo) else "",
                       f"{hi:.6f}" if np.isfinite(hi) else "", "Wald interval on the Fisher z scale"])
        wcsv.writerow(["corr_V_caz_avi", f"{r_v:.6f}", "", "", "", "separate volume deviates"])
        wcsv.writerow(["sigma_prop_caz", f"{sigma[0]:.6f}", "", "", "", "proportional residual error"])
        wcsv.writerow(["sigma_prop_avi", f"{sigma[1]:.6f}", "", "", "", "proportional residual error"])
        wcsv.writerow(["corr_residual_caz_avi", f"{c:.6f}", "", "", "",
                       "correlation of the 2 drugs' residual errors in the same sample"])
        wcsv.writerow(["OFV_full", f"{full.fun:.6f}", "", "", "", ""])
        wcsv.writerow(["OFV_null_no_crossdrug", f"{null.fun:.6f}", "", "", "", ""])
        wcsv.writerow(["dOFV", f"{d_ofv:.6f}", "", "", "", f"likelihood ratio, 1 df, p = {p_lrt:.4g}"])
    print()
    print(f"  wrote {path}")

    # Estimates and their asymptotic covariance on the estimation scale (log, log, Fisher z, log,
    # Fisher z), from the finite-difference Hessian of the OFV. crrt_parameter_uncertainty.py
    # draws parameter vectors from this distribution.
    cpath = os.path.join(OUT, "model1_parameter_covariance.csv")
    labels = ["log_CL_caz", "log_CL_avi", "log_V_caz", "log_V_avi", "log_omega_CL_caz", "log_omega_CL_avi",
              "log_omega_V_caz", "log_omega_V_avi", "atanh_corr_CL", "atanh_corr_V", "log_sigma_caz",
              "log_sigma_avi", "atanh_corr_residual"]
    with open(cpath, "w", newline="", encoding="utf-8") as fh:
        wcsv = csv.writer(fh, lineterminator="\n")
        wcsv.writerow(["parameter", "estimate"] + labels)
        for i, nm in enumerate(labels):
            row = [f"{v:.10g}" for v in cov[i]] if cov is not None else [""] * len(labels)
            wcsv.writerow([nm, f"{full.x[i]:.10f}"] + row)
    print(f"  wrote {cpath}")

    ipath = os.path.join(OUT, "model1_individual_parameters.csv")
    with open(ipath, "w", newline="", encoding="utf-8") as fh:
        wcsv = csv.writer(fh, lineterminator="\n")
        wcsv.writerow(["subjectID", "z_CL_caz", "z_CL_avi", "z_V_caz", "z_V_avi",
                       "CL_caz_L_h", "CL_avi_L_h", "V_caz_L", "V_avi_L"])
        for s, e in zip(subjects, etas):
            wcsv.writerow([s.sid] + [f"{x:.6f}" for x in e] + [
                f"{theta[0] * np.exp(w[0] * e[0]):.6f}", f"{theta[1] * np.exp(w[1] * e[1]):.6f}",
                f"{theta[2] * np.exp(w[2] * e[2]):.6f}", f"{theta[3] * np.exp(w[3] * e[3]):.6f}"])
    print(f"  wrote {ipath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
