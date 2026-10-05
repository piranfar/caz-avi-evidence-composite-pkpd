"""Simulation-based checks of the Model 1 estimator (primary model: 4 deviates plus the residual
correlation between the 2 drugs measured in the same sample).

PART A - ESTIMATOR RECOVERY
    Simulate replicate data sets from the fitted primary model at a known clearance correlation, on the
    observed design (same patients, sampling times and infusion durations), and refit each with the
    primary estimator. Report bias, root mean squared error and their Monte Carlo standard errors.
    This tests the estimator, not the biology: recovering the value that was put in is a check on the
    software and the design, not evidence about patients.

PART B - CALIBRATION OF THE LIKELIHOOD-RATIO TEST FOR RHO = 0.94
    In the data sets simulated at rho = 0.94, each replicate is refitted a second time with rho fixed at
    0.94. The share of replicates whose likelihood-ratio statistic exceeds the chi-square(1) critical
    value is the empirical type I error of the test used in the paper.

PART C - STRUCTURAL MISSPECIFICATION
    Simulate from a 2-compartment model and fit the 1-compartment model. With 5 to 7 samples per drug,
    starting at the trough, a distribution phase is unlikely to be identifiable, so this measures what
    the 1-compartment assumption costs if it is wrong.

Each replicate draws from its own random stream, seeded from (SEED, scenario, replicate), so replicates
are independent of run order and of the number of worker processes. Values are written at full precision
(6 decimals) and rounded once, when they are reported.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.stats import chi2

import joint_popk_nlme as M

sys.stdout.reconfigure(encoding="utf-8")

OUT = M.OUT
SEED = 20260811
CRIT = float(chi2.ppf(0.95, 1))
WORKERS = int(os.environ.get("SIM_WORKERS", max(1, min(14, (os.cpu_count() or 2) - 2))))
SCENARIOS = {                       # name: (scenario code for the seed, rho_true or None = estimate, 2-cmt)
    "correctly specified": (1, None, False),
    "correct, rho = 0.94": (2, 0.94, False),
    "two-compartment truth": (3, None, True),
}
CHECKPOINT = os.path.join(OUT, "model1_sbc_replicates.csv")
FIELDS = ["scenario", "replicate", "rho_true", "rho_estimated", "r_v_estimated", "c_estimated",
          "cl_caz", "cl_avi", "ofv", "ofv_rho_fixed_094", "lrt_rho_094"]


def fitted():
    """Primary-model estimates on the estimation scale, from model1_joint_popk_parameters.csv."""
    with open(os.path.join(OUT, "model1_joint_popk_parameters.csv"), encoding="utf-8") as fh:
        v = {r["parameter"]: float(r["estimate"]) for r in csv.DictReader(fh)}
    return np.array([np.log(v["CL ceftazidime (L/h)"]), np.log(v["CL avibactam (L/h)"]),
                     np.log(v["V ceftazidime (L)"]), np.log(v["V avibactam (L)"]),
                     np.log(v["omega_CL_caz"]), np.log(v["omega_CL_avi"]),
                     np.log(v["omega_V_caz"]), np.log(v["omega_V_avi"]),
                     np.arctanh(v["corr_CL_caz_avi"]), np.arctanh(v["corr_V_caz_avi"]),
                     np.log(v["sigma_prop_caz"]), np.log(v["sigma_prop_avi"]),
                     np.arctanh(v["corr_residual_caz_avi"])])


def log_pred_2cmt(subj, theta, z, w, q_frac=0.35, vp_frac=0.6):
    """2-compartment steady state by superposition of the biexponential impulse response.
    Peripheral volume and intercompartmental clearance are modest fractions of the central values,
    a plausible mild distribution phase rather than an extreme one."""
    cl = (theta[0] * np.exp(w[0] * z[0]), theta[1] * np.exp(w[1] * z[1]))
    v1 = (theta[2] * np.exp(w[2] * z[2]), theta[3] * np.exp(w[3] * z[3]))
    outs = []
    for j, a in enumerate(M.ANALYTES):
        CL, V1 = cl[j], v1[j]
        Q, V2 = q_frac * CL, vp_frac * V1
        k10, k12, k21 = CL / V1, Q / V1, Q / V2
        b = k10 + k12 + k21
        disc = np.sqrt(max(b * b - 4 * k10 * k21, 1e-12))
        alpha, beta = (b + disc) / 2, (b - disc) / 2
        A = (alpha - k21) / (V1 * (alpha - beta))
        B = (k21 - beta) / (V1 * (alpha - beta))
        rate = M.DOSE[a] / subj.t_inf
        t = np.asarray(subj.times[a], float)
        total = np.zeros_like(t)
        for i in range(60):
            s_ = t + i * M.TAU
            for amp, lam in ((A, alpha), (B, beta)):
                tin = np.minimum(s_, subj.t_inf)
                total += (rate * amp / lam) * (1.0 - np.exp(-lam * tin)) * np.exp(-lam * np.maximum(s_ - subj.t_inf, 0.0))
        outs.append(np.log(np.clip(total, 1e-10, None)))
    return np.concatenate(outs)


def simulate(subjects, x, rho, two_cmt, rng):
    """One data set on the observed design: correlated deviates, correlated residuals in each sample."""
    theta, om, w, _, _, sigma, c = M.unpack(x)
    om = om.copy()
    om[0, 1] = om[1, 0] = rho
    chol = np.linalg.cholesky(om)
    out = []
    for s in subjects:
        z = chol @ rng.standard_normal(4)
        lp = log_pred_2cmt(s, theta, z, w) if two_cmt else M.log_pred(s, theta, z, w)
        n = len(s.times["caz"])
        u1, u2 = rng.standard_normal(n), rng.standard_normal(n)
        e_caz = sigma[0] * u1
        e_avi = sigma[1] * (c * u1 + np.sqrt(1.0 - c * c) * u2)
        out.append(M.Subject(s.sid, s.t_inf, s.times, {"caz": lp[:n] + e_caz, "avi": lp[n:] + e_avi}))
    return out


def one(job):
    scenario, rep, x_hat = job
    code, rho_true, two_cmt = SCENARIOS[scenario]
    x_hat = np.asarray(x_hat, float)
    rho = float(np.tanh(x_hat[M.IDX_RHO])) if rho_true is None else rho_true
    rng = np.random.default_rng([SEED, code, rep])
    data = simulate(M.load(), x_hat, rho, two_cmt, rng)
    x0 = x_hat.copy()
    x0[M.IDX_RHO] = np.arctanh(rho)
    fr, _ = M.fit(data, M.build_omega4, M.N_OMEGA, "", x0, quiet=True, max_rounds=6)
    th, _, _, r_cl, r_v, _, c = M.unpack(fr.x)
    row = dict(scenario=scenario, replicate=rep, rho_true=rho, rho_estimated=r_cl, r_v_estimated=r_v,
               c_estimated=c, cl_caz=th[0], cl_avi=th[1], ofv=fr.fun, ofv_rho_fixed_094="", lrt_rho_094="")
    if rho_true == 0.94:
        f0, _ = M.fit(data, M.build_omega4, M.N_OMEGA, "", fr.x, quiet=True, max_rounds=6,
                      fixed={M.IDX_RHO: float(np.arctanh(0.94))})
        # The constrained fit cannot be better than the free one; if it is, the free fit stopped short.
        best = min(fr.fun, f0.fun)
        row.update(ofv=best, ofv_rho_fixed_094=f0.fun, lrt_rho_094=max(f0.fun - best, 0.0))
    return row


def load_done():
    """Stored replicates, one row per (scenario, replicate). Replicates are seeded by their key, so a
    key written twice (for example by an interrupted and a resumed run) must hold identical values."""
    if not os.path.exists(CHECKPOINT):
        return []
    with open(CHECKPOINT, encoding="utf-8") as fh:
        rr = list(csv.DictReader(fh))
    if not rr or list(rr[0]) != FIELDS:
        return []
    seen = {}
    for r in rr:
        k = (r["scenario"], int(r["replicate"]))
        if k in seen:
            assert seen[k]["rho_estimated"] == r["rho_estimated"], f"replicate {k} stored twice with different values"
            continue
        seen[k] = r
    return list(seen.values())


def append(row):
    new = not os.path.exists(CHECKPOINT)
    with open(CHECKPOINT, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        if new:
            w.writeheader()
        w.writerow({k: (f"{v:.6f}" if isinstance(v, (float, np.floating)) else v) for k, v in row.items()})


def summarise(rr, scenario):
    e = np.array([float(r["rho_estimated"]) for r in rr])
    t = float(rr[0]["rho_true"])
    n = len(e)
    out = dict(scenario=scenario, n_replicates=n, rho_true=t, mean_estimate=e.mean(), bias=e.mean() - t,
               bias_mcse=e.std(ddof=1) / np.sqrt(n), relative_bias_pct=100 * (e.mean() - t) / t,
               sd=e.std(ddof=1), rmse=float(np.sqrt(np.mean((e - t) ** 2))),
               p2_5=np.percentile(e, 2.5), p97_5=np.percentile(e, 97.5), pct_below_0_94=100 * np.mean(e < 0.94),
               lrt_type1_pct="", lrt_type1_mcse="")
    if scenario == "correct, rho = 0.94":
        lr = np.array([float(r["lrt_rho_094"]) for r in rr])
        p = float(np.mean(lr > CRIT))
        out.update(lrt_type1_pct=100 * p, lrt_type1_mcse=100 * np.sqrt(p * (1 - p) / n))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=50)
    args = ap.parse_args()
    x_hat = fitted()
    print("=" * 78)
    print("MODEL 1 simulation checks: recovery, LRT calibration, misspecification")
    print("=" * 78)
    print(f"  primary model, rho = {np.tanh(x_hat[M.IDX_RHO]):.4f}; {args.reps} replicates per scenario; "
          f"{WORKERS} workers")
    done = load_done()
    have = {(r["scenario"], int(r["replicate"])) for r in done}
    if not done and os.path.exists(CHECKPOINT):
        os.replace(CHECKPOINT, CHECKPOINT + ".previous_model")      # different layout: earlier model
    jobs = [(s, k, x_hat.tolist()) for s in SCENARIOS for k in range(1, args.reps + 1) if (s, k) not in have]
    print(f"  {len(have)} replicate(s) already stored, {len(jobs)} to run")
    if jobs:
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            futs = [ex.submit(one, j) for j in jobs]
            for i, f in enumerate(as_completed(futs), 1):
                row = f.result()
                append(row)
                if i % 10 == 0 or i == len(jobs):
                    print(f"    {i}/{len(jobs)} done", flush=True)
    rr = load_done()
    summary = [summarise([r for r in rr if r["scenario"] == s and int(r["replicate"]) <= args.reps], s)
               for s in SCENARIOS]
    a = next(s for s in summary if s["scenario"] == "correctly specified")
    b = next(s for s in summary if s["scenario"] == "two-compartment truth")
    ea = np.array([float(r["rho_estimated"]) for r in rr if r["scenario"] == "correctly specified"])
    eb = np.array([float(r["rho_estimated"]) for r in rr if r["scenario"] == "two-compartment truth"])
    diff, diff_se = eb.mean() - ea.mean(), float(np.sqrt(ea.var(ddof=1) / ea.size + eb.var(ddof=1) / eb.size))
    path = os.path.join(OUT, "model1_sbc_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]) + ["note"], lineterminator="\n")
        w.writeheader()
        for s in summary:
            note = (f"2-compartment minus correct mean estimate {diff:.6f} (SE {diff_se:.6f})"
                    if s["scenario"] == "two-compartment truth" else "")
            w.writerow({**{k: (f"{v:.6f}" if isinstance(v, (float, np.floating)) else v) for k, v in s.items()},
                        "note": note})
    print(f"\n  wrote {path}")
    for s in summary:
        extra = (f"  LRT type I error {s['lrt_type1_pct']:.1f}% (MCSE {s['lrt_type1_mcse']:.1f})"
                 if s["lrt_type1_pct"] != "" else "")
        print(f"  {s['scenario']:24} true {s['rho_true']:.4f}  mean {s['mean_estimate']:.4f}  "
              f"bias {s['bias']:+.4f} (MCSE {s['bias_mcse']:.4f}, {s['relative_bias_pct']:+.1f}%)  "
              f"SD {s['sd']:.4f}{extra}")
    print(f"  2-compartment minus correct: {diff:+.4f} (SE {diff_se:.4f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
