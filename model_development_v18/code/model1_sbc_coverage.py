"""Coverage of the profile-likelihood interval and the null distribution of the change in OFV at rho = 0.94,
from 500 simulated data sets per scenario (primary model, observed design).

Scenario "correctly specified": true rho = the estimate. Each replicate is fitted 3 times: free, with rho fixed
at the true value (coverage: the 95% profile-likelihood interval contains the true value when the change in OFV
at the true value is at most 3.841), and with rho fixed at 0.94 (how often these data would reject 0.94).
Scenario "correct, rho = 0.94": true rho = 0.94; fitted free and with rho fixed at 0.94, which gives the
coverage, the type I error and the empirical distribution of the change in OFV under the null.

Data sets are generated exactly as in model1_sbc.py (same seed and stream per scenario and replicate), so
replicates 1-50 must reproduce its stored free-fit estimates; main() checks this before summarising.
Values are written at full precision (6 decimals) and rounded once, when they are reported.
Outputs: outputs/model1_sbc_coverage_replicates.csv (checkpoint), outputs/model1_sbc_coverage_summary.csv
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.stats import chi2

import joint_popk_nlme as M
import model1_sbc as S

sys.stdout.reconfigure(encoding="utf-8")

CRIT = float(chi2.ppf(0.95, 1))
WORKERS = int(os.environ.get("SIM_WORKERS", "10"))
SCEN = ("correct, rho = 0.94", "correctly specified")
CHECKPOINT = os.path.join(M.OUT, "model1_sbc_coverage_replicates.csv")
FIELDS = ["scenario", "replicate", "rho_true", "rho_estimated", "c_estimated", "ofv", "ofv_rho_fixed_true",
          "lrt_rho_true", "ofv_rho_fixed_094", "lrt_rho_094", "seconds"]


def fit_fast(subjects, p0, fixed=None, max_rounds=6, tol=1e-4):
    """M.fit with the same objective and stopping rule, but each round starts with the quasi-Newton step
    (L-BFGS-B) and polishes with a short Nelder-Mead pass. Simulated data sets start at their true values,
    so the 2,500-evaluation exploratory simplex of M.fit is not needed; check_against_sbc() compares the
    result with the stored M.fit estimates of the same replicates."""
    from scipy.optimize import minimize
    fixed = fixed or {}
    free = [i for i in range(len(p0)) if i not in fixed]

    def expand(xf):
        x = np.array(p0, float)
        x[free] = xf
        for i, v in fixed.items():
            x[i] = v
        return x

    cache = {}

    def obj(xf):
        return M.ofv(expand(xf), subjects, M.build_omega4, M.N_OMEGA, cache)

    x = np.array([p0[i] for i in free], float)
    prev = np.inf
    for _ in range(max_rounds):
        r1 = minimize(obj, x, method="L-BFGS-B", options={"maxiter": 200, "ftol": 1e-12, "gtol": 1e-9, "eps": 1e-5})
        r2 = minimize(obj, r1.x, method="Nelder-Mead",
                      options={"maxiter": 600, "maxfev": 600, "xatol": 1e-6, "fatol": 1e-6, "adaptive": True})
        x, cur = (r2.x, r2.fun) if r2.fun < r1.fun else (r1.x, r1.fun)
        if prev - cur < tol:
            break
        prev = cur
    xx = expand(x)
    return xx, M.ofv(xx, subjects, M.build_omega4, M.N_OMEGA, {})


def one(job):
    scenario, rep, x_hat = job
    t0 = time.time()
    code, rho_true, two_cmt = S.SCENARIOS[scenario]
    x_hat = np.asarray(x_hat, float)
    rho = float(np.tanh(x_hat[M.IDX_RHO])) if rho_true is None else rho_true
    rng = np.random.default_rng([S.SEED, code, rep])
    data = S.simulate(M.load(), x_hat, rho, two_cmt, rng)
    x0 = x_hat.copy()
    x0[M.IDX_RHO] = np.arctanh(rho)
    xf, f_free = fit_fast(data, x0)
    _, _, _, r_cl, _, _, c = M.unpack(xf)
    _, f_true = fit_fast(data, xf, fixed={M.IDX_RHO: float(np.arctanh(rho))})
    best = min(f_free, f_true)                  # a constrained fit below the free one means the free fit stopped short
    is_094 = abs(rho - 0.94) < 1e-9
    return dict(scenario=scenario, replicate=rep, rho_true=rho, rho_estimated=r_cl, c_estimated=c, ofv=best,
                ofv_rho_fixed_true=f_true, lrt_rho_true=max(f_true - best, 0.0),
                ofv_rho_fixed_094=f_true if is_094 else "", lrt_rho_094=max(f_true - best, 0.0) if is_094 else "",
                seconds=time.time() - t0)


def load_done():
    if not os.path.exists(CHECKPOINT):
        return []
    with open(CHECKPOINT, encoding="utf-8") as fh:
        rr = list(csv.DictReader(fh))
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


def observed_lrt_094():
    with open(os.path.join(M.OUT, "model1_final_parameters.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["parameter"] == "dOFV_rho_0.94":
                return float(r["estimate"])
    raise KeyError("dOFV_rho_0.94")


def check_against_sbc(rr):
    """Replicates 1-50 must reproduce the free-fit estimates stored by model1_sbc.py."""
    with open(os.path.join(M.OUT, "model1_sbc_replicates.csv"), encoding="utf-8") as fh:
        ref = {(r["scenario"], int(r["replicate"])): float(r["rho_estimated"]) for r in csv.DictReader(fh)}
    worst, n = 0.0, 0
    for r in rr:
        k = (r["scenario"], int(r["replicate"]))
        if k in ref:
            worst = max(worst, abs(float(r["rho_estimated"]) - ref[k]))
            n += 1
    print(f"  {n} replicates also in model1_sbc_replicates.csv; largest difference in rho {worst:.2e}")
    return n, worst


def summarise(rr, obs):
    out = []
    for s in SCEN:
        sub = [r for r in rr if r["scenario"] == s]
        n = len(sub)
        e = np.array([float(r["rho_estimated"]) for r in sub])
        lt = np.array([float(r["lrt_rho_true"]) for r in sub])
        t = float(sub[0]["rho_true"])
        cov = float(np.mean(lt <= CRIT))
        row = dict(scenario=s, n_replicates=n, rho_true=t, mean_estimate=e.mean(), bias=e.mean() - t,
                   bias_mcse=e.std(ddof=1) / np.sqrt(n), coverage_95_profile_pct=100 * cov,
                   coverage_mcse=100 * np.sqrt(cov * (1 - cov) / n), reject_094_pct="", reject_094_mcse="",
                   dofv_094_p50="", dofv_094_p90="", dofv_094_p95="", dofv_094_p99="", dofv_094_max="",
                   n_at_or_above_observed="", observed_dofv_094=obs, empirical_p_observed="",
                   median_seconds=float(np.median([float(r["seconds"]) for r in sub])))
        if abs(t - 0.94) < 1e-9:                # under the null the constrained fit is the fit at 0.94
            rej = float(np.mean(lt > CRIT))
            k_obs = int(np.sum(lt >= obs))
            q = np.percentile(lt, [50, 90, 95, 99])
            row.update(reject_094_pct=100 * rej, reject_094_mcse=100 * np.sqrt(rej * (1 - rej) / n),
                       dofv_094_p50=q[0], dofv_094_p90=q[1], dofv_094_p95=q[2], dofv_094_p99=q[3],
                       dofv_094_max=float(lt.max()), n_at_or_above_observed=k_obs,
                       empirical_p_observed=(k_obs + 1) / (n + 1))
        out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=500)
    args = ap.parse_args()
    x_hat = S.fitted()
    done = load_done()
    have = {(r["scenario"], int(r["replicate"])) for r in done}
    jobs = [(s, k, x_hat.tolist()) for s in SCEN for k in range(1, args.reps + 1) if (s, k) not in have]   # null first
    print(f"rho estimate {np.tanh(x_hat[M.IDX_RHO]):.4f}; {args.reps} replicates per scenario; {len(have)} stored, "
          f"{len(jobs)} to run; {WORKERS} workers", flush=True)
    t0 = time.time()
    if jobs:
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            futs = [ex.submit(one, j) for j in jobs]
            for i, f in enumerate(as_completed(futs), 1):
                append(f.result())
                if i % 20 == 0 or i == len(jobs):
                    el = time.time() - t0
                    print(f"  {i}/{len(jobs)} done, {el / 60:.1f} min, about {el / i * (len(jobs) - i) / 60:.0f} min left",
                          flush=True)
    rr = [r for r in load_done() if int(r["replicate"]) <= args.reps]
    check_against_sbc(rr)
    summ = summarise(rr, observed_lrt_094())
    path = os.path.join(M.OUT, "model1_sbc_coverage_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0]), lineterminator="\n")
        w.writeheader()
        for s in summ:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, (float, np.floating)) else v) for k, v in s.items()})
    print(f"wrote {os.path.relpath(path, os.path.dirname(M.OUT))}")
    for s in summ:
        extra = "" if s["reject_094_pct"] == "" else (
            f"  reject 0.94 {s['reject_094_pct']:.1f}%  dOFV(0.94) p95 {s['dofv_094_p95']:.2f} p99 "
            f"{s['dofv_094_p99']:.2f} max {s['dofv_094_max']:.2f}  >= observed {s['n_at_or_above_observed']}")
        print(f"  {s['scenario']:22} n {s['n_replicates']}  coverage {s['coverage_95_profile_pct']:.1f}% "
              f"(MCSE {s['coverage_mcse']:.1f}){extra}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
