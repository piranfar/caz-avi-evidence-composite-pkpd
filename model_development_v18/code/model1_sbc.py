"""Simulation-based validation of the Model 1 estimator, and a misspecification check.

PART A — ESTIMATOR VALIDATION
    Simulate replicate datasets from the fitted Model 1 at a KNOWN correlation, using
    the observed design exactly (same patients, same sampling times, same infusion
    durations), then refit each with the real estimator. Measure bias, root mean
    squared error and interval coverage.

    This tests the ESTIMATOR, not the biology. Recovering the value that was put in
    is a check on the software; it is not evidence about patients, and nothing here
    may be described as validation of the model against data.

PART B — STRUCTURAL MISSPECIFICATION
    Simulate from a TWO-compartment model and fit the one-compartment model. Model 1
    assumes one compartment and that assumption was never tested, because with 5-7
    samples per analyte beginning at the trough a distribution phase is unlikely to
    be identifiable. This measures the bias that assumption costs if it is wrong.
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import joint_popk_nlme as M   # noqa: E402

OUT = M.OUT
SEED = 20260811

# fitted Model 1
THETA = np.array([2.5723, 3.2234, 19.9770, 27.0298])
W = np.array([0.2022, 0.1411, 0.2673, 0.1994])
RHO = 0.7031
SIGMA = np.array([0.0997, 0.0926])
P0 = np.array([np.log(2.57), np.log(3.22), np.log(20.0), np.log(27.0),
               np.log(0.20), np.log(0.14), np.log(0.27), np.log(0.20), 0.9,
               np.log(0.10), np.log(0.093)])


def simulate_dataset(subjects, rng, rho=RHO, theta=THETA, w=W, sigma=SIGMA,
                     two_compartment=False):
    """Replace the observed concentrations with simulated ones, same design."""
    om = np.eye(3)
    om[0, 1] = om[1, 0] = rho
    chol = np.linalg.cholesky(om)
    out = []
    for s in subjects:
        z = chol @ rng.standard_normal(3)
        if two_compartment:
            logc = _log_pred_2cmt(s, theta, z, w)
        else:
            logc = M.log_pred(s, theta, z, w)
        k = 0
        newlog = {}
        for j, a in enumerate(M.ANALYTES):
            n = len(s.times[a])
            newlog[a] = logc[k:k + n] + rng.normal(0, sigma[j], n)
            k += n
        out.append(M.Subject(s.sid, s.t_inf, s.times, newlog))
    return out


def _log_pred_2cmt(subj, theta, z, w, q_frac=0.35, vp_frac=0.6):
    """Two-compartment steady state by superposition of the biexponential impulse.

    Peripheral compartment sized as a fraction of the central one and inter-
    compartmental clearance as a fraction of elimination clearance; both are
    deliberately modest, so this represents a plausible mild distribution phase
    rather than an extreme one.
    """
    cl = (theta[0] * np.exp(w[0] * z[0]), theta[1] * np.exp(w[1] * z[1]))
    v1 = (theta[2] * np.exp(w[2] * z[2]), theta[3] * np.exp(w[3] * z[2]))
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
                rise = (rate * amp / lam) * (1.0 - np.exp(-lam * tin))
                total += rise * np.exp(-lam * np.maximum(s_ - subj.t_inf, 0.0))
        outs.append(np.log(np.clip(total, 1e-10, None)))
    return np.concatenate(outs)


CHECKPOINT = os.path.join(OUT, "model1_sbc_replicates.csv")
FIELDS = ["replicate", "scenario", "rho_true", "rho_estimated",
          "cl_caz", "cl_avi", "ofv"]


def _load_checkpoint():
    """Replicates already completed, so an interrupted run resumes rather than restarts.

    Each replicate costs about 80 seconds, so losing a part-finished run is expensive.
    The random stream is reseeded per scenario and advanced deterministically, so a
    resumed run reproduces exactly the replicates a single uninterrupted run would.
    """
    if not os.path.exists(CHECKPOINT):
        return []
    with open(CHECKPOINT, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _append(row):
    os.makedirs(OUT, exist_ok=True)
    new_file = not os.path.exists(CHECKPOINT)
    with open(CHECKPOINT, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        if new_file:
            w.writeheader()
        w.writerow(row)


def run(n_rep, two_compartment=False, rho_true=RHO, label=""):
    subjects = M.load()
    rng = np.random.default_rng(SEED + (999 if two_compartment else 0)
                                + (7 if rho_true != RHO else 0))
    done = [r for r in _load_checkpoint() if r["scenario"] == label]
    rows = [{"replicate": int(r["replicate"]), "scenario": r["scenario"],
             "rho_true": float(r["rho_true"]),
             "rho_estimated": float(r["rho_estimated"]),
             "cl_caz": float(r["cl_caz"]), "cl_avi": float(r["cl_avi"]),
             "ofv": float(r["ofv"])} for r in done]
    if rows:
        print(f"    [{label}] resuming: {len(rows)} replicate(s) already complete",
              flush=True)

    for r in range(n_rep):
        sim = simulate_dataset(subjects, rng, rho=rho_true,
                               two_compartment=two_compartment)
        if r < len(rows):
            continue                      # stream advanced; result already stored
        fit, _ = M.fit(sim, M.build_omega, 5, "", P0, quiet=True, max_rounds=6)
        rho_hat = float(np.tanh(fit.x[8]))
        theta_hat = np.exp(fit.x[:4])
        row = {"replicate": r + 1, "scenario": label,
               "rho_true": rho_true, "rho_estimated": round(rho_hat, 4),
               "cl_caz": round(float(theta_hat[0]), 4),
               "cl_avi": round(float(theta_hat[1]), 4),
               "ofv": round(float(fit.fun), 3)}
        rows.append(row)
        _append(row)
        print(f"    [{label}] {r+1:3}/{n_rep}  rho_hat {rho_hat:.4f}", flush=True)
    return rows


def summarise(rows, rho_true, label):
    e = np.array([r["rho_estimated"] for r in rows])
    return {
        "scenario": label, "n_replicates": len(e), "rho_true": rho_true,
        "mean_estimate": round(float(e.mean()), 4),
        "bias": round(float(e.mean() - rho_true), 4),
        "relative_bias_pct": round(100 * float((e.mean() - rho_true) / rho_true), 2),
        "sd": round(float(e.std(ddof=1)), 4),
        "rmse": round(float(np.sqrt(np.mean((e - rho_true) ** 2))), 4),
        "p2.5": round(float(np.percentile(e, 2.5)), 4),
        "p97.5": round(float(np.percentile(e, 97.5)), 4),
        "pct_below_0.94": round(100 * float(np.mean(e < 0.94)), 1),
    }


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=40)
    args = ap.parse_args()

    print("=" * 74)
    print("MODEL 1 — estimator validation and misspecification check")
    print("=" * 74)
    print(f"  {args.reps} replicates per scenario, observed design reused exactly\n")

    all_rows, summary = [], []
    print("  A. correctly specified (one compartment simulated, one fitted)")
    rows = run(args.reps, False, RHO, "correctly specified")
    all_rows += rows
    summary.append(summarise(rows, RHO, "correctly specified"))

    print("\n  B. misspecified (two compartments simulated, one fitted)")
    rows = run(args.reps, True, RHO, "two-compartment truth")
    all_rows += rows
    summary.append(summarise(rows, RHO, "two-compartment truth"))

    print("\n  C. correctly specified at the published correlation")
    rows = run(args.reps, False, 0.94, "correct, rho = 0.94")
    all_rows += rows
    summary.append(summarise(rows, 0.94, "correct, rho = 0.94"))

    os.makedirs(OUT, exist_ok=True)
    for data, name in ((summary, "model1_sbc_summary.csv"),):
        p = os.path.join(OUT, name)
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]), lineterminator="\n")
            w.writeheader(); w.writerows(data)
        print(f"\n  wrote {os.path.relpath(p, os.path.dirname(M.HERE))}")

    print("\n" + "=" * 74)
    for s in summary:
        print(f"  {s['scenario']:24} true {s['rho_true']:.3f}  "
              f"estimate {s['mean_estimate']:.4f}  bias {s['bias']:+.4f} "
              f"({s['relative_bias_pct']:+.1f}%)  SD {s['sd']:.4f}  "
              f"RMSE {s['rmse']:.4f}")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
