"""MODEL 1 finalisation — convergence, profile-likelihood interval, diagnostics,
visual predictive check, and the prespecified sensitivity analyses.

Everything here is deterministic. The only stochastic step is the visual predictive
check, which uses a fixed seed recorded in the output.

Outputs (all under ../outputs and ../figures):
    model1_final_parameters.csv        converged estimates
    model1_profile_likelihood.csv      OFV against the fixed correlation
    model1_diagnostics.csv             per-observation PRED, IPRED, CWRES
    model1_vpc.csv                     observed and simulated percentiles by time bin
    model1_sensitivity.csv             infusion duration and structural assumptions
    model1_gof.png, model1_vpc.png
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2

import joint_popk_nlme as M

sys.stdout.reconfigure(encoding="utf-8")

OUT = M.OUT
FIG = os.path.join(os.path.dirname(M.HERE), "figures")
VPC_SEED = 20260811
P0 = np.array([np.log(2.57), np.log(3.22), np.log(20.0), np.log(27.0),
               np.log(0.20), np.log(0.14), np.log(0.27), np.log(0.20), 0.9,
               np.log(0.10), np.log(0.093)])
IDX_RHO = 8
NAMES = ["CL_caz_L_h", "CL_avi_L_h", "V_caz_L", "V_avi_L",
         "omega_CL_caz", "omega_CL_avi", "omega_V_caz", "omega_V_avi",
         "corr_CL_caz_avi", "sigma_prop_caz", "sigma_prop_avi"]


def unpack(x):
    theta = np.exp(x[:4])
    om, w, r_cl, _ = M.build_omega(x[4:9])
    sigma = np.exp(x[9:11])
    return theta, om, w, r_cl, sigma


# --------------------------------------------------------------- diagnostics -

def diagnostics(subjects, x, cache):
    """Per-observation PRED, IPRED and conditional weighted residuals.

    CWRES follows Hooker et al. (2007): linearise around the conditional estimate
    eta*, so that
        E[y_i]   = f(eta*) - J eta*
        Var[y_i] = J Omega J' + Sigma
        CWRES_i  = chol(Var)^-1 (y_i - E[y_i])
    """
    theta, om, w, r_cl, sigma = unpack(x)
    rows = []
    for s in subjects:
        y = M.obs_vector(s)
        sig = M.sigma_vector(s, sigma)
        eta = cache[s.sid]
        ipred = M.log_pred(s, theta, eta, w)
        pred = M.log_pred(s, theta, np.zeros(M.N_ETA), w)

        J = np.empty((len(y), M.N_ETA))
        h = 1e-5
        for a in range(M.N_ETA):
            e = eta.copy()
            e[a] += h
            J[:, a] = (M.log_pred(s, theta, e, w) - ipred) / h

        mean = ipred - J @ eta
        var = J @ om @ J.T + np.diag(sig ** 2)
        L = np.linalg.cholesky(var + 1e-12 * np.eye(len(y)))
        cwres = np.linalg.solve(L, y - mean)

        k = 0
        for j, an in enumerate(M.ANALYTES):
            for t in s.times[an]:
                rows.append({
                    "subjectID": s.sid, "analyte": an, "time_h": t,
                    "dv_mg_l": float(np.exp(y[k])),
                    "pred_mg_l": float(np.exp(pred[k])),
                    "ipred_mg_l": float(np.exp(ipred[k])),
                    "cwres": float(cwres[k]),
                    "iwres": float((y[k] - ipred[k]) / sig[k]),
                })
                k += 1
    return rows


# ------------------------------------------------------ visual predictive check

def vpc(subjects, x, n_rep=1000, seed=VPC_SEED):
    """Simulate the observed design n_rep times and summarise by time bin."""
    theta, om, w, r_cl, sigma = unpack(x)
    rng = np.random.default_rng(seed)
    chol = np.linalg.cholesky(om)

    sim = {a: {} for a in M.ANALYTES}
    obs = {a: {} for a in M.ANALYTES}
    for s in subjects:
        for j, an in enumerate(M.ANALYTES):
            for t, c in zip(s.times[an], np.exp(s.logconc[an])):
                obs[an].setdefault(round(float(t), 3), []).append(float(c))

    for _ in range(n_rep):
        for s in subjects:
            z = chol @ rng.standard_normal(M.N_ETA)
            lp = M.log_pred(s, theta, z, w)
            k = 0
            for j, an in enumerate(M.ANALYTES):
                for t in s.times[an]:
                    c = float(np.exp(lp[k] + rng.normal(0.0, sigma[j])))
                    sim[an].setdefault(round(float(t), 3), []).append(c)
                    k += 1

    rows = []
    for an in M.ANALYTES:
        for t in sorted(obs[an]):
            o = np.array(obs[an][t])
            m = np.array(sim[an][t])
            rows.append({
                "analyte": an, "time_h": t, "n_obs": len(o),
                "obs_p5": np.percentile(o, 5), "obs_p50": np.percentile(o, 50),
                "obs_p95": np.percentile(o, 95),
                "sim_p5": np.percentile(m, 5), "sim_p50": np.percentile(m, 50),
                "sim_p95": np.percentile(m, 95),
                "pct_obs_within_sim_90": 100.0 * float(np.mean(
                    (o >= np.percentile(m, 5)) & (o <= np.percentile(m, 95)))),
            })
    return rows


# ------------------------------------------------------------------- figures -

def make_figures(diag, vpc_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(FIG, exist_ok=True)
    labels = {"caz": "Ceftazidime", "avi": "Avibactam"}
    colours = {"caz": "#1F4E85", "avi": "#C86438"}

    fig, ax = plt.subplots(2, 3, figsize=(13.5, 8))
    for r, an in enumerate(M.ANALYTES):
        d = [x for x in diag if x["analyte"] == an]
        dv = np.array([x["dv_mg_l"] for x in d])
        pr = np.array([x["pred_mg_l"] for x in d])
        ip = np.array([x["ipred_mg_l"] for x in d])
        cw = np.array([x["cwres"] for x in d])
        tm = np.array([x["time_h"] for x in d])

        for c, (xv, xl) in enumerate(((pr, "Population prediction (mg/L)"),
                                      (ip, "Individual prediction (mg/L)"))):
            a = ax[r, c]
            lim = [min(dv.min(), xv.min()) * 0.8, max(dv.max(), xv.max()) * 1.2]
            a.plot(lim, lim, color="0.4", lw=1)
            a.scatter(xv, dv, s=18, alpha=0.75, color=colours[an], edgecolor="none")
            a.set_xscale("log"); a.set_yscale("log")
            a.set_xlim(lim); a.set_ylim(lim)
            a.set_xlabel(xl); a.set_ylabel("Observed (mg/L)")
            a.set_title(f"{labels[an]}")

        a = ax[r, 2]
        a.axhline(0, color="0.4", lw=1)
        for h in (-2, 2):
            a.axhline(h, color="0.7", lw=0.8, ls="--")
        a.scatter(tm, cw, s=18, alpha=0.75, color=colours[an], edgecolor="none")
        a.set_xlabel("Time after dose (h)")
        a.set_ylabel("Conditional weighted residual")
        a.set_title(f"{labels[an]}  (SD {cw.std(ddof=1):.2f})")
    fig.suptitle("Model 1 goodness of fit — CRRT cohort, 21 patients", y=0.99)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "model1_gof.png"), dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for i, an in enumerate(M.ANALYTES):
        v = [x for x in vpc_rows if x["analyte"] == an]
        t = np.array([x["time_h"] for x in v])
        a = ax[i]
        a.fill_between(t, [x["sim_p5"] for x in v], [x["sim_p95"] for x in v],
                       color=colours[an], alpha=0.18,
                       label="simulated 5th-95th percentile")
        a.plot(t, [x["sim_p50"] for x in v], color=colours[an], lw=2,
               label="simulated median")
        a.plot(t, [x["obs_p50"] for x in v], "o--", color="0.15", ms=5, lw=1.2,
               label="observed median")
        a.plot(t, [x["obs_p5"] for x in v], ".", color="0.45", ms=6)
        a.plot(t, [x["obs_p95"] for x in v], ".", color="0.45", ms=6,
               label="observed 5th and 95th")
        a.set_xlabel("Time after dose (h)")
        a.set_ylabel("Concentration (mg/L)")
        a.set_title(labels[an])
        a.legend(fontsize=8, frameon=False)
    fig.suptitle("Model 1 visual predictive check — 1000 replicates of the observed design", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "model1_vpc.png"), dpi=300)
    plt.close(fig)
    print(f"  wrote {os.path.join(FIG, 'model1_gof.png')}")
    print(f"  wrote {os.path.join(FIG, 'model1_vpc.png')}")


def write_csv(rows, name):
    path = os.path.join(OUT, name)
    os.makedirs(OUT, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path} ({len(rows)} rows)")


# ---------------------------------------------------------------------- main -

def main():
    M.structural_self_check()
    subjects = M.load()
    nobs = sum(len(M.obs_vector(s)) for s in subjects)
    print("=" * 78)
    print("MODEL 1 FINALISATION")
    print("=" * 78)
    print(f"  {len(subjects)} subjects, {nobs} observations\n")

    # 1 -------------------------------------------------------- convergence --
    print("1. Convergence")
    full, cache = M.fit(subjects, M.build_omega, 5, "full model", P0)
    null, _ = M.fit(subjects, M.build_omega_diag, 4, "null (correlation fixed at 0)",
                    np.concatenate([full.x[:8], full.x[9:]]))
    theta, om, w, r_cl, sigma = unpack(full.x)
    d_ofv = null.fun - full.fun
    print(f"    dOFV against the null {d_ofv:.3f} on 1 df, p = {chi2.sf(d_ofv,1):.4g}")
    print(f"    correlation estimate {r_cl:.4f}")

    # 2 ------------------------------------------------- profile likelihood --
    print("\n2. Profile likelihood for the clearance correlation")
    print("   (the Hessian-based interval is replaced by this)")
    grid = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.703,
                     0.75, 0.80, 0.85, 0.88, 0.91, 0.94, 0.96])
    prof = []
    for r in grid:
        z = np.arctanh(np.clip(r, -0.999, 0.999))
        fr, _ = M.fit(subjects, M.build_omega, 5, f"rho={r:.3f}", full.x,
                      quiet=True, fixed={IDX_RHO: z}, max_rounds=8)
        prof.append({"rho": float(r), "ofv": float(fr.fun),
                     "delta_ofv": float(fr.fun - full.fun)})
        print(f"    rho {r:5.3f}   OFV {fr.fun:10.4f}   dOFV {fr.fun-full.fun:8.4f}")

    crit = chi2.ppf(0.95, 1)          # 3.841
    rhos = np.array([p["rho"] for p in prof])
    dofv = np.array([p["delta_ofv"] for p in prof])

    def crossing(lo_side):
        idx = np.where(rhos < r_cl)[0] if lo_side else np.where(rhos > r_cl)[0]
        if len(idx) == 0:
            return float("nan")
        rr, dd = rhos[idx], dofv[idx]
        order = np.argsort(rr)
        rr, dd = rr[order], dd[order]
        for k in range(len(rr) - 1):
            a, b = (k, k + 1) if lo_side else (len(rr) - 2 - k, len(rr) - 1 - k)
            if (dd[a] - crit) * (dd[b] - crit) <= 0 and dd[b] != dd[a]:
                return float(rr[a] + (crit - dd[a]) * (rr[b] - rr[a]) / (dd[b] - dd[a]))
        return float(rr[0] if lo_side else rr[-1])

    pl_lo, pl_hi = crossing(True), crossing(False)
    print(f"\n    profile-likelihood 95% interval: {pl_lo:.3f} to {pl_hi:.3f}")
    print(f"    excludes 0.94: {'YES' if pl_hi < 0.94 else 'NO'}")
    for p in prof:
        p["in_95_interval"] = "yes" if p["delta_ofv"] <= crit else "no"
    write_csv(prof, "model1_profile_likelihood.csv")

    # 3 -------------------------------------------------------- diagnostics --
    print("\n3. Goodness-of-fit diagnostics")
    diag = diagnostics(subjects, full.x, cache)
    for an in M.ANALYTES:
        cw = np.array([d["cwres"] for d in diag if d["analyte"] == an])
        out = 100.0 * float(np.mean(np.abs(cw) > 2))
        print(f"    {an}: CWRES mean {cw.mean():+.3f}, SD {cw.std(ddof=1):.3f}, "
              f"{out:.1f}% beyond +/-2")
    write_csv(diag, "model1_diagnostics.csv")

    # 4 --------------------------------------------- visual predictive check --
    print(f"\n4. Visual predictive check (1000 replicates, seed {VPC_SEED})")
    vrows = vpc(subjects, full.x)
    cov = np.mean([r["pct_obs_within_sim_90"] for r in vrows])
    print(f"    observations inside the simulated 5th-95th interval: {cov:.1f}% "
          f"(nominal 90%)")
    write_csv(vrows, "model1_vpc.csv")

    # 5 --------------------------------------------------------- sensitivity --
    print("\n5. Sensitivity analyses")
    sens = [{"analysis": "reference", "variant": "as fitted",
             "corr_CL": round(r_cl, 4), "ofv": round(full.fun, 3),
             "note": "infusion duration read as hours from the dataset"}]

    original = {s.sid: s.t_inf for s in subjects}
    for fixed_t in (1.0, 2.0, 3.0):
        alt = [M.Subject(s.sid, fixed_t, s.times, s.logconc) for s in subjects]
        fr, _ = M.fit(alt, M.build_omega, 5, "", full.x, quiet=True, max_rounds=8)
        _, _, _, rr, _ = unpack(fr.x)
        sens.append({"analysis": "infusion duration", "variant": f"all {fixed_t:.0f} h",
                     "corr_CL": round(rr, 4), "ofv": round(fr.fun, 3),
                     "note": "tests the categorical-versus-hours ambiguity in the source"})
        print(f"    infusion duration all {fixed_t:.0f} h -> corr {rr:.3f}")

    # residual error model: additive-on-log is the reference; test analyte-shared sigma
    def build_omega_shared(p):
        return M.build_omega(p)

    x_shared = full.x.copy()
    mean_sig = 0.5 * (x_shared[9] + x_shared[10])
    x_shared[9] = x_shared[10] = mean_sig
    fr, _ = M.fit(subjects, M.build_omega, 5, "", x_shared, quiet=True, max_rounds=8,
                  fixed={10: mean_sig})
    _, _, _, rr, _ = unpack(fr.x)
    sens.append({"analysis": "residual error", "variant": "shared across analytes",
                 "corr_CL": round(rr, 4), "ofv": round(fr.fun, 3),
                 "note": f"dOFV {fr.fun - full.fun:+.2f} on 1 df"})
    print(f"    shared residual error -> corr {rr:.3f}, dOFV {fr.fun-full.fun:+.2f}")

    # leave-one-subject-out influence on the correlation
    loo = []
    for s in subjects:
        rest = [t for t in subjects if t.sid != s.sid]
        fr, _ = M.fit(rest, M.build_omega, 5, "", full.x, quiet=True, max_rounds=6)
        _, _, _, rr, _ = unpack(fr.x)
        loo.append(rr)
    loo = np.array(loo)
    print(f"    leave-one-subject-out correlation: min {loo.min():.3f}, "
          f"median {np.median(loo):.3f}, max {loo.max():.3f}")
    sens.append({"analysis": "leave-one-subject-out", "variant": "range over 21 refits",
                 "corr_CL": f"{loo.min():.3f} to {loo.max():.3f}",
                 "ofv": "", "note": f"median {np.median(loo):.3f}; "
                                    f"most influential subject "
                                    f"{subjects[int(np.argmax(np.abs(loo - r_cl)))].sid}"})
    write_csv(sens, "model1_sensitivity.csv")

    # 6 ------------------------------------------------------------ figures --
    print("\n6. Figures")
    make_figures(diag, vrows)

    # 7 ------------------------------------------------- final parameter table
    print("\n7. Final parameter table")
    vals = list(theta) + list(w) + [r_cl] + list(sigma)
    rows = []
    for nm, v in zip(NAMES, vals):
        rows.append({"parameter": nm, "estimate": f"{v:.4f}",
                     "ci_low": f"{pl_lo:.4f}" if nm == "corr_CL_caz_avi" else "",
                     "ci_high": f"{pl_hi:.4f}" if nm == "corr_CL_caz_avi" else "",
                     "ci_method": "profile likelihood" if nm == "corr_CL_caz_avi" else "",
                     "note": ""})
    rows.append({"parameter": "OFV", "estimate": f"{full.fun:.4f}", "ci_low": "",
                 "ci_high": "", "ci_method": "",
                 "note": f"converged on tolerance: {full.converged}, {full.rounds} rounds"})
    rows.append({"parameter": "dOFV_vs_no_correlation", "estimate": f"{d_ofv:.4f}",
                 "ci_low": "", "ci_high": "", "ci_method": "",
                 "note": f"1 df, p = {chi2.sf(d_ofv,1):.4g}"})
    rows.append({"parameter": "vpc_coverage_pct", "estimate": f"{cov:.1f}", "ci_low": "",
                 "ci_high": "", "ci_method": "",
                 "note": f"nominal 90%, seed {VPC_SEED}"})
    write_csv(rows, "model1_final_parameters.csv")

    print("\n" + "=" * 78)
    print(f"  correlation {r_cl:.3f}, profile-likelihood 95% interval "
          f"{pl_lo:.3f} to {pl_hi:.3f}")
    print(f"  assumed value 0.94 is {'EXCLUDED' if pl_hi < 0.94 else 'NOT excluded'}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
