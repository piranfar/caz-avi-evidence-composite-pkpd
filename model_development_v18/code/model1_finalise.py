"""MODEL 1 finalisation: convergence, profile-likelihood interval, diagnostics, visual predictive
check, and the prespecified sensitivity analyses, for the primary model (4 deviates plus the
residual correlation between the 2 drugs measured in the same sample).

Everything here is deterministic. The only stochastic step is the visual predictive check,
which uses a fixed seed recorded in the output. Independent refits (profile points, sensitivity
analyses, leave-one-out) run in parallel processes; each starts from the same values with its own
cache, so the results do not depend on the order in which they finish.

Outputs (all under ../outputs and ../figures):
    model1_final_parameters.csv        converged estimates (6 decimals)
    model1_profile_likelihood.csv      OFV against the fixed clearance correlation
    model1_diagnostics.csv             per-observation PRED, IPRED, CWRES, IWRES
    model1_vpc.csv                     observed and simulated percentiles by time bin
    model1_sensitivity.csv             infusion duration, variance structure, residual error, leave-one-out
    model1_gof.png, model1_vpc.png
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import csv
import glob
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.stats import chi2

import joint_popk_nlme as M

sys.stdout.reconfigure(encoding="utf-8")

OUT = M.OUT
FIG = os.path.join(os.path.dirname(M.HERE), "figures")
VPC_SEED = 20260811
P0 = M.P0
N_OMEGA = M.N_OMEGA
IDX_RHO, IDX_RV, IDX_C = M.IDX_RHO, M.IDX_RV, M.IDX_C
CRIT = float(chi2.ppf(0.95, 1))          # 3.841
WORKERS = max(1, min(14, (os.cpu_count() or 2) - 2))
NAMES = ["CL_caz_L_h", "CL_avi_L_h", "V_caz_L", "V_avi_L",
         "omega_CL_caz", "omega_CL_avi", "omega_V_caz", "omega_V_avi",
         "corr_CL_caz_avi", "corr_V_caz_avi", "sigma_prop_caz", "sigma_prop_avi",
         "corr_residual_caz_avi"]


def unpack(x):
    return M.unpack(x)


# ------------------------------------------------------------ parallel refits -

def _refit(job):
    """One independent refit. job = (label, structure, x0, fixed, t_inf_override, drop_sid, max_rounds)."""
    label, structure, x0, fixed, t_inf, drop, rounds = job
    subjects = M.load()
    if t_inf is not None:
        subjects = [M.Subject(s.sid, t_inf, s.times, s.logconc) for s in subjects]
    if drop is not None:
        subjects = [s for s in subjects if s.sid != drop]
    builder, n_omega = {"4dev": (M.build_omega4, 6), "3dev": (M.build_omega, 5)}[structure]
    fr, _ = M.fit(subjects, builder, n_omega, label, np.asarray(x0, float), quiet=True,
                  max_rounds=rounds, fixed=fixed)
    return dict(label=label, ofv=float(fr.fun), rho=float(np.tanh(fr.x[8])), x=fr.x.tolist(),
                converged=bool(fr.converged))


# Optional resume cache. When MODEL1_REFIT_CACHE names a directory, each finished refit is stored
# there under a key built from the job, the estimation code and the data, so an interrupted run
# resumes instead of restarting. Results are identical with or without the cache.
CACHE_DIR = os.environ.get("MODEL1_REFIT_CACHE", "")


def _fingerprint():
    h = hashlib.sha256()
    for p in [os.path.join(M.HERE, "joint_popk_nlme.py")] + sorted(glob.glob(os.path.join(M.DATA, "*.csv"))):
        with open(p, "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()


def _key(job, fp):
    label, structure, x0, fixed, t_inf, drop, rounds = job
    spec = [label, structure, [round(float(v), 12) for v in x0],
            sorted((int(k), round(float(v), 12)) for k, v in (fixed or {}).items()), t_inf, drop, rounds, fp]
    return hashlib.sha256(json.dumps(spec).encode()).hexdigest()[:32]


def run_jobs(jobs):
    """Run independent refits in parallel; results come back in the order of `jobs`."""
    results = [None] * len(jobs)
    fp = _fingerprint() if CACHE_DIR else ""
    todo = []
    for i, job in enumerate(jobs):
        path = os.path.join(CACHE_DIR, _key(job, fp) + ".json") if CACHE_DIR else ""
        if path and os.path.exists(path):
            with open(path, encoding="utf-8") as fh:
                results[i] = json.load(fh)
        else:
            todo.append((i, path))
    if todo:
        with ProcessPoolExecutor(max_workers=min(WORKERS, len(todo))) as ex:
            futs = {ex.submit(_refit, jobs[i]): (i, path) for i, path in todo}
            for f in as_completed(futs):
                i, path = futs[f]
                results[i] = f.result()
                if path:
                    os.makedirs(CACHE_DIR, exist_ok=True)
                    with open(path, "w", encoding="utf-8") as fh:
                        json.dump(results[i], fh)
    return results


# --------------------------------------------------------------- diagnostics -

def residual_cov(s, sig, c):
    """Residual covariance of one subject: proportional SDs, correlation c within each sample."""
    n_caz = len(s.times["caz"])
    S = np.diag(sig ** 2)
    if c != 0.0:
        for i in range(n_caz):
            S[i, n_caz + i] = S[n_caz + i, i] = c * sig[i] * sig[n_caz + i]
    return S


def diagnostics(subjects, x, cache):
    """Per-observation PRED, IPRED, CWRES (Hooker et al. 2007) and IWRES.

    Linearise around the conditional estimate eta*:
        E[y_i] = f(eta*) - J eta*,  Var[y_i] = J Omega J' + Sigma_i,  CWRES_i = chol(Var)^-1 (y_i - E[y_i])
    Sigma_i carries the residual correlation between the 2 drugs in the same sample.
    """
    theta, om, w, r_cl, r_v, sigma, c = unpack(x)
    rows = []
    for s in subjects:
        y = M.obs_vector(s)
        sig = M.sigma_vector(s, sigma)
        eta = cache[s.sid]
        ipred = M.log_pred(s, theta, eta, w)
        pred = M.log_pred(s, theta, np.zeros(len(eta)), w)
        J = np.empty((len(y), len(eta)))
        h = 1e-5
        for a in range(len(eta)):
            e = eta.copy()
            e[a] += h
            J[:, a] = (M.log_pred(s, theta, e, w) - ipred) / h
        mean = ipred - J @ eta
        var = J @ om @ J.T + residual_cov(s, sig, c)
        L = np.linalg.cholesky(var + 1e-12 * np.eye(len(y)))
        cwres = np.linalg.solve(L, y - mean)
        k = 0
        for an in M.ANALYTES:
            for t in s.times[an]:
                rows.append({"subjectID": s.sid, "analyte": an, "time_h": t,
                             "dv_mg_l": float(np.exp(y[k])), "pred_mg_l": float(np.exp(pred[k])),
                             "ipred_mg_l": float(np.exp(ipred[k])), "cwres": float(cwres[k]),
                             "iwres": float((y[k] - ipred[k]) / sig[k])})
                k += 1
    return rows


# ------------------------------------------------------ visual predictive check

def vpc(subjects, x, n_rep=1000, seed=VPC_SEED):
    """Simulate the observed design n_rep times (correlated residuals) and summarise by time bin."""
    theta, om, w, r_cl, r_v, sigma, c = unpack(x)
    rng = np.random.default_rng(seed)
    chol = np.linalg.cholesky(om)
    sim = {a: {} for a in M.ANALYTES}
    obs = {a: {} for a in M.ANALYTES}
    for s in subjects:
        for an in M.ANALYTES:
            for t, conc in zip(s.times[an], np.exp(s.logconc[an])):
                obs[an].setdefault(round(float(t), 3), []).append(float(conc))
    for _ in range(n_rep):
        for s in subjects:
            z = chol @ rng.standard_normal(om.shape[0])
            lp = M.log_pred(s, theta, z, w)
            n = len(s.times["caz"])
            u1, u2 = rng.standard_normal(n), rng.standard_normal(n)
            eps = np.concatenate([sigma[0] * u1, sigma[1] * (c * u1 + np.sqrt(1.0 - c * c) * u2)])
            k = 0
            for an in M.ANALYTES:
                for t in s.times[an]:
                    sim[an].setdefault(round(float(t), 3), []).append(float(np.exp(lp[k] + eps[k])))
                    k += 1
    rows = []
    for an in M.ANALYTES:
        for t in sorted(obs[an]):
            o, m = np.array(obs[an][t]), np.array(sim[an][t])
            lo, hi = np.percentile(m, 5), np.percentile(m, 95)
            rows.append({"analyte": an, "time_h": t, "n_obs": len(o),
                         "n_obs_within_sim_90": int(np.sum((o >= lo) & (o <= hi))),
                         "obs_p5": np.percentile(o, 5), "obs_p50": np.percentile(o, 50),
                         "obs_p95": np.percentile(o, 95),
                         "sim_p5": lo, "sim_p50": np.percentile(m, 50), "sim_p95": hi,
                         "pct_obs_within_sim_90": 100.0 * float(np.mean((o >= lo) & (o <= hi)))})
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
        for col, (xv, xl) in enumerate(((pr, "Population prediction (mg/L)"),
                                        (ip, "Individual prediction (mg/L)"))):
            a = ax[r, col]
            lim = [min(dv.min(), xv.min()) * 0.8, max(dv.max(), xv.max()) * 1.2]
            a.plot(lim, lim, color="0.4", lw=1)
            a.scatter(xv, dv, s=18, alpha=0.75, color=colours[an], edgecolor="none")
            a.set_xscale("log"); a.set_yscale("log")
            a.set_xlim(lim); a.set_ylim(lim)
            a.set_xlabel(xl); a.set_ylabel("Observed (mg/L)")
            a.set_title(labels[an])
        a = ax[r, 2]
        a.axhline(0, color="0.4", lw=1)
        for h in (-2, 2):
            a.axhline(h, color="0.7", lw=0.8, ls="--")
        a.scatter(tm, cw, s=18, alpha=0.75, color=colours[an], edgecolor="none")
        a.set_xlabel("Time after dose (h)")
        a.set_ylabel("Conditional weighted residual")
        a.set_title(f"{labels[an]}  (SD {cw.std(ddof=1):.2f})")
    fig.suptitle("Model 1 goodness of fit, CRRT cohort, 21 patients", y=0.99)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "model1_gof.png"), dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for i, an in enumerate(M.ANALYTES):
        v = [x for x in vpc_rows if x["analyte"] == an]
        t = np.array([x["time_h"] for x in v])
        a = ax[i]
        a.fill_between(t, [x["sim_p5"] for x in v], [x["sim_p95"] for x in v],
                       color=colours[an], alpha=0.18, label="simulated 5th-95th percentile")
        a.plot(t, [x["sim_p50"] for x in v], color=colours[an], lw=2, label="simulated median")
        a.plot(t, [x["obs_p50"] for x in v], "o--", color="0.15", ms=5, lw=1.2, label="observed median")
        a.plot(t, [x["obs_p5"] for x in v], ".", color="0.45", ms=6)
        a.plot(t, [x["obs_p95"] for x in v], ".", color="0.45", ms=6, label="observed 5th and 95th")
        a.set_xlabel("Time after dose (h)")
        a.set_ylabel("Concentration (mg/L)")
        a.set_title(labels[an])
        a.legend(fontsize=8, frameon=False)
    fig.suptitle("Model 1 visual predictive check, 1000 replicates of the observed design", y=1.0)
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
        w.writerows({k: (f"{v:.6f}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()}
                    for r in rows)
    print(f"  wrote {path} ({len(rows)} rows)")


# ------------------------------------------------------------ profile likelihood

def profile(full):
    """Coarse grid, then a fine grid inside each bracket that contains a crossing of 3.84."""
    coarse = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.94, 0.96]
    r_hat = float(np.tanh(full.x[IDX_RHO]))
    jobs = [(f"rho={r:.3f}", "4dev", full.x, {IDX_RHO: float(np.arctanh(r))}, None, None, 8) for r in coarse]
    pts = {r: res["ofv"] - full.fun for r, res in zip(coarse, run_jobs(jobs))}
    pts[r_hat] = 0.0

    def bracket(lo_side):
        side = sorted(r for r in pts if (r <= r_hat if lo_side else r >= r_hat))
        pairs = list(zip(side[:-1], side[1:]))
        for a, b in (reversed(pairs) if lo_side else pairs):
            if (pts[a] - CRIT) * (pts[b] - CRIT) <= 0:
                return a, b
        return None

    fine = []
    for br in (bracket(True), bracket(False)):
        if br:
            fine += [float(v) for v in np.linspace(br[0], br[1], 8)[1:-1]]
    if fine:
        jobs = [(f"rho={r:.4f}", "4dev", full.x, {IDX_RHO: float(np.arctanh(r))}, None, None, 8) for r in fine]
        for r, res in zip(fine, run_jobs(jobs)):
            pts[r] = res["ofv"] - full.fun

    rr = np.array(sorted(pts))
    dd = np.array([pts[r] for r in rr])

    def crossing(lo_side):
        if lo_side:
            x_, y_ = rr[rr <= r_hat][::-1], dd[rr <= r_hat][::-1]       # walk down from the estimate
        else:
            x_, y_ = rr[rr >= r_hat], dd[rr >= r_hat]                    # walk up from the estimate
        for k in range(len(x_) - 1):
            if (y_[k] - CRIT) * (y_[k + 1] - CRIT) <= 0 and y_[k + 1] != y_[k]:
                return float(x_[k] + (CRIT - y_[k]) * (x_[k + 1] - x_[k]) / (y_[k + 1] - y_[k]))
        return float("nan")

    lo, hi = crossing(True), crossing(False)
    rows = [{"rho": float(r), "ofv": float(full.fun + d), "delta_ofv": float(d),
             "in_95_interval": "yes" if d <= CRIT else "no"} for r, d in zip(rr, dd)]
    return lo, hi, rows, pts


# ---------------------------------------------------------------------- main -

def main():
    M.structural_self_check()
    subjects = M.load()
    nobs = sum(len(M.obs_vector(s)) for s in subjects)
    print("=" * 78)
    print("MODEL 1 FINALISATION")
    print("=" * 78)
    print(f"  {len(subjects)} subjects, {nobs} observations; {WORKERS} worker processes\n")

    # 1 -------------------------------------------------------- convergence --
    print("1. Convergence")
    full, cache = M.fit(subjects, M.build_omega4, N_OMEGA, "full model", P0)
    theta, om, w, r_cl, r_v, sigma, c = unpack(full.x)
    alt = run_jobs([
        ("clearance correlation fixed at 0", "4dev", full.x, {IDX_RHO: 0.0}, None, None, 12),
        ("residual correlation fixed at 0", "4dev", full.x, {IDX_C: 0.0}, None, None, 12),
        ("shared volume deviate", "3dev", np.delete(full.x, IDX_RV), None, None, None, 12),
        ("3 deviates, no residual correlation", "3dev", np.delete(full.x, [IDX_RV, IDX_C]), None, None, None, 12),
    ])
    d_rho0 = alt[0]["ofv"] - full.fun
    d_c0 = alt[1]["ofv"] - full.fun
    d_shared = alt[2]["ofv"] - full.fun
    print(f"    OFV {full.fun:.4f}; rho {r_cl:.4f}; r_V {r_v:.4f}; residual correlation {c:.4f}")
    print(f"    dOFV with rho = 0: {d_rho0:.3f}; with c = 0: {d_c0:.3f}; with a shared volume deviate: {d_shared:.3f}")
    print(f"    3 deviates, no residual correlation: OFV {alt[3]['ofv']:.4f}, rho {alt[3]['rho']:.4f}")

    # 2 ------------------------------------------------- profile likelihood --
    print("\n2. Profile likelihood for the clearance correlation")
    pl_lo, pl_hi, prof, pts = profile(full)
    for p in prof:
        print(f"    rho {p['rho']:7.4f}   dOFV {p['delta_ofv']:8.4f}")
    d_094 = pts[0.94]
    print(f"\n    profile-likelihood 95% interval: {pl_lo:.4f} to {pl_hi:.4f}; dOFV at 0.94 {d_094:.3f}")
    print(f"    excludes 0.94: {'YES' if pl_hi < 0.94 else 'NO'}")
    write_csv(prof, "model1_profile_likelihood.csv")

    # 3 -------------------------------------------------------- diagnostics --
    print("\n3. Goodness-of-fit diagnostics")
    diag = diagnostics(subjects, full.x, cache)
    for an in M.ANALYTES:
        cw = np.array([d["cwres"] for d in diag if d["analyte"] == an])
        print(f"    {an}: CWRES mean {cw.mean():+.3f}, SD {cw.std(ddof=1):.3f}, "
              f"{100.0 * float(np.mean(np.abs(cw) > 2)):.1f}% beyond +/-2")
    write_csv(diag, "model1_diagnostics.csv")

    # 4 --------------------------------------------- visual predictive check --
    print(f"\n4. Visual predictive check (1000 replicates, seed {VPC_SEED})")
    vrows = vpc(subjects, full.x)
    n_in = sum(r["n_obs_within_sim_90"] for r in vrows)
    n_all = sum(r["n_obs"] for r in vrows)
    cov = 100.0 * n_in / n_all                      # per observation: the statistic the text reports
    cov_bins = float(np.mean([r["pct_obs_within_sim_90"] for r in vrows]))
    cov_an = {an: 100.0 * sum(r["n_obs_within_sim_90"] for r in vrows if r["analyte"] == an)
              / sum(r["n_obs"] for r in vrows if r["analyte"] == an) for an in M.ANALYTES}
    print(f"    observations inside the simulated 5th-95th interval: {n_in} of {n_all} = {cov:.1f}% "
          f"(nominal 90%); ceftazidime {cov_an['caz']:.1f}%, avibactam {cov_an['avi']:.1f}%; "
          f"unweighted mean over time bins {cov_bins:.1f}%")
    write_csv(vrows, "model1_vpc.csv")

    # 5 --------------------------------------------------------- sensitivity --
    print("\n5. Sensitivity analyses")
    mean_sig = 0.5 * (full.x[10] + full.x[11])
    x_sh = full.x.copy()
    x_sh[10] = x_sh[11] = mean_sig
    jobs = [(f"all {t:.0f} h", "4dev", full.x, None, t, None, 8) for t in (1.0, 2.0, 3.0)]
    jobs += [("shared residual SD", "4dev", x_sh, {11: mean_sig}, None, None, 8)]
    jobs += [(f"without patient {s.sid}", "4dev", full.x, None, None, s.sid, 6) for s in subjects]
    res = run_jobs(jobs)
    sens = [{"analysis": "reference", "variant": "as fitted", "corr_CL": round(r_cl, 6), "ofv": round(full.fun, 6),
             "note": "4 deviates; residual correlation between the 2 drugs in the same sample"}]
    for r, t in zip(res[:3], (1.0, 2.0, 3.0)):
        sens.append({"analysis": "infusion duration", "variant": f"all {t:.0f} h", "corr_CL": round(r["rho"], 6),
                     "ofv": round(r["ofv"], 6), "note": "tests the categorical-versus-hours ambiguity in the source"})
    sens.append({"analysis": "variance structure", "variant": "shared volume deviate",
                 "corr_CL": round(alt[2]["rho"], 6), "ofv": round(alt[2]["ofv"], 6),
                 "note": f"dOFV {d_shared:+.2f} on 1 df against the reference"})
    sens.append({"analysis": "variance structure", "variant": "3 deviates, no residual correlation",
                 "corr_CL": round(alt[3]["rho"], 6), "ofv": round(alt[3]["ofv"], 6),
                 "note": f"dOFV {alt[3]['ofv'] - full.fun:+.2f} on 2 df against the reference"})
    sens.append({"analysis": "residual error", "variant": "SD shared across drugs",
                 "corr_CL": round(res[3]["rho"], 6), "ofv": round(res[3]["ofv"], 6),
                 "note": f"dOFV {res[3]['ofv'] - full.fun:+.2f} on 1 df"})
    loo = np.array([r["rho"] for r in res[4:]])
    sens.append({"analysis": "leave-one-subject-out", "variant": "range over 21 refits",
                 "corr_CL": f"{loo.min():.6f} to {loo.max():.6f}", "ofv": "",
                 "note": f"median {np.median(loo):.6f}; most influential subject "
                         f"{subjects[int(np.argmax(np.abs(loo - r_cl)))].sid}"})
    for r in sens:
        print(f"    {r['analysis']:24} {r['variant']:40} {r['corr_CL']}")
    write_csv(sens, "model1_sensitivity.csv")

    # 6 ------------------------------------------------------------ figures --
    print("\n6. Figures")
    make_figures(diag, vrows)

    # 7 ------------------------------------------------- final parameter table
    print("\n7. Final parameter table")
    vals = list(theta) + list(w) + [r_cl, r_v] + list(sigma) + [c]
    rows = []
    for nm, v in zip(NAMES, vals):
        is_rho = nm == "corr_CL_caz_avi"
        rows.append({"parameter": nm, "estimate": f"{v:.6f}",
                     "ci_low": f"{pl_lo:.6f}" if is_rho else "", "ci_high": f"{pl_hi:.6f}" if is_rho else "",
                     "ci_method": "profile likelihood" if is_rho else "", "note": ""})
    extra = [("OFV", full.fun, f"converged on tolerance: {full.converged}, {full.rounds} rounds"),
             ("dOFV_vs_no_correlation", d_rho0,
              f"clearance correlation fixed at 0; 1 df, p = {chi2.sf(d_rho0, 1):.4g}"),
             ("dOFV_rho_0.94", d_094, f"clearance correlation fixed at 0.94; 1 df, p = {chi2.sf(d_094, 1):.4g}"),
             ("dOFV_vs_no_residual_correlation", d_c0,
              f"residual correlation fixed at 0; 1 df, p = {chi2.sf(d_c0, 1):.4g}"),
             ("dOFV_vs_shared_volume", d_shared, f"shared volume deviate; 1 df, p = {chi2.sf(d_shared, 1):.4g}"),
             ("vpc_coverage_pct", cov, f"{n_in} of {n_all} observations inside the simulated 90% interval; "
                                       f"nominal 90%, seed {VPC_SEED}"),
             ("vpc_coverage_caz_pct", cov_an["caz"], "ceftazidime observations"),
             ("vpc_coverage_avi_pct", cov_an["avi"], "avibactam observations"),
             ("vpc_coverage_mean_over_bins_pct", cov_bins, "unweighted mean over sampling-time bins")]
    for nm, v, note in extra:
        rows.append({"parameter": nm, "estimate": f"{v:.6f}", "ci_low": "", "ci_high": "", "ci_method": "",
                     "note": note})
    write_csv(rows, "model1_final_parameters.csv")

    print("\n" + "=" * 78)
    print(f"  clearance correlation {r_cl:.4f}, profile-likelihood 95% interval {pl_lo:.4f} to {pl_hi:.4f}")
    print(f"  0.94 is {'EXCLUDED' if pl_hi < 0.94 else 'NOT excluded'}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
