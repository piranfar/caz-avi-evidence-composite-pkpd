"""Further refits of Model 1: CRRT modality, infusion duration and rounding of the deposited values.

Modality: CL_j = theta_j x exp(beta_j x CVVH) for both drugs, 2 extra parameters.
Infusion duration: duration = f x the deposited category (1, 2 or 3 h), f estimated; and a visual predictive
check stratified by category, for the primary model and the model with f estimated.
Rounding: every deposited value is a multiple of 0.5 mg/L. 20 data sets with independent uniform noise of
+/-0.25 mg/L added to each value, each refitted with the primary model.
The covariance structure of the random effects is examined in model1_covariance_structure.py.

All fits use the Model 1 objective (joint_popk_nlme.laplace_subject) and the optimizer of joint_popk_nlme.fit
(alternating Nelder-Mead and L-BFGS-B until a round improves the OFV by less than 1e-4).
Values are written at full precision (6 decimals) and rounded once, when they are reported.
Outputs: outputs/model1_structural_refits.csv, outputs/model1_rounding_sensitivity.csv,
         outputs/model1_vpc_by_infusion.csv
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import csv
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

import joint_popk_nlme as M
import model1_sbc as S

OUT = M.OUT
WORKERS = int(os.environ.get("REFIT_WORKERS", "4"))
ROUND_SEED = 20261005
N_ROUND = 20
VPC_SEED = 20260811
D = 6


def modality():
    with open(os.path.join(M.DATA, "CRRT_parameters_data.csv"), encoding="utf-8-sig") as fh:
        m = {int(r["subjectID"]): int(r["crrt_modality"]) for r in csv.DictReader(fh)}
    assert sorted(m.values()).count(1) == 6 and sorted(m.values()).count(0) == 15, m   # 6 CVVH, 15 CVVHD
    return m


# ------------------------------------------------------------------ objectives --------------------------------
def ofv_ext(params, subjects, cache, kind, cvvh=None):
    """Primary model (13 parameters) plus extras at params[13:]: 'modality' beta_caz, beta_avi; 'infusion' log f."""
    theta = np.exp(params[:4])
    om, w, _, _ = M.build_omega4(params[4:10])
    sigma = np.exp(params[10:12])
    c = float(np.tanh(params[12]))
    extra = params[13:]
    if abs(c) >= 0.999:
        return 1e10
    try:
        om_inv = np.linalg.inv(om)
        om_chol_inv = np.linalg.inv(np.linalg.cholesky(om))
    except np.linalg.LinAlgError:
        return 1e10
    if kind == "infusion":
        f = float(np.exp(extra[0]))
        if not 0.05 < f * 3.0 < M.TAU - 0.05:
            return 1e10
    total = 0.0
    for s in subjects:
        th, subj = theta, s
        if kind == "modality" and cvvh[s.sid]:
            th = theta.copy()
            th[0] *= np.exp(extra[0])
            th[1] *= np.exp(extra[1])
        if kind == "infusion":
            subj = M.Subject(s.sid, s.t_inf * f, s.times, s.logconc)
        val, z = M.laplace_subject(subj, th, w, om, om_inv, om_chol_inv, sigma, cache.get(s.sid, np.zeros(4)), c)
        cache[s.sid] = z
        total += val
    return total if np.isfinite(total) else 1e10


def fit_obj(objective, p0, max_rounds=12, tol=1e-4):
    """The optimizer of joint_popk_nlme.fit, for an arbitrary objective(x, cache)."""
    x = np.array(p0, float)
    prev, rounds, cache = np.inf, 0, {}
    for rounds in range(1, max_rounds + 1):
        budget = 2500 if rounds == 1 else 600
        r1 = minimize(objective, x, args=(cache,), method="Nelder-Mead",
                      options={"maxiter": budget, "maxfev": budget, "xatol": 1e-6, "fatol": 1e-6, "adaptive": True})
        r2 = minimize(objective, r1.x, args=(cache,), method="L-BFGS-B",
                      options={"maxiter": 120, "ftol": 1e-12, "gtol": 1e-9, "eps": 1e-5})
        x, cur = (r2.x, r2.fun) if r2.fun < r1.fun else (r1.x, r1.fun)
        if prev - cur < tol:
            break
        prev = cur
    final = objective(x, {})
    return x, final, rounds, (prev - final) < tol or rounds < max_rounds


# ------------------------------------------------------------------ jobs --------------------------------------
def job(spec):
    t0 = time.time()
    kind, arg = spec
    subjects = M.load()
    x_hat = S.fitted()
    if kind == "modality":
        cv = modality()
        x, f, rounds, conv = fit_obj(lambda p, cache: ofv_ext(p, subjects, cache, "modality", cv), np.r_[x_hat, 0.0, 0.0])
        _, _, _, r_cl, r_v, _, c = M.unpack(x[:13])
        extra = dict(cl_caz_cvvh_vs_cvvhd=float(np.exp(x[13])), cl_avi_cvvh_vs_cvvhd=float(np.exp(x[14])))
        return dict(model="modality on both clearances", n_par=15, df=2, ofv=f, rho=r_cl, r_v=r_v, c=c,
                    extra=extra, rounds=rounds, converged=conv, seconds=time.time() - t0)
    if kind == "infusion":
        x, f, rounds, conv = fit_obj(lambda p, cache: ofv_ext(p, subjects, cache, "infusion"), np.r_[x_hat, 0.0])
        _, _, _, r_cl, r_v, _, c = M.unpack(x[:13])
        return dict(model="infusion duration = f x category", n_par=14, df=1, ofv=f, rho=r_cl, r_v=r_v, c=c,
                    extra=dict(f=float(np.exp(x[13]))), rounds=rounds, converged=conv, seconds=time.time() - t0,
                    x=x.tolist())
    if kind == "rounding":
        rng = np.random.default_rng([ROUND_SEED, arg])
        data = []
        for s in subjects:
            lc = {}
            for a in M.ANALYTES:
                conc = np.exp(s.logconc[a])
                assert np.allclose(conc * 2, np.round(conc * 2), atol=1e-9)        # multiples of 0.5 mg/L
                lc[a] = np.log(conc + rng.uniform(-0.25, 0.25, conc.size))
            data.append(M.Subject(s.sid, s.t_inf, s.times, lc))
        fr, _ = M.fit(data, M.build_omega4, M.N_OMEGA, "", x_hat, quiet=True)
        _, _, _, r_cl, r_v, sigma, c = M.unpack(fr.x)
        return dict(dataset=arg, rho=r_cl, r_v=r_v, c=c, sigma_caz=float(sigma[0]), sigma_avi=float(sigma[1]),
                    ofv=fr.fun, rounds=fr.rounds, converged=fr.converged, seconds=time.time() - t0)
    raise ValueError(kind)


# ------------------------------------------------------------------ VPC by infusion category --------------------
def vpc_by_category(x13, f=1.0, label="primary", n_rep=1000, seed=VPC_SEED):
    subjects = M.load()
    sim_subj = [M.Subject(s.sid, s.t_inf * f, s.times, s.logconc) for s in subjects]
    theta, om, w, _, _, sigma, c = M.unpack(np.asarray(x13, float))
    chol = np.linalg.cholesky(om)
    rng = np.random.default_rng(seed)
    sims = {s.sid: [] for s in subjects}
    for _ in range(n_rep):
        for s in sim_subj:
            z = chol @ rng.standard_normal(4)
            lp = M.log_pred(s, theta, z, w)
            n = len(s.times["caz"])
            u1, u2 = rng.standard_normal(n), rng.standard_normal(n)
            sims[s.sid].append(np.r_[lp[:n] + sigma[0] * u1, lp[n:] + sigma[1] * (c * u1 + np.sqrt(1 - c * c) * u2)])
    rows = []
    for cat in sorted({s.t_inf for s in subjects}):
        grp = [s for s in subjects if s.t_inf == cat]
        for j, a in enumerate(M.ANALYTES):
            for t in sorted({float(t) for s in grp for t in s.times[a]}):
                obs, sim = [], []
                for s in grp:
                    n = len(s.times[a])
                    hit = np.where(np.isclose(s.times[a], t))[0]
                    for i in hit:
                        obs.append(np.exp(s.logconc[a][i]))
                        sim.extend(np.exp(np.array(sims[s.sid])[:, j * n + i]))
                if not obs:
                    continue
                o, sm = np.array(obs), np.array(sim)
                lo, md, hi = np.percentile(sm, [5, 50, 95])
                rows.append(dict(model=label, infusion_category_h=cat, analyte=a, time_h=t, n_obs=len(o),
                                 obs_p50=float(np.median(o)), sim_p5=float(lo), sim_p50=float(md), sim_p95=float(hi),
                                 n_obs_within_sim_90=int(np.sum((o >= lo) & (o <= hi))),
                                 pct_obs_within_sim_90=100.0 * float(np.mean((o >= lo) & (o <= hi)))))
    return rows


def write(name, rows):
    keys = list(rows[0])
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.{D}f}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    print(f"  wrote outputs/{name} ({len(rows)} rows)", flush=True)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ref_ofv = M.ofv(S.fitted(), M.load(), M.build_omega4, M.N_OMEGA, {})
    print(f"primary model OFV at the stored estimates {ref_ofv:.4f}", flush=True)
    specs = [("modality", None), ("infusion", None)] + [("rounding", k) for k in range(1, N_ROUND + 1)]
    t0 = time.time()
    refits, rounding = [], []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for spec, res in zip(specs, ex.map(job, specs)):
            print(f"  {spec[0]}{'' if spec[1] is None else ' ' + str(spec[1])} done after {(time.time() - t0) / 60:.1f} min: "
                  f"rho {res['rho']:.4f}, OFV {res['ofv']:.4f}", flush=True)
            (rounding if spec[0] == "rounding" else refits).append(res)
    rows = []
    for r in refits:
        d = ref_ofv - r["ofv"]
        rows.append(dict(model=r["model"], n_par=r["n_par"], df=r["df"], ofv=r["ofv"], dofv_vs_primary=d,
                         p_value=float(chi2.sf(max(d, 0.0), r["df"])), rho=r["rho"], r_v=r["r_v"], c=r["c"],
                         extra="; ".join(f"{k} {v:.6f}" for k, v in r["extra"].items()), rounds=r["rounds"],
                         converged=r["converged"], seconds=r["seconds"]))
    write("model1_structural_refits.csv", rows)
    write("model1_rounding_sensitivity.csv", rounding)
    rho = np.array([r["rho"] for r in rounding])
    cc = np.array([r["c"] for r in rounding])
    print(f"  rounding: rho mean {rho.mean():.4f}, SD {rho.std(ddof=1):.4f}, range {rho.min():.4f}-{rho.max():.4f}; "
          f"c mean {cc.mean():.4f}, SD {cc.std(ddof=1):.4f}", flush=True)
    lfit = next(r for r in refits if r["model"].startswith("infusion"))
    vrows = vpc_by_category(S.fitted(), 1.0, "primary") + \
        vpc_by_category(np.asarray(lfit["x"][:13]), lfit["extra"]["f"], "infusion duration estimated")
    write("model1_vpc_by_infusion.csv", vrows)
    print(f"total {(time.time() - t0) / 60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
