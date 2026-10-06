"""Covariance structure of the 4 random effects.

Structures (all with the residual correlation between the 2 drugs in the same sample):
  primary  rho (clearance-clearance) and r_V (volume-volume)                        13 parameters
  Ka       + clearance-volume correlation within each drug                          15 parameters
  Kb       all 6 correlations, parametrised by C-vine partial correlations          17 parameters
           (variable order CL ceftazidime, CL avibactam, V ceftazidime, V avibactam), so every
           parameter vector gives a positive-definite matrix and the first partial correlation is rho.
Each structure is fitted from several starting points with the Model 1 objective and optimizer
(joint_popk_nlme.fit); the best fit is kept. The profile likelihood of rho is then computed on a grid,
every other parameter re-estimated, and the 95% interval read where the change in OFV crosses 3.841.
Values are written at full precision (6 decimals) and rounded once, when they are reported.
Outputs: outputs/model1_covariance_structure.csv, outputs/model1_covariance_profile.csv
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
from scipy.stats import chi2

import joint_popk_nlme as M
import model1_sbc as S

WORKERS = int(os.environ.get("REFIT_WORKERS", "8"))
CRIT = float(chi2.ppf(0.95, 1))
GRID = (0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.94)
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
NAMES = {(0, 1): "rho_cl_caz_cl_avi", (2, 3): "r_v_caz_v_avi", (0, 2): "r_cl_caz_v_caz", (1, 3): "r_cl_avi_v_avi",
         (0, 3): "r_cl_caz_v_avi", (1, 2): "r_cl_avi_v_caz"}
IDX_RHO = 8
D = 6


# ------------------------------------------------------------------ structures ---------------------------------
def cvine(pc):
    """Correlation matrix from C-vine partial correlations pc[(i, j)] (partial on variables 0..i-1)."""
    R = np.eye(4)
    for i in range(4):
        for j in range(i + 1, 4):
            p = pc[(i, j)]
            for k in range(i - 1, -1, -1):
                p = p * np.sqrt((1 - pc[(k, i)] ** 2) * (1 - pc[(k, j)] ** 2)) + pc[(k, i)] * pc[(k, j)]
            R[i, j] = R[j, i] = p
    return R


def to_cvine(R):
    """Inverse of cvine(): partial correlation of (i, j) given variables 0..i-1."""
    pc = {}
    for i, j in PAIRS:
        s = list(range(i))
        a = [i, j]
        if s:
            cond = R[np.ix_(a, a)] - R[np.ix_(a, s)] @ np.linalg.solve(R[np.ix_(s, s)], R[np.ix_(s, a)])
        else:
            cond = R[np.ix_(a, a)]
        pc[(i, j)] = cond[0, 1] / np.sqrt(cond[0, 0] * cond[1, 1])
    return pc


def build_kb(p):
    """[log w1..w4, z01, z02, z03, z12, z13, z23] with z = atanh(partial correlation)."""
    w = np.exp(p[:4])
    pc = {pair: float(np.tanh(p[4 + k])) for k, pair in enumerate(PAIRS)}
    R = cvine(pc)
    return R, w, float(R[0, 1]), float(R[2, 3])


def build_ka(p):
    """[log w1..w4, z_cl, z_v, z_clv_caz, z_clv_avi]; non-positive-definite values are rejected by M.ofv."""
    w = np.exp(p[:4])
    R = np.eye(4)
    for k, (i, j) in enumerate(((0, 1), (2, 3), (0, 2), (1, 3))):
        R[i, j] = R[j, i] = np.tanh(p[4 + k])
    return R, w, float(R[0, 1]), float(R[2, 3])


STRUCT = {"primary": (M.build_omega4, 6), "Ka": (build_ka, 8), "Kb": (build_kb, 10)}


def start(structure, R):
    """Full parameter vector for a structure from the primary estimates and a target correlation matrix."""
    x = S.fitted()
    head, tail = x[:8], x[10:]
    if structure == "primary":
        return x
    if structure == "Ka":
        z = [np.arctanh(R[0, 1]), np.arctanh(R[2, 3]), np.arctanh(R[0, 2]), np.arctanh(R[1, 3])]
    else:
        pc = to_cvine(R)
        z = [np.arctanh(pc[pair]) for pair in PAIRS]
    return np.r_[head, z, tail]


def corr_of(structure, x):
    b, n_om = STRUCT[structure]
    R, w, _, _ = b(x[4:4 + n_om])
    return R, w


def fit_job(job):
    t0 = time.time()
    structure, label, x0, fixed_rho = job
    b, n_om = STRUCT[structure]
    x0 = np.asarray(x0, float)
    fixed = None
    if fixed_rho is not None:
        x0 = x0.copy()
        x0[IDX_RHO] = np.arctanh(fixed_rho)
        fixed = {IDX_RHO: float(np.arctanh(fixed_rho))}
    fr, _ = M.fit(M.load(), b, n_om, "", x0, quiet=True, fixed=fixed)
    return dict(structure=structure, start=label, fixed_rho=fixed_rho, x=fr.x.tolist(), ofv=fr.fun,
                rounds=fr.rounds, converged=fr.converged, seconds=time.time() - t0)


def run(jobs):
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        return list(ex.map(fit_job, jobs))


def write(name, rows):
    with open(os.path.join(M.OUT, name), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.{D}f}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    print(f"  wrote outputs/{name} ({len(rows)} rows)", flush=True)


def crossing(grid_rows, best_rho, side):
    """Linear interpolation of the 3.841 crossing of the change in OFV on one side of the estimate."""
    pts = sorted((r["rho"], r["delta_ofv"]) for r in grid_rows)
    pts = [p for p in pts if (p[0] <= best_rho if side == "low" else p[0] >= best_rho)]
    pts = pts if side == "high" else pts[::-1]
    prev = (best_rho, 0.0)
    for rho, d in pts:
        if d >= CRIT:
            (r0, d0), (r1, d1) = prev, (rho, d)
            return r0 + (CRIT - d0) * (r1 - r0) / (d1 - d0)
        prev = (rho, d)
    return float("nan")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    R_primary = np.eye(4)
    _, om0, _, r_cl0, r_v0, _, _ = M.unpack(S.fitted())
    R_primary[0, 1] = R_primary[1, 0] = r_cl0
    R_primary[2, 3] = R_primary[3, 2] = r_v0
    R_ka = R_primary.copy()                        # moderate within-drug terms; Ka has zero cross terms, so the
    for (i, j), v in {(0, 2): 0.25, (1, 3): 0.25}.items():   # Bayesian medians (0.87, 0.78) are not admissible
        R_ka[i, j] = R_ka[j, i] = v
    R_cross = R_primary.copy()                     # Bayesian medians of the full fit (model1_bayes_numpyro.py), cross terms implied
    for (i, j), v in {(0, 1): 0.45, (2, 3): 0.45, (0, 2): 0.87, (1, 3): 0.78}.items():
        R_cross[i, j] = R_cross[j, i] = v
    R_cross[0, 3] = R_cross[3, 0] = 0.87 * 0.45
    R_cross[1, 2] = R_cross[2, 1] = 0.78 * 0.45
    for R in (R_ka, R_cross):
        assert np.linalg.eigvalsh(R).min() > 0, np.linalg.eigvalsh(R)
    pc = {pair: np.tanh(v) for pair, v in zip(PAIRS, np.random.default_rng(1).normal(size=6))}
    assert np.linalg.eigvalsh(cvine(pc)).min() > 0
    assert np.allclose(cvine(to_cvine(R_cross)), R_cross)
    t0 = time.time()
    stage1 = [("Ka", "primary values", start("Ka", R_primary).tolist(), None),
              ("Ka", "within-drug 0.25", start("Ka", R_ka).tolist(), None),
              ("Kb", "primary values", start("Kb", R_primary).tolist(), None),
              ("Kb", "Bayesian medians, implied cross terms", start("Kb", R_cross).tolist(), None)]
    fits = run(stage1)
    print(f"stage 1 done in {(time.time() - t0) / 60:.1f} min", flush=True)
    best = {}
    for f in fits:
        print(f"  {f['structure']} from {f['start']}: OFV {f['ofv']:.4f}, rho {corr_of(f['structure'], np.array(f['x']))[0][0, 1]:.4f}",
              flush=True)
        if f["structure"] not in best or f["ofv"] < best[f["structure"]]["ofv"]:
            best[f["structure"]] = f
    # Kb nests Ka: if a Kb start did worse than the best Ka fit, restart Kb from the Ka solution
    if best["Kb"]["ofv"] > best["Ka"]["ofv"] + 1e-3:
        Rka, _ = corr_of("Ka", np.array(best["Ka"]["x"]))
        xa = np.array(best["Ka"]["x"])
        pcs = to_cvine(Rka)
        x0 = np.r_[xa[:8], [np.arctanh(pcs[p]) for p in PAIRS], xa[12:]]
        f = run([("Kb", "Ka solution", x0.tolist(), None)])[0]
        fits.append(f)
        print(f"  Kb from the Ka solution: OFV {f['ofv']:.4f}", flush=True)
        if f["ofv"] < best["Kb"]["ofv"]:
            best["Kb"] = f
    ref = M.ofv(S.fitted(), M.load(), M.build_omega4, M.N_OMEGA, {})
    jobs = [(s, "profile", best[s]["x"], g) for s in ("Ka", "Kb") for g in GRID]
    prof = run(jobs)
    print(f"stage 2 done in {(time.time() - t0) / 60:.1f} min", flush=True)
    prow, srow = [], []
    for s in ("Ka", "Kb"):
        b_ofv = min([best[s]["ofv"]] + [p["ofv"] for p in prof if p["structure"] == s])
        if b_ofv < best[s]["ofv"] - 1e-3:
            print(f"  WARNING {s}: a profile point is below the free fit by {best[s]['ofv'] - b_ofv:.4f}", flush=True)
        rho_hat = corr_of(s, np.array(best[s]["x"]))[0][0, 1]
        g_rows = []
        for p in prof:
            if p["structure"] == s:
                g_rows.append(dict(structure=s, rho=p["fixed_rho"], ofv=p["ofv"], delta_ofv=p["ofv"] - b_ofv))
        prow += sorted(g_rows, key=lambda r: r["rho"])
        R, w = corr_of(s, np.array(best[s]["x"]))
        x = np.array(best[s]["x"])
        n_om = STRUCT[s][1]
        sig = np.exp(x[4 + n_om:6 + n_om])
        c = float(np.tanh(x[6 + n_om]))
        d094 = next(r["delta_ofv"] for r in g_rows if abs(r["rho"] - 0.94) < 1e-9)
        d0 = next(r["delta_ofv"] for r in g_rows if r["rho"] == 0.0)
        row = dict(structure=s, n_par=len(x), ofv=best[s]["ofv"], dofv_vs_primary=ref - best[s]["ofv"],
                   df_vs_primary=len(x) - 13, p_vs_primary=float(chi2.sf(max(ref - best[s]["ofv"], 0), len(x) - 13)),
                   aic=best[s]["ofv"] + 2 * len(x), rho=rho_hat, rho_ci_low=crossing(g_rows, rho_hat, "low"),
                   rho_ci_high=crossing(g_rows, rho_hat, "high"), dofv_rho_0=d0, dofv_rho_094=d094,
                   cl_caz=float(np.exp(x[0])), cl_avi=float(np.exp(x[1])), v_caz=float(np.exp(x[2])),
                   v_avi=float(np.exp(x[3])), omega_cl_caz=w[0], omega_cl_avi=w[1], omega_v_caz=w[2], omega_v_avi=w[3],
                   sigma_caz=sig[0], sigma_avi=sig[1], c=c, min_eigenvalue=float(np.linalg.eigvalsh(R).min()),
                   start=best[s]["start"], rounds=best[s]["rounds"], converged=best[s]["converged"])
        for pair, nm in NAMES.items():
            row[nm] = float(R[pair])
        srow.append(row)
    srow.insert(0, dict({k: "" for k in srow[0]}, structure="primary", n_par=13, ofv=ref, dofv_vs_primary=0.0,
                        df_vs_primary=0, aic=ref + 26, rho=r_cl0, r_v_caz_v_avi=r_v0))
    write("model1_covariance_structure.csv", srow)
    write("model1_covariance_profile.csv", prow)
    for r in srow:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if v != ""}, flush=True)
    print(f"total {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    raise SystemExit(main())
