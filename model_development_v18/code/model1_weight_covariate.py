"""Deposited weight category as a covariate on clearance and volume.

The structure with all 6 random-effect correlations (Kb, model1_covariance_structure.py) fits better than the primary
structure, but its clearance-volume correlations (0.93 and 0.88) leave the correlation matrix close to singular. Body
size moves clearance and volume together, so the deposited weight category w (0 to 3) may carry that shared variation:
    CL_d = theta_CL,d exp(b_CL (w - 1.5)),   V_d = theta_V,d exp(b_V (w - 1.5)),   b shared by the 2 drugs.
Fits: weight + primary structure (15 parameters) and weight + Kb structure (19 parameters, several starting points, the
best kept). Same objective (joint_popk_nlme.laplace_subject) and optimizer (alternating Nelder-Mead and L-BFGS-B) as
Model 1. The change in OFV is against the same structure without weight (2 df).

Values are written at full precision (6 decimals) and rounded once, when they are reported.
Output: outputs/model1_weight_covariate.csv
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
import model1_covariance_structure as CS
import model1_sbc as S
import model1_structural_refits as SR

sys.stdout.reconfigure(encoding="utf-8")
OUT = M.OUT
WORKERS = int(os.environ.get("REFIT_WORKERS", "4"))
D = 6


def weight_category():
    with open(os.path.join(M.DATA, "Demographic_data.csv"), encoding="utf-8-sig") as fh:
        return {int(r["subjectID"]): int(r["weight_cat"]) for r in csv.DictReader(fh)}


def objective(params, subjects, cache, structure, wcat):
    builder, n_om = CS.STRUCT[structure]
    theta = np.exp(params[:4])
    om, w, _, _ = builder(params[4:4 + n_om])
    sigma = np.exp(params[4 + n_om:6 + n_om])
    c = float(np.tanh(params[6 + n_om]))
    b_cl, b_v = params[7 + n_om:9 + n_om]
    if abs(c) >= 0.999 or abs(b_cl) > 3 or abs(b_v) > 3:
        return 1e10
    try:
        om_inv = np.linalg.inv(om)
        om_chol_inv = np.linalg.inv(np.linalg.cholesky(om))
    except np.linalg.LinAlgError:
        return 1e10
    total = 0.0
    for s in subjects:
        dw = wcat[s.sid] - 1.5
        th = theta * np.exp(np.array([b_cl, b_cl, b_v, b_v]) * dw)
        val, z = M.laplace_subject(s, th, w, om, om_inv, om_chol_inv, sigma, cache.get(s.sid, np.zeros(4)), c)
        cache[s.sid] = z
        total += val
    return total if np.isfinite(total) else 1e10


def kb_start(label):
    """Kb starting vectors: from the reported Kb estimates, or from the primary estimates with Kb's correlations."""
    with open(os.path.join(OUT, "model1_covariance_structure.csv"), encoding="utf-8") as fh:
        kb = next(r for r in csv.DictReader(fh) if r["structure"] == "Kb")
    R = np.eye(4)
    for (i, j), name in CS.NAMES.items():
        R[i, j] = R[j, i] = float(kb[name])
    if label == "Kb estimates":
        pc = CS.to_cvine(R)
        z = [np.arctanh(np.clip(pc[p], -0.999, 0.999)) for p in CS.PAIRS]
        head = np.log([float(kb[k]) for k in ("cl_caz", "cl_avi", "v_caz", "v_avi")] +
                      [float(kb[k]) for k in ("omega_cl_caz", "omega_cl_avi", "omega_v_caz", "omega_v_avi")])
        tail = np.r_[np.log([float(kb["sigma_caz"]), float(kb["sigma_avi"])]), np.arctanh(float(kb["c"]))]
        return np.r_[head, z, tail]
    if label == "primary estimates, Kb correlations":
        return CS.start("Kb", R)
    x = S.fitted()
    return CS.start("Kb", np.array(M.build_omega4(x[4:10])[0]))      # primary correlations, other pairs 0


def job(spec):
    t0 = time.time()
    structure, label = spec
    subjects, wcat = M.load(), weight_category()
    x0 = S.fitted() if structure == "primary" else kb_start(label)
    x, f, rounds, conv = SR.fit_obj(lambda p, cache: objective(p, subjects, cache, structure, wcat), np.r_[x0, 0.0, 0.0])
    return dict(structure=structure, start=label, x=x.tolist(), ofv=f, rounds=rounds, converged=conv,
                seconds=time.time() - t0)


def main():
    t0 = time.time()
    wcat = weight_category()
    counts = {k: sum(1 for v in wcat.values() if v == k) for k in sorted(set(wcat.values()))}
    print(f"weight categories: {counts}")
    specs = [("primary", "primary estimates"), ("Kb", "Kb estimates"), ("Kb", "primary estimates, Kb correlations"),
             ("Kb", "primary estimates, primary correlations")]
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        res = list(ex.map(job, specs))
    base = {}
    with open(os.path.join(OUT, "model1_covariance_structure.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            base[r["structure"]] = float(r["ofv"])
    with open(os.path.join(OUT, "model1_final_parameters.csv"), encoding="utf-8") as fh:
        base["primary"] = float(next(r for r in csv.DictReader(fh) if r["parameter"] == "OFV")["estimate"])
    rows = []
    for r in res:
        builder, n_om = CS.STRUCT[r["structure"]]
        x = np.array(r["x"])
        R, w, rho, r_v = builder(x[4:4 + n_om])
        eig = float(np.linalg.eigvalsh(R).min())
        d_ofv = base[r["structure"]] - r["ofv"]
        row = dict(structure=f"{r['structure']} + weight category", start=r["start"], n_par=len(x), ofv=r["ofv"],
                   ofv_without_weight=base[r["structure"]], dofv_vs_without_weight=d_ofv, df=2,
                   p_value=float(chi2.sf(d_ofv, 2)) if d_ofv > 0 else 1.0,
                   b_cl=float(x[7 + n_om]), b_v=float(x[8 + n_om]),
                   cl_per_category_pct=100 * float(np.expm1(x[7 + n_om])), v_per_category_pct=100 * float(np.expm1(x[8 + n_om])),
                   rho=rho, r_v=r_v, r_cl_caz_v_caz=float(R[0, 2]), r_cl_avi_v_avi=float(R[1, 3]),
                   r_cl_caz_v_avi=float(R[0, 3]), r_cl_avi_v_caz=float(R[1, 2]), min_eigenvalue=eig,
                   omega_cl_caz=float(w[0]), omega_cl_avi=float(w[1]), omega_v_caz=float(w[2]), omega_v_avi=float(w[3]),
                   c=float(np.tanh(x[6 + n_om])), rounds=r["rounds"], converged=r["converged"], seconds=r["seconds"])
        rows.append(row)
        print(f"  {row['structure']:28s} [{r['start']}]: OFV {r['ofv']:.4f} (dOFV {d_ofv:+.2f}); b_CL {row['b_cl']:+.3f}, "
              f"b_V {row['b_v']:+.3f}; rho {rho:.3f}; CL-V {row['r_cl_caz_v_caz']:.3f}/{row['r_cl_avi_v_avi']:.3f}; "
              f"min eig {eig:.2e}")
    path = os.path.join(OUT, "model1_weight_covariate.csv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{v:.{D}f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote outputs/model1_weight_covariate.csv; total {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
