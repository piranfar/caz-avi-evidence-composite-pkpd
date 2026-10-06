"""Calibrated profile-likelihood interval for rho, and the profile of the infusion scaling factor.

1. Calibrated interval. In the 500 datasets simulated at the estimate (model1_sbc_coverage.py), the nominal 95%
   profile interval (cutoff 3.841) contained the true rho in 90.6%. The 95th percentile c* of the OFV change at the
   true rho in those datasets is used as the cutoff instead; the 2 crossings are found by root finding on the profile
   (rho fixed, every other parameter re-estimated with the Model 1 fit), bracketed by the stored profile grid.
2. Infusion scaling factor. The deposited infusion duration is a category (completed within 1, 2 or 3 hours); the
   duration was f x category with f estimated at 1.000 (model1_structural_refits.py). Here the OFV is computed with f
   fixed at 0.8, 0.9, 1.0, 1.1 and 1.2, every other parameter re-estimated, to show the shape of the profile.

Values are written at full precision (6 decimals) and rounded once, when they are reported.
Outputs: outputs/model1_calibrated_interval.csv, outputs/model1_infusion_profile.csv
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
from scipy.optimize import brentq

import joint_popk_nlme as M
import model1_sbc as S
import model1_structural_refits as SR

sys.stdout.reconfigure(encoding="utf-8")
OUT = M.OUT
IDX_RHO = 8
F_GRID = (0.8, 0.9, 1.0, 1.1, 1.2)
D = 6


def rows(name):
    with open(os.path.join(OUT, name), encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def ofv_min():
    return float(next(r for r in rows("model1_final_parameters.csv") if r["parameter"] == "OFV")["estimate"])


def profile_ofv(rho, x0):
    fr, _ = M.fit(M.load(), M.build_omega4, M.N_OMEGA, "", np.r_[x0[:IDX_RHO], np.arctanh(rho), x0[IDX_RHO + 1:]],
                  quiet=True, fixed={IDX_RHO: float(np.arctanh(rho))})
    return fr.fun, fr.x


def crossing(side):
    """Root of OFV(rho) - OFV_min - c* on one side, bracketed by the stored profile grid."""
    t0 = time.time()
    cstar = cutoff()[0]
    grid = sorted((float(r["rho"]), float(r["delta_ofv"])) for r in rows("model1_profile_likelihood.csv"))
    best = ofv_min()
    if side == "lower":
        pts = [g for g in grid if g[0] < 0.5884]
        a = max(g for g in pts if g[1] > cstar)
        b = min(g for g in pts if g[1] <= cstar and g[0] > a[0])
    else:
        pts = [g for g in grid if g[0] > 0.5884]
        a = max(g for g in pts if g[1] <= cstar)
        b = min(g for g in pts if g[1] > cstar)
    evals = []
    state = {"x": S.fitted()}

    def g(rho):
        f, x = profile_ofv(rho, state["x"])
        state["x"] = x
        evals.append((rho, f - best))
        return f - best - cstar

    root = brentq(g, a[0], b[0], xtol=5e-4)
    return dict(side=side, cutoff=cstar, rho=root, bracket_low=a[0], bracket_high=b[0], n_fits=len(evals),
                evaluations="; ".join(f"{r:.6f}:{d:.6f}" for r, d in evals), seconds=time.time() - t0)


def cutoff():
    v = np.array([float(r["lrt_rho_true"]) for r in rows("model1_sbc_coverage_replicates.csv")
                  if r["scenario"] == "correctly specified"])
    assert v.size == 500, v.size
    return float(np.percentile(v, 95)), float(np.mean(v > 3.841459)), v.size


def infusion(f):
    t0 = time.time()
    subjects = M.load()
    lf = float(np.log(f))
    x, fun, rounds, conv = SR.fit_obj(lambda p, cache: SR.ofv_ext(np.r_[p, lf], subjects, cache, "infusion"), S.fitted())
    _, _, _, r_cl, r_v, _, c = M.unpack(x)
    return dict(f=f, ofv=fun, rho=r_cl, r_v=r_v, c=c, rounds=rounds, converged=conv, seconds=time.time() - t0)


def task(spec):
    kind, arg = spec
    return crossing(arg) if kind == "crossing" else infusion(arg)


def write(name, recs):
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]), lineterminator="\n")
        w.writeheader()
        for r in recs:
            w.writerow({k: (f"{v:.{D}f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"  wrote outputs/{name} ({len(recs)} rows)")


def main():
    t0 = time.time()
    cstar, above, n = cutoff()
    print(f"cutoff: 95th percentile of the OFV change at the true rho = {cstar:.4f} ({n} datasets; "
          f"{100 * above:.1f}% above 3.841)")
    specs = [("crossing", "lower"), ("crossing", "upper")] + [("infusion", f) for f in F_GRID]
    with ProcessPoolExecutor(max_workers=int(os.environ.get("REFIT_WORKERS", "7"))) as ex:
        res = list(ex.map(task, specs))
    cr, inf = res[:2], res[2:]
    best = ofv_min()
    for r in inf:
        r["dofv_vs_primary"] = r["ofv"] - best
    write("model1_calibrated_interval.csv", cr)
    write("model1_infusion_profile.csv", inf)
    print(f"calibrated interval {cr[0]['rho']:.4f} to {cr[1]['rho']:.4f} (cutoff {cstar:.4f})")
    for r in inf:
        print(f"  f {r['f']:.1f}: OFV {r['ofv']:.4f} (change {r['dofv_vs_primary']:+.4f}), rho {r['rho']:.4f}")
    print(f"total {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
