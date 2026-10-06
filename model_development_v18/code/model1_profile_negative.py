"""Profile likelihood of rho below 0.

rho is estimated as tanh(z), so it is free between -1 and 1; the stored profile (model1_finalise.py) covers 0.05 to
0.96 and the fit with rho fixed at 0. Here rho is fixed at -0.25, -0.5 and -0.75, every other parameter re-estimated
with the Model 1 fit (joint_popk_nlme.fit) from the primary estimates.

Values are written at full precision (6 decimals) and rounded once, when they are reported.
Output: outputs/model1_profile_negative.csv
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

import joint_popk_nlme as M
import model1_sbc as S

sys.stdout.reconfigure(encoding="utf-8")
IDX_RHO = 8
GRID = (-0.25, -0.5, -0.75)
D = 6


def fit_at(rho):
    t0 = time.time()
    x0 = S.fitted().copy()
    x0[IDX_RHO] = np.arctanh(rho)
    fr, _ = M.fit(M.load(), M.build_omega4, M.N_OMEGA, "", x0, quiet=True, fixed={IDX_RHO: float(np.arctanh(rho))})
    return dict(rho=rho, ofv=fr.fun, rounds=fr.rounds, converged=fr.converged, seconds=time.time() - t0)


def main():
    t0 = time.time()
    with open(os.path.join(M.OUT, "model1_final_parameters.csv"), encoding="utf-8") as fh:
        best = float(next(r for r in csv.DictReader(fh) if r["parameter"] == "OFV")["estimate"])
    with ProcessPoolExecutor(max_workers=int(os.environ.get("REFIT_WORKERS", "3"))) as ex:
        rows = list(ex.map(fit_at, GRID))
    for r in rows:
        r["delta_ofv"] = r["ofv"] - best
    path = os.path.join(M.OUT, "model1_profile_negative.csv")
    keys = ["rho", "ofv", "delta_ofv", "rounds", "converged", "seconds"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, lineterminator="\n")
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{r[k]:.{D}f}" if isinstance(r[k], float) else r[k]) for k in keys})
    for r in rows:
        print(f"  rho {r['rho']:+.2f}: delta OFV {r['delta_ofv']:.4f} ({r['rounds']} rounds, {r['seconds']:.0f} s)")
    print(f"wrote outputs/model1_profile_negative.csv; total {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
