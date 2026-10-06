"""Primary model refitted without the patients whose infusion was completed within 1 hour.

The 3 patients with 1-hour infusions lay below the prediction at every sampling time in the visual predictive check by
infusion category (model1_structural_refits.py). The leave-one-patient-out refits remove them one at a time; here all 3
are removed together and the primary model is refitted from the primary estimates with the Model 1 objective and
optimizer (joint_popk_nlme.fit). The change in OFV against the 21-patient fit is not a test: the data differ.
Within the 18 patients, rho is then fixed at 0.94 and at 0 (every other parameter re-estimated from the subset
estimates), and the OFV increases over the subset fit are reported.

Values are written at full precision (6 decimals) and rounded once, when they are reported.
Output: outputs/model1_without_1h_infusions.csv
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
D = 6


def fixed_fit(args):
    """The 18 patients with rho fixed (index 8 of the parameter vector is atanh(rho))."""
    rho, x0 = args
    kept = [s for s in M.load() if s.t_inf != 1.0]
    x0 = np.array(x0, float)
    x0[8] = np.arctanh(rho)
    fr, _ = M.fit(kept, M.build_omega4, M.N_OMEGA, "", x0, quiet=True, fixed={8: float(np.arctanh(rho))})
    return rho, fr.fun


def main():
    t0 = time.time()
    subjects = M.load()
    kept = [s for s in subjects if s.t_inf != 1.0]
    dropped = sorted(s.sid for s in subjects if s.t_inf == 1.0)
    assert len(dropped) == 3 and len(kept) == 18, dropped
    fr, _ = M.fit(kept, M.build_omega4, M.N_OMEGA, "without 1-hour infusions", S.fitted(), quiet=True)
    theta, _, w, r_cl, r_v, sigma, c = M.unpack(fr.x)
    with ProcessPoolExecutor(max_workers=2) as ex:
        fixed = dict(ex.map(fixed_fit, [(0.94, fr.x.tolist()), (0.0, fr.x.tolist())]))
    row = dict(analysis="primary model without the patients with 1-hour infusions",
               patients_excluded=" ".join(str(s) for s in dropped), n_patients=len(kept),
               n_obs=sum(2 * len(s.times["caz"]) for s in kept), ofv=fr.fun, rho=r_cl, r_v=r_v, c=c,
               cl_caz=float(theta[0]), cl_avi=float(theta[1]), v_caz=float(theta[2]), v_avi=float(theta[3]),
               omega_cl_caz=float(w[0]), omega_cl_avi=float(w[1]), dofv_rho_094=fixed[0.94] - fr.fun,
               dofv_rho_0=fixed[0.0] - fr.fun, rounds=fr.rounds, converged=fr.converged, seconds=time.time() - t0)
    path = os.path.join(M.OUT, "model1_without_1h_infusions.csv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(row), lineterminator="\n")
        wr.writeheader()
        wr.writerow({k: (f"{v:.{D}f}" if isinstance(v, float) else v) for k, v in row.items()})
    print(f"without patients {row['patients_excluded']}: rho {r_cl:.4f}, c {c:.4f}, r_v {r_v:.4f}, OFV {fr.fun:.4f} "
          f"({row['n_patients']} patients, {row['n_obs']} observations; {row['seconds'] / 60:.1f} min)")
    print("wrote outputs/model1_without_1h_infusions.csv")


if __name__ == "__main__":
    main()
