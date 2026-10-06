"""Avibactam clearance variability left unexplained by ceftazidime clearance, with intervals.

Here: CV of omega_avi x sqrt(1 - rho^2) and rho^2 (the share of the variance of log avibactam clearance explained by
log ceftazidime clearance) over the 500 parameter draws of crrt_parameter_uncertainty.py: same asymptotic covariance
(outputs/model1_parameter_covariance.csv), same random stream (SEED + 3) and order, so the draws are identical;
2.5th to 97.5th percentiles.
Without renal replacement therapy (Cojutti et al): rho 0.94 with a relative standard error of 23.8% gives a Wald 95%
interval of 0.50 to 1.38, capped at 1; the unexplained CV and rho^2 are computed at the ends of that interval with the
reported avibactam clearance variability (CV 76.91%).

Values are written at full precision (6 decimals) and rounded once, when they are reported.
Output: outputs/crrt_unexplained_variability.csv
"""
import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import crrt_circuit_clearance as CC  # noqa: E402
import crrt_parameter_uncertainty as PU  # noqa: E402
import crrt_virtual_tdm as sim  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
D = 6


def cv_pct(w):
    return 100 * math.sqrt(math.exp(w * w) - 1)


def main():
    names, x_hat, cov = PU.load_covariance()
    _, om0, rho0, _, _ = PU.unpack(x_hat)
    fp = {r["parameter"]: float(r["estimate"]) for r in CC.rows(os.path.join(CC.OUT, "model1_final_parameters.csv"))}
    assert abs(om0["cl_avi"] - fp["omega_CL_avi"]) < 1e-5 and abs(rho0 - fp["corr_CL_caz_avi"]) < 1e-5, (om0, rho0)
    L = np.linalg.cholesky(cov)
    rng = np.random.default_rng(sim.SEED + 3)
    draws = [x_hat + L @ rng.standard_normal(x_hat.size) for _ in range(PU.N_DRAWS)]
    ref = CC.rows(os.path.join(CC.OUT, "crrt_parameter_uncertainty_draws.csv"))
    for d, r in zip(draws, ref):                                  # the same draws as crrt_parameter_uncertainty.py
        _, _, rho, r_v, _ = PU.unpack(d)
        assert abs(rho - float(r["rho"])) < 1e-6 and abs(r_v - float(r["r_v"])) < 1e-6, r["draw"]
    recs = []

    def add(population, quantity, est, lo, hi, note):
        recs.append(dict(population=population, quantity=quantity, estimate=float(est), ci_low=float(lo),
                         ci_high=float(hi), note=note))

    _, om_hat, rho_hat, _, _ = PU.unpack(x_hat)
    u = np.array([cv_pct(PU.unpack(d)[1]["cl_avi"] * math.sqrt(1 - PU.unpack(d)[2] ** 2)) for d in draws])
    tot = np.array([cv_pct(PU.unpack(d)[1]["cl_avi"]) for d in draws])
    r2 = np.array([PU.unpack(d)[2] ** 2 for d in draws])
    note = f"{len(draws)} draws from the asymptotic covariance (as crrt_parameter_uncertainty.py); 2.5th-97.5th percentiles"
    add("CRRT (Model 1)", "unexplained_cv_cl_avi_pct", cv_pct(om_hat["cl_avi"] * math.sqrt(1 - rho_hat ** 2)),
        *np.percentile(u, [2.5, 97.5]), note)
    add("CRRT (Model 1)", "cv_cl_avi_pct", cv_pct(om_hat["cl_avi"]), *np.percentile(tot, [2.5, 97.5]), note)
    add("CRRT (Model 1)", "variance_explained_rho2", rho_hat ** 2, *np.percentile(r2, [2.5, 97.5]), note)

    coj = {(r["parameter"], r["drug"]): r for r in CC.rows(CC.COJ)}
    rho_c = float(coj[("CL_CAZ_CL_AVI", "both")]["value"])
    rse = float(coj[("CL_CAZ_CL_AVI", "both")]["rse_pct"]) / 100
    cv_avi = float(coj[("Omega_CL_CV", "AVI")]["value"]) / 100
    w_c = math.sqrt(math.log(1 + cv_avi ** 2))
    lo_r, hi_r = rho_c - 1.96 * rse * rho_c, rho_c + 1.96 * rse * rho_c
    hi_cap = min(hi_r, 1.0)
    note_c = (f"Cojutti et al: rho {rho_c} (relative standard error {100 * rse:.1f}%), Wald 95% interval {lo_r:.6f} to "
              f"{hi_r:.6f}, capped at 1; avibactam clearance CV {100 * cv_avi:.2f}% held at its estimate")
    add("without RRT (Cojutti et al)", "unexplained_cv_cl_avi_pct", cv_pct(w_c * math.sqrt(1 - rho_c ** 2)),
        cv_pct(w_c * math.sqrt(1 - hi_cap ** 2)), cv_pct(w_c * math.sqrt(1 - lo_r ** 2)), note_c)
    add("without RRT (Cojutti et al)", "variance_explained_rho2", rho_c ** 2, lo_r ** 2, hi_cap ** 2, note_c)
    add("without RRT (Cojutti et al)", "rho", rho_c, lo_r, hi_cap, note_c)

    path = os.path.join(CC.OUT, "crrt_unexplained_variability.csv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]), lineterminator="\n")
        w.writeheader()
        for r in recs:
            w.writerow({k: (f"{v:.{D}f}" if isinstance(v, float) else v) for k, v in r.items()})
    for r in recs:
        print(f"  {r['population']:30s} {r['quantity']:28s} {r['estimate']:.4f} ({r['ci_low']:.4f} to {r['ci_high']:.4f})")
    print("  wrote outputs/crrt_unexplained_variability.csv")


if __name__ == "__main__":
    main()
