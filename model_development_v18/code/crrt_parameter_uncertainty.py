"""Sensitivity analyses of the CRRT virtual population: parameter uncertainty and the avibactam
unbound fraction.

1. Parameter uncertainty. Parameter vectors are drawn from the asymptotic normal distribution of the
   Model 1 estimates on the estimation scale (log for typical values, between-patient SDs and residual
   SDs; Fisher z for correlations), using the covariance written by joint_popk_nlme.py
   (outputs/model1_parameter_covariance.csv). Each draw defines a virtual population of 20,000 patients
   per regimen, and its own ceftazidime-only classifier (0.5 rule, built from the same draw). Reported:
   median and 2.5th to 97.5th percentiles across draws of the percentage below the avibactam target and
   the percentage wrongly reassured. If the covariance matrix is not positive definite, nothing is
   simulated and the output says so.

2. Avibactam unbound fraction. O'Jeanson et al (Int J Antimicrob Agents 2025) measured a median avibactam
   unbound fraction of 0.73 (range 0.55-0.92) in 4 patients on CVVHDF, against 0.92 in the primary
   analysis. At 0.73: the percentage below target at each regimen, and classification at the 2 reduced
   regimens with the classifier recalibrated to 0.73 and with the classifier left at 0.92.

Random streams: draws from SEED + 3, one stream per draw from (SEED, 3, draw), and the unbound-fraction
block from SEED + 4, so no value in crrt_attainment.csv or crrt_classifier.csv changes.
Values are written at full precision (6 decimals) and rounded once, when they are reported.

Outputs: outputs/crrt_parameter_uncertainty_draws.csv, outputs/crrt_parameter_uncertainty.csv,
         outputs/crrt_unbound_fraction_sensitivity.csv
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import csv
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import crrt_virtual_tdm as sim

sys.stdout.reconfigure(encoding="utf-8")
OUT = sim.OUT
N_DRAWS = 500
N_PATIENTS = 20_000
N_PRIOR = 100_000
FU_AVI_CRRT = 0.73
TARGETS = [r for r in sim.REGIMENS + sim.EXTRA_REGIMENS if r[0] in ("2.5 g q8h", "1.25 g q8h")]
WORKERS = int(os.environ.get("SIM_WORKERS", max(1, min(14, (os.cpu_count() or 2) - 2))))


def load_covariance():
    with open(os.path.join(OUT, "model1_parameter_covariance.csv"), encoding="utf-8") as fh:
        rr = list(csv.DictReader(fh))
    names = [r["parameter"] for r in rr]
    x = np.array([float(r["estimate"]) for r in rr])
    cov = np.array([[float(r[n]) if r[n] != "" else np.nan for n in names] for r in rr])
    return names, x, cov


def unpack(x):
    th = dict(cl_caz=np.exp(x[0]), cl_avi=np.exp(x[1]), v_caz=np.exp(x[2]), v_avi=np.exp(x[3]))
    om = dict(cl_caz=np.exp(x[4]), cl_avi=np.exp(x[5]), v_caz=np.exp(x[6]), v_avi=np.exp(x[7]))
    return th, om, float(np.tanh(x[8])), float(np.tanh(x[9])), float(np.exp(x[10]))


def one_draw(job):
    k, x = job
    th, om, rho, r_v, sig = unpack(np.asarray(x))
    rng = np.random.default_rng([sim.SEED, 3, k])
    row = dict(draw=k, rho=rho, r_v=r_v, cl_caz=th["cl_caz"], cl_avi=th["cl_avi"])
    for lab, cz, av, tau in TARGETS:
        p = sim.population(rng, N_PATIENTS, rho, 1.0, th=th, om=om, r_v=r_v)
        cmin_caz, cmin_avi, _ = sim.exposures(p, cz, av, tau)
        truth = sim.FU_AVI * cmin_avi >= sim.AVI_CT
        obs = np.log(cmin_caz) + sig * rng.standard_normal(N_PATIENTS)
        grid = np.linspace(obs.min() - 0.05, obs.max() + 0.05, 200)
        curve = sim.posterior_curve(rng, rho, cz, av, tau, 1.0, grid, n_prior=N_PRIOR, sig_caz=sig,
                                    th=th, om=om, r_v=r_v)
        pred = np.interp(obs, grid, curve) >= 0.5
        key = lab.replace(" ", "_").replace(".", "p")
        row[f"below_target_pct_{key}"] = 100 * float(np.mean(~truth))
        row[f"wrongly_reassured_pct_{key}"] = 100 * float(np.mean(pred & ~truth))
    return row


def classify_fu(rng, rho, cz, av, tau, fu_true, fu_model):
    p = sim.population(rng, sim.N, rho)
    cmin_caz, cmin_avi, _ = sim.exposures(p, cz, av, tau)
    truth = fu_true * cmin_avi >= sim.AVI_CT
    obs = np.log(cmin_caz) + sim.SIG_CAZ * rng.standard_normal(sim.N)
    grid = np.linspace(obs.min() - 0.05, obs.max() + 0.05, 300)
    curve = sim.posterior_curve(rng, rho, cz, av, tau, 1.0, grid, fu_avi=fu_model)
    return sim.summary(truth, np.interp(obs, grid, curve) >= 0.5)


def write(path, rows_):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_[0]), lineterminator="\n")
        w.writeheader()
        w.writerows({k: (f"{float(v):.6f}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()}
                    for r in rows_)
    print(f"  wrote {path} ({len(rows_)} rows)")


def main():
    # 1. Parameter uncertainty ------------------------------------------------------------------
    names, x_hat, cov = load_covariance()
    summ = []
    try:
        if not np.all(np.isfinite(cov)):
            raise np.linalg.LinAlgError("covariance has missing values")
        L = np.linalg.cholesky(cov)
        rng = np.random.default_rng(sim.SEED + 3)
        draws = [x_hat + L @ rng.standard_normal(x_hat.size) for _ in range(N_DRAWS)]
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            res = list(ex.map(one_draw, [(k, d.tolist()) for k, d in enumerate(draws, 1)], chunksize=4))
        write(os.path.join(OUT, "crrt_parameter_uncertainty_draws.csv"), res)
        sd = np.sqrt(np.diag(cov))
        print("  standard errors on the estimation scale: " + ", ".join(f"{n} {s:.4f}" for n, s in zip(names, sd)))
        for key in [k for k in res[0] if k.endswith(("q8h",))]:
            v = np.array([r[key] for r in res])
            summ.append(dict(quantity=key, n_draws=len(v), fixed_parameter_note="see crrt_classifier.csv",
                             median=float(np.median(v)), p2_5=float(np.percentile(v, 2.5)),
                             p97_5=float(np.percentile(v, 97.5)), mean=float(v.mean()), status="done"))
        rho_d = np.array([r["rho"] for r in res])
        summ.append(dict(quantity="rho_of_draws", n_draws=len(rho_d), fixed_parameter_note="",
                         median=float(np.median(rho_d)), p2_5=float(np.percentile(rho_d, 2.5)),
                         p97_5=float(np.percentile(rho_d, 97.5)), mean=float(rho_d.mean()), status="done"))
    except np.linalg.LinAlgError as e:
        summ.append(dict(quantity="all", n_draws=0, fixed_parameter_note="", median="", p2_5="", p97_5="",
                         mean="", status=f"not done: covariance not positive definite ({e})"))
    write(os.path.join(OUT, "crrt_parameter_uncertainty.csv"), summ)
    for s in summ:
        print("   ", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()})

    # 2. Avibactam unbound fraction ---------------------------------------------------------------
    rng = np.random.default_rng(sim.SEED + 4)
    fu_rows = []
    for lab, cz, av, tau in sim.REGIMENS + sim.EXTRA_REGIMENS:
        p = sim.population(rng, sim.N, sim.RHO_EST)
        _, cmin_avi, ft1 = sim.exposures(p, cz, av, tau, fu_avi=FU_AVI_CRRT)
        fu_rows.append(dict(regimen=lab, analysis="attainment", fu_avi_true=FU_AVI_CRRT, fu_avi_classifier="",
                            below_target_pct=100 * float(np.mean(FU_AVI_CRRT * cmin_avi < sim.AVI_CT)),
                            registrational_attainment_pct=100 * float(np.mean(ft1 >= 0.5)),
                            false_reassurance=""))
    for lab, cz, av, tau in sim.REGIMENS + sim.EXTRA_REGIMENS:
        if lab not in ("1.25 g q8h", "2.5 g q12h"):
            continue
        for fu_model, kind in ((FU_AVI_CRRT, "classification, recalibrated"),
                               (sim.FU_AVI, "classification, not recalibrated")):
            r = classify_fu(rng, sim.RHO_EST, cz, av, tau, FU_AVI_CRRT, fu_model)
            fu_rows.append(dict(regimen=lab, analysis=kind, fu_avi_true=FU_AVI_CRRT, fu_avi_classifier=fu_model,
                                below_target_pct=r["missed_by_population_prior"], registrational_attainment_pct="",
                                false_reassurance=r["false_reassurance"]))
    write(os.path.join(OUT, "crrt_unbound_fraction_sensitivity.csv"), fu_rows)
    for r in fu_rows:
        print("   ", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
