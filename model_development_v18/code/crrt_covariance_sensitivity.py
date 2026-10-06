"""Monitoring results under the full covariance structure of the random effects.

The primary simulations (crrt_virtual_tdm.py) use the primary model, in which only the 2 clearances and the 2
volumes are correlated. Estimating all 6 correlations (model1_covariance_structure.py) improved the fit, and a
trough depends on clearance and volume together, so the main monitoring results are recomputed here with the
full model: typical values, between-patient SDs, residual SD and the full 4 x 4 correlation matrix.
The counterfactual rho = 0.94 replaces the first C-vine partial correlation (clearance-clearance) and keeps the
other partial correlations, so the matrix stays positive definite.

For each regimen and correlation: percentage below the 4-mg/L free avibactam trough target, registrational
attainment (free avibactam above 1 mg/L for 50% of the interval), percentage wrongly reassured by the 0.5 rule,
percentage flagged to identify half of the patients below target, area under the ROC curve of the measured
ceftazidime trough, and the 90% range of the free avibactam trough overall and given a measured ceftazidime
trough near the median. The primary model is run through the same function as a reference.
Random stream SEED + 6. Output: outputs/crrt_covariance_sensitivity.csv
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import crrt_virtual_tdm as sim  # noqa: E402
import model1_covariance_structure as K  # noqa: E402

OUT = sim.OUT
N, N_PRIOR = sim.N, 400_000
SEED = sim.SEED + 6
D = 6


def models():
    rr = {r["structure"]: r for r in csv.DictReader(open(os.path.join(OUT, "model1_covariance_structure.csv"), encoding="utf-8"))}
    R0 = np.eye(4)
    R0[0, 1] = R0[1, 0] = sim.RHO_EST
    R0[2, 3] = R0[3, 2] = sim.R_V
    out = {"primary": (dict(sim.TH), dict(sim.OM), sim.SIG_CAZ, R0)}
    r = rr["Kb"]
    th = {k: float(r[k]) for k in ("cl_caz", "cl_avi", "v_caz", "v_avi")}
    om = {k: float(r[f"omega_{k}"]) for k in ("cl_caz", "cl_avi", "v_caz", "v_avi")}
    R = np.eye(4)
    for pair, nm in K.NAMES.items():
        R[pair] = R[pair[::-1]] = float(r[nm])
    assert np.linalg.eigvalsh(R).min() > 0
    out["full"] = (th, om, float(r["sigma_caz"]), R)
    return out


def with_rho(R, rho):
    pc = K.to_cvine(R)
    pc[(0, 1)] = rho
    return K.cvine(pc)


def population(rng, n, R, th, om):
    z = rng.standard_normal((n, 4)) @ np.linalg.cholesky(R).T
    return dict(cl_caz=th["cl_caz"] * np.exp(om["cl_caz"] * z[:, 0]), cl_avi=th["cl_avi"] * np.exp(om["cl_avi"] * z[:, 1]),
                v_caz=th["v_caz"] * np.exp(om["v_caz"] * z[:, 2]), v_avi=th["v_avi"] * np.exp(om["v_avi"] * z[:, 3]))


def evaluate(rng, R, th, om, sig, reg):
    lab, cz, av, tau = reg
    p = population(rng, N, R, th, om)
    c_caz, c_avi, ft1 = sim.exposures(p, cz, av, tau)
    f_avi = sim.FU_AVI * c_avi
    below = f_avi < sim.AVI_CT
    obs = np.log(c_caz) + sig * rng.standard_normal(N)
    q = population(rng, N_PRIOR, R, th, om)
    qc, qa = sim.css(tau, q["cl_caz"], q["v_caz"], cz, tau), sim.css(tau, q["cl_avi"], q["v_avi"], av, tau)
    att = (sim.FU_AVI * qa >= sim.AVI_CT).astype(float)
    lq = np.log(qc)
    grid = np.linspace(obs.min() - 0.05, obs.max() + 0.05, 300)
    curve = np.empty(grid.size)
    for i, g in enumerate(grid):
        u = (g - lq) / sig
        w = np.exp(-0.5 * (u * u - np.min(u * u)))
        curve[i] = (w * att).sum() / w.sum()
    prob = np.interp(obs, grid, curve)
    n_below = int(below.sum())
    row = dict(regimen=lab, below_target_pct=100 * below.mean(), registrational_attainment_pct=100 * float(np.mean(ft1 >= 0.5)),
               avi_fcmin_p5=float(np.percentile(f_avi, 5)), avi_fcmin_p95=float(np.percentile(f_avi, 95)),
               wrongly_reassured_pct=100 * float(np.mean((prob >= 0.5) & below)), flagged_pct_for_half="", auroc_trough="")
    if n_below >= 20:
        order = np.argsort(prob, kind="stable")
        k50 = int(np.searchsorted(np.cumsum(below[order]), 0.5 * n_below))
        row["flagged_pct_for_half"] = 100 * (k50 + 1) / N
        ranks = np.argsort(np.argsort(np.r_[obs[below], obs[~below]])) + 1
        row["auroc_trough"] = 1 - (ranks[:n_below].sum() - n_below * (n_below + 1) / 2) / (n_below * (N - n_below))
    mid = (obs >= np.percentile(obs, 47.5)) & (obs <= np.percentile(obs, 52.5))
    row["avi_fcmin_p5_given_median_caz"], row["avi_fcmin_p95_given_median_caz"] = (float(v) for v in np.percentile(f_avi[mid], [5, 95]))
    return row


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    rows = []
    for name, (th, om, sig, R) in models().items():
        for reg in sim.REGIMENS + sim.EXTRA_REGIMENS:
            for scen, Rs in (("estimate", R), ("rho 0.94", with_rho(R, 0.94))):
                r = evaluate(np.random.default_rng([SEED, len(rows)]), Rs, th, om, sig, reg)
                rows.append(dict(model=name, scenario=scen, rho=float(Rs[0, 1]), **r))
                print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    path = os.path.join(OUT, "crrt_covariance_sensitivity.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.{D}f}" if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote", os.path.relpath(path, os.path.dirname(OUT)))


if __name__ == "__main__":
    main()
