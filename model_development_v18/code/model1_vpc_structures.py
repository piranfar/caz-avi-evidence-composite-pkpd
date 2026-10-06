"""VPC of the primary and the full-covariance model on the observed design: which predicts the observed spread?

Same simulation for both: 1,000 replicates of the 21 patients, random effects from each model's correlation
matrix, correlated residual errors as in Model 1. Reports the percentage of observations inside the simulated
90% interval (overall, by drug, and at the 2 trough samples, 0 h and 8 h), and the observed against the
simulated between-patient SD of the log trough. Seed 20260811 (the VPC seed of model1_finalise.py).
Output: outputs/model1_vpc_structures.csv
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import joint_popk_nlme as M  # noqa: E402
import model1_covariance_structure as K  # noqa: E402
import model1_sbc as S  # noqa: E402

N_REP = 1000
SEED = 20260811


def models():
    th, om, w, r_cl, r_v, sig, c = M.unpack(S.fitted())
    out = {"primary": (th, w, om, sig, c)}
    r = {x["structure"]: x for x in csv.DictReader(open(os.path.join(M.OUT, "model1_covariance_structure.csv"), encoding="utf-8"))}["Kb"]
    R = np.eye(4)
    for pair, nm in K.NAMES.items():
        R[pair] = R[pair[::-1]] = float(r[nm])
    out["full"] = (np.array([float(r[k]) for k in ("cl_caz", "cl_avi", "v_caz", "v_avi")]),
                   np.array([float(r[f"omega_{k}"]) for k in ("cl_caz", "cl_avi", "v_caz", "v_avi")]), R,
                   np.array([float(r["sigma_caz"]), float(r["sigma_avi"])]), float(r["c"]))
    return out


def simulate(subjects, th, w, R, sig, c, rng):
    L = np.linalg.cholesky(R + 1e-12 * np.eye(4))
    sims = {s.sid: [] for s in subjects}
    for _ in range(N_REP):
        for s in subjects:
            z = L @ rng.standard_normal(4)
            lp = M.log_pred(s, th, z, w)
            n = len(s.times["caz"])
            u1, u2 = rng.standard_normal(n), rng.standard_normal(n)
            sims[s.sid].append(np.r_[lp[:n] + sig[0] * u1, lp[n:] + sig[1] * (c * u1 + np.sqrt(1 - c * c) * u2)])
    return {k: np.array(v) for k, v in sims.items()}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    subjects = M.load()
    rows = []
    for name, (th, w, R, sig, c) in models().items():
        sims = simulate(subjects, th, w, R, sig, c, np.random.default_rng(SEED))
        inside = {"caz": [], "avi": []}
        trough = {"caz": [], "avi": []}
        sd_sim = {"caz": [], "avi": []}
        for j, a in enumerate(M.ANALYTES):
            times = sorted({float(t) for s in subjects for t in s.times[a]})
            for t in times:
                obs, sim = [], []
                for s in subjects:
                    n = len(s.times[a])
                    for i in np.where(np.isclose(s.times[a], t))[0]:
                        obs.append(s.logconc[a][i])
                        sim.append(sims[s.sid][:, j * n + i])
                obs = np.array(obs)
                sim_all = np.concatenate(sim)
                lo, hi = np.percentile(sim_all, [5, 95])
                ok = (obs >= lo) & (obs <= hi)
                inside[a] += ok.tolist()
                if t in (0.0, 8.0):
                    trough[a] += ok.tolist()
                    if t == 0.0:
                        sim_mat = np.stack(sim, axis=1)                 # replicate x patient
                        sd_sim[a] = [float(np.median(sim_mat.std(axis=1, ddof=1))),
                                     float(np.percentile(sim_mat.std(axis=1, ddof=1), 5)),
                                     float(np.percentile(sim_mat.std(axis=1, ddof=1), 95)), float(obs.std(ddof=1))]
        for a in M.ANALYTES:
            rows.append(dict(model=name, analyte=a, pct_obs_within_sim_90=100 * np.mean(inside[a]), n_obs=len(inside[a]),
                             pct_trough_within_sim_90=100 * np.mean(trough[a]), n_trough=len(trough[a]),
                             sd_log_trough_0h_observed=sd_sim[a][3], sd_log_trough_0h_sim_median=sd_sim[a][0],
                             sd_log_trough_0h_sim_p5=sd_sim[a][1], sd_log_trough_0h_sim_p95=sd_sim[a][2]))
        allin = inside["caz"] + inside["avi"]
        rows.append(dict(model=name, analyte="both", pct_obs_within_sim_90=100 * np.mean(allin), n_obs=len(allin),
                         pct_trough_within_sim_90=100 * np.mean(trough["caz"] + trough["avi"]), n_trough=len(trough["caz"]) * 2,
                         sd_log_trough_0h_observed="", sd_log_trough_0h_sim_median="", sd_log_trough_0h_sim_p5="",
                         sd_log_trough_0h_sim_p95=""))
    path = os.path.join(M.OUT, "model1_vpc_structures.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    for r in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
