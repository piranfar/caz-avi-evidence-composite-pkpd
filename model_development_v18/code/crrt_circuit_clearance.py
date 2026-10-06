"""Circuit and other clearance in the CRRT cohort, and what a ceftazidime trough says about avibactam.

1. Circuit and non-circuit clearance. Li et al measured pre- and post-filter plasma concentrations of both drugs.
   Their protocol kept blood flow at 160 mL/min and replacement fluid or dialysate at 2 L/h (post-dilution CVVH,
   n = 6; CVVHD, n = 15). With their formulas (CVVHD: Qout = Qp - Quf; post-dilution CVVH: Qout = Qp - Quf - QR;
   Qp = blood flow x (1 - haematocrit)), the circuit clearance of each drug over the sampled interval is
   CL_CRRT = Qp - Qout x AUCpost / AUCpre. Net ultrafiltration is deposited only as a category, so it is set to
   0 L/h and varied (0.1 and 0.2 L/h) in a sensitivity analysis. Total clearance is the Model 1 empirical Bayes
   estimate (non-compartmental dose/AUC as a check); non-circuit clearance = total - circuit.
2. What a ceftazidime trough says about avibactam. The between-patient variability of avibactam clearance left
   unexplained by ceftazidime clearance, CV of omega_avi x sqrt(1 - rho^2), here and in the source model without
   renal replacement therapy (Cojutti et al); and, in the virtual CRRT population, the 90% prediction interval of
   the free avibactam trough given one measured ceftazidime trough (model residual error), at the fitted
   correlation, its profile-likelihood bounds and 0.94. Random stream SEED + 5.
3. Uncertainty of every Model 1 parameter from the asymptotic covariance written by joint_popk_nlme.py.

Values are written at full precision (6 decimals) and rounded once, when they are reported.
Outputs: outputs/crrt_extracorporeal_components.csv, outputs/crrt_extracorporeal_summary.csv,
         outputs/crrt_conditional_avibactam.csv, outputs/model1_parameter_ci.csv
"""
import csv
import math
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import crrt_virtual_tdm as sim  # noqa: E402

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "outputs")
DRY = os.path.join(ROOT, "data_external", "dryad_Li2025_CRRT")
COJ = os.path.join(ROOT, "data_external", "Cojutti2024_ANCHOR_PopPK", "Cojutti2024_PopPK_parameters.csv")
QB_L_H = 160 * 60 / 1000          # blood flow, 160 mL/min (Li et al, Methods)
QR_L_H = 2.0                      # replacement fluid or dialysate flow, 2 L/h (Li et al, Methods)
QUF_GRID = (0.0, 0.1, 0.2)        # net ultrafiltration, L/h (category only in the deposit)
DOSE = {"caz": 2000.0, "avi": 500.0}
LI_REPORTED = {("CVVHD", "caz"): (33.98, 5.62), ("CVVHD", "avi"): (35.53, 7.79),
               ("CVVH", "caz"): (40.80, 2.14), ("CVVH", "avi"): (43.19, 2.94)}   # mL/min, mean and SD
D = 6


def rows(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write(name, recs):
    path = os.path.join(OUT, name)
    keys = list(recs[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        for r in recs:
            w.writerow({k: (f"{v:.{D}f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"  wrote {os.path.relpath(path, ROOT)} ({len(recs)} rows)")


def auc(t, c):
    o = np.argsort(t)
    return float(np.trapezoid(np.asarray(c)[o], np.asarray(t)[o]))


def fisher_ci(r, n):
    z, se = math.atanh(r), 1 / math.sqrt(n - 3)
    return math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)


# --------------------------------------------------------------------------- 1. circuit components
def components():
    print("1. Circuit and non-circuit clearance")
    conc = {}
    for drug, fname, pre, post in (("caz", "Ceftazidime_concentration.csv", "caz_pre", "caz_post"),
                                   ("avi", "Avibactam_concentration.csv", "avi_pre", "avi_post")):
        for r in rows(os.path.join(DRY, fname)):
            conc.setdefault((int(r["subjectID"]), drug), []).append((float(r["Time_h"]), float(r[pre]), float(r[post])))
    hct = {int(r["subjectID"]): float(r["HCT"]) for r in rows(os.path.join(DRY, "Demographic_data.csv"))}
    urine = {int(r["subjectID"]): float(r["urine_24h_mL"]) for r in rows(os.path.join(DRY, "Demographic_data.csv"))}
    mod = {int(r["subjectID"]): int(r["crrt_modality"]) for r in rows(os.path.join(DRY, "CRRT_parameters_data.csv"))}
    counts = {k: sum(1 for v in mod.values() if v == k) for k in set(mod.values())}
    assert counts == {0: 15, 1: 6}, counts           # Li et al: 15 CVVHD, 6 post-dilution CVVH
    label = {0: "CVVHD", 1: "CVVH"}
    ebe = {int(r["subjectID"]): {"caz": float(r["CL_caz_L_h"]), "avi": float(r["CL_avi_L_h"])}
           for r in rows(os.path.join(OUT, "model1_individual_parameters.csv"))}
    recs = []
    for sid in sorted(mod):
        for quf in QUF_GRID:
            rec = dict(subjectID=sid, modality=label[mod[sid]], hct=hct[sid], urine_24h_mL=urine[sid], quf_l_h=quf)
            qp = QB_L_H * (1 - hct[sid])
            qout = qp - quf - (QR_L_H if mod[sid] == 1 else 0.0)
            for drug in ("caz", "avi"):
                t, cpre, cpost = map(np.array, zip(*conc[(sid, drug)]))
                a_pre, a_post = auc(t, cpre), auc(t, cpost)
                cl_crrt = qp - qout * a_post / a_pre
                rec[f"cl_crrt_{drug}_l_h"] = cl_crrt
                rec[f"cl_total_ebe_{drug}_l_h"] = ebe[sid][drug]
                rec[f"cl_total_nca_{drug}_l_h"] = DOSE[drug] / a_pre
                rec[f"cl_noncircuit_{drug}_l_h"] = ebe[sid][drug] - cl_crrt
                rec[f"circuit_share_{drug}"] = cl_crrt / ebe[sid][drug]
            recs.append(rec)
    write("crrt_extracorporeal_components.csv", recs)

    summ = []

    def add(name, value, n="", note=""):
        summ.append(dict(quantity=name, value=float(value), n=n, note=note))

    for quf in QUF_GRID:
        sub = [r for r in recs if r["quf_l_h"] == quf]
        tag = f"net ultrafiltration {quf:.1f} L/h"
        for (m, drug), (mean_li, sd_li) in LI_REPORTED.items():
            v = np.array([r[f"cl_crrt_{drug}_l_h"] for r in sub if r["modality"] == m]) * 1000 / 60
            add(f"cl_crrt_mL_min_{m}_{drug}", v.mean(), len(v),
                f"{tag}; SD {v.std(ddof=1):.6f}; Li et al {mean_li} (SD {sd_li}); ratio {v.mean() / mean_li:.6f}")
        for drug in ("caz", "avi"):
            s = np.array([r[f"circuit_share_{drug}"] for r in sub])
            add(f"circuit_share_median_{drug}", float(np.median(s)), len(s), f"{tag}; range {s.min():.6f}-{s.max():.6f}")
            for part in ("cl_crrt", "cl_noncircuit", "cl_total_ebe"):
                v = np.array([r[f"{part}_{drug}_l_h"] for r in sub])
                add(f"{part}_sd_l_h_{drug}", float(v.std(ddof=1)), len(v), f"{tag}; mean {v.mean():.6f}")
        for part in ("cl_crrt", "cl_noncircuit", "cl_total_ebe", "cl_total_nca"):
            a = np.array([r[f"{part}_caz_l_h"] for r in sub])
            b = np.array([r[f"{part}_avi_l_h"] for r in sub])
            r_p = float(np.corrcoef(a, b)[0, 1])
            lo, hi = fisher_ci(r_p, len(a))
            r_s = float(stats.spearmanr(a, b).statistic)
            add(f"corr_{part}_caz_avi", r_p, len(a), f"{tag}; Pearson, Fisher 95% interval {lo:.6f} to {hi:.6f}; Spearman {r_s:.6f}")
        for drug in ("caz", "avi"):
            for part in ("cl_crrt", "cl_noncircuit"):
                a = np.array([r[f"{part}_{drug}_l_h"] for r in sub if r["modality"] == "CVVH"])
                b = np.array([r[f"{part}_{drug}_l_h"] for r in sub if r["modality"] == "CVVHD"])
                p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
                add(f"modality_{part}_{drug}_p", p, len(a) + len(b),
                    f"{tag}; CVVH median {np.median(a):.6f}, CVVHD median {np.median(b):.6f} L/h; Mann-Whitney")
            nc = np.array([r[f"cl_noncircuit_{drug}_l_h"] for r in sub])
            u = np.array([r["urine_24h_mL"] for r in sub])
            add(f"spearman_urine_vs_noncircuit_{drug}", float(stats.spearmanr(u, nc).statistic), len(u), tag)
    write("crrt_extracorporeal_summary.csv", summ)
    return recs, summ


# --------------------------------------------------------------------------- 2. conditional avibactam
def conditional():
    print("2. Avibactam given a ceftazidime trough")
    cv = lambda w: 100 * math.sqrt(math.exp(w * w) - 1)  # noqa: E731
    coj = {(r["parameter"], r["drug"]): r for r in rows(COJ)}
    cv_avi_coj = float(coj[("Omega_CL_CV", "AVI")]["value"]) / 100
    rho_coj = float(coj[("CL_CAZ_CL_AVI", "both")]["value"])
    w_coj = math.sqrt(math.log(1 + cv_avi_coj ** 2))
    recs = []
    scen = [("estimate", sim.RHO_EST), ("lower bound", sim.RHO_LO), ("upper bound", sim.RHO_HI), ("0.94", 0.94)]
    for name, rho in scen:
        recs.append(dict(population="CRRT (Model 1)", regimen="", scenario=name, rho=float(rho), band="",
                         quantity="unexplained_cv_cl_avi_pct", value=cv(sim.OM["cl_avi"] * math.sqrt(1 - rho ** 2)),
                         total=cv(sim.OM["cl_avi"])))
    recs.append(dict(population="without RRT (Cojutti et al)", regimen="", scenario="source model", rho=rho_coj, band="",
                     quantity="unexplained_cv_cl_avi_pct", value=cv(w_coj * math.sqrt(1 - rho_coj ** 2)), total=100 * cv_avi_coj))
    rng = np.random.default_rng(sim.SEED + 5)
    n = 400_000
    for lab, cz, av, tau in (sim.REGIMENS[0], sim.REGIMENS[1]):
        for name, rho in scen:
            p = sim.population(rng, n, rho)
            cmin_caz, cmin_avi, _ = sim.exposures(p, cz, av, tau)
            f_avi = sim.FU_AVI * cmin_avi
            obs = np.log(cmin_caz) + sim.SIG_CAZ * rng.standard_normal(n)
            q5, q95 = np.percentile(f_avi, [5, 95])
            recs.append(dict(population="CRRT (Model 1)", regimen=lab, scenario=name, rho=float(rho), band="all patients",
                             quantity="free_avi_trough_p5_p95_mg_L", value=float(q5), total=float(q95)))
            for band, (a, b) in (("measured ceftazidime trough, 47.5th-52.5th percentile", (47.5, 52.5)),
                                 ("measured ceftazidime trough, 7.5th-12.5th percentile", (7.5, 12.5))):
                lo, hi = np.percentile(obs, [a, b])
                sel = f_avi[(obs >= lo) & (obs <= hi)]
                c5, c95 = np.percentile(sel, [5, 95])
                recs.append(dict(population="CRRT (Model 1)", regimen=lab, scenario=name, rho=float(rho), band=band,
                                 quantity="free_avi_trough_p5_p95_mg_L", value=float(c5), total=float(c95)))
    write("crrt_conditional_avibactam.csv", recs)
    return recs


# --------------------------------------------------------------------------- 3. parameter uncertainty
def parameter_ci():
    print("3. Uncertainty of the Model 1 parameters")
    cov = rows(os.path.join(OUT, "model1_parameter_covariance.csv"))
    names = [r["parameter"] for r in cov]
    m = np.array([[float(r[k]) for k in names] for r in cov])
    assert np.allclose(m, m.T, atol=1e-10)
    recs = []
    for i, r in enumerate(cov):
        e, se = float(r["estimate"]), math.sqrt(m[i, i])
        nm = r["parameter"]
        if nm.startswith("atanh_"):
            est, lo, hi, scale = math.tanh(e), math.tanh(e - 1.96 * se), math.tanh(e + 1.96 * se), "Fisher z"
            rse = float("nan")
        else:
            est, lo, hi, scale = math.exp(e), math.exp(e - 1.96 * se), math.exp(e + 1.96 * se), "log"
            rse = 100 * se
        rec = dict(parameter=nm.replace("log_", "").replace("atanh_", ""), estimation_scale=scale, estimate=est,
                   se_estimation_scale=se, rse_pct=rse, wald_ci_low=lo, wald_ci_high=hi)
        if nm.startswith("log_omega"):
            cvf = lambda w: 100 * math.sqrt(math.exp(w * w) - 1)  # noqa: E731
            rec.update(cv_pct=cvf(est), cv_ci_low=cvf(lo), cv_ci_high=cvf(hi))
        else:
            rec.update(cv_pct=float("nan"), cv_ci_low=float("nan"), cv_ci_high=float("nan"))
        recs.append(rec)
    eig = np.linalg.eigvalsh(m)
    print(f"  covariance: {len(names)} parameters, smallest eigenvalue {eig.min():.3e}")
    write("model1_parameter_ci.csv", recs)
    return recs


if __name__ == "__main__":
    components()
    conditional()
    parameter_ci()
