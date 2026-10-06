"""Measurement error in circuit and other clearance.

Circuit clearance comes from the ratio of postfilter to prefilter AUC (crrt_circuit_clearance.py), so errors in single
concentrations reach it directly. Parametric bootstrap, 2,000 replicates, random stream SEED + 7: every deposited
concentration gets uniform noise of +/-0.25 mg/L (the deposit is rounded to 0.5 mg/L) and a log-normal error with
Model 1's residual SD for that drug; the 2 drugs measured in the same tube share the residual correlation c, and the
prefilter and postfilter tubes are independent. Each replicate recomputes circuit clearance (Li et al's formula, net
ultrafiltration 0 L/h) and non-compartmental total clearance (dose / AUC over the sampled interval); other clearance =
total - circuit, and also empirical Bayes total (Model 1, held fixed) - circuit, the definition used in
crrt_circuit_clearance.py. Model 1's residual SD includes model misfit as well as assay error, so it bounds the measurement
error from above; half that SD is the sensitivity analysis.
A second bootstrap resamples patients and then measurement errors (2,000 replicates, same stream) for the
correlations and CVs of the components across patients.

Values are written at full precision (6 decimals) and rounded once, when they are reported.
Outputs: outputs/crrt_circuit_uncertainty.csv (per patient), outputs/crrt_circuit_uncertainty_summary.csv
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import crrt_circuit_clearance as CC  # noqa: E402
import crrt_virtual_tdm as sim  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
N_REP = 2000
D = 6


def load():
    conc = {}
    for drug, fname, pre, post in (("caz", "Ceftazidime_concentration.csv", "caz_pre", "caz_post"),
                                   ("avi", "Avibactam_concentration.csv", "avi_pre", "avi_post")):
        for r in CC.rows(os.path.join(CC.DRY, fname)):
            conc.setdefault(int(r["subjectID"]), {}).setdefault(drug, []).append(
                (float(r["Time_h"]), float(r[pre]), float(r[post])))
    hct = {int(r["subjectID"]): float(r["HCT"]) for r in CC.rows(os.path.join(CC.DRY, "Demographic_data.csv"))}
    mod = {int(r["subjectID"]): int(r["crrt_modality"]) for r in CC.rows(os.path.join(CC.DRY, "CRRT_parameters_data.csv"))}
    pats = []
    for sid in sorted(conc):
        caz, avi = sorted(conc[sid]["caz"]), sorted(conc[sid]["avi"])
        t = np.array([x[0] for x in caz])
        assert np.array_equal(t, [x[0] for x in avi]), sid                  # both drugs from the same tubes
        pats.append(dict(sid=sid, t=t, qp=CC.QB_L_H * (1 - hct[sid]),
                         qout=CC.QB_L_H * (1 - hct[sid]) - (CC.QR_L_H if mod[sid] == 1 else 0.0),
                         pre=np.array([[x[1] for x in caz], [x[1] for x in avi]]),
                         post=np.array([[x[2] for x in caz], [x[2] for x in avi]])))
    fp = {r["parameter"]: float(r["estimate"]) for r in CC.rows(os.path.join(CC.OUT, "model1_final_parameters.csv"))}
    return pats, (fp["sigma_prop_caz"], fp["sigma_prop_avi"]), fp["corr_residual_caz_avi"]


def components(p, pre, post):
    """Circuit, non-compartmental total and other clearance for both drugs (arrays of 2)."""
    dose = np.array([CC.DOSE["caz"], CC.DOSE["avi"]])
    a_pre = np.trapezoid(pre, p["t"], axis=-1)
    a_post = np.trapezoid(post, p["t"], axis=-1)
    circ = p["qp"] - p["qout"] * a_post / a_pre
    tot = dose / a_pre
    return circ, tot, tot - circ


def perturb(rng, x, sig, c, scale):
    """Rounding noise, then a correlated log-normal error per tube (rows: drugs; columns: tubes)."""
    n = x.shape[1]
    z = rng.standard_normal((2, n))
    e = np.vstack([z[0], c * z[0] + np.sqrt(1 - c * c) * z[1]]) * (scale * np.array(sig))[:, None]
    return np.clip(x + rng.uniform(-0.25, 0.25, x.shape), 0.05, None) * np.exp(e)


def main():
    pats, sig, c = load()
    n = len(pats)
    rows, summ = [], []
    base = np.array([[*components(p, p["pre"], p["post"])] for p in pats])        # [patient, part, drug]
    ref = {int(r["subjectID"]): r for r in CC.rows(os.path.join(CC.OUT, "crrt_extracorporeal_components.csv"))
           if float(r["quf_l_h"]) == 0}
    for i, p in enumerate(pats):                        # same formula as crrt_circuit_clearance.py
        for k, drug in enumerate(("caz", "avi")):
            assert abs(base[i, 0, k] - float(ref[p["sid"]][f"cl_crrt_{drug}_l_h"])) < 1e-5, (p["sid"], drug)
            assert abs(base[i, 1, k] - float(ref[p["sid"]][f"cl_total_nca_{drug}_l_h"])) < 1e-5, (p["sid"], drug)
    names = ("circuit", "total_nca", "other", "other_ebe")
    ebe_all = np.array([[float(ref[p["sid"]][f"cl_total_ebe_{d}_l_h"]) for d in ("caz", "avi")] for p in pats])
    base = np.concatenate([base, (ebe_all - base[:, 0, :])[:, None, :]], axis=1)   # other = EBE total - circuit
    reps = {}
    for scale in (1.0, 0.5):
        rng = np.random.default_rng([sim.SEED, 7, int(scale * 10)])
        out = np.empty((N_REP, n, 4, 2))
        for b in range(N_REP):
            for i, p in enumerate(pats):
                out[b, i, :3] = components(p, perturb(rng, p["pre"], sig, c, scale), perturb(rng, p["post"], sig, c, scale))
        out[:, :, 3, :] = ebe_all[None, :, :] - out[:, :, 0, :]
        reps[scale] = out
    full = reps[1.0]
    for i, p in enumerate(pats):
        for k, drug in enumerate(("caz", "avi")):
            rec = dict(subjectID=p["sid"], drug=drug)
            for j, nm in enumerate(names):
                v = full[:, i, j, k]
                rec[f"{nm}_l_h"] = float(base[i, j, k])
                rec[f"{nm}_sd"] = float(v.std(ddof=1))
                rec[f"{nm}_p5"], rec[f"{nm}_p95"] = (float(q) for q in np.percentile(v, [5, 95]))
            share = full[:, i, 0, k] / full[:, i, 1, k]
            rec["circuit_share"] = float(base[i, 0, k] / base[i, 1, k])
            rec["circuit_share_p5"], rec["circuit_share_p95"] = (float(q) for q in np.percentile(share, [5, 95]))
            rows.append(rec)

    def add(q, value, note=""):
        summ.append(dict(quantity=q, value=float(value), note=note))

    for scale, out in reps.items():
        tag = "residual SD" if scale == 1.0 else "half the residual SD"
        for k, drug in enumerate(("caz", "avi")):
            for j, nm in enumerate(names):
                err_var = out[:, :, j, k].var(axis=0, ddof=1)
                obs_var = base[:, j, k].var(ddof=1)
                add(f"median_patient_sd_{nm}_{drug}", np.median(np.sqrt(err_var)), tag)
                add(f"between_patient_sd_{nm}_{drug}", np.sqrt(obs_var), tag)
                add(f"error_share_of_variance_{nm}_{drug}", err_var.mean() / obs_var, tag)
            share = out[:, :, 0, k] / out[:, :, 1, k]
            add(f"replicates_with_circuit_share_above_1_{drug}", float(np.mean((share > 1).any(axis=1))), tag)
            add(f"median_circuit_share_{drug}", float(np.median(base[:, 0, k] / base[:, 1, k])),
                f"{tag}; 95% interval of the median {np.percentile(np.median(share, axis=1), 2.5):.6f} to "
                f"{np.percentile(np.median(share, axis=1), 97.5):.6f}")
            ebe = np.array([float(ref[p["sid"]][f"cl_total_ebe_{drug}_l_h"]) for p in pats])
            se = out[:, :, 0, k] / ebe[None, :]
            add(f"median_circuit_share_ebe_{drug}", float(np.median(base[:, 0, k] / ebe)),
                f"{tag}; empirical Bayes total as denominator (as in crrt_extracorporeal_components.csv); 95% interval of "
                f"the median {np.percentile(np.median(se, axis=1), 2.5):.6f} to {np.percentile(np.median(se, axis=1), 97.5):.6f}")
        # patients and errors resampled together
        rng = np.random.default_rng([sim.SEED, 7, 100 + int(scale * 10)])
        stats = {f"corr_{nm}": [] for nm in names} | {f"cv_{nm}_{d}": [] for nm in names for d in ("caz", "avi")}
        for _ in range(N_REP):
            idx = rng.integers(0, n, n)
            v = out[rng.integers(0, N_REP), idx]                     # [patient, part, drug]
            for j, nm in enumerate(names):
                stats[f"corr_{nm}"].append(np.corrcoef(v[:, j, 0], v[:, j, 1])[0, 1])
                for k, d in enumerate(("caz", "avi")):
                    stats[f"cv_{nm}_{d}"].append(v[:, j, k].std(ddof=1) / v[:, j, k].mean())
        for key, vals in stats.items():
            vals = np.array(vals)
            if key.startswith("corr_"):
                j = names.index(key[5:])
                point = float(np.corrcoef(base[:, j, 0], base[:, j, 1])[0, 1])
            else:
                nm, d = key[3:].rsplit("_", 1)
                j, k = names.index(nm), ("caz", "avi").index(d)
                point = float(base[:, j, k].std(ddof=1) / base[:, j, k].mean())
            lo, hi = np.percentile(vals, [2.5, 97.5])
            add(key, point, f"{tag}; patients and errors resampled: 95% interval {lo:.6f} to {hi:.6f}")
        diff = np.array(stats["corr_circuit"]) - np.array(stats["corr_other"])
        lo, hi = np.percentile(diff, [2.5, 97.5])
        add("corr_circuit_minus_corr_other", float(np.median(diff)), f"{tag}; 95% interval {lo:.6f} to {hi:.6f}")
        diff = np.array(stats["corr_circuit"]) - np.array(stats["corr_other_ebe"])
        lo, hi = np.percentile(diff, [2.5, 97.5])
        add("corr_circuit_minus_corr_other_ebe", float(np.median(diff)), f"{tag}; 95% interval {lo:.6f} to {hi:.6f}")
        for d in ("caz", "avi"):
            dd = np.array(stats[f"cv_circuit_{d}"]) - np.array(stats[f"cv_total_nca_{d}"])
            lo, hi = np.percentile(dd, [2.5, 97.5])
            add(f"cv_circuit_minus_cv_total_{d}", float(np.median(dd)), f"{tag}; 95% interval {lo:.6f} to {hi:.6f}")

    for name, recs in (("crrt_circuit_uncertainty.csv", rows), ("crrt_circuit_uncertainty_summary.csv", summ)):
        path = os.path.join(CC.OUT, name)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(recs[0]), lineterminator="\n")
            w.writeheader()
            for r in recs:
                w.writerow({k: (f"{v:.{D}f}" if isinstance(v, float) else v) for k, v in r.items()})
        print(f"  wrote outputs/{name} ({len(recs)} rows)")
    for r in summ:
        if r["note"].startswith("residual SD"):
            print(f"  {r['quantity']:48s} {r['value']:.4f}  {r['note'][13:]}")


if __name__ == "__main__":
    main()
