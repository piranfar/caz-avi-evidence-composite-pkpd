"""Supplementary analyses for the CRRT virtual-population analysis.

1. Decision thresholds. For the 2 regimens at which avibactam approaches its target (1.25 g every
   8 h and 2.5 g every 12 h), the probability of avibactam attainment given the measured
   ceftazidime trough is recomputed for the same virtual patients as outputs/crrt_classifier.csv.
   The random streams of crrt_virtual_tdm.main() are replayed call by call, so the 0.5 threshold
   reproduces crrt_classifier.csv exactly (checked below). Reported per scenario: the area under
   the receiver operating characteristic curve (AUROC) of the measured trough and of the
   probability for separating patients above and below the avibactam target; detection and
   flagging at thresholds of 0.5, 0.8, 0.9 and 0.95; and the full trade-off curve.

2. Checks on the deposited cohort (Dryad doi:10.5061/dryad.fxpnvx16s):
   - deposited sampling times against the published schedule for each infusion-duration category
     (1 h: 0,1,2,3,4,6,8 h; 2 h: 0,2,3,4,6,8 h; 3 h: 0,3,4,6,8 h);
   - total ceftazidime above 104 mg/L in any sample, in the trough, and as the interval average
     (linear trapezoid over the 8-h steady-state interval);
   - the individual avibactam-to-ceftazidime clearance ratio (Model 1 empirical Bayes estimates)
     against urine output and serum creatinine (Spearman) and between the 2 deposited CRRT
     modality codes (two-sided Mann-Whitney test).

3. Second CRRT cohort: correlation of log clearances in the 17 occasions of Gatti et al
   (J Crit Care 2023;76:154301, Table 2; 8 patients on CVVHDF, continuous infusion), with a
   patient-cluster bootstrap interval.

4. Model-free check in the deposited cohort: correlation of log non-compartmental clearances.

Values are written at full precision (6 decimals) and rounded once, when they are reported.
Outputs: outputs/crrt_decision_threshold.csv, outputs/crrt_triage_curve.csv,
         outputs/crrt_equal_detection.csv, outputs/crrt_data_checks.csv
"""
from __future__ import annotations

import csv
import math
import os
import sys

import numpy as np
from scipy import stats

import crrt_virtual_tdm as sim

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "outputs")
DRY = os.path.join(ROOT, "data_external", "dryad_Li2025_CRRT")

THRESHOLDS = (0.5, 0.8, 0.9, 0.95)
TARGET_REGIMENS = ("1.25 g q8h", "2.5 g q12h")
OMEGA = 1.0
DECIMALS = 6
BOOT = 20_000
SCHEDULE = {"1": [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0], "2": [0.0, 2.0, 3.0, 4.0, 6.0, 8.0],
            "3": [0.0, 3.0, 4.0, 6.0, 8.0]}


def rows(path):
    with open(path, encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def classify_detail(rng, rho_true, rho_model, caz_mg, avi_mg, tau, om_scale=1.0):
    """sim.classify, returning per-patient truth, posterior and measured log trough.
    Draws from rng in the same order, so the patients are the ones behind crrt_classifier.csv."""
    p = sim.population(rng, sim.N, rho_true, om_scale)
    cmin_caz, cmin_avi, _ = sim.exposures(p, caz_mg, avi_mg, tau)
    truth = sim.FU_AVI * cmin_avi >= sim.AVI_CT
    obs = np.log(cmin_caz) + sim.SIG_CAZ * rng.standard_normal(sim.N)
    grid = np.linspace(obs.min() - 0.05, obs.max() + 0.05, 300)
    curve = sim.posterior_curve(rng, rho_model, caz_mg, avi_mg, tau, om_scale, grid)
    return truth, np.interp(obs, grid, curve), obs


def replay(seed, regimens):
    """Replay crrt_virtual_tdm.main() for one random stream; keep the omega x1.0 classifier draws."""
    rng = np.random.default_rng(seed)
    for om in (1.0, 1.5, 2.0):
        for lab, cz, av, tau in regimens:
            sim.attainment_row(rng, om, lab, cz, av, tau)
    kept = {}
    for lab, cz, av, tau in regimens:          # omega 1.0 is the first classifier block
        for name, rt, rm in sim.SCENARIOS:
            truth, post, obs = classify_detail(rng, rt, rm, cz, av, tau, OMEGA)
            if lab in TARGET_REGIMENS:
                kept[(lab, name)] = (rt, rm, truth, post, obs)
    return kept


def auroc(score, positive):
    """Probability that a randomly chosen positive scores above a randomly chosen negative (ties count 1/2)."""
    u = stats.mannwhitneyu(score[positive], score[~positive], alternative="two-sided").statistic
    return float(u) / (positive.sum() * (~positive).sum())


def threshold_rows(kept, published):
    out, curve = [], []
    grid = np.unique(np.concatenate([np.linspace(0.0, 1.0, 201), 1.0 - np.logspace(-4, -1, 61)]))
    for (lab, name), (rt, rm, truth, post, obs) in kept.items():
        n_nan = int(np.isnan(post).sum())
        assert n_nan == 0, (lab, name, n_nan)
        fail = ~truth
        # The 0.5 rule must reproduce crrt_classifier.csv exactly (same patients, same rule).
        pred = post >= 0.5
        ref = published[(lab, name)]
        assert abs(100 * np.mean(pred & fail) - float(ref["false_reassurance"])) < 1e-6, (lab, name)
        assert abs(100 * fail.mean() - float(ref["missed_by_population_prior"])) < 1e-6, (lab, name)
        auc_post = auroc(post, truth)
        auc_trough = auroc(obs, truth)
        # The probability rises with the measured trough, so each threshold is a trough cut-off.
        order = np.argsort(obs)
        mono_viol = int(np.sum(np.diff(post[order]) < -1e-9))
        for t in THRESHOLDS:
            flag = ~(post >= t)
            cut = float(np.exp(obs[~flag].min())) if (~flag).any() else float("nan")
            n_out = int(np.sum(flag & (obs >= np.log(cut)))) if (~flag).any() else 0
            out.append(dict(
                omega_scale=OMEGA, regimen=lab, scenario=name, rho_true=rt, rho_model=rm,
                below_target_pct=100 * fail.mean(), auroc_trough=auc_trough, auroc_probability=auc_post,
                threshold=t, total_caz_trough_cutoff_mg_L=cut, patients_off_cutoff=n_out,
                detected_pct_of_below_target=100 * (flag & fail).sum() / fail.sum(),
                attainers_flagged_pct=100 * (flag & truth).sum() / truth.sum(),
                flagged_pct_of_all=100 * flag.mean(),
                wrongly_reassured_pct_of_all=100 * (~flag & fail).mean(),
                monotonicity_violations=mono_viol))
        for t in grid:
            flag = ~(post >= t)
            curve.append(dict(regimen=lab, scenario=name, threshold=float(t),
                              flagged_pct_of_all=100 * flag.mean(),
                              detected_pct_of_below_target=100 * (flag & fail).sum() / fail.sum()))
    return out, curve


def nca_correlation(caz_rows):
    """Model-free check: CL = dose / AUC(0-8 h) by linear trapezoid on pre-filter concentrations."""
    avi_rows = rows(os.path.join(DRY, "Avibactam_concentration.csv"))
    dose = {"caz": 2000.0, "avi": 500.0}
    cl = {}
    for drug, rr, col in (("caz", caz_rows, "caz_pre"), ("avi", avi_rows, "avi_pre")):
        by = {}
        for r in rr:
            by.setdefault(r["subjectID"], []).append((float(r["Time_h"]), float(r[col])))
        for sid, v in by.items():
            v.sort()
            t = np.array([a for a, _ in v]); c = np.array([b for _, b in v])
            assert t[0] == 0.0 and t[-1] == 8.0, (drug, sid)
            cl.setdefault(sid, {})[drug] = dose[drug] / float(np.sum(np.diff(t) * (c[1:] + c[:-1]) / 2))
    ids = sorted(cl, key=int)
    lc = np.log([cl[s]["caz"] for s in ids]); la = np.log([cl[s]["avi"] for s in ids])
    r = float(np.corrcoef(lc, la)[0, 1])
    z, se = np.arctanh(r), 1.0 / np.sqrt(len(ids) - 3)
    return dict(check="nca_corr_log_cl", value=r, n=len(ids),
                note=f"non-compartmental clearances (dose / AUC 0-8 h, linear trapezoid); Fisher 95% interval "
                     f"{np.tanh(z - 1.959964 * se):.6f} to {np.tanh(z + 1.959964 * se):.6f}")


def data_checks():
    caz = rows(os.path.join(DRY, "Ceftazidime_concentration.csv"))
    med = {r["subjectID"]: r["infusion_duration_cat"] for r in rows(os.path.join(DRY, "Medication_Information.csv"))}
    demo = {r["subjectID"]: r for r in rows(os.path.join(DRY, "Demographic_data.csv"))}
    crrt = {r["subjectID"]: r for r in rows(os.path.join(DRY, "CRRT_parameters_data.csv"))}
    ip = {r["subjectID"]: r for r in rows(os.path.join(OUT, "model1_individual_parameters.csv"))}

    by = {}
    for r in caz:
        by.setdefault(r["subjectID"], []).append((float(r["Time_h"]), float(r["caz_pre"])))
    ids = sorted(by, key=int)
    n = len(ids)
    match = [sid for sid in ids if sorted(t for t, _ in by[sid]) == SCHEDULE[med[sid]]]
    differ = [sid for sid in ids if sid not in match]
    with_1h = sum(1 for sid in ids if any(t == 1.0 for t, _ in by[sid]))

    any104 = trough104 = avg104 = 0
    cavg = []
    for sid in ids:
        v = sorted(by[sid])
        t = np.array([a for a, _ in v]); c = np.array([b for _, b in v])
        assert t[0] == 0.0 and t[-1] == 8.0, sid
        any104 += bool(c.max() > sim.TOX)
        trough104 += bool(c[0] > sim.TOX)
        a = float(np.sum(np.diff(t) * (c[1:] + c[:-1]) / 2)) / 8.0
        cavg.append(a)
        avg104 += bool(a > sim.TOX)

    ratio = np.array([math.log(float(ip[s]["CL_avi_L_h"]) / float(ip[s]["CL_caz_L_h"])) for s in ids])
    uo = np.array([float(demo[s]["urine_24h_mL"]) for s in ids])
    scr = np.array([float(demo[s]["SCr_umol_L"]) for s in ids])
    mod = np.array([crrt[s]["crrt_modality"] for s in ids])
    r_uo, p_uo = stats.spearmanr(uo, ratio)
    r_scr, p_scr = stats.spearmanr(scr, ratio)
    g0, g1 = ratio[mod == "0"], ratio[mod == "1"]
    mw = stats.mannwhitneyu(g0, g1, alternative="two-sided")

    def spearman_ci(r):
        """Approximate 95% interval, Fisher z with the Fieller variance 1.06 / (n - 3)."""
        z, se = np.arctanh(r), math.sqrt(1.06 / (n - 3))
        return f"95% CI {np.tanh(z - 1.959964 * se):.6f} to {np.tanh(z + 1.959964 * se):.6f} (Fieller)"

    # Observed troughs (pre-dose, 0-h samples) and the lowest avibactam concentration.
    avi = rows(os.path.join(DRY, "Avibactam_concentration.csv"))
    caz0 = np.array([c for sid in ids for t, c in by[sid] if t == 0.0])
    avi0 = np.array([float(r["avi_pre"]) for r in avi if float(r["Time_h"]) == 0.0])
    avi_all = np.array([float(r["avi_pre"]) for r in avi])
    every = np.concatenate([avi_all, [float(r["caz_pre"]) for r in caz]])
    n_half = int(np.sum(np.isclose(every * 2.0, np.round(every * 2.0))))
    li_min = 7.89     # Li et al, AAC 2026: "the lowest Cmin of CAZ and AVI were 37.60 mg/L and 7.89 mg/L"

    return [
        dict(check="sampling_schedule_matches_published", value=len(match), n=n,
             note="patients whose deposited sampling times equal the published schedule for their "
                  "infusion-duration category; differing: " + ", ".join(
                      f"patient {s} (category {med[s]}, times {sorted(t for t, _ in by[s])})" for s in differ)),
        dict(check="patients_with_1h_sample", value=with_1h, n=n, note="samples at 1 h after the start of infusion"),
        dict(check="caz_any_sample_gt104_patients", value=any104, n=n, note="total ceftazidime, any sampling time"),
        dict(check="caz_trough_gt104_patients", value=trough104, n=n, note="total ceftazidime at 0 h"),
        dict(check="caz_interval_average_gt104_patients", value=avg104, n=n,
             note="total ceftazidime, linear-trapezoid average over 0-8 h"),
        dict(check="caz_interval_average_median_mg_L", value=float(np.median(cavg)), n=n,
             note=f"range {min(cavg):.6f}-{max(cavg):.6f}"),
        dict(check="cl_ratio_avi_caz_median", value=float(np.exp(np.median(ratio))), n=n,
             note=f"range {np.exp(ratio.min()):.6f}-{np.exp(ratio.max()):.6f}; empirical Bayes clearances"),
        dict(check="spearman_urine_output_vs_log_cl_ratio", value=float(r_uo), n=n,
             note=f"P = {p_uo:.6f}; {spearman_ci(r_uo)}"),
        dict(check="spearman_serum_creatinine_vs_log_cl_ratio", value=float(r_scr), n=n,
             note=f"P = {p_scr:.6f}; {spearman_ci(r_scr)}"),
        dict(check="mannwhitney_modality_code_cl_ratio_p", value=float(mw.pvalue), n=n,
             note=f"code 0: n = {g0.size}, median ratio {np.exp(np.median(g0)):.6f}; "
                  f"code 1: n = {g1.size}, median ratio {np.exp(np.median(g1)):.6f}"),
        dict(check="urine_output_24h_median_mL", value=float(np.median(uo)), n=n,
             note=f"below 100 mL: {int(np.sum(uo < 100))}; at most 170 mL: {int(np.sum(uo <= 170))}"),
        dict(check="caz_trough_0h_median_mg_L", value=float(np.median(caz0)), n=int(caz0.size),
             note=f"total ceftazidime, pre-dose (0 h) sample; range {caz0.min():.6f}-{caz0.max():.6f}"),
        dict(check="avi_trough_0h_median_mg_L", value=float(np.median(avi0)), n=int(avi0.size),
             note=f"total avibactam, pre-dose (0 h) sample; range {avi0.min():.6f}-{avi0.max():.6f}"),
        dict(check="avi_min_any_sample_mg_L", value=float(avi_all.min()), n=int(avi_all.size),
             note=f"total; deposited values are rounded; Li et al report a lowest avibactam Cmin of {li_min} mg/L, "
                  f"free {sim.FU_AVI * li_min:.6f} mg/L = {sim.FU_AVI * li_min / sim.AVI_CT:.6f} times the "
                  f"{sim.AVI_CT:g}-mg/L target"),
        dict(check="deposited_values_multiple_of_0_5_mg_L", value=n_half, n=int(every.size),
             note="concentrations of both drugs that are exact multiples of 0.5 mg/L"),
    ]


def second_cohort():
    """Gatti et al, J Crit Care 2023;76:154301, Table 2: 8 patients on CVVHDF receiving continuous
    infusion, 17 occasions with both clearances (data_external/Gatti2023_individual_patient_data.csv).
    Pearson correlation of log clearances across occasions; 95% percentile interval from a bootstrap
    that resamples patients, because occasions within a patient are not independent."""
    g = [r for r in rows(os.path.join(ROOT, "data_external", "Gatti2023_individual_patient_data.csv"))
         if r["record_type"] == "tdm"]
    pid = np.array([r["id_case"] for r in g])
    lc = np.log([float(r["caz_cl_l_h"]) for r in g])
    la = np.log([float(r["avi_cl_l_h"]) for r in g])
    r_occ = float(np.corrcoef(lc, la)[0, 1])
    pts = sorted(set(pid), key=int)
    idx = {p: np.where(pid == p)[0] for p in pts}
    rng = np.random.default_rng(sim.SEED + 2)
    boot = []
    for _ in range(BOOT):
        take = np.concatenate([idx[p] for p in rng.choice(pts, len(pts), replace=True)])
        if np.ptp(lc[take]) > 0 and np.ptp(la[take]) > 0:
            boot.append(np.corrcoef(lc[take], la[take])[0, 1])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    mc = np.array([lc[idx[p]].mean() for p in pts]); ma = np.array([la[idx[p]].mean() for p in pts])
    r_pat = float(np.corrcoef(mc, ma)[0, 1])
    z, se = np.arctanh(r_pat), 1.0 / np.sqrt(len(pts) - 3)
    pat_lo, pat_hi = np.tanh(z - 1.959964 * se), np.tanh(z + 1.959964 * se)
    multi = [p for p in pts if len(idx[p]) >= 2]
    wc = np.concatenate([lc[idx[p]] - lc[idx[p]].mean() for p in multi])
    wa = np.concatenate([la[idx[p]] - la[idx[p]].mean() for p in multi])
    r_within = float(np.corrcoef(wc, wa)[0, 1])
    # Within-patient correlation = partial correlation given patient (len(multi) - 1 indicators):
    # Fisher z with SE 1 / sqrt(n - q - 3), q = len(multi) - 1.
    se_w = 1.0 / np.sqrt(wc.size - (len(multi) - 1) - 3)
    w_lo, w_hi = np.tanh(np.arctanh(r_within) - 1.959964 * se_w), np.tanh(np.arctanh(r_within) + 1.959964 * se_w)
    loo = [float(np.corrcoef(np.delete(lc, i), np.delete(la, i))[0, 1]) for i in range(len(g))]
    lpo = [float(np.corrcoef(np.delete(lc, idx[p]), np.delete(la, idx[p]))[0, 1]) for p in pts]
    # Within patients, both clearances against total effluent flow (log scale), and the within-patient
    # correlation of the 2 clearances after adjusting for effluent flow (partial correlation).
    le = np.log([float(r["total_effluent_ml_h"]) for r in g])
    we = np.concatenate([le[idx[p]] - le[idx[p]].mean() for p in multi])
    r_ce, r_ae = float(np.corrcoef(wc, we)[0, 1]), float(np.corrcoef(wa, we)[0, 1])
    r_part = (r_within - r_ce * r_ae) / math.sqrt((1 - r_ce ** 2) * (1 - r_ae ** 2))
    # Is between-patient variance detectable at all? One-way ANOVA of each drug's log clearance by patient.
    f_c = stats.f_oneway(*[lc[idx[p]] for p in pts])
    f_a = stats.f_oneway(*[la[idx[p]] for p in pts])
    return [
        dict(check="gatti2023_occasion_corr_log_cl", value=r_occ, n=len(g),
             note=f"{len(pts)} patients; patient-cluster bootstrap 95% interval {lo:.6f} to {hi:.6f} "
                  f"({len(boot)} of {BOOT} resamples usable, seed SEED + 2); "
                  f"excludes 0.94: {'yes' if hi < 0.94 else 'no'}"),
        dict(check="gatti2023_within_patient_corr_log_cl", value=r_within, n=int(wc.size),
             note=f"{len(multi)} patients with 2 or more occasions; log clearances centered on each patient's mean; "
                  f"Fisher 95% interval for a partial correlation {w_lo:.6f} to {w_hi:.6f}"),
        dict(check="gatti2023_patient_mean_corr_log_cl", value=r_pat, n=len(pts),
             note=f"correlation of per-patient mean log clearances; Fisher 95% interval {pat_lo:.6f} to {pat_hi:.6f}"),
        dict(check="gatti2023_leave_one_occasion_out_range", value=min(loo), n=len(g),
             note=f"minimum shown; maximum {max(loo):.6f}; dropped occasion at the minimum: "
                  f"patient {pid[int(np.argmin(loo))]}"),
        dict(check="gatti2023_leave_one_patient_out_range", value=min(lpo), n=len(pts),
             note=f"minimum shown; maximum {max(lpo):.6f}"),
        dict(check="gatti2023_within_patient_corr_cl_caz_vs_effluent", value=r_ce, n=int(wc.size),
             note="log ceftazidime clearance against log total effluent flow, both centered on each patient's mean"),
        dict(check="gatti2023_within_patient_corr_cl_avi_vs_effluent", value=r_ae, n=int(wc.size),
             note="log avibactam clearance against log total effluent flow, both centered on each patient's mean"),
        dict(check="gatti2023_within_patient_partial_corr_given_effluent", value=r_part, n=int(wc.size),
             note="within-patient correlation of the 2 log clearances after adjusting for log total effluent flow"),
        dict(check="gatti2023_between_patient_anova_caz_p", value=float(f_c.pvalue), n=len(g),
             note=f"one-way ANOVA of log ceftazidime clearance by patient, F = {f_c.statistic:.6f}"),
        dict(check="gatti2023_between_patient_anova_avi_p", value=float(f_a.pvalue), n=len(g),
             note=f"one-way ANOVA of log avibactam clearance by patient, F = {f_a.statistic:.6f}"),
    ]


DETECTION_LEVELS = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0)


def equal_detection(curve, thr):
    """Patients flagged (% of all) at each correlation scenario to identify a given percentage of the
    patients below the avibactam target; linear interpolation along each scenario's trade-off curve.
    Rows: the detection that threshold 0.9 gives at the estimated correlation, then fixed detection
    levels of 50% to 95%. ratio_estimate_to_nonrrt compares the estimated correlation with 0.94."""
    out = []
    for lab in TARGET_REGIMENS:
        cur = {}
        for name, _, _ in sim.SCENARIOS:
            pts = sorted({(r["detected_pct_of_below_target"], r["flagged_pct_of_all"])
                          for r in curve if r["regimen"] == lab and r["scenario"] == name})
            cur[name] = (np.array([p[0] for p in pts]), np.array([p[1] for p in pts]))
        ref = next(r for r in thr if r["regimen"] == lab and r["scenario"] == "estimate" and r["threshold"] == 0.9)
        levels = [("estimate, threshold 0.9", ref["detected_pct_of_below_target"])] + \
                 [("fixed detection level", v) for v in DETECTION_LEVELS]
        for kind, det in levels:
            row = dict(regimen=lab, reference=kind, detected_pct_of_below_target=det)
            for name, _, _ in sim.SCENARIOS:
                d, f = cur[name]
                row["flagged_pct_" + name.split(":")[0].replace(" ", "_").replace("-", "")] = float(np.interp(det, d, f))
            row["ratio_estimate_to_nonrrt"] = row["flagged_pct_estimate"] / row["flagged_pct_nonRRT_value"]
            out.append(row)
    return out


def seed_replication(n_seeds=5):
    """Monte Carlo precision at 1.25 g every 8 h: the estimate and non-RRT scenarios repeated in further
    random streams (SEED + 10 onwards), each with its own virtual patients, measurements and posterior.
    Per stream: below target, wrongly reassured (0.5 rule), detection and flagging at threshold 0.9, and
    the flagging ratio (estimate / 0.94) needed to identify 50% to 95% of patients below target, from the
    exact ranking of the measured trough (the classifier is a cut-off on the trough)."""
    lab, cz, av, tau = next(r for r in sim.REGIMENS if r[0] == "1.25 g q8h")
    out = []
    for k in range(n_seeds):
        seed = sim.SEED + 10 + k
        rng = np.random.default_rng(seed)
        need = {}
        for name, rt, rm in sim.SCENARIOS:
            if name not in ("estimate", "non-RRT value"):
                continue
            truth, post, obs = classify_detail(rng, rt, rm, cz, av, tau, OMEGA)
            fail = ~truth
            cum = np.cumsum(fail[np.argsort(obs)])
            need[name] = np.array([100.0 * (np.searchsorted(cum, lvl / 100.0 * fail.sum()) + 1) / fail.size
                                   for lvl in DETECTION_LEVELS])
            flag = post < 0.9
            out.append(dict(seed=seed, regimen=lab, scenario=name, below_target_pct=100 * fail.mean(),
                            wrongly_reassured_pct=100 * np.mean((post >= 0.5) & fail),
                            detected_pct_at_0_9=100 * (flag & fail).sum() / fail.sum(),
                            flagged_pct_at_0_9=100 * flag.mean(), flag_ratio_min=float("nan"),
                            flag_ratio_max=float("nan")))
        ratio = need["estimate"] / need["non-RRT value"]
        out.append(dict(seed=seed, regimen=lab, scenario="ratio, estimate to non-RRT value",
                        below_target_pct=float("nan"), wrongly_reassured_pct=float("nan"),
                        detected_pct_at_0_9=float("nan"), flagged_pct_at_0_9=float("nan"),
                        flag_ratio_min=float(ratio.min()), flag_ratio_max=float(ratio.max())))
    return out


def write(path, rows_):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_[0]))
        w.writeheader()
        w.writerows({k: (round(float(v), DECIMALS) if isinstance(v, (float, np.floating)) else v)
                     for k, v in r.items()} for r in rows_)


def main():
    published = {(r["regimen"], r["scenario"]): r for r in rows(os.path.join(OUT, "crrt_classifier.csv"))
                 if float(r["omega_scale"]) == OMEGA}
    kept = replay(sim.SEED, sim.REGIMENS)
    kept.update(replay(sim.EXTRA_SEED, sim.EXTRA_REGIMENS))
    thr, curve = threshold_rows(kept, published)
    write(os.path.join(OUT, "crrt_decision_threshold.csv"), thr)
    write(os.path.join(OUT, "crrt_triage_curve.csv"), curve)
    write(os.path.join(OUT, "crrt_equal_detection.csv"), equal_detection(curve, thr))
    chk = data_checks() + [nca_correlation(rows(os.path.join(DRY, "Ceftazidime_concentration.csv")))] + second_cohort()
    write(os.path.join(OUT, "crrt_data_checks.csv"), chk)
    seeds = seed_replication()
    write(os.path.join(OUT, "crrt_seed_replication.csv"), seeds)
    for r in thr:
        if r["scenario"] in ("estimate", "non-RRT value"):
            print({k: (round(float(v), 3) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    for r in chk:
        print(r)
    for r in seeds:
        print({k: (round(float(v), 3) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
