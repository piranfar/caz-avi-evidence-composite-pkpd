"""Two multi-panel greyscale figures for the CRRT analysis.

Figure 1 (joint model): A observed vs individual prediction; B CWRES vs time;
                        C profile likelihood of rho; D circuit and non-circuit clearance of each patient,
                          ceftazidime against avibactam (outputs/crrt_extracorporeal_components.csv).
Figure 2 (virtual CRRT population): A free avibactam trough by regimen with observed troughs;
                        B joint attainment vs MIC; C percentage wrongly reassured, at the estimated
                          correlation and at 0.94;
                        D patients below the avibactam target identified vs patients flagged,
                          as the decision threshold varies (outputs/crrt_triage_curve.csv).
Writes PDF (vector) and TIFF (1200 dpi) to figures/.
"""
from __future__ import annotations

import csv
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import crrt_virtual_tdm as sim

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "outputs")
FIG = os.path.join(ROOT, "figures")
DRY = os.path.join(ROOT, "data_external", "dryad_Li2025_CRRT")
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6})
DARK, MID, LIGHT = "#000000", "#707070", "#bdbdbd"


def rows(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def label(ax, s):
    ax.text(-0.18, 1.04, s, transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom")


def save(fig, name):
    """Vector PDF, and a 1,200-dpi LZW-compressed RGB TIFF (no alpha channel)."""
    from PIL import Image
    fig.savefig(os.path.join(FIG, f"{name}.pdf"))
    tif = os.path.join(FIG, f"{name}.tif")
    fig.savefig(tif, dpi=1200, pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tif) as im:
        rgb = im.convert("RGB")
    rgb.save(tif, compression="tiff_lzw", dpi=(1200, 1200))


def figure1():
    dg = rows(os.path.join(OUT, "model1_diagnostics.csv"))
    fig, axs = plt.subplots(2, 2, figsize=(6.8, 5.6))
    a, b, c, d = axs.ravel()
    style = {"caz": dict(marker="o", mfc=DARK, mec=DARK, label="Ceftazidime"),
             "avi": dict(marker="^", mfc="white", mec=MID, label="Avibactam")}
    for an, st in style.items():
        r = [x for x in dg if x["analyte"] == an]
        a.plot([float(x["ipred_mg_l"]) for x in r], [float(x["dv_mg_l"]) for x in r], ls="none",
               ms=3.2, mew=0.6, **st)
        b.plot([float(x["time_h"]) + (0.12 if an == "avi" else -0.12) for x in r],
               [float(x["cwres"]) for x in r], ls="none", ms=3.2, mew=0.6, **st)
    a.set_xscale("log"); a.set_yscale("log")
    lim = (5, 300)
    a.plot(lim, lim, color=MID, lw=0.7); a.set_xlim(lim); a.set_ylim(lim)
    a.set_xlabel("Individual prediction (mg/L)"); a.set_ylabel("Observed concentration (mg/L)")
    a.legend(frameon=False, loc="upper left")
    for y, ls in ((0, "-"), (2, "--"), (-2, "--")):
        b.axhline(y, color=MID, lw=0.7, ls=ls)
    b.set_xlabel("Time after dose (h)"); b.set_ylabel("Conditional weighted residual")
    b.set_ylim(-3.5, 3.5)

    pl = rows(os.path.join(OUT, "model1_profile_likelihood.csv"))
    r_ = np.array([float(x["rho"]) for x in pl]); dofv = np.array([float(x["delta_ofv"]) for x in pl])
    top = float(np.ceil(1.1 * dofv[r_ <= 0.96].max() / 2.0) * 2.0)
    c.axvspan(sim.RHO_LO, sim.RHO_HI, color=LIGHT, alpha=0.5, lw=0)
    c.plot(r_, dofv, "o-", color=DARK, ms=3, lw=0.9, label="Primary model")
    kb = sorted((float(x["rho"]), float(x["delta_ofv"])) for x in rows(os.path.join(OUT, "model1_covariance_profile.csv"))
                if x["structure"] == "Kb")
    c.plot([p[0] for p in kb], [p[1] for p in kb], "s--", color=MID, ms=2.8, lw=0.8, mfc="white",
           label="All 6 correlations")
    c.legend(frameon=False, loc="upper left", fontsize=6.5)
    c.axhline(3.84, color=MID, ls="--", lw=0.7)
    c.axvline(0.94, color=DARK, ls=":", lw=0.9)
    c.text(0.925, 0.92 * top, "0.94", ha="right", fontsize=7)
    c.set_xlim(0, 1); c.set_ylim(0, top)
    c.set_xlabel("Clearance correlation (ρ)"); c.set_ylabel("ΔOFV")

    comp = [x for x in rows(os.path.join(OUT, "crrt_extracorporeal_components.csv")) if float(x["quf_l_h"]) == 0.0]
    summ = {x["quantity"]: x for x in rows(os.path.join(OUT, "crrt_extracorporeal_summary.csv"))
            if x["note"].startswith("net ultrafiltration 0.0 L/h")}
    pts = {k: (np.array([float(x[f"{k}_caz_l_h"]) for x in comp]), np.array([float(x[f"{k}_avi_l_h"]) for x in comp]))
           for k in ("cl_crrt", "cl_noncircuit")}
    r_c = float(summ["corr_cl_crrt_caz_avi"]["value"])
    r_n = float(summ["corr_cl_noncircuit_caz_avi"]["value"])
    assert abs(r_c - np.corrcoef(*pts["cl_crrt"])[0, 1]) < 1e-6 and abs(r_n - np.corrcoef(*pts["cl_noncircuit"])[0, 1]) < 1e-6
    assert len(comp) == 21 and min(pts["cl_noncircuit"][0].min(), pts["cl_noncircuit"][1].min()) > 0
    d.plot(*pts["cl_crrt"], "o", mfc=DARK, mec=DARK, ms=4, label=f"Circuit (r = {r_c:.2f})")
    d.plot(*pts["cl_noncircuit"], "^", mfc="white", mec=MID, ms=4, mew=0.8, label=f"Outside the circuit (r = {r_n:.2f})")
    lo, hi = 0.0, 3.0
    d.plot((lo, hi), (lo, hi), color=MID, lw=0.7, ls="--")
    d.set_xlim(lo, hi); d.set_ylim(lo, hi)
    d.set_xlabel("Ceftazidime clearance (L/h)"); d.set_ylabel("Avibactam clearance (L/h)")
    d.legend(frameon=False, loc="upper left", fontsize=6.5)
    for ax, s in zip((a, b, c, d), "ABCD"):
        label(ax, s)
    fig.tight_layout()
    save(fig, "TDM_Figure1_joint_model")


ORDER = ["2.5 g q8h", "2.5 g q12h", "1.25 g q8h", "0.94 g q12h"]   # by daily dose


def simulate_regimens():
    """Virtual populations at omega x1.0, drawn exactly as in crrt_virtual_tdm.main():
    the original regimens from the SEED stream in REGIMENS order, the added regimen from EXTRA_SEED.
    Returns the SEED stream (used afterwards only for jitter) and the exposures by regimen."""
    rng = np.random.default_rng(sim.SEED)
    out = {}
    for lab_, cz, av, tau in sim.REGIMENS:
        out[lab_] = sim.exposures(sim.population(rng, sim.N, sim.RHO_EST), cz, av, tau)
    rng2 = np.random.default_rng(sim.EXTRA_SEED)
    for lab_, cz, av, tau in sim.EXTRA_REGIMENS:
        out[lab_] = sim.exposures(sim.population(rng2, sim.N, sim.RHO_EST), cz, av, tau)
    att = {x["regimen"]: x for x in rows(os.path.join(OUT, "crrt_attainment.csv")) if float(x["omega_scale"]) == 1.0}
    for lab_, (cmin_caz, cmin_avi, _) in out.items():
        # The plotted populations must be the ones behind the reported numbers.
        assert abs(np.median(sim.FU_AVI * cmin_avi) - float(att[lab_]["avi_fcmin_median"])) < 1e-5, lab_
        assert abs(100 * np.mean(sim.FU_AVI * cmin_avi >= sim.AVI_CT) - float(att[lab_]["avi_attain_fcmin4"])) < 1e-5, lab_
    return rng, out


def figure2():
    rng, sims = simulate_regimens()
    fig = plt.figure(figsize=(6.8, 5.9))
    gs = fig.add_gridspec(2, 2)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 0]); d = fig.add_subplot(gs[1, 1])

    data = []
    mics = np.array([0.5, 1, 2, 4, 8, 16, 32])
    joint = {}
    for lab_ in ORDER:
        cmin_caz, cmin_avi, _ = sims[lab_]
        f_avi = sim.FU_AVI * cmin_avi
        data.append(f_avi)
        joint[lab_] = [100 * np.mean((sim.FU_CAZ * cmin_caz >= sim.CAZ_K * m) & (f_avi >= sim.AVI_CT))
                       for m in mics]
    bp = a.boxplot(data, widths=0.45, whis=(5, 95), showfliers=False, patch_artist=True)
    for box in bp["boxes"]:
        box.set(facecolor=LIGHT, edgecolor=DARK, linewidth=0.7)
    for k in ("whiskers", "caps", "medians"):
        for ln in bp[k]:
            ln.set(color=DARK, linewidth=0.7)
    obs = [sim.FU_AVI * float(x["avi_pre"]) for x in rows(os.path.join(DRY, "Avibactam_concentration.csv"))
           if float(x["Time_h"]) == 0]
    a.plot(1.38 + 0.03 * rng.standard_normal(len(obs)), obs, "o", mfc="white", mec=DARK, ms=3.2, mew=0.6,
           label="Observed (n = 21)")
    a.axhline(4, color=DARK, ls="--", lw=0.7); a.axhline(1, color=MID, ls=":", lw=0.8)
    a.text(4.45, 4.3, "4 mg/L", ha="right", fontsize=7); a.text(4.45, 1.07, "1 mg/L", ha="right", fontsize=7, color=MID)
    a.set_yscale("log"); a.set_ylim(0.3, 40)
    a.set_xticks([1, 2, 3, 4], [s.replace(" q", "\nq") for s in ORDER])
    a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    a.set_ylabel("Free avibactam trough (mg/L)")
    a.legend(frameon=False, loc="lower left", fontsize=7)

    for lab_, ls, mk in zip(ORDER, ("-", "-.", "--", ":"), ("o", "D", "s", "^")):
        b.plot(mics, joint[lab_], ls=ls, marker=mk, color=DARK, ms=3.5, lw=0.9, mfc="white", label=lab_)
    b.axhline(90, color=LIGHT, lw=0.8, ls="-")
    b.set_xscale("log", base=2); b.set_xticks(mics, [f"{m:g}" for m in mics])
    b.set_ylim(-3, 103); b.set_xlabel("Ceftazidime-avibactam MIC (mg/L)"); b.set_ylabel("Joint attainment (%)")
    b.legend(frameon=False, loc=(0.02, 0.2), fontsize=7)

    cls = rows(os.path.join(OUT, "crrt_classifier.csv"))
    def val(reg, scen, key):
        return next(float(x[key]) for x in cls if x["regimen"] == reg and x["scenario"] == scen
                    and float(x["omega_scale"]) == 1.0)
    crrt = ORDER[:3]
    groups = [g.replace(" q", "\nq") for g in crrt]
    nomeas = [val(g, "estimate", "missed_by_population_prior") for g in crrt]
    est = [val(g, "estimate", "false_reassurance") for g in crrt]
    r94 = [val(g, "non-RRT value", "false_reassurance") for g in crrt]
    # Range across the 95% CI of rho: min and max over the estimate and both bounds (not always monotone).
    span = [[val(g, s, "false_reassurance") for s in ("estimate", "upper bound", "lower bound")] for g in crrt]
    lo = [min(v) for v in span]
    hi = [max(v) for v in span]
    x = np.arange(len(crrt))
    w = 0.26
    c.bar(x - w, nomeas, w, color=LIGHT, edgecolor=DARK, lw=0.6, label="No measurement")
    c.bar(x, est, w, color=DARK, label=f"Ceftazidime trough, ρ = {sim.RHO_EST:.2f}")
    c.bar(x + w, r94, w, color="white", edgecolor=DARK, lw=0.6, hatch="////", label="Ceftazidime trough, ρ = 0.94")
    c.errorbar(x, est, yerr=[np.array(est) - lo, np.array(hi) - est], fmt="none", ecolor=MID, capsize=2.5, lw=0.8)
    for xs, vals, tops in ((x - w, nomeas, nomeas), (x, est, hi), (x + w, r94, r94)):
        for xi, v, t in zip(xs, vals, tops):
            c.text(xi, t + 0.25, "<0.1" if v < 0.1 else f"{v:.1f}", ha="center", fontsize=6)
    c.set_xticks(x, groups, fontsize=7); c.set_ylim(0, 12.5)
    c.set_ylabel("Wrongly reassured (% of patients)")
    c.legend(frameon=False, loc="upper left", fontsize=6.3)

    # D: what a stricter decision threshold buys, at the estimated and at the non-RRT correlation.
    tc = rows(os.path.join(OUT, "crrt_triage_curve.csv"))
    thr = rows(os.path.join(OUT, "crrt_decision_threshold.csv"))
    d.plot([0, 100], [0, 100], color=LIGHT, lw=0.8, ls="-")
    for reg, col in (("1.25 g q8h", DARK), ("2.5 g q12h", MID)):
        for scen, ls, rho_lab in (("estimate", "-", f"{sim.RHO_EST:.2f}"), ("non-RRT value", "--", "0.94")):
            pts = sorted((float(r["flagged_pct_of_all"]), float(r["detected_pct_of_below_target"]))
                         for r in tc if r["regimen"] == reg and r["scenario"] == scen)
            d.plot([p[0] for p in pts], [p[1] for p in pts], color=col, ls=ls, lw=1.0,
                   label=f"{reg}, ρ = {rho_lab}")
            for t, mk in (("0.5", "o"), ("0.9", "s")):
                r = next(r for r in thr if r["regimen"] == reg and r["scenario"] == scen and r["threshold"] == t)
                d.plot(float(r["flagged_pct_of_all"]), float(r["detected_pct_of_below_target"]), mk,
                       mfc="white", mec=col, ms=3.6, mew=0.8)
    d.plot([], [], "o", mfc="white", mec=DARK, ms=3.6, mew=0.8, label="Threshold 0.5")
    d.plot([], [], "s", mfc="white", mec=DARK, ms=3.6, mew=0.8, label="Threshold 0.9")
    d.set_xlim(0, 40); d.set_ylim(0, 100)
    d.set_xlabel("Patients flagged (% of all)")
    d.set_ylabel("Patients below avibactam target\nidentified (%)")
    d.legend(frameon=False, loc="lower right", fontsize=6.2)
    for ax, s in zip((a, b, c, d), "ABCD"):
        label(ax, s)
    fig.tight_layout()
    save(fig, "TDM_Figure2_virtual_crrt")


if __name__ == "__main__":
    figure1()
    figure2()
