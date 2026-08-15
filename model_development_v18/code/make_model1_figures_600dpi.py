"""Publication renders of the joint-model diagnostics, at 600 dpi.

WHY THIS EXISTS SEPARATELY FROM model1_finalise.py

`model1_finalise.py` writes `model1_gof.png` and `model1_vpc.png` at 300 dpi with an in-figure
title, which is right for a working diagnostic and wrong for a journal. Elsevier wants 600 dpi
for combination art, no title inside the image (the caption carries it), and a vector copy for
production. Re-running the fit to change a font size would also be wasteful and would put the
estimates at risk of drifting.

So this script re-renders from the frozen CSV outputs of that fit — `model1_diagnostics.csv` and
`model1_vpc.csv` — and touches nothing else. The numbers cannot change, because no model is
refitted here.

Writes Figure 7 (goodness of fit) and Figure 8 (visual predictive check) as 600 dpi PNG plus
vector PDF, into model_development_v18/figures/.
"""

from __future__ import annotations

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "outputs")
FIG = os.path.join(HERE, "..", "figures")

INK, MUTED = "#1A1A1A", "#6B6B6B"
COLOURS = {"caz": "#1F4E85", "avi": "#C86438"}
LABELS = {"caz": "Ceftazidime", "avi": "Avibactam"}
ANALYTES = ("caz", "avi")


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": MUTED,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "figure.dpi": 600,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def read(name):
    with open(os.path.join(OUT, name), encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  wrote {name}.png and {name}.pdf")


def log_ticks(lim, want=(3, 6)):
    """Round tick values inside `lim`, using the coarsest ladder that still gives enough ticks.

    Greedy thinning of a fine ladder anchors on whichever value happens to fall lowest in the
    range and can strand the axis on odd labels, so choose the ladder instead of filtering one.
    """
    for mantissas in ((1, 2, 5), (1, 1.5, 2, 3, 5, 7), (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8)):
        ticks = sorted(m * 10 ** d for d in range(-1, 4) for m in mantissas
                       if lim[0] <= m * 10 ** d <= lim[1])
        if len(ticks) >= want[0]:
            return [int(t) if t == int(t) else t for t in ticks[:want[1]]]
    return [int(t) if t == int(t) else t for t in ticks]


def figure7(diag):
    """Observed vs population and individual predictions, and CWRES against time."""
    # 190 mm is Elsevier's double-column width; 7.48 in is the same in inches.
    fig, ax = plt.subplots(2, 3, figsize=(7.48, 4.7))
    panel = iter("ABCDEF")
    for r, an in enumerate(ANALYTES):
        d = [x for x in diag if x["analyte"] == an]
        dv = np.array([float(x["dv_mg_l"]) for x in d])
        pr = np.array([float(x["pred_mg_l"]) for x in d])
        ip = np.array([float(x["ipred_mg_l"]) for x in d])
        cw = np.array([float(x["cwres"]) for x in d])
        tm = np.array([float(x["time_h"]) for x in d])

        for c, (xv, xl) in enumerate(((pr, "Population prediction (mg/L)"),
                                      (ip, "Individual prediction (mg/L)"))):
            a = ax[r, c]
            lim = [min(dv.min(), xv.min()) * 0.8, max(dv.max(), xv.max()) * 1.2]
            a.plot(lim, lim, color="0.45", lw=0.8, zorder=1)
            a.scatter(xv, dv, s=9, alpha=0.7, color=COLOURS[an], edgecolor="none", zorder=2)
            a.set_xscale("log"); a.set_yscale("log")
            a.set_xlim(lim); a.set_ylim(lim)
            # Default log ticks collide on a narrow axis ("3 x 10^1 4 x 10^1"); pick round
            # concentrations inside the range and label them as plain numbers.
            ticks = log_ticks(lim)
            for setter, fmt in ((a.set_xticks, a.set_xticklabels), (a.set_yticks, a.set_yticklabels)):
                setter(ticks); fmt([str(t) for t in ticks])
            a.minorticks_off()
            a.set_xlabel(xl); a.set_ylabel("Observed (mg/L)")
            a.set_title(f"({next(panel)}) {LABELS[an]}", loc="left")

        a = ax[r, 2]
        a.axhline(0, color="0.45", lw=0.8)
        for h in (-2, 2):
            a.axhline(h, color="0.75", lw=0.6, ls="--")
        a.scatter(tm, cw, s=9, alpha=0.7, color=COLOURS[an], edgecolor="none")
        a.set_xlabel("Time after dose (h)")
        a.set_ylabel("Conditional weighted residual")
        a.set_title(f"({next(panel)}) {LABELS[an]}, SD {cw.std(ddof=1):.2f}", loc="left")
    fig.tight_layout(pad=0.4)
    save(fig, "Figure7_model1_goodness_of_fit")


def figure8(vpc):
    """Observed versus simulated concentration-time profiles."""
    fig, ax = plt.subplots(1, 2, figsize=(7.48, 2.9))
    for i, an in enumerate(ANALYTES):
        v = [x for x in vpc if x["analyte"] == an]
        t = np.array([float(x["time_h"]) for x in v])
        g = lambda k: [float(x[k]) for x in v]  # noqa: E731
        a = ax[i]
        a.fill_between(t, g("sim_p5"), g("sim_p95"), color=COLOURS[an], alpha=0.18, lw=0,
                       label="simulated 5th–95th percentile")
        a.plot(t, g("sim_p50"), color=COLOURS[an], lw=1.4, label="simulated median")
        a.plot(t, g("obs_p50"), "o--", color="0.15", ms=3, lw=0.9, label="observed median")
        a.plot(t, g("obs_p5"), ".", color="0.45", ms=4)
        a.plot(t, g("obs_p95"), ".", color="0.45", ms=4, label="observed 5th and 95th")
        a.set_xlabel("Time after dose (h)")
        a.set_ylabel("Concentration (mg/L)")
        a.set_title(f"({'AB'[i]}) {LABELS[an]}", loc="left")
        a.legend(frameon=False, loc="upper right", handlelength=1.6)
    fig.tight_layout(pad=0.4)
    save(fig, "Figure8_model1_visual_predictive_check")


def main() -> int:
    style()
    diag, vpc = read("model1_diagnostics.csv"), read("model1_vpc.csv")
    print(f"diagnostics rows {len(diag)}  |  vpc rows {len(vpc)}")
    for an in ANALYTES:
        n = sum(1 for x in diag if x["analyte"] == an)
        cw = np.array([float(x["cwres"]) for x in diag if x["analyte"] == an])
        print(f"  {LABELS[an]:12s} n={n:4d}  CWRES mean {cw.mean():+.3f}  SD {cw.std(ddof=1):.3f}")
    cov = [float(x["pct_obs_within_sim_90"]) for x in vpc if x.get("pct_obs_within_sim_90")]
    if cov:
        print(f"  VPC coverage of the nominal 90% interval: {np.mean(cov):.1f}%")
    figure7(diag)
    figure8(vpc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
