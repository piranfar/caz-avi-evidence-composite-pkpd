"""Figures for the repositioned manuscript.

Two new panels replace the therapeutic-window figure, which reported an
identity rather than a simulation result.

  fig_individualised_dose.png
    (a) the daily dose that reaches target, by renal class across the MIC range,
        against the licensed maximum
    (b) the proportion of each class whose required dose lies within the
        licensed range

  fig_second_assay.png
    misclassification of avibactam status when it is predicted from a measured
    ceftazidime concentration, against the assumed clearance correlation

Usage:
    python make_v9_figures.py
"""

from __future__ import annotations

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator

from cazavi_analyses import DEFAULT_OUT as OUT

FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
CLASSES = ["0–30", "31–60", "61–90", "91–120", "121–150"]
COLOURS = ["#1b3a6b", "#2f6f8f", "#3e8e7e", "#c1663c", "#8a6ea8"]
LICENSED = 10.0


def read(name):
    with open(os.path.join(OUT, name)) as fh:
        return list(csv.DictReader(fh))


def individualised():
    grid = read("prescriptive_decision_grid.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

    ax = axes[0]
    for cls, col in zip(CLASSES, COLOURS):
        rows = sorted((r for r in grid if r["ekfc_class"] == cls),
                      key=lambda r: float(r["mic_mg_l"]))
        x = [float(r["mic_mg_l"]) for r in rows]
        y = [float(r["median_placing_dose_g_day"]) for r in rows]
        ax.plot(x, y, color=col, lw=2.0, marker="o", ms=3.5, label=cls)
    ax.axhline(LICENSED, color="#b03030", ls="--", lw=1.2)
    ax.text(0.07, LICENSED * 1.08, "licensed maximum", fontsize=8, color="#b03030")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xticks([0.0625, 0.25, 1, 4, 16, 64])
    ax.set_xticklabels(["0.06", "0.25", "1", "4", "16", "64"])
    ax.set_yticks([0.5, 1, 2, 5, 10, 20, 50])
    ax.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "50"])
    ax.set_xlabel("MIC (mg/L)")
    ax.set_ylabel("Daily dose reaching target (g/day, product)")
    ax.set_title("(a) Median dose required, by renal function", loc="left", fontsize=11)
    ax.legend(title="EKFC class", frameon=False, fontsize=8, title_fontsize=8,
              loc="upper left")
    ax.grid(alpha=0.25, lw=0.6)

    ax = axes[1]
    width = 0.15
    mics = [4.0, 8.0, 16.0]
    idx = np.arange(len(mics))
    for k, (cls, col) in enumerate(zip(CLASSES, COLOURS)):
        vals = []
        for mic in mics:
            r = next(x for x in grid
                     if x["ekfc_class"] == cls and float(x["mic_mg_l"]) == mic)
            vals.append(float(r["placing_dose_within_licensed_pct"]))
        ax.bar(idx + (k - 2) * width, vals, width, color=col, label=cls)
    ax.set_xticks(idx)
    ax.set_xticklabels([f"MIC {m:g}" for m in mics])
    ax.set_ylabel("Required dose within licensed range (%)")
    ax.set_ylim(0, 105)
    ax.set_title("(b) Who the licensed range can reach", loc="left", fontsize=11)
    ax.legend(title="EKFC class", frameon=False, fontsize=8, title_fontsize=8,
              ncol=2, loc="upper right")
    ax.grid(alpha=0.25, lw=0.6, axis="y")

    fig.tight_layout()
    path = os.path.join(FIGDIR, "fig_individualised_dose.png")
    fig.savefig(path, dpi=300)
    print(f"  wrote {os.path.basename(path)}")


def second_assay():
    rows = sorted(read("critique_e_second_assay.csv"), key=lambda r: float(r["rho"]))
    rho = [float(r["rho"]) for r in rows]
    miss = [float(r["avibactam_status_misclassified_pct"]) for r in rows]
    cv = [float(r["avi_caz_ratio_cv_pct"]) for r in rows]

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(rho, miss, color="#1b3a6b", lw=2.2, marker="o", ms=6)
    for x, y in zip(rho, miss):
        ax.annotate(f"{y:.1f}%", (x, y), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8.5, color="#1b3a6b")
    ax.axvline(0.94, color="#b03030", ls=":", lw=1.3)
    ax.text(0.915, max(miss) * 0.55, "correlation reported\nin the source model",
            fontsize=8, color="#b03030", ha="right")
    ax2 = ax.twinx()
    ax2.plot(rho, cv, color="#c1663c", lw=1.6, ls="--", marker="s", ms=4.5)
    ax2.set_ylabel("Between-subject CV of the avibactam:ceftazidime ratio (%)",
                   color="#c1663c", fontsize=9)
    ax2.tick_params(axis="y", labelcolor="#c1663c")
    ax.set_xlabel("Assumed correlation between component clearances (ρ)")
    ax.set_ylabel("Avibactam status misclassified (%)", color="#1b3a6b")
    ax.tick_params(axis="y", labelcolor="#1b3a6b")
    ax.set_title("Value of a second assay, against the assumed clearance correlation",
                 loc="left", fontsize=10.5)
    ax.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    path = os.path.join(FIGDIR, "fig_second_assay.png")
    fig.savefig(path, dpi=300)
    print(f"  wrote {os.path.basename(path)}")


if __name__ == "__main__":
    os.makedirs(FIGDIR, exist_ok=True)
    individualised()
    second_assay()
