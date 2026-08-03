"""Figure: the same simulation under four independently published PK models.

Three panels, left to right:
  (a) the clearance models themselves, so the reader sees how far apart they are
  (b) joint attainment across the MIC range under each model, for the regimen
      selected in each renal class
  (c) the dose that reaches 90% attainment at the breakpoint versus the dose at
      which the exposure ceiling is crossed, per model and renal class

Usage:
    python make_structural_figure.py
"""

from __future__ import annotations

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator

from cazavi_analyses import DEFAULT_OUT as OUT, SELECTED_REGIMENS
from reproduce_primary_run import REGIMENS
from structural_uncertainty import CLASS_ORDER, MODELS

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
# in the repository, figures sit at the top level; standalone they sit beside the script
FIGDIR = (os.path.join(_REPO, "figures")
          if os.path.isdir(os.path.join(_REPO, "data", "processed"))
          else os.path.join(_HERE, "figures"))
os.makedirs(FIGDIR, exist_ok=True)
COLOURS = {"M1_Cojutti_2024": "#1b3a6b", "M2_Chen_2025": "#3e8e7e",
           "M3_Registrational": "#c1663c", "M4_Bensman_2017": "#8a6ea8"}
SHORT = {"M1_Cojutti_2024": "Cojutti 2024 (primary)", "M2_Chen_2025": "Chen 2025",
         "M3_Registrational": "Registrational", "M4_Bensman_2017": "Bensman 2017"}


def read(name):
    with open(os.path.join(OUT, name)) as fh:
        return list(csv.DictReader(fh))


def main():
    pta = [r for r in read("structural_uncertainty_pta.csv") if r["bsv"] == "structural"]
    cross = [r for r in read("structural_escalation_crossing.csv")
             if r["bsv"] == "structural"]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0))

    # (a) the clearance models
    ax = axes[0]
    renal = np.linspace(5, 150, 300)
    for m in MODELS:
        lo = {"M4_Bensman_2017": 61.0}.get(m.key, 5.0)
        x = renal[renal >= lo]
        ax.plot(x, m.cl_caz(x), color=COLOURS[m.key], lw=2.2, label=SHORT[m.key])
    ax.set_xlabel("Renal function (mL/min/1.73 m$^2$)")
    ax.set_ylabel("Typical ceftazidime clearance (L/h)")
    ax.set_title("(a) Four published clearance models", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.grid(alpha=0.25, lw=0.6)

    # (b) joint attainment across the MIC range
    ax = axes[1]
    for m in MODELS:
        xs, ys = [], []
        rows = [r for r in pta if r["model"] == m.key
                and r["regimen"] in SELECTED_REGIMENS]
        mics = sorted({float(r["mic_mg_l"]) for r in rows})
        for mic in mics:
            vals = [float(r["joint_pta_pct"]) for r in rows
                    if float(r["mic_mg_l"]) == mic]
            xs.append(mic)
            ys.append(np.mean(vals))
        ax.plot(xs, ys, color=COLOURS[m.key], lw=2.2, marker="o", ms=3.5,
                label=SHORT[m.key])
    ax.axhline(90, color="#888", ls="--", lw=1.0)
    ax.axvline(8, color="#b03030", ls=":", lw=1.2)
    ax.text(8.4, 5, "clinical\nbreakpoint", fontsize=7.5, color="#b03030")
    ax.text(20, 91.5, "90% target", fontsize=7.5, color="#666")
    ax.set_xscale("log", base=2)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticks([0.0625, 0.25, 1, 4, 16, 64])
    ax.set_xticklabels(["0.06", "0.25", "1", "4", "16", "64"])
    ax.set_xlabel("MIC (mg/L)")
    ax.set_ylabel("Joint target attainment (%)")
    ax.set_ylim(0, 103)
    ax.set_title("(b) Attainment, averaged over the selected regimens",
                 loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    ax.grid(alpha=0.25, lw=0.6)

    # (c) which limit is reached first
    ax = axes[2]
    width = 0.19
    idx = np.arange(len(CLASS_ORDER))
    for k, m in enumerate(MODELS):
        eff, saf, pos = [], [], []
        for j, cls in enumerate(CLASS_ORDER):
            hit = [r for r in cross if r["model"] == m.key and r["ekfc_class"] == cls]
            if not hit:
                continue
            r = hit[0]
            try:
                eff.append(float(r["daily_dose_for_90pct_pta_g"]))
                saf.append(float(r["daily_dose_at_15pct_exceedance_g"]))
                pos.append(j + (k - 1.5) * width)
            except ValueError:
                continue
        ax.bar(pos, eff, width, color=COLOURS[m.key], alpha=0.30,
               label=SHORT[m.key] if k == 0 else None)
        ax.bar(pos, saf, width, color=COLOURS[m.key])
    ax.set_xticks(idx)
    ax.set_xticklabels(CLASS_ORDER, fontsize=8.5)
    ax.set_xlabel("EKFC class (mL/min/1.73 m$^2$)")
    ax.set_ylabel("Daily dose, product (g/day)")
    ax.set_title("(c) Exposure ceiling is reached first, in every model",
                 loc="left", fontsize=11)
    solid = plt.Rectangle((0, 0), 1, 1, fc="#444")
    faint = plt.Rectangle((0, 0), 1, 1, fc="#444", alpha=0.30)
    ax.legend([solid, faint],
              ["dose at the 15% exposure ceiling",
               "dose needed for 90% attainment at MIC 8"],
              frameon=False, fontsize=8.5, loc="upper left")
    ax.grid(alpha=0.25, lw=0.6, axis="y")

    fig.tight_layout()
    path = os.path.join(FIGDIR, "fig_structural_uncertainty.png")
    fig.savefig(path, dpi=300)
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
