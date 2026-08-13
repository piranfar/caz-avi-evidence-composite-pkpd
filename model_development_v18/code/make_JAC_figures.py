"""Publication figures for the JAC draft.

Three figures, each generated from data rather than redrawn by hand:

  Figure 1  induced attainment correlation against clearance correlation, with the
            Frechet-Hoeffding bound. THE central figure -- it is the paper's argument.
  Figure 2  accuracy gain from measuring rather than inferring, against clearance
            correlation, at four assay imprecisions. Shows the gain never reaches zero.
  Figure 3  triage: accuracy against assay budget, uncertainty rule vs ratio rule.

DESIGN NOTES

Print, not screen. The house data-viz guidance assumes an interactive HTML surface with
hover layers and a dark mode; none of that applies to a journal figure, and both are
deliberately absent. What does carry over, and is applied:

  * form chosen by the data's job, not by habit
  * categorical hues assigned in fixed order (#2a78d6 blue, #eb6834 orange), validated for
    colour-vision deficiency -- worst adjacent pair dE 24.7 protan, 33.6 normal, all six
    checks pass
  * ordered variables get a SEQUENTIAL single-hue ramp, not categorical hues: assay
    imprecision in Figure 2 is ordered, so it is encoded light-to-dark in one hue
  * one y-axis per panel, never two
  * thin marks, recessive axes, no chartjunk, selective direct labels
  * identity never by colour alone -- every series also carries a distinct line style and
    marker, so the figures survive greyscale printing, which journals still do

Output: 600 dpi TIFF/PNG for review plus vector PDF for production, which is what OUP asks
for. Sized to 89 mm single-column.

Writes only into model_development_v18/.
"""

from __future__ import annotations

import csv
import os
import sys
from dataclasses import replace

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import model2_engine as E          # noqa: E402
import model2_hujam as H           # noqa: E402
import model2_monitoring as M      # noqa: E402
import reproduce_primary_run as P  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "manuscript_JAC", "figures")
DATA = os.path.join(HERE, "..", "outputs")

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, MUTED = "#0b0b0b", "#52514e"
MM = 1.0 / 25.4
SINGLE_COL = 89 * MM          # OUP single column

# Sequential ramp for the ORDERED assay-imprecision variable (light -> dark, one hue).
SEQ = ["#a8c8ee", "#6ea3e3", "#2a78d6", "#14508f"]


def house_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7.5,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
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


def save(fig, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT_DIR, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  wrote {name}.pdf and {name}.png")


def read_csv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# --------------------------------------------------------------------------------------
# Figure 1 -- the bound
# --------------------------------------------------------------------------------------
def figure1():
    rows = read_csv(os.path.join(DATA, "dispute_boundary_fresan_gatti.csv"))
    rho = np.array([float(r["clearance_rho"]) for r in rows])
    phi = np.array([float(r["induced_attainment_phi"]) for r in rows])
    bound = np.array([float(r["frechet_bound_phi"]) for r in rows])

    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.78))

    # Drawn per rho from the source data rather than as a single mean line, so each point
    # is compared with its own bound. Since model2_dispute_boundary.py adopted common
    # random numbers the prevalences no longer drift across the grid, so this line is flat
    # to within 0.006 -- as it must be, because rho changes the joint distribution and
    # neither margin. An earlier version used the mean bound and put two points above a
    # line labelled "unattainable", which was exactly the wrong impression.
    ax.fill_between(rho, bound, 1.0, color=MUTED, alpha=0.08, lw=0, zorder=0)
    ax.plot(rho, bound, color=ORANGE, lw=1.4, ls=(0, (5, 2)), zorder=2)
    ax.plot(rho, phi, color=BLUE, lw=1.6, marker="o", ms=3.6,
            mfc="white", mew=1.2, mec=BLUE, zorder=3, clip_on=False)

    # identity line as the reference the reader expects and does not get
    ax.plot([0.28, 1.0], [0.28, 1.0], color=MUTED, lw=0.7, ls=":", zorder=1)

    ax.text(0.315, 0.95, "unattainable for any\njoint distribution",
            fontsize=6.2, color=MUTED, va="top", linespacing=1.35)
    ax.text(0.315, 0.565, "Fréchet–Hoeffding bound",
            fontsize=6.4, color=ORANGE, ha="left", va="bottom")
    ax.text(0.615, 0.745, "if attainment tracked\nclearance 1:1",
            fontsize=6.2, color=MUTED, ha="center", va="center",
            rotation=36, rotation_mode="anchor", linespacing=1.3)

    # the headline point
    ax.annotate(f"$\\rho$ = 0.99 reaches\n$\\varphi$ = {phi[-1]:.2f}, its ceiling",
                xy=(rho[-1], phi[-1]), xytext=(0.90, 0.245),
                fontsize=6.4, color=INK, ha="center", linespacing=1.3,
                arrowprops=dict(arrowstyle="-", lw=0.7, color=MUTED,
                                shrinkA=2, shrinkB=3))

    ax.set_xlabel("Between-subject clearance correlation, $\\rho$")
    ax.set_ylabel("Attainment correlation, $\\varphi$")
    ax.set_xlim(0.28, 1.0)
    ax.set_ylim(0, 1.0)
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    save(fig, "Figure1_attainment_correlation_bound")


# --------------------------------------------------------------------------------------
# Figure 2 -- the gain never closes
# --------------------------------------------------------------------------------------
def figure2(n_per_class=2500, n_draws=120):
    """Gain from measuring, across rho, at four assay imprecisions.

    COMMON RANDOM NUMBERS, as in model2_dispute_boundary.py and for the same reason. The
    expected ordering is monotone -- a noisier assay can only reduce the advantage of
    measuring -- so curves that cross are reporting Monte Carlo noise as if it were signal.
    Re-seeding identically at every (CV, rho) point gives the whole grid the same parameter
    draws and the same assay-noise realisations, leaving only the effects being plotted.
    """
    grid = [0.30, 0.45, 0.60, 0.703, 0.80, 0.90, 0.94, 0.97, 0.99]
    cvs = [0.0, 0.10, 0.20, 0.30]
    pop = E.draw_population(n_per_class, H.MASTER_SEED)

    curves = {}
    for cv in cvs:
        gains = []
        for rho in grid:
            rng = np.random.default_rng(H.MASTER_SEED + 909)
            vals = []
            for _ in range(n_draws):
                pr = E.draw_parameters(rng)
                pr = replace(pr, rho=rho, avi_target=P.AVI_CT)
                m = M.classify(pop, pr, rng, assay_cv_caz=cv, assay_cv_avi=cv)
                vals.append(m["accuracy_measure"] - m["accuracy_infer"])
            gains.append(float(np.median(vals)))
        curves[cv] = gains
        print(f"    CV {cv:.0%}: gain {min(gains):.2f} to {max(gains):.2f} pp")

    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.78))
    ax.axhline(0, color=INK, lw=0.8, zorder=2)
    ax.text(0.30, 0.42, "inference catches up", fontsize=6.2, color=MUTED,
            ha="left", va="bottom")

    styles = ["-", (0, (5, 1.5)), (0, (3, 1.2)), (0, (1.2, 1.2))]
    marks = ["o", "s", "^", "D"]
    for (cv, g), col, ls, mk in zip(curves.items(), SEQ, styles, marks):
        ax.plot(grid, g, color=col, lw=1.5, ls=ls, marker=mk, ms=3.2,
                mfc="white", mew=1.0, mec=col, label=f"{cv:.0%}", clip_on=False)

    ax.set_xlabel("Between-subject clearance correlation, $\\rho$")
    ax.set_ylabel("Accuracy gain from measuring (pp)")
    ax.set_xlim(0.28, 1.0)
    ax.set_ylim(bottom=-1.0)
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    leg = ax.legend(title="Assay CV", frameon=False, loc="upper right",
                    handlelength=2.4, borderpad=0.2, labelspacing=0.35)
    leg.get_title().set_fontsize(6.5)
    save(fig, "Figure2_gain_from_measuring")


# --------------------------------------------------------------------------------------
# Figure 3 -- triage
# --------------------------------------------------------------------------------------
def figure3():
    """Two panels, because the comparison is policy-dependent and saying so is the honest
    reading. The ratio-based rule was derived in a cohort largely on a FIXED dose; under
    this project's renally-adjusted grid it has far less to discriminate on. Plotting only
    one panel would either flatter or handicap it depending which was chosen.
    """
    all_rows = read_csv(os.path.join(DATA, "model2_triage_vs_gatti.csv"))
    panels = [("renally-adjusted", "Renally-adjusted dosing"),
              ("FIXED", "Fixed 2.5 g q8h dosing")]

    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL * 2.0, SINGLE_COL * 0.80),
                             sharey=True)
    for ax, (key, title), tag in zip(axes, panels, ["(a)", "(b)"]):
        rows = [r for r in all_rows
                if r["dosing_policy"].startswith(key)
                and r["rho_scenario"].startswith("published")]
        x = np.array([float(r["pct_measured"]) for r in rows])
        r2 = np.array([float(r["accuracy_R2_rule_pct"]) for r in rows])
        ga = np.array([float(r["accuracy_Gatti_rule_pct"]) for r in rows])
        o = np.argsort(x)
        x, r2, ga = x[o], r2[o], ga[o]

        ax.fill_between(x, ga, r2, where=(r2 >= ga), color=BLUE, alpha=0.10,
                        lw=0, zorder=1, interpolate=True)
        ax.plot(x, r2, color=BLUE, lw=1.6, ls="-", marker="o", ms=3.2,
                mfc="white", mew=1.1, mec=BLUE, label="Uncertainty-based rule", zorder=3)
        ax.plot(x, ga, color=ORANGE, lw=1.6, ls=(0, (5, 2)), marker="s", ms=3.2,
                mfc="white", mew=1.1, mec=ORANGE, label="Ratio-based rule", zorder=3)
        ax.set_title(f"{tag}  {title}", fontsize=7, loc="left", color=INK, pad=4)
        ax.set_xlabel("Patients assayed for avibactam (%)")
        ax.set_xlim(0, 100)
        ax.xaxis.set_major_locator(MultipleLocator(25))

    k = 3   # the 12.5% budget
    rows0 = [r for r in all_rows if r["dosing_policy"].startswith("renally-adjusted")
             and r["rho_scenario"].startswith("published")]
    rows0.sort(key=lambda r: float(r["pct_measured"]))
    axes[0].annotate(f"{float(rows0[k]['accuracy_R2_rule_pct']):.1f}% at a "
                     f"{float(rows0[k]['pct_measured']):.1f}% budget",
                     xy=(float(rows0[k]["pct_measured"]),
                         float(rows0[k]["accuracy_R2_rule_pct"])),
                     xytext=(34, 93.4), fontsize=6.4, color=INK, ha="left",
                     arrowprops=dict(arrowstyle="-", lw=0.7, color=MUTED,
                                     shrinkA=0, shrinkB=3))

    axes[0].set_ylabel("Correct classification (%)")
    axes[1].legend(frameon=False, loc="lower right", handlelength=2.4,
                   borderpad=0.2, labelspacing=0.35)
    fig.subplots_adjust(wspace=0.08)
    save(fig, "Figure3_triage_budget")


def main():
    house_style()
    print("Figure 1 — the bound")
    figure1()
    print("Figure 2 — gain from measuring")
    figure2()
    print("Figure 3 — triage")
    figure3()
    print(f"\nAll figures in {os.path.relpath(OUT_DIR, HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
