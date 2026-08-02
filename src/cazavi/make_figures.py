"""Manuscript figures, built from the verified analysis outputs.

Reads only `outputs/*.csv` produced by `cazavi_analyses.py`, so every figure
traces back to a run that has been checked against the frozen RC1 tables.
"""

from __future__ import annotations

import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))


def _paths():
    """Locate the analysis outputs and the folder figures belong in."""
    repo = os.path.dirname(os.path.dirname(HERE))
    if os.path.isdir(os.path.join(repo, "data", "processed")):
        return os.path.join(repo, "data", "processed"), os.path.join(repo, "figures")
    return os.path.join(HERE, "outputs"), os.path.join(HERE, "figures")


OUT, FIG = _paths()
os.makedirs(FIG, exist_ok=True)

CAZ, AVI, JOINT, WARN = "#2E6FB7", "#C1443C", "#1B1B1B", "#E08A1E"
SELECTED = [("R1", "0–30"), ("R8", "31–60"), ("R10", "61–90"),
            ("R12", "91–120"), ("R13", "121–150")]
BREAKPOINT = 8.0     # EUCAST clinical breakpoint for ceftazidime-avibactam


def read(name):
    with open(os.path.join(OUT, name)) as fh:
        return list(csv.DictReader(fh))


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  {os.path.relpath(path, os.path.dirname(os.path.dirname(HERE)))}")


def fig_pta_vs_mic():
    """Where each component becomes limiting, across the whole MIC range."""
    rows = read("primary_pta_results.csv")
    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by[r["regimen"]]["mic"].append(float(r["mic_mg_l"]))
        for k, f in (("caz", "caz_pta_pct"), ("avi", "avi_attainment_pct"),
                     ("joint", "joint_pta_pct")):
            by[r["regimen"]][k].append(float(r[f]))

    fig, axes = plt.subplots(1, 5, figsize=(19, 4.1), sharey=True)
    for ax, (reg, cls) in zip(axes, SELECTED):
        d = by[reg]
        ax.plot(d["mic"], d["caz"], "o-", color=CAZ, lw=2, ms=4, label="Ceftazidime")
        ax.plot(d["mic"], d["avi"], "s-", color=AVI, lw=2, ms=4, label="Avibactam")
        ax.plot(d["mic"], d["joint"], "^-", color=JOINT, lw=2.4, ms=5, label="Joint")
        ax.axvline(BREAKPOINT, color="#888", ls="--", lw=1.2)
        ax.axhline(90, color="#999", ls=":", lw=1)
        ax.set_xscale("log", base=2)
        ax.set_xlim(0.0625, 64)
        ax.set_ylim(0, 102)
        ax.set_title(f"EKFC {cls}\n{reg}", fontsize=10)
        ax.set_xlabel("MIC (mg/L)", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=.25, lw=.5)
    axes[0].set_ylabel("Target attainment (%)", fontsize=10)
    axes[0].legend(fontsize=8, loc="lower left")
    axes[2].text(BREAKPOINT * 1.15, 4, "EUCAST\nbreakpoint", fontsize=8, color="#555")
    fig.suptitle("Avibactam caps attainment below the breakpoint; ceftazidime fails at it",
                 fontsize=11.5, y=1.03)
    save(fig, "fig_pta_vs_mic.png")


def fig_oat_tornado():
    """Which assumption moves the answer most."""
    rows = [r for r in read("gsa_sensitivity_ranking.csv") if r["scenario"] != "BASE"]
    rows.sort(key=lambda r: float(r["max_abs_delta_joint_cfr_pp"]))
    labels = [r["scenario"] for r in rows]
    vals = [float(r["max_abs_delta_joint_cfr_pp"]) for r in rows]
    policy = {"AVI_CT_2", "AVI_CT_6", "CAZ_TARGET_2", "CAZ_TARGET_6", "TOX_80", "TOX_128"}
    cols = [AVI if l in policy else CAZ for l in labels]

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    ax.barh(range(len(vals)), vals, color=cols, height=.68)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Maximum |Δ joint CFR| versus baseline (percentage points)", fontsize=10)
    for i, v in enumerate(vals):
        ax.text(v + .25, i, f"{v:.2f}", va="center", fontsize=8.5)
    ax.set_xlim(0, max(vals) * 1.16)
    ax.grid(axis="x", alpha=.25, lw=.5)
    ax.set_axisbelow(True)
    ax.legend(handles=[Patch(color=AVI, label="Analyst-specified target or threshold"),
                       Patch(color=CAZ, label="Pharmacokinetic parameter")],
              fontsize=8.5, loc="lower right")
    ax.set_title("The avibactam critical concentration outweighs every PK parameter",
                 fontsize=11.5)
    save(fig, "fig_oat_tornado.png")


def fig_dose_response():
    """Joint CFR against daily dose, with the toxicity ceiling overlaid."""
    rows = [r for r in read("cfr_all_distributions.csv")
            if r["distribution_id"] == "LEE2022_KPC_KP"]
    by_class = defaultdict(list)
    for r in rows:
        daily = float(r["dose_g"]) * 24.0 / float(r["interval_h"])
        by_class[r["ekfc_class"]].append(
            (daily, float(r["joint_cfr_pct"]), float(r["toxicity_pct"]), r["regimen"]))

    order = ["0–30", "31–60", "61–90", "91–120", "121–150"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
    cmap = plt.get_cmap("viridis")
    for i, cls in enumerate(order):
        pts = sorted(by_class.get(cls, []))
        if not pts:
            continue
        c = cmap(i / max(len(order) - 1, 1))
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=c, lw=2, ms=6,
                label=f"EKFC {cls}")
        ax2.plot([p[0] for p in pts], [p[2] for p in pts], "o-", color=c, lw=2, ms=6,
                 label=f"EKFC {cls}")
        for d, j, t, reg in pts:
            ax.annotate(reg, (d, j), textcoords="offset points", xytext=(0, 7),
                        fontsize=7.5, ha="center", color="#444")

    ax.axhline(80, color="#999", ls=":", lw=1.2)
    ax.text(10.3, 80.6, "80% CFR", fontsize=8, color="#666")
    ax.set_xlabel("Total daily product dose (g)", fontsize=10)
    ax.set_ylabel("Joint CFR (%), KPC-KP distribution", fontsize=10)
    ax.set_title("Efficacy rises with daily dose", fontsize=11)

    ax2.axhline(15, color=AVI, ls="--", lw=1.4)
    ax2.text(10.3, 15.6, "15% ceiling", fontsize=8, color=AVI)
    ax2.set_xlabel("Total daily product dose (g)", fontsize=10)
    ax2.set_ylabel("Exposure-screen exceedance (%)", fontsize=10)
    ax2.set_title("...and so does the ceftazidime exposure screen", fontsize=11)

    for a in (ax, ax2):
        a.grid(alpha=.25, lw=.5)
        a.tick_params(labelsize=9)
        a.legend(fontsize=8)
    fig.tight_layout()
    save(fig, "fig_dose_response.png")


def fig_cfr_distributions():
    """Joint CFR for the selected regimens across all four MIC distributions."""
    rows = read("cfr_all_distributions.csv")
    keep = {r for r, _ in SELECTED}
    dists = sorted({r["distribution_id"] for r in rows})
    lookup = {(r["regimen"], r["distribution_id"]): float(r["joint_cfr_pct"])
              for r in rows if r["regimen"] in keep}

    fig, ax = plt.subplots(figsize=(11, 5))
    width = 0.8 / len(dists)
    cmap = plt.get_cmap("cividis")
    for j, did in enumerate(dists):
        xs = [i + j * width - 0.4 + width / 2 for i in range(len(SELECTED))]
        ys = [lookup.get((reg, did), 0) for reg, _ in SELECTED]
        ax.bar(xs, ys, width=width, label=did.replace("_", " "),
               color=cmap(j / max(len(dists) - 1, 1)))
    ax.axhline(80, color="#999", ls=":", lw=1.2)
    ax.axhline(90, color="#999", ls="--", lw=1.2)
    ax.set_xticks(range(len(SELECTED)))
    ax.set_xticklabels([f"{reg}\nEKFC {cls}" for reg, cls in SELECTED], fontsize=9)
    ax.set_ylabel("Joint CFR (%)", fontsize=10)
    ax.set_ylim(0, 100)
    ax.grid(axis="y", alpha=.25, lw=.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, ncol=2)
    ax.set_title("Joint CFR is stable across MIC distributions but never reaches 90%",
                 fontsize=11.5)
    fig.tight_layout()
    save(fig, "fig_cfr_distributions.png")


def fig_toxicity_gate():
    """Regimens that lose permissibility under a modest clearance error."""
    rows = read("gsa_toxicity_gate_crossings.csv")
    if not rows:
        print("  (no gate crossings)")
        return
    rows.sort(key=lambda r: float(r["scenario_toxicity_pct"]))
    labels = [f"{r['regimen']}  ({r['scenario']})" for r in rows]
    base = [float(r["base_toxicity_pct"]) for r in rows]
    scen = [float(r["scenario_toxicity_pct"]) for r in rows]

    fig, ax = plt.subplots(figsize=(9.5, 0.55 * len(rows) + 2.2))
    y = range(len(rows))
    for i, (b, s) in enumerate(zip(base, scen)):
        ax.plot([b, s], [i, i], color="#BBB", lw=2, zorder=1)
    ax.scatter(base, y, color=CAZ, s=55, zorder=2, label="Baseline")
    ax.scatter(scen, y, color=WARN, s=55, zorder=2, label="Under scenario")
    ax.axvline(15, color=AVI, ls="--", lw=1.5)
    ax.text(15.2, -0.7, "15% permissibility ceiling", fontsize=8.5, color=AVI)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Probability of exceeding the ceftazidime exposure screen (%)", fontsize=10)
    ax.grid(axis="x", alpha=.25, lw=.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8.5, loc="lower right")
    ax.set_title("A 20% error in ceftazidime clearance makes the selected\nhigh-dose regimens impermissible",
                 fontsize=11.5)
    fig.tight_layout()
    save(fig, "fig_toxicity_gate.png")


if __name__ == "__main__":
    print("figures:")
    fig_pta_vs_mic()
    fig_oat_tornado()
    fig_dose_response()
    fig_cfr_distributions()
    fig_toxicity_gate()
