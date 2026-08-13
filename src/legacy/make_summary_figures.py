"""Regenerate summary figures for the CAZ-AVI evidence-composite pilot.

This script uses the processed CSV summaries included in data/processed/.
It regenerates the manuscript-style summary figures for selected regimen
outputs, calibration status, convergence, CFR, PSA robustness, and uncertainty ranking.

Usage:
    python src/make_summary_figures.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "legacy_rc1"
FIGS = ROOT / "figures" / "regenerated"
FIGS.mkdir(parents=True, exist_ok=True)


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def save_current(name: str) -> None:
    plt.tight_layout()
    plt.savefig(FIGS / name, dpi=300, bbox_inches="tight")
    plt.close()


def plot_selected_regimens() -> None:
    rows = read_csv("selected_regimen_outputs.csv")
    labels = [f"{r['ekfc_class']}\n{r['selected_regimen']}" for r in rows]
    joint = [as_float(r, "joint_pta_mic4_percent") for r in rows]
    tox = [as_float(r, "toxicity_percent") for r in rows]
    x = range(len(rows))
    width = 0.38
    plt.figure(figsize=(7.5, 4.5))
    plt.bar([i - width / 2 for i in x], joint, width, label="Joint PTA at MIC 4 (%)")
    plt.bar([i + width / 2 for i in x], tox, width, label="Toxicity (%)")
    plt.xticks(list(x), labels)
    plt.ylabel("Percent")
    plt.title("Selected regimen performance in the primary simulation")
    plt.legend()
    save_current("fig3_selected_regimens_regenerated.png")


def plot_calibration() -> None:
    rows = read_csv("calibration_summary.csv")
    labels = [r["endpoint"] for r in rows]
    passes = [int(r["pass_rows"]) for r in rows]
    reviews = [int(r["review_rows"]) for r in rows]
    fails = [int(r["fail_rows"]) for r in rows]
    x = range(len(rows))
    plt.figure(figsize=(6.5, 4.5))
    plt.bar(x, passes, label="PASS")
    plt.bar(x, reviews, bottom=passes, label="REVIEW")
    bottom_fail = [p + rv for p, rv in zip(passes, reviews)]
    plt.bar(x, fails, bottom=bottom_fail, label="FAIL")
    plt.xticks(list(x), labels)
    plt.ylabel("Rows")
    plt.title("Calibration status across 72 published-comparison rows")
    plt.legend()
    save_current("fig4_calibration_regenerated.png")


def plot_convergence() -> None:
    rows = read_csv("convergence_summary.csv")
    xs = [int(r["virtual_subjects"]) for r in rows]
    mean_delta = [as_float(r, "mean_abs_delta_joint_pta_pp") for r in rows]
    max_delta = [as_float(r, "max_abs_delta_pp") for r in rows]
    plt.figure(figsize=(7.5, 4.5))
    plt.plot(xs, mean_delta, marker="o", label="Mean |Delta joint PTA|")
    plt.plot(xs, max_delta, marker="o", label="Maximum |Delta joint PTA|")
    plt.xscale("log")
    plt.xlabel("Virtual subjects")
    plt.ylabel("Absolute difference vs. 100,000 reference (percentage points)")
    plt.title("Monte Carlo convergence for joint PTA")
    plt.legend()
    save_current("fig5_convergence_regenerated.png")


def plot_cfr() -> None:
    rows = read_csv("cfr_selected_outputs.csv")
    by_label: dict[str, dict[str, float]] = {}
    for r in rows:
        label = f"{r['ekfc_class']} {r['regimen']}"
        by_label.setdefault(label, {})[r["mic_distribution"]] = as_float(r, "joint_cfr_percent")
    labels = list(by_label.keys())
    kpc = [by_label[l].get("LEE2022_KPC_KP", 0) for l in labels]
    oxa = [by_label[l].get("LEE2022_OXA_KP", 0) for l in labels]
    x = range(len(labels))
    width = 0.38
    plt.figure(figsize=(7.8, 4.5))
    plt.bar([i - width / 2 for i in x], kpc, width, label="Lee KPC-KP")
    plt.bar([i + width / 2 for i in x], oxa, width, label="Lee OXA-like KP")
    plt.axhline(80, linestyle="--", linewidth=1)
    plt.axhline(90, linestyle="--", linewidth=1)
    plt.xticks(list(x), labels, rotation=25, ha="right")
    plt.ylabel("Joint CFR (%)")
    plt.title("MIC-weighted joint CFR for selected regimens")
    plt.legend()
    save_current("fig6_cfr_regenerated.png")


def plot_psa_robustness() -> None:
    rows = read_csv("psa_robustness_summary.csv")
    labels = [r["classification"] for r in rows]
    vals = [int(r["regimen_distribution_pairs"]) for r in rows]
    plt.figure(figsize=(6.5, 4.5))
    bars = plt.bar(labels, vals)
    for bar, val in zip(bars, vals):
        plt.text(bar.get_x() + bar.get_width() / 2, val + 0.2, str(val), ha="center")
    plt.ylabel("Regimen-distribution pairs")
    plt.title("PSA baseline robustness classification")
    save_current("fig7_psa_robustness_regenerated.png")


def plot_uncertainty_ranking() -> None:
    rows = read_csv("psa_parameter_ranking.csv")
    labels = [r["description"] for r in rows]
    vals = [as_float(r, "mean_abs_rank_correlation") for r in rows]
    labels = labels[::-1]
    vals = vals[::-1]
    plt.figure(figsize=(7.5, 5.5))
    plt.barh(labels, vals)
    plt.xlabel("Mean absolute rank correlation")
    plt.title("Dominant uncertainty drivers for baseline joint CFR")
    save_current("fig8_uncertainty_drivers_regenerated.png")


def main() -> None:
    plot_selected_regimens()
    plot_calibration()
    plot_convergence()
    plot_cfr()
    plot_psa_robustness()
    plot_uncertainty_ranking()
    print(f"Regenerated figures written to {FIGS}")


if __name__ == "__main__":
    main()
