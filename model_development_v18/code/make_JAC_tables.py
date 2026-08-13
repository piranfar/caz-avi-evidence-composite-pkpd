"""The three tables for the JAC draft, built from analysis outputs rather than typed.

  Table 1  clearance correlation -> induced attainment correlation, against the
           Frechet-Hoeffding bound and the fraction of it reached
  Table 2  accuracy of inference and of measurement, and the gain, by correlation
           scenario and assay imprecision
  Table 3  EVPPI ranking across uncertain inputs, for both decisions

Each is written as CSV (for the submission system) and as Markdown (to paste into the
draft). Nothing is hand-entered: every cell traces to a file in outputs/, so the tables
cannot drift from the analysis the way typed ones do.

Writes only into model_development_v18/.
"""

from __future__ import annotations

import csv
import os
import sys
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "outputs")
OUT = os.path.join(HERE, "..", "manuscript_JAC", "tables")


def read(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write(name, header, rows, note=""):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"{name}.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    with open(os.path.join(OUT, f"{name}.md"), "w", encoding="utf-8") as fh:
        fh.write("| " + " | ".join(header) + " |\n")
        fh.write("|" + "|".join("---" for _ in header) + "|\n")
        for r in rows:
            fh.write("| " + " | ".join(str(c) for c in r) + " |\n")
        if note:
            fh.write("\n" + note + "\n")
    print(f"  wrote {name}.csv and {name}.md  ({len(rows)} rows)")


# --------------------------------------------------------------------------------------
def table1():
    src = read("dispute_boundary_fresan_gatti.csv")
    rows = []
    for r in src:
        rows.append([
            f"{float(r['clearance_rho']):.3f}".rstrip("0").rstrip("."),
            f"{float(r['induced_attainment_phi']):.3f}",
            f"{float(r['frechet_bound_phi']):.3f}",
            f"{float(r['pct_of_bound_reached']):.0f}%",
            "no" if r["inference_catches_up"] == "NO" else "yes",
        ])
    note = ("Attainment correlation is the phi coefficient between the binary ceftazidime and "
            "avibactam attainment indicators, evaluated at the EUCAST clinical breakpoint of "
            "8 mg/L. The bound is the Frechet-Hoeffding limit for the two attainment "
            "prevalences and holds for any joint distribution with those margins. Common random "
            "numbers across the grid, so the bound is flat by construction: rho changes the "
            "joint distribution and neither margin.")
    write("Table1_attainment_correlation_bound",
          ["Clearance correlation, rho", "Induced attainment correlation, phi",
           "Frechet-Hoeffding bound", "Fraction of bound reached",
           "Inference catches up?"], rows, note)


# --------------------------------------------------------------------------------------
def table2():
    src = read("model2_monitoring_decision.csv")
    label = {"C1_cojutti": "0.94 (published)",
             "C2_model1": "0.703 (Model 1)",
             "C3_agnostic": "agnostic, U(0.38, 0.98)"}
    order = ["C1_cojutti", "C2_model1", "C3_agnostic"]
    rows = []
    for key in order:
        for r in sorted((x for x in src if x["rho_scenario"] == key
                         and x["target_scenario"] == "T2_point_4"),
                        key=lambda x: int(x["assay_cv_caz_pct"])):
            rows.append([
                label[key],
                f"{int(r['assay_cv_caz_pct'])}%",
                f"{float(r['accuracy_infer_median']):.1f}",
                f"{float(r['accuracy_measure_median']):.1f}",
                f"{float(r['accuracy_gain_median_pp']):.1f}",
                f"{float(r['accuracy_gain_p2.5_pp']):.1f} to {float(r['accuracy_gain_p97.5_pp']):.1f}",
                f"{float(r['p_gain_positive_pct']):.1f}%",
            ])
    note = ("Avibactam target fixed at 4 mg/L. Accuracy is correct classification of avibactam "
            "target attainment. Interval is the 2.5th to 97.5th percentile across parameter "
            "draws. The final column is the proportion of draws in which measuring was better, "
            "not a P value.")
    write("Table2_inference_vs_measurement",
          ["Correlation scenario", "Assay CV", "Infer (%)", "Measure (%)",
           "Gain (pp)", "95% interval (pp)", "P(gain > 0)"], rows, note)


# --------------------------------------------------------------------------------------
def table3():
    src = read("model2_evppi.csv")
    agg = defaultdict(list)
    for r in src:
        if r["rho_scenario"] == "C1_cojutti" and r["target_scenario"] == "T4_uniform":
            agg[r["parameter"]].append(float(r["evppi_pp"]))
    means = {k: sum(v) / len(v) for k, v in agg.items()}
    pretty = {"target": "Avibactam target", "avi_target": "Avibactam target", "cl0_caz": "Ceftazidime clearance",
              "omega_avi": "Avibactam between-subject variability",
              "cl0_avi": "Avibactam clearance", "fu_avi": "Avibactam unbound fraction",
              "rho": "Clearance correlation"}
    ranked = sorted(means.items(), key=lambda kv: -kv[1])
    base = means.get("rho", 0.0)
    rows = []
    for name, val in ranked:
        rel = f"{val / base:.0f}x" if base > 0 else "n/a"
        rows.append([pretty.get(name, name), f"{val:.3f}", rel])
    note = ("Decision A (which regimen), averaged over renal-function classes, correlation "
            "scenario C1 and target scenario T4. For decision B (measure or infer), EVPPI on "
            "the clearance correlation is 0.0000 pp in every scenario examined -- these are two "
            "different decisions and the two rho figures are not the same number. Units are "
            "percentage points of attainment, not currency or lives.")
    write("Table3_evppi_ranking",
          ["Uncertain input", "EVPPI (pp)", "Relative to correlation"], rows, note)


def main():
    print("Building JAC tables")
    table1()
    table2()
    table3()
    print(f"\nAll tables in {os.path.relpath(OUT, HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
