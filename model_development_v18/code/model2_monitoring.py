"""MODEL 2, decision problem B — is a separate avibactam assay worth running?

WHY THIS EXISTS
    Decision problem A (which regimen, in `model2_hujam.py`) turned out to be almost
    insensitive to the ceftazidime-avibactam clearance correlation. That is not a bug
    and it is worth stating: rho shifts the LEVEL of joint attainment but affects every
    regimen within a renal class in the same direction, so it barely changes the
    RANKING. The first run made this explicit — misselection and EVPI were nearly
    identical under rho = 0.94 and rho = 0.703.

    Where rho is decisive is the MONITORING decision: whether avibactam must be
    measured, or can be inferred from a ceftazidime concentration. That decision is
    the clinically actionable one, and it is the one Model 1 informs.

THE TWO STRATEGIES
    INFER    measure ceftazidime only; predict avibactam attainment from the
             conditional distribution of the avibactam random effect given the
             observed ceftazidime concentration. Accuracy depends on rho.
    MEASURE  assay avibactam directly. Accuracy is limited only by assay error.

THE OUTPUT
    The ACCURACY GAIN from measuring, in percentage points of correctly classified
    patients, with its uncertainty. Reported as a break-even quantity: measuring is
    worth it when its cost, expressed in the same units, is below the gain. No
    monetary cost is invented, because none is defensible.

    Expected value of perfect information on rho is then the value of knowing rho
    BEFORE choosing a monitoring strategy. If the gain is positive across the whole
    plausible range of rho, the EVPI is zero — which is itself the answer: measure
    avibactam regardless of which value of rho turns out to be right.
"""
from __future__ import annotations

import csv
import os
import sys
from dataclasses import replace
from math import erf, sqrt

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import model2_engine as E          # noqa: E402
import model2_hujam as H           # noqa: E402
import reproduce_primary_run as P  # noqa: E402

OUT = H.OUT
SELECTED = E.SELECTED_REGIMENS


def _phi(x):
    return 0.5 * (1.0 + np.vectorize(erf)(x / sqrt(2.0)))


def classify(pop, pr, rng, assay_cv_caz=0.0, assay_cv_avi=0.0,
             regimens=SELECTED):
    """Accuracy of inferring avibactam attainment, and of measuring it directly.

    INFER: from the observed ceftazidime concentration, back out the implied
    ceftazidime random effect, then use the conditional normal distribution of the
    avibactam random effect given it. Predict attainment when that conditional
    probability is at least one half. This is the same classifier the manuscript
    uses; it is reimplemented here so it can be run inside the uncertainty loop.

    MEASURE: compare the assayed avibactam concentration with the target. Errors
    arise only from assay noise near the threshold.
    """
    tp = fp = tn = fn = 0
    correct_measure = total = 0
    for reg in regimens:
        cls, dose_g, interval_h = P.REGIMENS[reg]
        ekfc, z = pop[cls]
        cl_caz, cl_avi = E.clearances(ekfc, z, pr)
        rate_caz = dose_g * 1000.0 * P.CAZ_FRACTION / interval_h
        rate_avi = dose_g * 1000.0 * P.AVI_FRACTION / interval_h
        css_caz, css_avi = rate_caz / cl_caz, rate_avi / cl_avi

        truth = css_avi * pr.fu_avi >= pr.avi_target

        # --- strategy INFER -------------------------------------------------
        obs_caz = css_caz if assay_cv_caz <= 0 else css_caz * np.exp(
            rng.normal(0.0, np.sqrt(np.log(1 + assay_cv_caz ** 2)), css_caz.size))
        typ_caz = pr.cl0_caz * (ekfc / P.EKFC_REF) ** pr.exp_caz
        eta_caz_hat = np.log((rate_caz / obs_caz) / typ_caz)
        mu = pr.rho * (pr.omega_avi / pr.omega_caz) * eta_caz_hat
        sd = pr.omega_avi * np.sqrt(max(1.0 - pr.rho ** 2, 1e-12))
        typ_avi = pr.cl0_avi * (ekfc / P.EKFC_REF) ** pr.exp_avi
        thresh = np.log((rate_avi * pr.fu_avi / pr.avi_target) / typ_avi)
        pred = _phi((thresh - mu) / sd) >= 0.5

        tp += int(np.sum(pred & truth)); fp += int(np.sum(pred & ~truth))
        tn += int(np.sum(~pred & ~truth)); fn += int(np.sum(~pred & truth))

        # --- strategy MEASURE -----------------------------------------------
        obs_avi = css_avi if assay_cv_avi <= 0 else css_avi * np.exp(
            rng.normal(0.0, np.sqrt(np.log(1 + assay_cv_avi ** 2)), css_avi.size))
        correct_measure += int(np.sum((obs_avi * pr.fu_avi >= pr.avi_target) == truth))
        total += truth.size

    n = tp + fp + tn + fn
    return {
        "accuracy_infer": 100.0 * (tp + tn) / n,
        "accuracy_measure": 100.0 * correct_measure / total,
        "ppv": 100.0 * tp / max(tp + fp, 1),
        "npv": 100.0 * tn / max(tn + fn, 1),
        "specificity": 100.0 * tn / max(tn + fp, 1),
        "sensitivity": 100.0 * tp / max(tp + fn, 1),
        "false_reassurance": 100.0 * fp / n,
        "prevalence_attaining": 100.0 * (tp + fn) / n,
    }


def run(n_draws, n_per_class, rho_scenario, target_scenario,
        assay_cv_caz, assay_cv_avi, seed=H.MASTER_SEED):
    pop = E.draw_population(n_per_class, seed)
    rng = np.random.default_rng(seed + 104729)
    rows = []
    for _ in range(n_draws):
        pr = E.draw_parameters(rng)
        pr = replace(pr, rho=H.sample_rho(rng, rho_scenario),
                     avi_target=H.sample_target(rng, target_scenario))
        m = classify(pop, pr, rng, assay_cv_caz, assay_cv_avi)
        m["rho"] = pr.rho
        m["target"] = pr.avi_target
        m["gain"] = m["accuracy_measure"] - m["accuracy_infer"]
        rows.append(m)
    return rows


def summarise(rows, rho_scenario, target_scenario, cv_caz, cv_avi):
    g = np.array([r["gain"] for r in rows])
    return {
        "rho_scenario": rho_scenario, "target_scenario": target_scenario,
        "assay_cv_caz_pct": round(100 * cv_caz), "assay_cv_avi_pct": round(100 * cv_avi),
        "accuracy_infer_median": round(float(np.median([r["accuracy_infer"] for r in rows])), 2),
        "accuracy_measure_median": round(float(np.median([r["accuracy_measure"] for r in rows])), 2),
        "accuracy_gain_median_pp": round(float(np.median(g)), 2),
        "accuracy_gain_p2.5_pp": round(float(np.percentile(g, 2.5)), 2),
        "accuracy_gain_p97.5_pp": round(float(np.percentile(g, 97.5)), 2),
        "p_gain_positive_pct": round(100.0 * float(np.mean(g > 0)), 1),
        "npv_infer_median": round(float(np.median([r["npv"] for r in rows])), 2),
        "specificity_infer_median": round(float(np.median([r["specificity"] for r in rows])), 2),
        "false_reassurance_median_pct": round(
            float(np.median([r["false_reassurance"] for r in rows])), 2),
    }


def evpi_monitoring(rows):
    """Value of knowing rho before choosing a monitoring strategy.

    Two strategies, so the perfect-information decision is simply to take whichever
    is better in each draw. If one strategy dominates in every draw, EVPI is zero:
    the decision does not depend on the unknown.
    """
    infer = np.array([r["accuracy_infer"] for r in rows])
    meas = np.array([r["accuracy_measure"] for r in rows])
    with_pi = np.maximum(infer, meas).mean()
    without_pi = max(infer.mean(), meas.mean())
    return float(with_pi - without_pi), float(np.mean(meas > infer))


def evppi_rho(rows, degree=4):
    """Partial EVPI for the correlation alone (Strong & Oakley regression)."""
    x = np.array([r["rho"] for r in rows])
    if x.std() < 1e-12:
        return 0.0
    xs = (x - x.mean()) / x.std()
    infer = np.array([r["accuracy_infer"] for r in rows])
    meas = np.array([r["accuracy_measure"] for r in rows])
    gi = np.polyval(np.polyfit(xs, infer, degree), xs)
    gm = np.polyval(np.polyfit(xs, meas, degree), xs)
    return float(np.maximum(gi, gm).mean() - max(infer.mean(), meas.mean()))


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    n_draws = 200 if args.quick else 1000
    n_per_class = 2000 if args.quick else 4000

    print("=" * 78)
    print("MODEL 2, decision problem B — is a separate avibactam assay worth running?")
    print("=" * 78)
    print(f"  outer draws {n_draws}, virtual subjects {n_per_class*5:,}, "
          f"selected regimens {', '.join(SELECTED)}\n")

    summary, voi = [], []
    combos = [("C1_cojutti", "T2_point_4"), ("C2_model1", "T2_point_4"),
              ("C3_agnostic", "T2_point_4"), ("C1_cojutti", "T4_uniform"),
              ("C2_model1", "T4_uniform"), ("C3_agnostic", "T4_uniform")]
    assays = [(0.0, 0.0), (0.20, 0.20)]

    hdr = (f"{'rho':13}{'target':14}{'assayCV':>8}{'infer%':>9}{'measure%':>10}"
           f"{'gain pp':>9}{'P(gain>0)':>11}{'EVPPI rho':>11}")
    print(hdr)
    print("-" * len(hdr))
    for rs, ts in combos:
        for cvc, cva in assays:
            rows = run(n_draws, n_per_class, rs, ts, cvc, cva)
            s = summarise(rows, rs, ts, cvc, cva)
            ev, p_meas = evpi_monitoring(rows)
            ep = evppi_rho(rows)
            s["evpi_monitoring_pp"] = round(ev, 4)
            s["evppi_rho_pp"] = round(max(ep, 0.0), 4)
            s["p_measure_better_pct"] = round(100 * p_meas, 1)
            summary.append(s)
            print(f"{rs:13}{ts:14}{100*cvc:>7.0f}%{s['accuracy_infer_median']:>9.1f}"
                  f"{s['accuracy_measure_median']:>10.1f}"
                  f"{s['accuracy_gain_median_pp']:>9.2f}"
                  f"{s['p_gain_positive_pct']:>10.1f}%{max(ep,0.0):>11.4f}")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "model2_monitoring_decision.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
    print(f"\n  wrote {os.path.relpath(path, os.path.dirname(E.HERE))} "
          f"({len(summary)} rows)")

    print("\n  READING")
    print("  The accuracy gain is the break-even quantity: measuring avibactam is")
    print("  worth it when the cost of the assay, expressed in percentage points of")
    print("  correctly classified patients, falls below this gain. No monetary cost")
    print("  is assumed, because none would be defensible.")
    print("  An EVPPI for rho near zero means the monitoring decision does not depend")
    print("  on resolving rho — the same strategy wins across its plausible range.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
