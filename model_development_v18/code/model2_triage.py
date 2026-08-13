"""WHICH patients need the avibactam assay? — a triage monitoring rule.

THE PROBLEM WITH "MEASURE AVIBACTAM"
    Model 2 shows measuring beats inferring at every plausible correlation. But a
    hospital hearing "measure both components in every patient" may simply decline on
    cost, and then nothing changes. The useful question is narrower:

        in WHICH patients does the second assay actually change the answer?

THE RULE
    Inference fails near the decision boundary. For a patient whose ceftazidime
    concentration is far above or far below the level that implies avibactam
    attainment, the conditional probability of attainment is near 1 or near 0 and the
    inference is safe. For a patient near the boundary it is close to one half and the
    inference is worthless. So:

        measure avibactam only when the conditional probability of attainment given
        the observed ceftazidime concentration falls between t and 1 - t;
        otherwise infer.

    Sweeping t from 0 to 0.5 traces the whole trade-off between assays performed and
    accuracy achieved, from "measure nobody" to "measure everybody".

    The window is also reported in CEFTAZIDIME CONCENTRATION units, per renal class,
    because that is the form a clinician can act on: measure avibactam when the
    ceftazidime result falls between these two numbers.

NOT A DOSING RECOMMENDATION. This is a monitoring-strategy analysis under a model,
and the accuracy it reports is accuracy against the model's own definition of
attainment, not against a clinical outcome.
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
N_PER_CLASS = 20000
SELECTED = E.SELECTED_REGIMENS


def _phi(x):
    return 0.5 * (1.0 + np.vectorize(erf)(x / sqrt(2.0)))


def per_patient(pop, pr, regimen, rng, assay_cv_caz=0.0):
    """Truth, conditional probability of attainment, and the observed ceftazidime."""
    cls, dose_g, interval_h = P.REGIMENS[regimen]
    ekfc, z = pop[cls]
    cl_caz, cl_avi = E.clearances(ekfc, z, pr)
    rate_caz = dose_g * 1000.0 * P.CAZ_FRACTION / interval_h
    rate_avi = dose_g * 1000.0 * P.AVI_FRACTION / interval_h
    css_caz, css_avi = rate_caz / cl_caz, rate_avi / cl_avi

    truth = css_avi * pr.fu_avi >= pr.avi_target
    obs_caz = css_caz if assay_cv_caz <= 0 else css_caz * np.exp(
        rng.normal(0.0, np.sqrt(np.log(1 + assay_cv_caz ** 2)), css_caz.size))

    typ_caz = pr.cl0_caz * (ekfc / P.EKFC_REF) ** pr.exp_caz
    eta_caz_hat = np.log((rate_caz / obs_caz) / typ_caz)
    mu = pr.rho * (pr.omega_avi / pr.omega_caz) * eta_caz_hat
    sd = pr.omega_avi * np.sqrt(max(1.0 - pr.rho ** 2, 1e-12))
    typ_avi = pr.cl0_avi * (ekfc / P.EKFC_REF) ** pr.exp_avi
    thresh = np.log((rate_avi * pr.fu_avi / pr.avi_target) / typ_avi)
    p_attain = _phi((thresh - mu) / sd)
    return truth, p_attain, obs_caz


def triage_curve(pop, pr, rng, thresholds, assay_cv_caz=0.0, assay_cv_avi=0.0,
                 regimens=SELECTED):
    rows = []
    for t in thresholds:
        measured = correct = total = 0
        for reg in regimens:
            truth, p, _ = per_patient(pop, pr, reg, rng, assay_cv_caz)
            need = (p > t) & (p < 1 - t)              # uncertain -> measure
            inferred_ok = ((p >= 0.5) == truth) & ~need
            if assay_cv_avi > 0:
                # measurement is imperfect near the threshold too
                cls, dose_g, interval_h = P.REGIMENS[reg]
                ekfc, z = pop[cls]
                _, cl_avi = E.clearances(ekfc, z, pr)
                css_avi = dose_g * 1000.0 * P.AVI_FRACTION / interval_h / cl_avi
                obs = css_avi * np.exp(rng.normal(
                    0.0, np.sqrt(np.log(1 + assay_cv_avi ** 2)), css_avi.size))
                measured_ok = ((obs * pr.fu_avi >= pr.avi_target) == truth) & need
            else:
                measured_ok = need
            measured += int(need.sum())
            correct += int(inferred_ok.sum() + measured_ok.sum())
            total += truth.size
        rows.append({
            "probability_band": round(float(t), 3),
            "pct_patients_measured": round(100.0 * measured / total, 2),
            "accuracy_pct": round(100.0 * correct / total, 2),
        })
    return rows


def concentration_window(pop, pr, regimen, t):
    """The triage band expressed as a ceftazidime concentration range."""
    cls, dose_g, interval_h = P.REGIMENS[regimen]
    ekfc, z = pop[cls]
    rate_caz = dose_g * 1000.0 * P.CAZ_FRACTION / interval_h
    rate_avi = dose_g * 1000.0 * P.AVI_FRACTION / interval_h
    # invert the conditional probability at the class-median renal function
    med = float(np.median(ekfc))
    typ_caz = pr.cl0_caz * (med / P.EKFC_REF) ** pr.exp_caz
    typ_avi = pr.cl0_avi * (med / P.EKFC_REF) ** pr.exp_avi
    thresh = np.log((rate_avi * pr.fu_avi / pr.avi_target) / typ_avi)
    sd = pr.omega_avi * np.sqrt(max(1.0 - pr.rho ** 2, 1e-12))
    slope = pr.rho * (pr.omega_avi / pr.omega_caz)

    from scipy.stats import norm
    if t <= 0.0 or t >= 0.5:
        return float("nan"), float("nan")     # band 0 measures everyone; no window
    out = []
    for target_p in (1 - t, t):
        # p = Phi((thresh - slope*eta)/sd)  ->  eta = (thresh - sd*z_p)/slope
        eta = (thresh - sd * norm.ppf(target_p)) / slope
        cl = typ_caz * np.exp(eta)
        out.append(rate_caz / cl)
    lo, hi = sorted(out)
    return lo, hi


def main():
    pop = E.draw_population(N_PER_CLASS, H.MASTER_SEED)
    rng = np.random.default_rng(H.MASTER_SEED + 31)
    thresholds = np.round(np.arange(0.0, 0.501, 0.025), 3)

    print("=" * 78)
    print("TRIAGE MONITORING RULE — which patients need the avibactam assay?")
    print("=" * 78)

    rows = []
    for rho_label, rho in (("published 0.94", 0.94), ("Model 1 0.703", 0.703)):
        pr = replace(E.BASE, rho=rho, avi_target=4.0)
        # both assays carry the same imprecision; giving ceftazidime a perfect
        # assay while avibactam has a 20% CV would flatter the inference strategy
        curve = triage_curve(pop, pr, rng, thresholds, 0.20, 0.20)
        for c in curve:
            c["rho_scenario"] = rho_label
            c["rho"] = rho
        rows += curve

        # band 0 measures every patient whose probability is not exactly 0 or 1;
        # band 0.5 measures nobody. The curve therefore runs from ALL to NONE.
        all_acc = curve[0]["accuracy_pct"]
        none_acc = curve[-1]["accuracy_pct"]
        gain = all_acc - none_acc
        print(f"\n  rho = {rho} ({rho_label}), avibactam assay CV 20%")
        print(f"    infer in everyone : accuracy {none_acc:.2f}%")
        print(f"    measure everyone  : accuracy {all_acc:.2f}%   "
              f"(gain {gain:.2f} pp, for an assay in every patient)")
        print(f"\n    {'band':>6} {'% measured':>11} {'accuracy':>10} "
              f"{'% of the gain':>15} {'gain per 100 assays':>21}")
        for c in curve:
            if c["probability_band"] in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5):
                share = (c["accuracy_pct"] - none_acc) / gain * 100 if gain else 0.0
                m = c["pct_patients_measured"]
                eff = ((c["accuracy_pct"] - none_acc) / m * 100) if m > 0.05 else float("nan")
                eff_s = "—" if np.isnan(eff) else f"{eff:.2f}"
                print(f"    {c['probability_band']:>6} "
                      f"{c['pct_patients_measured']:>10.1f}% "
                      f"{c['accuracy_pct']:>9.2f}% {share:>14.1f}% {eff_s:>20}")

        peak = max(curve, key=lambda c: c["accuracy_pct"])
        if peak["accuracy_pct"] > all_acc + 0.01:
            print(f"\n    NOTE: selective measurement BEATS measuring everyone. "
                  f"Peak accuracy {peak['accuracy_pct']:.2f}% at band "
                  f"{peak['probability_band']} ({peak['pct_patients_measured']:.1f}% "
                  f"measured), against {all_acc:.2f}% when everyone is measured.")
            print("    Assaying a patient whose inference was already confident can")
            print("    only introduce assay error, so measuring everyone is not the")
            print("    optimum.")

        # the operating point: smallest fraction measured that recovers 90% of the gain
        # scan from fewest assays upward and stop at the FIRST band that reaches
        # 90% of the benefit, which is the cheapest operating point that does
        best = None
        for c in reversed(curve):
            if gain > 0 and (c["accuracy_pct"] - none_acc) / gain >= 0.90:
                best = c
                break
        if best:
            print(f"\n    90% of the benefit is reached by measuring "
                  f"{best['pct_patients_measured']:.1f}% of patients "
                  f"(band {best['probability_band']})")
            print("    ceftazidime concentration window per renal class at that band:")
            for reg in SELECTED:
                lo, hi = concentration_window(pop, pr, reg, best["probability_band"])
                cls = P.REGIMENS[reg][0]
                if np.isnan(lo):
                    print(f"      {cls:9} ({reg:3}) band 0 — measure everyone")
                else:
                    print(f"      {cls:9} ({reg:3}) measure avibactam when total "
                          f"ceftazidime Css is {lo:6.1f} to {hi:6.1f} mg/L")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "model2_triage_curve.csv")
    fields = ["rho_scenario", "rho", "probability_band", "pct_patients_measured",
              "accuracy_pct"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {os.path.relpath(path, os.path.dirname(E.HERE))} "
          f"({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
