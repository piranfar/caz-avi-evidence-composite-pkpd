"""R2's model-derived triage rule versus Gatti 2024's clinically-derived one.

THE TWO RULES ANSWER THE SAME QUESTION BY DIFFERENT CRITERIA
    Both ask: which patients need the avibactam assay? They disagree on how to tell.

    R2 (model-derived, `model2_triage.py`)
        Measure when the conditional probability of avibactam attainment given the
        observed ceftazidime concentration is near one half -- i.e. near the decision
        boundary, where the inference is worthless. Confident either way -> infer.

    Gatti, Viale & Pea 2024 (clinically derived, J Antimicrob Chemother 79:195-9,
    doi:10.1093/jac/dkad367, archived at data_external/Gatti2024_ratio_one_leg/)
        Measure when renal function is high. Their ROC analysis found CrCL
        > 75 mL/min/1.73 m2 identifies ceftazidime-to-avibactam ratios > 5:1
        (AUC 0.694) and > 78 identifies ratios > 6:1 (AUC 0.694), i.e. the patients
        in whom extrapolating avibactam from the fixed 4:1 vial ratio fails.

WHY THEY DISAGREE -- AND WHAT THE RUN ACTUALLY SHOWED
    The prediction written here before running was that Gatti's high-ratio patients
    would sit at a conditional attainment probability near ZERO (confidently
    non-attaining), where R2 says no assay is needed. THAT PREDICTION WAS WRONG IN
    DIRECTION, and the real answer is more interesting.

    Under this project's renally-adjusted dosing grid the dose escalates steeply with
    renal function (1.25 g/day in the lowest EKFC class up to 10 g/day in the highest).
    So the high-EKFC patients Gatti's rule selects have their ratio distorted -- median
    CAZ:AVI 5.14:1 versus 4.04:1 in the rest, exactly as Gatti found -- but their
    attainment is ALREADY PROTECTED by the larger dose: median P(attain) = 0.999, and
    inference is already correct in 91.6% of them versus 92.7% of everyone else. Almost
    no discrimination.

    Hold the dose FIXED instead, which is the policy Gatti's own cohort was largely on
    (85% started full dose), and the rule comes alive: inference is correct in 88.7% of
    the selected versus 97.9% of the rest at rho = 0.94, and 81.4% versus 97.1% at
    rho = 0.703. A 9-16 point separation where there had been ~1.

    So neither rule is wrong. GATTI'S RULE ENCODES A DOSING POLICY -- "high renal
    clearance means low avibactam means risk" -- which is true when the dose does not
    compensate and largely stops being true once it does. That is a transferability
    finding about clinically-derived triage rules, not a defect in theirs.

WHAT THIS SCRIPT DOES
    Puts both rules on the same population, the same model, the same assay
    imprecision, and the same budget axis (% of patients measured), and asks which
    delivers more correct attainment calls per assay spent. It also reports the
    agreement between them and the direction of their disagreement.

    Gatti's rule is operationalised on EKFC because that is the renal-function measure
    this project's model carries. Gatti used CKD-EPI creatinine clearance. The two are
    different equations on a similar scale; the cut-off is therefore approximate and
    the comparison is of RULE SHAPE, not of a specific mL/min value. This is stated
    again in the output.

NOT A DOSING OR MONITORING RECOMMENDATION. Accuracy here is accuracy against the
model's own definition of attainment, not against a clinical outcome.

Run:  python model2_triage_vs_gatti.py
"""
from __future__ import annotations

import csv
import os
import sys
from dataclasses import replace

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import model2_engine as E          # noqa: E402
import model2_hujam as H           # noqa: E402
import model2_triage as T          # noqa: E402
import reproduce_primary_run as P  # noqa: E402

OUT = H.OUT
N_PER_CLASS = 20000
SELECTED = E.SELECTED_REGIMENS
ASSAY_CV = 0.20            # same imprecision on both assays, as in model2_triage.py

# Gatti 2024 ROC cut-offs, in mL/min/1.73 m^2 (their CKD-EPI CLcr; applied here to EKFC)
GATTI_CUTOFFS = {"ratio>5:1 (Gatti AUC 0.694)": 75.0,
                 "ratio>6:1 (Gatti AUC 0.694)": 78.0}


def _rates(dose_g, interval_h):
    return (dose_g * 1000.0 * P.CAZ_FRACTION / interval_h,
            dose_g * 1000.0 * P.AVI_FRACTION / interval_h)


def patient_level(pop, pr, rng, fixed_dose=None):
    """Per-patient truth, inference probability, renal function and CAZ:AVI ratio.

    fixed_dose=None  -> the renally-adjusted policy (SELECTED regimens, one per class),
                        which is what model2_triage.py and the manuscript's grid use.
    fixed_dose=(g,h) -> the SAME dose in every renal class, which is the policy Gatti's
                        cohort was largely on (85% started full dose 2.5 g q8h).
    """
    truth, p_attain, ekfc_all, ratio = [], [], [], []
    for reg in SELECTED:
        cls, dose_g, interval_h = P.REGIMENS[reg]
        if fixed_dose is not None:
            dose_g, interval_h = fixed_dose
        ekfc, z = pop[cls]
        cl_caz, cl_avi = E.clearances(ekfc, z, pr)
        rate_caz, rate_avi = _rates(dose_g, interval_h)
        css_caz, css_avi = rate_caz / cl_caz, rate_avi / cl_avi

        tr = css_avi * pr.fu_avi >= pr.avi_target
        obs_caz = css_caz * np.exp(rng.normal(
            0.0, np.sqrt(np.log(1 + ASSAY_CV ** 2)), css_caz.size))
        typ_caz = pr.cl0_caz * (ekfc / P.EKFC_REF) ** pr.exp_caz
        eta_caz_hat = np.log((rate_caz / obs_caz) / typ_caz)
        mu = pr.rho * (pr.omega_avi / pr.omega_caz) * eta_caz_hat
        sd = pr.omega_avi * np.sqrt(max(1.0 - pr.rho ** 2, 1e-12))
        typ_avi = pr.cl0_avi * (ekfc / P.EKFC_REF) ** pr.exp_avi
        thresh = np.log((rate_avi * pr.fu_avi / pr.avi_target) / typ_avi)
        p = T._phi((thresh - mu) / sd)

        truth.append(tr); p_attain.append(p)
        ekfc_all.append(ekfc); ratio.append(css_caz / css_avi)
    return (np.concatenate(truth), np.concatenate(p_attain),
            np.concatenate(ekfc_all), np.concatenate(ratio))


def measured_correct_exact(pop, pr, rng, fixed_dose=None):
    """Exact per-patient 'measurement got it right' flag, matching model2_triage."""
    flags = []
    for reg in SELECTED:
        cls, dose_g, interval_h = P.REGIMENS[reg]
        if fixed_dose is not None:
            dose_g, interval_h = fixed_dose
        ekfc, z = pop[cls]
        _, cl_avi = E.clearances(ekfc, z, pr)
        css_avi = dose_g * 1000.0 * P.AVI_FRACTION / interval_h / cl_avi
        truth = css_avi * pr.fu_avi >= pr.avi_target
        obs = css_avi * np.exp(rng.normal(
            0.0, np.sqrt(np.log(1 + ASSAY_CV ** 2)), css_avi.size))
        flags.append((obs * pr.fu_avi >= pr.avi_target) == truth)
    return np.concatenate(flags)


def curve_for_rule(truth, p_attain, order_stat, meas_ok, budgets, descending=True):
    """Accuracy vs % measured, measuring patients ranked by `order_stat`."""
    inferred_ok = (p_attain >= 0.5) == truth
    idx = np.argsort(-order_stat if descending else order_stat)
    n = truth.size
    rows = []
    for b in budgets:
        k = int(round(n * b / 100.0))
        measure = np.zeros(n, dtype=bool)
        measure[idx[:k]] = True
        correct = np.where(measure, meas_ok, inferred_ok)
        rows.append((b, 100.0 * correct.mean()))
    return rows


def main():
    pop = E.draw_population(N_PER_CLASS, H.MASTER_SEED)
    budgets = [0, 5, 10, 12.5, 20, 25, 30, 40, 50, 60, 80, 100]

    print("=" * 78)
    print("R2's MODEL-DERIVED TRIAGE RULE  vs  GATTI 2024's CLINICAL ROC RULE")
    print("=" * 78)
    print("  Gatti's cut-offs are CKD-EPI CLcr; applied here to EKFC, which is what")
    print("  this model carries. Different equations, similar scale -- the comparison")
    print("  is of rule SHAPE and ranking, not of a specific mL/min value.")

    rows = []
    policies = [("renally-adjusted grid (this project)", None),
                ("FIXED 2.5 g q8h in every class (Gatti's cohort)", (2.5, 8))]
    for pol_label, fixed in policies:
      print(f"\n{'='*78}\n  DOSING POLICY: {pol_label}\n{'='*78}")
      for rho_label, rho in (("published 0.94", 0.94), ("Model 1 0.703", 0.703)):
        pr = replace(E.BASE, rho=rho, avi_target=4.0)
        rng = np.random.default_rng(H.MASTER_SEED + 31)
        truth, p_attain, ekfc, ratio = patient_level(pop, pr, rng, fixed)
        rng2 = np.random.default_rng(H.MASTER_SEED + 31)
        meas_ok = measured_correct_exact(pop, pr, rng2, fixed)

        inferred_ok = (p_attain >= 0.5) == truth
        base = 100.0 * inferred_ok.mean()
        allm = 100.0 * meas_ok.mean()

        print(f"\n{'-'*78}\n  rho = {rho} ({rho_label}), assay CV {ASSAY_CV:.0%} on both")
        print(f"    attainment prevalence: {100.0*truth.mean():.1f}% of patients attain")
        print(f"    infer in everyone: {base:.2f}%    measure everyone: {allm:.2f}%"
              f"    (gain {allm-base:+.2f} pp)")

        # --- rule statistics ------------------------------------------------
        uncertainty = -np.abs(p_attain - 0.5)         # higher = nearer the boundary
        print(f"\n    Where do the two rules point?")
        for name, cut in GATTI_CUTOFFS.items():
            sel = ekfc > cut
            pct = 100.0 * sel.mean()
            print(f"      Gatti {name}: EKFC>{cut:.0f} selects {pct:5.1f}% of patients")
            print(f"        their median P(attain) = {np.median(p_attain[sel]):.3f} "
                  f"vs {np.median(p_attain[~sel]):.3f} in the rest")
            print(f"        their median CAZ:AVI   = {np.median(ratio[sel]):.2f}:1 "
                  f"vs {np.median(ratio[~sel]):.2f}:1")
            print(f"        inference is ALREADY correct in "
                  f"{100.0*inferred_ok[sel].mean():.1f}% of them "
                  f"(vs {100.0*inferred_ok[~sel].mean():.1f}% of the rest)")

        # --- budget-matched comparison ---------------------------------------
        print(f"\n    Accuracy at a matched assay budget:")
        print(f"      {'% measured':>11} {'R2 rule':>10} {'Gatti rule':>12} {'difference':>12}")
        r2 = dict(curve_for_rule(truth, p_attain, uncertainty, meas_ok, budgets))
        gt = dict(curve_for_rule(truth, p_attain, ekfc, meas_ok, budgets))
        for b in budgets:
            d = r2[b] - gt[b]
            print(f"      {b:>10.1f}% {r2[b]:>9.2f}% {gt[b]:>11.2f}% {d:>+11.2f} pp")
            # The dosing policy MUST be recorded. Without it the file stacks two policy
            # blocks with identical rho_scenario and pct_measured values, so the rows are
            # indistinguishable and anything plotted from it silently overlays them.
            rows.append({"dosing_policy": pol_label,
                         "rho_scenario": rho_label, "rho": rho,
                         "pct_measured": b,
                         "accuracy_R2_rule_pct": round(r2[b], 3),
                         "accuracy_Gatti_rule_pct": round(gt[b], 3),
                         "difference_pp": round(d, 3)})

        # --- agreement at a common operating point ---------------------------
        k = int(round(truth.size * 0.125))
        r2_sel = np.zeros(truth.size, bool); r2_sel[np.argsort(-uncertainty)[:k]] = True
        gt_sel = np.zeros(truth.size, bool); gt_sel[np.argsort(-ekfc)[:k]] = True
        both = (r2_sel & gt_sel).sum()
        jac = both / max((r2_sel | gt_sel).sum(), 1)
        print(f"\n    At a 12.5% budget the two rules pick the SAME patient in "
              f"{100.0*both/k:.1f}% of cases (Jaccard {jac:.3f}).")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "model2_triage_vs_gatti.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\n  wrote {os.path.relpath(path, os.path.dirname(E.HERE))} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
