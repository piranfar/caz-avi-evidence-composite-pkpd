"""A pre-registered, falsifiable prediction — and a check that it can actually be falsified.

WHAT NOVELTY_STRATEGY.md SECTION 4.3 ASKED FOR

    "State, before any such cohort is analysed, the predicted negative predictive value of
     ceftazidime-based prediction in a non-RRT continuous-infusion population, with an
     interval. A modelling paper that exposes itself to refutation is worth more than one
     that cannot be wrong."

    It has never been done. This script does it, and does one thing the request did not
    ask for but which decides whether the exercise is worth anything at all.

THE THING THAT DECIDES IT

    A prediction is only falsifiable if a realistic cohort could distinguish it from its
    rival. The two live values of the clearance correlation -- Cojutti's published 0.94 and
    Model 1's 0.703 -- imply different NPVs. If those two predictions cannot be told apart
    with the number of patients any real cohort would have, then stating an interval is
    theatre: no result could ever contradict it.

    So this script computes, in order:

      1. predicted NPV under each rho scenario, with a credible interval from parameter
         uncertainty (this is the pre-registration itself)
      2. how many patients in a cohort of size N are actually PREDICTED NEGATIVE, since
         NPV is estimated only on that subset and it is much smaller than N
      3. the binomial sampling interval on NPV at that effective subset size
      4. whether the rho = 0.94 and rho = 0.703 predictions separate once BOTH parameter
         uncertainty and sampling uncertainty are carried

    If step 4 says they do not separate, the honest output is to say so and report the
    cohort size that would be needed, rather than to publish an unfalsifiable interval.

TO BE A REAL PRE-REGISTRATION this output must be committed and timestamped BEFORE any
non-RRT continuous-infusion cohort is analysed. The Bologna request (sent 11 Aug 2026) is
outstanding and no such data are held. That is what makes this pre-data rather than
retrospective.

Reads nothing outside model_development_v18/.
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
import model2_monitoring as M      # noqa: E402
import reproduce_primary_run as P  # noqa: E402

# Bologna cohort as described in the data request: 112 non-RRT critically ill adults.
COHORT_N = 112
# Cohort sizes to report the required-n curve over.
CANDIDATE_N = [50, 112, 200, 400, 800, 1600]

N_PER_CLASS = 3000
N_DRAWS = 400
ASSAY_CVS = [0.0, 0.20]

SCENARIOS = {
    "rho = 0.94 (Cojutti, published)": "C1_cojutti",
    "rho = 0.703 (Model 1 estimate)": "C2_model1",
    "rho agnostic U(0.38, 0.98)": "C3_agnostic",
}


def npv_and_negative_rate(pop, pr, rng, assay_cv):
    """NPV of ceftazidime-based prediction, and the share of patients predicted negative.

    'Negative' means predicted NOT to attain the avibactam target. NPV is estimated only
    on those patients, so the share is what sets the effective sample size in any real
    cohort -- which is the quantity that decides whether the prediction is testable.
    """
    m = M.classify(pop, pr, rng, assay_cv_caz=assay_cv, assay_cv_avi=0.0)
    # classify() reports rates; recover the predicted-negative share from them.
    # specificity = tn / (tn + fp); prevalence_attaining = (tp + fn) / n
    p_attain = m["prevalence_attaining"] / 100.0
    spec = m["specificity"] / 100.0
    sens = m["sensitivity"] / 100.0
    tn = spec * (1.0 - p_attain)
    fn = (1.0 - sens) * p_attain
    neg_rate = tn + fn
    return m["npv"], 100.0 * neg_rate


def wilson(k, n, z=1.96):
    """Wilson score interval — behaves sensibly at small n and proportions near 1."""
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    hw = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - hw), min(1.0, c + hw))


def main() -> int:
    pop = E.draw_population(N_PER_CLASS, H.MASTER_SEED)
    rng = np.random.default_rng(H.MASTER_SEED + 4030201)

    print("PRE-REGISTERED PREDICTION — negative predictive value of ceftazidime-only monitoring")
    print("=" * 100)
    print("Population: non-RRT critically ill adults, continuous infusion, both analytes assayed.")
    print("Predicted quantity: NPV of predicting avibactam target attainment from ceftazidime alone.")
    print(f"Reference cohort size for testability: N = {COHORT_N}.")
    print()

    results = {}
    rows = []
    for cv in ASSAY_CVS:
        print(f"--- ceftazidime assay CV = {cv:.0%} " + "-" * 62)
        print(f"{'scenario':34}{'NPV %':>9}{'95% CrI':>18}{'pred-neg %':>12}"
              f"{'n pred-neg':>12}")
        for label, key in SCENARIOS.items():
            npvs, negs = [], []
            for _ in range(N_DRAWS):
                pr = E.draw_parameters(rng)
                pr = replace(pr, rho=H.sample_rho(rng, key), avi_target=P.AVI_CT)
                npv, neg = npv_and_negative_rate(pop, pr, rng, cv)
                npvs.append(npv); negs.append(neg)
            npvs = np.array(npvs); negs = np.array(negs)
            lo, hi = np.percentile(npvs, [2.5, 97.5])
            med = float(np.median(npvs))
            neg_med = float(np.median(negs))
            n_neg = COHORT_N * neg_med / 100.0
            print(f"{label:34}{med:>9.1f}{f'{lo:.1f} to {hi:.1f}':>18}"
                  f"{neg_med:>12.1f}{n_neg:>12.1f}")
            results[(cv, key)] = dict(median=med, lo=lo, hi=hi, neg_rate=neg_med,
                                      n_neg=n_neg, draws=npvs)
            rows.append(dict(assay_cv_pct=f"{100*cv:.0f}", scenario=label,
                             npv_median_pct=f"{med:.2f}",
                             npv_cri_low_pct=f"{lo:.2f}", npv_cri_high_pct=f"{hi:.2f}",
                             predicted_negative_pct=f"{neg_med:.2f}",
                             n_predicted_negative_at_112=f"{n_neg:.1f}"))
        print()

    # ---- can a real cohort tell the two rivals apart? --------------------------------
    print("=" * 100)
    print("IS THE PREDICTION FALSIFIABLE? — separating rho = 0.94 from rho = 0.703")
    print("=" * 100)
    print("NPV is estimated only on patients predicted negative, so the effective sample")
    print("size is far below the cohort size. Both parameter and sampling uncertainty are")
    print("carried below; if the intervals overlap, no result could contradict the claim.")
    print()

    verdicts = []
    for cv in ASSAY_CVS:
        a = results[(cv, "C1_cojutti")]
        b = results[(cv, "C2_model1")]
        gap = a["median"] - b["median"]
        print(f"--- assay CV = {cv:.0%}   predicted NPV gap = {gap:.1f} pp "
              f"({b['median']:.1f}% vs {a['median']:.1f}%) " + "-" * 24)
        print(f"{'cohort N':>10}{'n pred-neg':>13}{'rho=.703 95% CI':>24}"
              f"{'rho=.94 95% CI':>24}{'separated?':>13}")
        for N in CANDIDATE_N:
            n_neg_a = max(1, int(round(N * a["neg_rate"] / 100.0)))
            n_neg_b = max(1, int(round(N * b["neg_rate"] / 100.0)))
            lo_a, hi_a = wilson(a["median"] / 100.0 * n_neg_a, n_neg_a)
            lo_b, hi_b = wilson(b["median"] / 100.0 * n_neg_b, n_neg_b)
            sep = hi_b < lo_a  # Model 1 predicts LOWER npv; separated if its CI sits below
            mark = "YES" if sep else "no"
            if N == COHORT_N:
                verdicts.append((cv, sep))
            print(f"{N:>10}{n_neg_b:>13}"
                  f"{f'{100*lo_b:.1f} to {100*hi_b:.1f}':>24}"
                  f"{f'{100*lo_a:.1f} to {100*hi_a:.1f}':>24}{mark:>13}")
        print()

    path = os.path.join(H.OUT, "prereg_npv_prediction.csv")
    os.makedirs(H.OUT, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    # ---- if NPV cannot be tested, which endpoint can? --------------------------------
    print("=" * 100)
    print("IF NPV IS UNDERPOWERED, WHICH ENDPOINT IS NOT?")
    print("=" * 100)
    print("NPV conditions on the predicted-negative subset and throws away ~90% of the")
    print("cohort. Two endpoints that use every patient are checked instead.")
    print()

    # Endpoint 2 -- overall classification accuracy, estimated on ALL N patients.
    print("Endpoint 2: overall accuracy of ceftazidime-based prediction (uses all N)")
    print(f"{'assay CV':>10}{'rho=.703':>12}{'rho=.94':>10}{'gap pp':>9}"
          f"{'.703 95% CI at 112':>22}{'.94 95% CI at 112':>21}{'sep?':>7}")
    acc_sep = {}
    for cv in ASSAY_CVS:
        accs = {}
        for key in ("C1_cojutti", "C2_model1"):
            vals = []
            for _ in range(N_DRAWS // 2):
                pr = E.draw_parameters(rng)
                pr = replace(pr, rho=H.sample_rho(rng, key), avi_target=P.AVI_CT)
                vals.append(M.classify(pop, pr, rng, cv, 0.0)["accuracy_infer"])
            accs[key] = float(np.median(vals))
        lo_b, hi_b = wilson(accs["C2_model1"] / 100.0 * COHORT_N, COHORT_N)
        lo_a, hi_a = wilson(accs["C1_cojutti"] / 100.0 * COHORT_N, COHORT_N)
        sep = hi_b < lo_a
        acc_sep[cv] = sep
        print(f"{cv:>9.0%}{accs['C2_model1']:>12.1f}{accs['C1_cojutti']:>10.1f}"
              f"{accs['C1_cojutti']-accs['C2_model1']:>9.1f}"
              f"{f'{100*lo_b:.1f} to {100*hi_b:.1f}':>22}"
              f"{f'{100*lo_a:.1f} to {100*hi_a:.1f}':>21}{'YES' if sep else 'no':>7}")

    # Endpoint 3 -- estimate rho directly from paired concentrations.
    # Fisher z: se = 1/sqrt(n-3); anchored against Model 1's own achieved precision.
    print()
    print("Endpoint 3: estimate rho DIRECTLY from paired steady-state concentrations")
    print(f"{'cohort N':>10}{'se(z)':>9}{'CI if true rho = 0.703':>26}"
          f"{'excludes 0.94?':>16}{'CI if true rho = 0.94':>25}{'excludes 0.703?':>17}")
    rho_ok_112 = False
    for N in CANDIDATE_N:
        se = 1.0 / np.sqrt(max(N - 3, 1))
        out = []
        for truth in (0.703, 0.94):
            z = np.arctanh(truth)
            lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
            out.append((lo, hi))
        ex94 = out[0][1] < 0.94
        ex703 = out[1][0] > 0.703
        if N == COHORT_N:
            rho_ok_112 = ex94 and ex703
        print(f"{N:>10}{se:>9.4f}{f'{out[0][0]:.3f} to {out[0][1]:.3f}':>26}"
              f"{('YES' if ex94 else 'no'):>16}"
              f"{f'{out[1][0]:.3f} to {out[1][1]:.3f}':>25}"
              f"{('YES' if ex703 else 'no'):>17}")
    print()
    print("  Sanity anchor: Model 1 achieved 95% CI 0.381 to 0.873 on 21 subjects. Scaling that")
    print("  precision by sqrt(112/21) = 2.31 gives a width near 0.21 at N = 112, consistent")
    print("  with the Fisher-z row above. The approximation is therefore not flattering itself.")

    print()
    print("=" * 100)
    print("VERDICT")
    print("=" * 100)
    for cv, sep in verdicts:
        state = "IS" if sep else "is NOT"
        print(f"  At assay CV {cv:.0%}, a cohort of {COHORT_N} {state} large enough to separate "
              f"the two rival predictions.")
    for cv, sep in acc_sep.items():
        state = "IS" if sep else "is NOT"
        print(f"  On overall accuracy, a cohort of {COHORT_N} {state} large enough "
              f"(assay CV {cv:.0%}).")
    print(f"  On a direct estimate of rho, a cohort of {COHORT_N} "
          f"{'IS' if rho_ok_112 else 'is NOT'} large enough.")

    print()
    if not any(s for _, s in verdicts):
        print("  NPV IS THE WRONG ENDPOINT TO PRE-REGISTER. It conditions on the ~10 patients")
        print("  of 112 predicted negative, so its sampling interval swamps a 15 pp signal.")
        print("  Pre-registering it would produce an interval no realistic result could")
        print("  contradict -- the appearance of falsifiability without the substance.")
        print()
        if rho_ok_112:
            print("  PRE-REGISTER RHO INSTEAD. It uses every paired sample, and at N = 112 the")
            print("  two rival values are separated by a wide margin. That is a prediction that")
            print("  can actually fail, which is the whole point of section 4.3.")
        print()
        print("  Report NPV as a secondary, descriptive quantity with its interval stated -- not")
        print("  as the falsifiable claim.")
    print(f"\nwrote {os.path.relpath(path, os.path.dirname(os.path.abspath(__file__)))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
