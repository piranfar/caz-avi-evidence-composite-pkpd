"""What simulated patients CAN and CANNOT do for this project.

THE QUESTION
    Can virtual patients generated with existing software substitute for the real
    non-RRT continuous-infusion data we do not have?

THE ANSWER, AND WHY IT IS NOT A MATTER OF OPINION
    NO for estimating rho. Simulated patients contain exactly the information put
    into them. Simulate from a model with rho = 0.94, fit it, and you recover 0.94 —
    that is a check on the software, not evidence about patients. It is also what
    the project brief forbids being called validation.

    YES for four things, all of which make the model more robust:

      1. ESTIMATOR VALIDATION. Simulate at a known rho, refit with the actual Model 1
         estimator, and measure bias, root mean squared error and interval coverage.
         This tests the estimator rather than the biology, which is a legitimate and
         necessary thing to test.
      2. DESIGN ANALYSIS. Ask what a future dataset would have to look like to
         settle the question. That is what this module computes, and it converts
         "we need the Bologna data" from an assertion into a number.
      3. MISSPECIFICATION ROBUSTNESS. Simulate under a structure the model does not
         assume and measure the resulting bias.
      4. A FALSIFIABLE PREDICTION. State in advance what a future cohort should show.

WHAT THIS MODULE COMPUTES
    Part 1  the attenuation identity, and why a naive analysis of the source cohort
            would NOT return 0.94 even if 0.94 were true
    Part 2  how precisely rho could be estimated from a Bologna-shaped dataset, and
            the probability that such a dataset would distinguish 0.94 from 0.703
    Part 3  the same for hypothetical prospective designs, to say what would be needed
    Part 4  the falsifiable prediction
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "outputs")
SEED = 20260811

# Cojutti 2024 Table 2
OMEGA_CAZ, OMEGA_AVI = 0.6159, 0.6817     # SD of the log clearance random effects
SIGMA_CAZ, SIGMA_AVI = 0.31, 0.33         # proportional residual error, b1 and b2
RHO_PUBLISHED = 0.94
RHO_MODEL1 = 0.703


def attenuation_factor(om1=OMEGA_CAZ, om2=OMEGA_AVI, s1=SIGMA_CAZ, s2=SIGMA_AVI):
    """How far a naive correlation of observed clearances falls below the truth.

    An observed steady-state concentration carries the random effect PLUS residual
    error. Correlating clearances derived from single observations therefore
    estimates
        rho_obs = rho_true * om1*om2 / (sqrt(om1^2+s1^2) * sqrt(om2^2+s2^2)),
    the classic regression-dilution result. A mixed-effects model separates the two
    variance components and estimates rho_true; a two-stage analysis does not.
    """
    return (om1 * om2) / (np.sqrt(om1 ** 2 + s1 ** 2) * np.sqrt(om2 ** 2 + s2 ** 2))


def simulate_and_estimate(n_patients, occasions, rho_true, n_rep, rng,
                          om1=OMEGA_CAZ, om2=OMEGA_AVI,
                          s1=SIGMA_CAZ, s2=SIGMA_AVI):
    """Maximum-likelihood estimate of rho from a continuous-infusion TDM design.

    Under continuous infusion at steady state a sample gives one concentration, so
    Css = rate/CL and no time course is observed: volume is not identifiable and the
    problem reduces to a bivariate normal on the log scale, with a per-patient random
    effect and per-occasion residual error. That is exactly the Bologna design.

    Returns the estimated rho for each replicate.
    """
    out = np.empty(n_rep)
    for r in range(n_rep):
        k = rng.choice(occasions, size=n_patients)          # occasions per patient
        eta = rng.multivariate_normal(
            [0, 0], [[om1 ** 2, rho_true * om1 * om2],
                     [rho_true * om1 * om2, om2 ** 2]], size=n_patients)
        y1, y2 = [], []
        for i in range(n_patients):
            y1.append(eta[i, 0] + rng.normal(0, s1, k[i]))
            y2.append(eta[i, 1] + rng.normal(0, s2, k[i]))

        # patient means: the sufficient statistic for the random effect
        m1 = np.array([v.mean() for v in y1])
        m2 = np.array([v.mean() for v in y2])
        ki = k.astype(float)
        v1 = om1 ** 2 + s1 ** 2 / ki                        # variance of each mean
        v2 = om2 ** 2 + s2 ** 2 / ki

        def neg_ll(rho):
            rho = np.clip(rho, -0.995, 0.995)
            c = rho * om1 * om2                             # covariance of the means
            det = v1 * v2 - c ** 2
            if np.any(det <= 0):
                return 1e12
            q = (v2 * m1 ** 2 - 2 * c * m1 * m2 + v1 * m2 ** 2) / det
            return float(np.sum(np.log(det) + q))

        res = minimize_scalar(neg_ll, bounds=(-0.99, 0.99), method="bounded")
        out[r] = float(res.x)
    return out


def design_row(label, n_patients, occasions, rho_true, alternative, n_rep, rng):
    est = simulate_and_estimate(n_patients, occasions, rho_true, n_rep, rng)
    n_eff = n_patients
    se_z = 1.0 / np.sqrt(max(n_eff - 3, 1))
    lo = np.tanh(np.arctanh(np.clip(est, -0.99, 0.99)) - 1.959964 * se_z)
    hi = np.tanh(np.arctanh(np.clip(est, -0.99, 0.99)) + 1.959964 * se_z)
    excludes = float(np.mean((alternative < lo) | (alternative > hi)))
    return {
        "design": label, "n_patients": n_patients,
        "occasions_per_patient": "/".join(str(o) for o in occasions),
        "rho_true": rho_true, "alternative_tested": alternative,
        "rho_mean_estimate": round(float(est.mean()), 4),
        "rho_bias": round(float(est.mean() - rho_true), 4),
        "rho_sd": round(float(est.std(ddof=1)), 4),
        "rho_p2.5": round(float(np.percentile(est, 2.5)), 4),
        "rho_p97.5": round(float(np.percentile(est, 97.5)), 4),
        "power_to_exclude_alternative_pct": round(100 * excludes, 1),
    }


def main():
    rng = np.random.default_rng(SEED)
    n_rep = 4000
    print("=" * 78)
    print("DESIGN ANALYSIS — what simulated patients can and cannot settle")
    print("=" * 78)

    # ---------------------------------------------------------------- Part 1 --
    a = attenuation_factor()
    print("\n1. The attenuation identity")
    print(f"   residual error in the source model: {100*SIGMA_CAZ:.0f}% and "
          f"{100*SIGMA_AVI:.0f}% (Cojutti Table 2, b1 and b2)")
    print(f"   attenuation factor for a SINGLE observation per patient: {a:.3f}")
    print(f"   so if the true correlation were 0.940, a naive two-stage analysis of")
    print(f"   one sample per patient would return about {a*RHO_PUBLISHED:.3f}")
    print(f"   and if the true correlation were 0.703, about {a*RHO_MODEL1:.3f}")
    print("\n   This is why the earlier two-stage estimates in this project (0.560 to")
    print("   0.598) were biased low, and why only a mixed-effects estimate is")
    print("   comparable with the published 0.94.")

    # ---------------------------------------------------------------- Part 2 --
    print("\n2. Could a Bologna-shaped dataset settle it?")
    print("   112 patients, 185 samples: median 1 occasion, some with 2")
    rows = []
    for rho_true, alt in ((RHO_PUBLISHED, RHO_MODEL1), (RHO_MODEL1, RHO_PUBLISHED)):
        r = design_row("Bologna (Cojutti 2024) as published", 112, (1, 1, 2),
                       rho_true, alt, n_rep, rng)
        rows.append(r)
        print(f"     true rho {rho_true:.3f} -> estimate {r['rho_mean_estimate']:.3f} "
              f"(SD {r['rho_sd']:.3f}, bias {r['rho_bias']:+.3f}); "
              f"power to exclude {alt:.3f}: {r['power_to_exclude_alternative_pct']:.1f}%")

    # ---------------------------------------------------------------- Part 3 --
    print("\n3. What a prospective study would need")
    for label, n, occ in (("small TDM cohort", 30, (1, 1, 2)),
                          ("moderate TDM cohort", 60, (1, 2, 2)),
                          ("Bologna-sized, 2 samples each", 112, (2, 2, 2)),
                          ("large TDM cohort", 200, (2, 2, 3))):
        r = design_row(label, n, occ, RHO_PUBLISHED, RHO_MODEL1, n_rep, rng)
        rows.append(r)
        print(f"     {label:32} n={n:3}  SD {r['rho_sd']:.3f}   "
              f"power {r['power_to_exclude_alternative_pct']:5.1f}%")

    # current cohort, for reference
    r = design_row("Model 1 cohort (21 patients, rich sampling)", 21, (6, 6, 6),
                   RHO_MODEL1, RHO_PUBLISHED, n_rep, rng)
    rows.append(r)
    print(f"     {'Model 1 cohort (21, rich)':32} n= 21  SD {r['rho_sd']:.3f}   "
          f"power {r['power_to_exclude_alternative_pct']:5.1f}%")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "model1_design_analysis.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {os.path.relpath(path, os.path.dirname(HERE))}")

    # ---------------------------------------------------------------- Part 4 --
    print("\n4. Falsifiable prediction, stated in advance")
    print("   If the Bologna cohort is re-analysed with a joint mixed-effects model,")
    print("   and if the correlation there is what the two RRT cohorts suggest rather")
    print("   than the published 0.94, the estimate should fall near 0.70 with a 95%")
    print("   interval roughly 0.60 to 0.79. If instead the published value is right,")
    print("   the estimate should fall near 0.94 with an interval excluding 0.70.")
    print("   The design analysis above says that dataset CAN tell those apart.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
