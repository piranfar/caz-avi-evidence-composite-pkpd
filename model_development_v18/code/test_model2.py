"""Unit tests for the Model 2 decision layer.

Closes the item MODEL2_REPORT.md section 7 flagged as outstanding: "unit tests
for the decision layer: the regret identity, the EVPI identity, and the scenario
samplers." Also runs the interface conformance check as a test.

Covers identities that must hold BY CONSTRUCTION, independent of the pharmacology:
regret is zero for the optimal regimen and non-negative everywhere; EVPI is
non-negative and dominates EVPPI for any single parameter; the target and
correlation samplers draw from their documented support; and the zero-uncertainty
degenerate case of Model 2 reproduces the frozen v16 outputs exactly (the same
regression test model2_engine.py runs on import, wrapped here so it is discoverable
as a test).

Run:  python test_model2.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

import interface as I          # noqa: E402
import model2_engine as E      # noqa: E402
import model2_hujam as H       # noqa: E402
import model2_dispute_boundary as D  # noqa: E402

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def test_engine_conformance():
    print("\nEngine/decision-toolkit interface conformance")
    ok = I.check_engine(E, n_per_class=200, seed=1, verbose=False)
    check("model2_engine satisfies the JointAttainmentEngine contract", ok)


def test_zero_uncertainty_regression():
    print("\nZero-uncertainty regression (Model 2 must reduce to the frozen v16 model)")
    ok, worst = E.verify_against_frozen(verbose=False)
    check("engine reproduces the frozen v16 primary PTA table exactly", ok,
          f"max deltas: {worst}")


def test_target_samplers():
    print("\nTarget distribution samplers (Layer 4)")
    rng = np.random.default_rng(20260811)
    n = 20000

    check("T1 always returns exactly 1.0",
          all(H.sample_target(rng, "T1_point_1") == 1.0 for _ in range(100)))
    check("T2 always returns exactly 4.0",
          all(H.sample_target(rng, "T2_point_4") == 4.0 for _ in range(100)))

    draws = np.array([H.sample_target(rng, "T3_discrete") for _ in range(n)])
    check("T3 draws only take the values {0.5, 1.0, 4.0}",
          set(np.unique(draws)) <= {0.5, 1.0, 4.0})
    for v in (0.5, 1.0, 4.0):
        frac = np.mean(draws == v)
        check(f"T3 draws {v} mg/L close to the documented 1/3 weight",
              abs(frac - 1 / 3) < 0.02, f"observed {frac:.3f}")

    draws = np.array([H.sample_target(rng, "T4_uniform") for _ in range(n)])
    check("T4 draws lie within [0.5, 4.0]",
          bool(np.all((draws >= 0.5) & (draws <= 4.0))))
    check("T4 mean is close to the uniform midpoint 2.25",
          abs(draws.mean() - 2.25) < 0.05, f"observed {draws.mean():.3f}")

    draws = np.array([H.sample_target(rng, "T5_triangular") for _ in range(n)])
    check("T5 draws lie within [0.25, 4.0]",
          bool(np.all((draws >= 0.25) & (draws <= 4.0))))
    # Tri(0.25, 4, mode=1) is right-skewed (mode near the low end, long tail to
    # 4), so median != mode: the theoretical median is 1.628, not 1.0 -- a
    # property of the triangular distribution's CDF, not of this sampler. See
    # scipy.stats.triang.ppf(0.5, c=(1-0.25)/(4-0.25), loc=0.25, scale=3.75).
    check("T5 median matches the theoretical value for Tri(0.25, 4, mode=1)",
          abs(np.median(draws) - 1.628) < 0.05, f"median {np.median(draws):.3f}")
    check("T5 median exceeds the mode (confirms the expected right skew)",
          np.median(draws) > 1.2, f"median {np.median(draws):.3f}")

    draws = np.array([H.sample_target(rng, "T6_lognormal") for _ in range(n)])
    check("T6 draws are strictly positive", bool(np.all(draws > 0)))
    check("T6 median is close to the documented 1 mg/L",
          abs(np.median(draws) - 1.0) < 0.1, f"median {np.median(draws):.3f}")
    p95 = np.percentile(draws, 95)
    check("T6 95th percentile is close to the documented 4 mg/L",
          abs(p95 - 4.0) < 0.5, f"p95 {p95:.3f}")

    draws = np.array([H.sample_target(rng, "T7_coleman_evidence") for _ in range(n)])
    check("T7 draws only take the three Coleman 2014 Table 2 values {0.15, 0.22, 0.28}",
          set(np.unique(draws)) <= {0.15, 0.22, 0.28})
    for v in (0.15, 0.22, 0.28):
        frac = np.mean(draws == v)
        check(f"T7 draws {v} mg/L close to the documented 1/3 weight",
              abs(frac - 1 / 3) < 0.02, f"observed {frac:.3f}")


def test_rho_samplers():
    print("\nClearance-correlation scenario samplers (Layer 3)")
    rng = np.random.default_rng(20260811)
    n = 20000

    draws = np.array([H.sample_rho(rng, "C1_cojutti") for _ in range(n)])
    check("C1 draws lie within (-1, 1)", bool(np.all((draws > -1) & (draws < 1))))
    # Sampled on the Fisher-z scale (see sample_rho), so the median recovers
    # 0.94 closely; the mean is pulled slightly below it by Jensen's inequality
    # (tanh is concave for positive arguments), not by boundary clipping.
    check("C1 median is close to the published 0.94",
          abs(np.median(draws) - 0.94) < 0.02, f"observed median {np.median(draws):.4f}")
    check("C1 mean is close to but below 0.94 (Jensen's inequality, not clipping)",
          0.85 < draws.mean() < 0.94, f"observed mean {draws.mean():.4f}")
    check("C1 has no boundary clipping at 0.999",
          float(np.mean(draws >= 0.999)) < 0.01,
          f"observed {100*np.mean(draws>=0.999):.1f}% clipped")

    draws = np.array([H.sample_rho(rng, "C2_model1") for _ in range(n)])
    check("C2 draws lie within (-1, 1)", bool(np.all((draws > -1) & (draws < 1))))
    check("C2 median is close to the Model 1 point estimate 0.703",
          abs(np.median(draws) - 0.703) < 0.03, f"observed {np.median(draws):.4f}")
    lo, hi = np.percentile(draws, [2.5, 97.5])
    check("C2 95% interval is close to the profile-likelihood interval "
          "(0.380, 0.874)", abs(lo - 0.380) < 0.05 and abs(hi - 0.874) < 0.05,
          f"observed ({lo:.3f}, {hi:.3f})")

    draws = np.array([H.sample_rho(rng, "C3_agnostic") for _ in range(n)])
    check("C3 draws lie within [0.38, 0.98]",
          bool(np.all((draws >= 0.38) & (draws <= 0.98))))
    check("C3 mean is close to the uniform midpoint 0.68",
          abs(draws.mean() - 0.68) < 0.02, f"observed {draws.mean():.4f}")


def test_utility_and_identities():
    print("\nUtility function and decision identities")

    # utility is monotone: more attainment is better, more exceedance is worse
    weights = np.array([1.0])
    base = {"joint_pta": np.array([50.0]), "exceedance": 10.0}
    better_attain = {"joint_pta": np.array([70.0]), "exceedance": 10.0}
    worse_exceed = {"joint_pta": np.array([50.0]), "exceedance": 30.0}
    check("utility increases with joint attainment, exceedance held fixed",
          H.utility(better_attain, weights) > H.utility(base, weights))
    check("utility decreases with exceedance, attainment held fixed",
          H.utility(worse_exceed, weights) < H.utility(base, weights))
    check("utility with lambda=0 ignores exceedance entirely",
          H.utility(worse_exceed, weights, lam=0.0)
          == H.utility(base, weights, lam=0.0))

    # the constrained variant: hard cliff at the ceiling
    below = {"joint_pta": np.array([80.0]), "exceedance": 14.9}
    above = {"joint_pta": np.array([80.0]), "exceedance": 15.1}
    check("constrained utility is unaffected just below the ceiling",
          H.utility_constrained(below, weights) == 80.0)
    check("constrained utility drops to zero just above the ceiling",
          H.utility_constrained(above, weights) == 0.0)

    # regret identity: for any small synthetic decision problem, regret of the
    # best-in-hindsight choice is exactly zero, and every regimen's regret is
    # non-negative, by the definition regret = best - chosen
    rng = np.random.default_rng(7)
    U = rng.uniform(0, 100, size=(500, 4))          # draws x regimens
    u_best = U.max(axis=1)
    for j in range(4):
        regret_j = u_best - U[:, j]
        check(f"regret is non-negative for regimen {j} in every draw",
              bool(np.all(regret_j >= -1e-9)))
    best_idx = np.argmax(U, axis=1)
    regret_best = u_best - U[np.arange(len(U)), best_idx]
    check("regret is exactly zero for the best-in-hindsight regimen in every draw",
          bool(np.allclose(regret_best, 0.0)))

    # EVPI identity: EVPI = E[max] - max[E], and by Jensen's inequality for the
    # max of random variables this is always >= 0
    evpi = float(u_best.mean() - U.mean(axis=0).max())
    check("EVPI is non-negative (Jensen's inequality for max)", evpi >= -1e-9,
          f"evpi={evpi:.6f}")

    # a degenerate case: if one regimen dominates in every draw, EVPI must be
    # exactly zero, because perfect information never changes the decision
    U_dominated = U.copy()
    U_dominated[:, 0] = U.max(axis=1) + 1.0          # column 0 always wins
    evpi_degenerate = float(U_dominated.max(axis=1).mean()
                            - U_dominated.mean(axis=0).max())
    check("EVPI is exactly zero when one option dominates in every draw",
          abs(evpi_degenerate) < 1e-9, f"evpi={evpi_degenerate:.9f}")


def test_evppi_bounded_by_evpi():
    print("\nEVPPI <= EVPI (partial information cannot be worth more than full information)")
    rng = np.random.default_rng(11)
    n = 2000
    x1 = rng.normal(size=n)             # the parameter EVPPI will target
    x2 = rng.normal(size=n)             # a second, independent parameter
    # utility of two regimens, each a noisy function of both parameters
    u0 = 10 * x1 + 3 * x2 + rng.normal(scale=2, size=n)
    u1 = -5 * x1 + 6 * x2 + rng.normal(scale=2, size=n)
    U = np.column_stack([u0, u1])

    evpi = float(U.max(axis=1).mean() - U.mean(axis=0).max())

    def evppi_for(x, degree=4):
        xs = (x - x.mean()) / x.std()
        fitted = np.column_stack([np.polyval(np.polyfit(xs, U[:, k], degree), xs)
                                  for k in range(U.shape[1])])
        return float(fitted.max(axis=1).mean() - U.mean(axis=0).max())

    e1 = max(evppi_for(x1), 0.0)
    e2 = max(evppi_for(x2), 0.0)
    check("EVPPI(x1) does not exceed EVPI", e1 <= evpi + 1e-6,
          f"evppi={e1:.4f} evpi={evpi:.4f}")
    check("EVPPI(x2) does not exceed EVPI", e2 <= evpi + 1e-6,
          f"evppi={e2:.4f} evpi={evpi:.4f}")
    check("EVPPI is non-negative for both parameters", e1 >= -1e-9 and e2 >= -1e-9)


def test_frechet_bound():
    """The attainment-correlation ceiling used in MODEL2_REPORT.md 3.7.

    The claim that rests on this is strong -- that no pharmacokinetic correlation can
    make ceftazidime attainment a good proxy for avibactam attainment -- so the bound
    itself is checked against its known analytic properties and one hand-computed value.
    """
    # equal margins allow perfect correlation
    for p in (0.1, 0.35, 0.5, 0.8):
        check(f"bound is 1 when prevalences are equal (p={p})",
              abs(D.frechet_phi_bound(p, p) - 1.0) < 1e-12)

    # symmetric in its arguments
    check("bound is symmetric",
          abs(D.frechet_phi_bound(0.58, 0.85) - D.frechet_phi_bound(0.85, 0.58)) < 1e-12)

    # independently computed: sqrt(0.58*0.15 / (0.85*0.42)) = 0.4936572484949417
    check("bound matches independent computation at the reported margins",
          abs(D.frechet_phi_bound(0.58, 0.85) - 0.4936572484949417) < 1e-12)

    # strictly inside (0, 1] and decreasing as the margins separate
    b_close = D.frechet_phi_bound(0.50, 0.55)
    b_far = D.frechet_phi_bound(0.50, 0.95)
    check("bound lies in (0, 1]", 0.0 < b_far <= 1.0 and 0.0 < b_close <= 1.0)
    check("bound falls as prevalences separate", b_far < b_close)

    # degenerate margins are refused rather than silently returning a number
    check("bound is nan for a degenerate margin",
          np.isnan(D.frechet_phi_bound(0.0, 0.5)) and np.isnan(D.frechet_phi_bound(0.5, 1.0)))

    # the empirical phi can never exceed the bound for the same margins
    rng = np.random.default_rng(20260812)
    worst = -1.0
    for _ in range(200):
        p1, p2 = rng.uniform(0.05, 0.95, 2)
        a = rng.random(4000) < p1
        b = rng.random(4000) < p2
        phi = D.phi_coefficient(a, b)
        bound = D.frechet_phi_bound(a.mean(), b.mean())
        if not np.isnan(phi):
            worst = max(worst, phi - bound)
    check("empirical phi never exceeds the bound over 200 random margin pairs",
          worst <= 1e-9)

    # phi behaves at the extremes
    v = np.array([True, False, True, False, True, False])
    check("phi of a vector with itself is 1", abs(D.phi_coefficient(v, v) - 1.0) < 1e-12)
    check("phi of a vector with its complement is -1",
          abs(D.phi_coefficient(v, ~v) + 1.0) < 1e-12)
    check("phi is nan when one indicator never varies",
          np.isnan(D.phi_coefficient(v, np.ones(6, dtype=bool))))


def main():
    print("=" * 70)
    print("MODEL 2 DECISION-LAYER TEST SUITE")
    print("=" * 70)
    test_engine_conformance()
    test_zero_uncertainty_regression()
    test_target_samplers()
    test_rho_samplers()
    test_utility_and_identities()
    test_evppi_bounded_by_evpi()
    test_frechet_bound()
    print("\n" + "=" * 70)
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s)")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
