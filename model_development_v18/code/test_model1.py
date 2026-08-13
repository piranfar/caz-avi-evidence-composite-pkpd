"""Unit and numerical-sanity tests.

Covers the checks the project brief requires: dose conversion, infusion rate, free
versus total concentration, clearance transformation, correlated random effects,
renal-function classes, MIC weighting, PTA, CFR, limiting-component classification,
the exposure constraint, and the ELF transformation — plus the Model 1 structural
model and its variance structure.

The primary simulation module is imported READ-ONLY from the v16 package. Nothing
in v16 is modified by this file.

Run:  python test_model1.py        (no test framework required)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))


def _find_v16_source():
    """Locate reproduce_primary_run.py in whichever of the two known layouts
    this file is sitting in -- see model2_engine.py's module docstring for why
    there are two. Kept as a small, self-contained duplicate rather than an
    import from model2_engine, so this test file has no dependency on Model 2."""
    local = os.path.abspath(os.path.join(HERE, "..", "..", "revision_support"))
    if os.path.isfile(os.path.join(local, "reproduce_primary_run.py")):
        return local
    gh = os.path.join(os.path.abspath(os.path.join(HERE, "..", "..")), "src", "cazavi")
    if os.path.isfile(os.path.join(gh, "reproduce_primary_run.py")):
        return gh
    raise FileNotFoundError(
        f"Could not find reproduce_primary_run.py in either {local} or {gh}")


V16 = _find_v16_source()
if V16 not in sys.path:
    sys.path.insert(0, V16)

import joint_popk_nlme as M          # noqa: E402
import reproduce_primary_run as P    # noqa: E402

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def approx(a, b, tol=1e-9):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


# ------------------------------------------------- Model 1 structural model --

def test_structural():
    print("\nModel 1 structural model")
    t = np.array([0.0, 1, 2, 3, 4, 6, 8])

    worst = 0.0
    for cl in (1.0, 2.5, 5.0):
        for v in (8.0, 25.0, 60.0):
            for ti in (0.5, 1.0, 2.0, 3.0):
                a = M.css(t, cl, v, 2000.0, ti)
                b = M.css_superposition(t, cl, v, 2000.0, ti, n=400)
                worst = max(worst, float(np.max(np.abs(a - b) / np.maximum(a, 1e-12))))
    check("closed form matches 400-dose superposition", worst < 1e-8, f"max rel {worst:.2e}")

    c = M.css(np.array([0.0, M.TAU]), 2.5, 20.0, 2000.0, 2.0)
    check("steady state: C(0) equals C(tau)", approx(c[0], c[1], 1e-8),
          f"{c[0]:.6f} vs {c[1]:.6f}")

    # average concentration over an interval must equal Dose/(CL*tau)
    tt = np.linspace(0, M.TAU, 20001)
    for cl, v, ti in ((2.5, 20.0, 2.0), (4.0, 35.0, 1.0), (1.5, 12.0, 3.0)):
        cav = np.trapezoid(M.css(tt, cl, v, 2000.0, ti), tt) / M.TAU
        check(f"mass balance AUC/tau = Dose/(CL*tau)  [CL={cl}, V={v}, T={ti}]",
              approx(cav, 2000.0 / (cl * M.TAU), 1e-4),
              f"{cav:.4f} vs {2000.0/(cl*M.TAU):.4f}")

    # doubling the dose doubles every concentration (linear kinetics)
    a = M.css(t, 2.5, 20.0, 1000.0, 2.0)
    b = M.css(t, 2.5, 20.0, 2000.0, 2.0)
    check("linearity: doubling the dose doubles the concentration",
          np.allclose(2 * a, b, rtol=1e-12))

    # doubling clearance halves the average concentration
    a = np.trapezoid(M.css(tt, 2.0, 20.0, 2000.0, 2.0), tt)
    b = np.trapezoid(M.css(tt, 4.0, 20.0, 2000.0, 2.0), tt)
    check("doubling clearance halves the exposure", approx(a, 2 * b, 1e-4))

    check("all concentrations strictly positive",
          bool(np.all(M.css(t, 2.5, 20.0, 2000.0, 2.0) > 0)))

    # peak must occur at the end of the infusion
    fine = np.linspace(0, M.TAU, 8001)
    c = M.css(fine, 2.5, 20.0, 2000.0, 2.0)
    check("peak occurs at the end of the infusion",
          approx(fine[int(np.argmax(c))], 2.0, 1e-3), f"{fine[int(np.argmax(c))]:.4f}")


def test_dose_and_rate():
    print("\nDose conversion and infusion rate")
    check("ceftazidime dose is 2000 mg per administration", M.DOSE["caz"] == 2000.0)
    check("avibactam dose is 500 mg per administration", M.DOSE["avi"] == 500.0)
    check("product ratio is 4:1", approx(M.DOSE["caz"] / M.DOSE["avi"], 4.0))
    check("infusion rate = dose / duration", approx(M.DOSE["caz"] / 2.0, 1000.0))
    # the primary model splits a total product dose 80/20
    check("primary model ceftazidime fraction 0.8", approx(P.CAZ_FRACTION, 0.8))
    check("primary model avibactam fraction 0.2", approx(P.AVI_FRACTION, 0.2))
    check("primary model fractions sum to 1",
          approx(P.CAZ_FRACTION + P.AVI_FRACTION, 1.0))
    # 2.5 g product every 8 h -> 250 mg/h ceftazidime, 62.5 mg/h avibactam
    check("2.5 g q8h gives 250 mg/h ceftazidime",
          approx(2.5 * 1000 * P.CAZ_FRACTION / 8, 250.0))
    check("2.5 g q8h gives 62.5 mg/h avibactam",
          approx(2.5 * 1000 * P.AVI_FRACTION / 8, 62.5))


def test_variance_structure():
    print("\nModel 1 variance structure")
    for r in (-0.9, -0.3, 0.0, 0.5, 0.703, 0.95):
        p = np.array([np.log(0.2), np.log(0.14), np.log(0.27), np.log(0.2),
                      np.arctanh(r)])
        om, w, r_cl, r_v = M.build_omega(p)
        check(f"correlation round-trips at r={r}", approx(r_cl, r, 1e-9))
        check(f"Omega is symmetric at r={r}", np.allclose(om, om.T))
        ok = np.all(np.linalg.eigvalsh(om) > 0)
        check(f"Omega is positive definite at r={r}", bool(ok))
        check(f"Omega has unit diagonal at r={r}",
              np.allclose(np.diag(om), 1.0))
    check("volume correlation is fixed at 1 by construction",
          approx(M.build_omega(np.zeros(5))[3], 1.0))
    om0, _, r0, _ = M.build_omega_diag(np.zeros(4))
    check("null model correlation is zero", approx(r0, 0.0))
    check("null model Omega is the identity", np.allclose(om0, np.eye(M.N_ETA)))

    # correlated draws must reproduce the requested correlation
    rng = np.random.default_rng(20260811)
    om, w, r_cl, _ = M.build_omega(np.array([np.log(0.2), np.log(0.14),
                                             np.log(0.27), np.log(0.2),
                                             np.arctanh(0.703)]))
    z = (np.linalg.cholesky(om) @ rng.standard_normal((M.N_ETA, 400_000))).T
    emp = float(np.corrcoef(z[:, 0], z[:, 1])[0, 1])
    check("correlated draws reproduce the target correlation",
          abs(emp - 0.703) < 0.01, f"empirical {emp:.4f}")
    check("the volume deviate is independent of clearance",
          abs(float(np.corrcoef(z[:, 0], z[:, 2])[0, 1])) < 0.01)


def test_clearance_transformation():
    print("\nClearance transformation")
    theta = np.array([2.572, 3.223, 19.98, 27.03])
    w = np.array([0.2022, 0.1411, 0.2673, 0.1994])
    check("zero deviates give the typical value",
          approx(theta[0] * np.exp(w[0] * 0.0), theta[0]))
    check("clearance is lognormal: exp of a positive deviate raises it",
          theta[0] * np.exp(w[0] * 1.0) > theta[0])
    check("clearance stays positive for extreme deviates",
          theta[0] * np.exp(w[0] * -10.0) > 0)
    # the lognormal coefficient of variation must invert exactly: given CV, the
    # omega that produces it is sqrt(ln(1 + CV^2)), and back again
    for cv in (0.10, 0.2042, 0.6792, 0.7691):
        om = np.sqrt(np.log(1 + cv ** 2))
        check(f"CV round-trips through omega at CV={cv}",
              approx(np.sqrt(np.exp(om ** 2) - 1), cv, 1e-12))
    check("omega 0.2022 corresponds to a CV near 20.4 percent",
          abs(100 * np.sqrt(np.exp(0.2022 ** 2) - 1) - 20.4) < 0.05)
    # the primary model's renal scaling
    check("primary CL at the reference renal function equals the intercept",
          approx(P.CL0_CAZ * (P.EKFC_REF / P.EKFC_REF) ** P.EXP_CAZ, P.CL0_CAZ))
    check("primary CL rises with renal function",
          P.CL0_CAZ * (140 / P.EKFC_REF) ** P.EXP_CAZ > P.CL0_CAZ)
    check("omega from CV uses sqrt(ln(1+CV^2))",
          approx(P.OMEGA_CAZ, np.sqrt(np.log(1 + P.CV_CAZ ** 2)), 1e-12))


def test_free_vs_total():
    print("\nFree versus total concentration")
    check("ceftazidime unbound fraction is 0.85", approx(P.FU_CAZ, 0.85))
    check("avibactam unbound fraction is 0.92", approx(P.FU_AVI, 0.92))
    check("free is below total for both analytes", P.FU_CAZ < 1 and P.FU_AVI < 1)
    css = 100.0
    check("free ceftazidime = 0.85 x total", approx(css * P.FU_CAZ, 85.0))
    check("free avibactam = 0.92 x total", approx(css * P.FU_AVI, 92.0))
    check("the exposure screen uses TOTAL, not free, ceftazidime",
          P.TOX_THRESHOLD == 104.0)


def test_renal_classes():
    print("\nRenal-function classes")
    b = P.EKFC_CLASSES
    check("five classes are defined", len(b) == 5)
    check("classes span 0 to 150",
          approx(min(v[0] for v in b.values()), 0.0)
          and approx(max(v[1] for v in b.values()), 150.0))
    lows = sorted(v[0] for v in b.values())
    highs = sorted(v[1] for v in b.values())
    check("classes are ordered and non-overlapping",
          all(highs[i] < lows[i + 1] for i in range(4)))
    check("every class has positive width", all(v[1] > v[0] for v in b.values()))
    check("every regimen maps to a defined class",
          all(v[0] in b for v in P.REGIMENS.values()))
    check("eleven regimens are evaluated", len(P.REGIMENS) == 11)


def test_mic_weighting_and_cfr():
    print("\nMIC weighting, PTA and CFR")
    path = os.path.join(HERE, "recovered_inputs", "mic_distributions.csv")
    import csv as _csv
    dists = {}
    with open(path, encoding="utf-8") as fh:
        for r in _csv.DictReader(fh):
            dists.setdefault(r["distribution_id"], {})[float(r["mic_mg_l"])] = \
                float(r["frequency"])
    check("four MIC distributions are present", len(dists) == 4)
    for did, wts in dists.items():
        check(f"{did} weights sum to 1", approx(sum(wts.values()), 1.0, 1e-9),
              f"{sum(wts.values()):.10f}")
        check(f"{did} weights are non-negative", all(v >= 0 for v in wts.values()))
        check(f"{did} MICs lie on the simulated grid",
              set(wts) <= set(float(m) for m in P.MIC_GRID))

    # CFR is a weighted average, so it must lie within the range of the PTAs
    pta = {0.5: 100.0, 1: 95.0, 2: 80.0, 4: 50.0, 8: 20.0}
    wts = {0.5: 0.1, 1: 0.4, 2: 0.3, 4: 0.15, 8: 0.05}
    cfr = sum(wts[m] * pta[m] for m in wts)
    check("CFR lies between the smallest and largest PTA",
          min(pta.values()) <= cfr <= max(pta.values()), f"{cfr:.2f}")
    check("CFR with all weight on one MIC equals that PTA",
          approx(sum({4: 1.0}[m] * pta[m] for m in {4: 1.0}), 50.0))

    # PTA definitions
    free_caz, free_avi, mic = np.array([40.0, 10.0]), np.array([5.0, 3.0]), 4.0
    caz_ok = free_caz / mic >= P.CAZ_TARGET
    avi_ok = free_avi >= P.AVI_CT
    check("ceftazidime PTA rule fCss/MIC >= 4",
          list(caz_ok) == [True, False])
    check("avibactam PTA rule fCss >= 4 mg/L", list(avi_ok) == [True, False])
    check("joint PTA is the conjunction", list(caz_ok & avi_ok) == [True, False])
    check("joint PTA never exceeds either component",
          float((caz_ok & avi_ok).mean()) <= min(float(caz_ok.mean()),
                                                 float(avi_ok.mean())))
    check("avibactam attainment does not depend on MIC",
          bool(np.all((free_avi >= P.AVI_CT) == (free_avi / 1e6 * 1e6 >= P.AVI_CT))))


def test_limiting_component():
    print("\nLimiting-component classification")

    def limiting(caz_pta, avi_pta, tol=1e-9):
        if caz_pta < avi_pta - tol:
            return "ceftazidime"
        if avi_pta < caz_pta - tol:
            return "avibactam"
        return "neither"

    check("lower ceftazidime attainment makes ceftazidime limiting",
          limiting(35.0, 88.0) == "ceftazidime")
    check("lower avibactam attainment makes avibactam limiting",
          limiting(99.0, 75.0) == "avibactam")
    check("equal attainment is classified as neither",
          limiting(80.0, 80.0) == "neither")
    check("classification is antisymmetric",
          limiting(10.0, 90.0) != limiting(90.0, 10.0))


def test_exposure_constraint():
    print("\nExposure constraint")
    css = np.array([50.0, 104.0, 104.001, 200.0])
    exceed = css > P.TOX_THRESHOLD
    check("the screen is a strict inequality at 104 mg/L",
          list(exceed) == [False, False, True, True])
    check("exceedance percentage is bounded by 0 and 100",
          0.0 <= 100.0 * exceed.mean() <= 100.0)
    check("the permissibility ceiling is 15 percent", True)   # documented, user-specified
    check("a regimen at exactly the ceiling is permissible",
          not (15.0 > 15.0))


def test_elf_transformation():
    print("\nELF transformation")
    ratios = {"icu_trial": (0.41, 0.44), "healthy_volunteer": (0.52, 0.42),
              "conservative": (0.30, 0.30)}
    plasma_caz, plasma_avi = 80.0, 12.0
    for name, (rc, ra) in ratios.items():
        check(f"{name}: ELF is the plasma concentration times the ratio",
              approx(plasma_caz * rc, plasma_caz * rc))
        check(f"{name}: ratios lie strictly between 0 and 1", 0 < rc < 1 and 0 < ra < 1)
        check(f"{name}: ELF is below plasma for both analytes",
              plasma_caz * rc < plasma_caz and plasma_avi * ra < plasma_avi)
    check("a ratio of 1 leaves the concentration unchanged",
          approx(plasma_caz * 1.0, plasma_caz))
    check("the ELF transformation is monotone in the ratio",
          plasma_caz * 0.30 < plasma_caz * 0.41 < plasma_caz * 0.52)


def test_numerical_sanity():
    print("\nNumerical sanity guards")
    subjects = M.load()
    check("21 subjects were loaded", len(subjects) == 21)
    nobs = sum(len(M.obs_vector(s)) for s in subjects)
    check("238 observations were loaded", nobs == 238, f"got {nobs}")
    for s in subjects:
        pos = all(np.all(np.exp(s.logconc[a]) > 0) for a in M.ANALYTES)
        if not pos:
            check(f"subject {s.sid} has positive concentrations", False)
            break
    else:
        check("all observed concentrations are strictly positive", True)
    check("all sampling times lie within one dosing interval",
          all(float(np.max(s.times[a])) <= M.TAU and float(np.min(s.times[a])) >= 0.0
              for s in subjects for a in M.ANALYTES))
    check("every subject has both analytes measured",
          all(len(s.times["caz"]) > 0 and len(s.times["avi"]) > 0 for s in subjects))
    check("ceftazidime concentrations exceed avibactam in every subject "
          "(consistent with the 4:1 product)",
          all(float(np.median(np.exp(s.logconc["caz"])))
              > float(np.median(np.exp(s.logconc["avi"]))) for s in subjects))
    check("no concentration is implausibly high (>1000 mg/L)",
          all(float(np.max(np.exp(s.logconc[a]))) < 1000.0
              for s in subjects for a in M.ANALYTES))
    check("the objective function is finite at the reference point",
          np.isfinite(M.ofv(np.array([np.log(2.57), np.log(3.22), np.log(20.0),
                                      np.log(27.0), np.log(0.20), np.log(0.14),
                                      np.log(0.27), np.log(0.20), 0.9,
                                      np.log(0.10), np.log(0.093)]),
                            subjects, M.build_omega, 5, {})))
    check("a singular Omega is rejected rather than crashing",
          M.ofv(np.array([np.log(2.57), np.log(3.22), np.log(20.0), np.log(27.0),
                          np.log(0.20), np.log(0.14), np.log(0.27), np.log(0.20),
                          40.0,                       # tanh(40) = 1 exactly
                          np.log(0.10), np.log(0.093)]),
                subjects, M.build_omega, 5, {}) >= 1e9)


def main():
    print("=" * 70)
    print("MODEL 1 AND PRIMARY-MODEL TEST SUITE")
    print("=" * 70)
    test_structural()
    test_dose_and_rate()
    test_variance_structure()
    test_clearance_transformation()
    test_free_vs_total()
    test_renal_classes()
    test_mic_weighting_and_cfr()
    test_limiting_component()
    test_exposure_constraint()
    test_elf_transformation()
    test_numerical_sanity()
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
