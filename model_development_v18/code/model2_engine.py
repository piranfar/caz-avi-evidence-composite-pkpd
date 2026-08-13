"""MODEL 2 engine — Hierarchical Uncertainty-integrated Joint Attainment Model (HU-JAM).

This module is the simulation core only. The decision-analytic outputs live in
`model2_hujam.py`; the value-of-information analysis in `model2_voi.py`.

WHAT THIS IS
    A layer over the primary simulation, not a replacement for it. The primary model
    treats every input as a known constant and perturbs them one at a time. This
    engine treats them as random variables carrying their PUBLISHED uncertainty and
    propagates them jointly to the decision.

WHAT IT IS NOT
    Not a pharmacokinetic model. It adds no pharmacokinetic knowledge and makes no
    pharmacokinetic claim.

REGRESSION TEST — the reason this module can be trusted
    With every uncertainty set to zero, rho at 0.94 and the avibactam target at
    4 mg/L, this engine must reproduce the frozen v16 output `primary_pta_results.csv`
    EXACTLY. `verify_against_frozen()` performs that check. If it fails, the wrapper
    is wrong, not the model.

    That requires reproducing the original random-number call order exactly: per
    renal class, uniform renal function first, then the correlated random effects,
    class by class. The draw is deliberately written to match `cazavi_analyses.py`
    in the v16 package rather than to be tidier.

TWO LAYOUTS THIS MODULE HAS TO FIND THE PRIMARY MODEL IN
    This package was developed inside the full local project directory, where the
    primary model lives at `revision_support/reproduce_primary_run.py` and the frozen
    reference at `revision_support/outputs/primary_pta_results.csv` -- both exact
    zero-tolerance matches, verified repeatedly throughout this project.

    It is ALSO distributed inside `caz-avi-evidence-composite-pkpd` on GitHub, whose
    layout is different: the same model file lives at `src/cazavi/reproduce_primary_run.py`
    (confirmed byte-identical to the local copy). That repository ships TWO copies of
    the primary PTA table: `data/reference/primary_pta_results.csv`, which was checked
    against a fresh run of its own bundled code and found to differ by up to 0.355
    percentage points (a pre-existing staleness in that one file, not something
    introduced here, and not the file that repo's own README points readers to for its
    figures/tables); and `data/processed/primary_pta_results.csv`, the results-facing
    copy, which matches a fresh run of the same code to floating-point exactness (0.0 pp
    on every metric, every regimen). `_locate_primary_model()` therefore points the
    GitHub layout at `data/processed/`, and `verify_against_frozen()` uses the same
    exact-zero tolerance in both layouts -- there is no fudge factor left. The stale
    `data/reference/` copy is a separate, low-severity cleanup item on that repository's
    own `main` branch; it is not read by anything in this package.

PROVENANCE OF THE UNCERTAINTY
    Every relative standard error below is published in Cojutti 2024 Table 2 and is
    ALREADY in the manuscript's source. The primary model discards it. Nothing here
    is invented; where a distributional form had to be chosen, the choice is stated
    in `UNCERTAINTY_NOTES` and is a scenario assumption, not a finding.
"""
from __future__ import annotations

import csv
import os
import sys
from dataclasses import dataclass, replace

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))


def _locate_primary_model():
    """Find reproduce_primary_run.py and its frozen reference, in whichever of the
    two known layouts (see module docstring) this file happens to be sitting in.

    Returns (source_dir, frozen_pta_path, tolerance_pp, schema). `schema` is
    "lower" for the revision_support/cazavi_analyses.py column names this module
    itself uses natively, or "title" for reproduce_primary_run.py's own
    Title-Case column names, so `verify_against_frozen` can read whichever it
    actually finds without guessing.
    """
    local_src = os.path.abspath(os.path.join(HERE, "..", "..", "revision_support"))
    local_pta = os.path.join(local_src, "outputs", "primary_pta_results.csv")
    if os.path.isfile(os.path.join(local_src, "reproduce_primary_run.py")):
        return local_src, local_pta, 0.0, "lower"

    repo_root = os.path.abspath(os.path.join(HERE, "..", ".."))
    gh_src = os.path.join(repo_root, "src", "cazavi")
    gh_pta = os.path.join(repo_root, "data", "processed", "primary_pta_results.csv")
    if os.path.isfile(os.path.join(gh_src, "reproduce_primary_run.py")):
        return gh_src, gh_pta, 0.0, "lower"

    raise FileNotFoundError(
        "Could not find reproduce_primary_run.py in either known layout. Looked in:\n"
        f"  {local_src}\n  {gh_src}\n"
        "This package depends on the primary model from the parent repository; "
        "see SOFTWARE.md.")


V16, FROZEN_PTA, _FROZEN_TOLERANCE, _FROZEN_SCHEMA = _locate_primary_model()
if V16 not in sys.path:
    sys.path.insert(0, V16)

import reproduce_primary_run as P   # noqa: E402  read-only import from the v16 package

MIC_INPUT = os.path.join(HERE, "recovered_inputs", "mic_distributions.csv")

SELECTED_REGIMENS = ("R1", "R8", "R10", "R12", "R13")
CLASS_REGIMENS = {}
for _r, (_c, _d, _i) in P.REGIMENS.items():
    CLASS_REGIMENS.setdefault(_c, []).append(_r)

# Renal-function weights for population-level CFR (Supplementary Table S3c).
CLASS_WEIGHTS = {"0–30": 0.15, "31–60": 0.159524, "61–90": 0.178571,
                 "91–120": 0.309202, "121–150": 0.202703}

UNCERTAINTY_NOTES = {
    "cl0_caz": "lognormal; RSE 6.36% (Cojutti 2024 Table 2)",
    "exp_caz": "normal; RSE 14.0% (Cojutti 2024 Table 2)",
    "cl0_avi": "lognormal; RSE 7.4% (Cojutti 2024 Table 2)",
    "exp_avi": "normal; RSE 12.3% (Cojutti 2024 Table 2)",
    "omega_caz": "lognormal on omega; RSE 33.3% reported on the variability term. "
                 "Whether the reported RSE applies to omega or to omega-squared is "
                 "not stated in the source; omega is assumed. SCENARIO ASSUMPTION.",
    "omega_avi": "lognormal on omega; RSE 30.1%. Same assumption as above.",
    "fu_caz": "uniform(0.80, 0.90); no uncertainty is published, so the existing "
              "PSA scenario bounds are reused. SCENARIO ASSUMPTION, not a standard error.",
    "fu_avi": "uniform(0.87, 0.97); as above.",
    "rho": "scenario set, see model2_hujam.RHO_SCENARIOS. Never pooled across populations.",
    "tox_threshold": "fixed at 104 mg/L; the exposure screen is an operational "
                     "ceiling, not a validated toxicity threshold.",
}


# ------------------------------------------------------------------ parameters

@dataclass(frozen=True)
class Params:
    """One draw of the uncertain parameters. Defaults are the point estimates."""
    cl0_caz: float = P.CL0_CAZ
    exp_caz: float = P.EXP_CAZ
    cl0_avi: float = P.CL0_AVI
    exp_avi: float = P.EXP_AVI
    omega_caz: float = P.OMEGA_CAZ
    omega_avi: float = P.OMEGA_AVI
    rho: float = P.RHO
    fu_caz: float = P.FU_CAZ
    fu_avi: float = P.FU_AVI
    avi_target: float = P.AVI_CT
    caz_target: float = P.CAZ_TARGET
    tox_threshold: float = P.TOX_THRESHOLD


BASE = Params()

RSE = {"cl0_caz": 0.0636, "exp_caz": 0.140, "cl0_avi": 0.074, "exp_avi": 0.123,
       "omega_caz": 0.333, "omega_avi": 0.301}


def draw_parameters(rng, base=BASE, scale=1.0, fu_uncertainty=True, tau_between=0.0):
    """Sample the pharmacokinetic parameters from their published uncertainty.

    `scale` multiplies every standard error, so scale=0 returns the point estimates
    exactly. That is what makes the regression test possible.

    `tau_between` is LAYER 2: the between-study standard deviation on the log
    clearance scale, applied as a single study-level effect SHARED by both analytes,
    because a model that predicts high clearance predicts it for both. Default 0.0
    leaves the engine's verified behaviour untouched. See `model2_heterogeneity.py`
    for the estimate and for why it is a scenario rather than a result.
    """
    if scale == 0.0 and tau_between == 0.0:
        return base
    study = np.exp(rng.normal(0.0, tau_between)) if tau_between > 0 else 1.0
    if scale == 0.0:
        return replace(base, cl0_caz=base.cl0_caz * study,
                       cl0_avi=base.cl0_avi * study)
    return replace(
        base,
        cl0_caz=base.cl0_caz * study * np.exp(rng.normal(0, RSE["cl0_caz"] * scale)),
        exp_caz=base.exp_caz + rng.normal(0, base.exp_caz * RSE["exp_caz"] * scale),
        cl0_avi=base.cl0_avi * study * np.exp(rng.normal(0, RSE["cl0_avi"] * scale)),
        exp_avi=base.exp_avi + rng.normal(0, base.exp_avi * RSE["exp_avi"] * scale),
        omega_caz=base.omega_caz * np.exp(rng.normal(0, RSE["omega_caz"] * scale)),
        omega_avi=base.omega_avi * np.exp(rng.normal(0, RSE["omega_avi"] * scale)),
        fu_caz=rng.uniform(0.80, 0.90) if fu_uncertainty else base.fu_caz,
        fu_avi=rng.uniform(0.87, 0.97) if fu_uncertainty else base.fu_avi,
    )


# ------------------------------------------------------------------ population

def draw_population(n_per_class, seed):
    """Renal function and whitened standard-normal deviates, per renal class.

    Deliberately mirrors `cazavi_analyses.draw_population` in the v16 package,
    including the order of random-number calls, so that the regression test against
    the frozen outputs is exact. Do not "tidy" this function.
    """
    rng = np.random.default_rng(seed)
    base_l = _cholesky(P.OMEGA_CAZ, P.OMEGA_AVI, P.RHO)
    cov = np.array([[P.OMEGA_CAZ ** 2, P.RHO * P.OMEGA_CAZ * P.OMEGA_AVI],
                    [P.RHO * P.OMEGA_CAZ * P.OMEGA_AVI, P.OMEGA_AVI ** 2]])
    population = {}
    for cls, (lo, hi) in P.EKFC_CLASSES.items():
        ekfc = rng.uniform(lo, hi, n_per_class)
        eta = rng.multivariate_normal(np.zeros(2), cov, size=n_per_class)
        z = np.linalg.solve(base_l, eta.T).T
        population[cls] = (ekfc, z)
    return population


def _cholesky(omega_caz, omega_avi, rho):
    rho = float(np.clip(rho, -0.999, 0.999))
    return np.array([[omega_caz, 0.0],
                     [rho * omega_avi, omega_avi * np.sqrt(1.0 - rho ** 2)]])


def clearances(ekfc, z, pr: Params):
    l = _cholesky(pr.omega_caz, pr.omega_avi, pr.rho)
    eta = z @ l.T
    cl_caz = pr.cl0_caz * (ekfc / P.EKFC_REF) ** pr.exp_caz * np.exp(eta[:, 0])
    cl_avi = pr.cl0_avi * (ekfc / P.EKFC_REF) ** pr.exp_avi * np.exp(eta[:, 1])
    return cl_caz, cl_avi


# ------------------------------------------------------------------ evaluation

def evaluate(population, pr: Params = BASE, regimens=None):
    """Attainment for every regimen at every MIC, under one parameter draw.

    Returns a dict keyed by regimen holding arrays over the MIC grid, plus the
    MIC-independent avibactam attainment and the exposure-screen exceedance.
    """
    regimens = regimens or list(P.REGIMENS)
    mics = np.asarray(P.MIC_GRID, float)
    out = {}
    for reg in regimens:
        cls, dose_g, interval_h = P.REGIMENS[reg]
        ekfc, z = population[cls]
        cl_caz, cl_avi = clearances(ekfc, z, pr)
        css_caz = dose_g * 1000.0 * P.CAZ_FRACTION / interval_h / cl_caz
        css_avi = dose_g * 1000.0 * P.AVI_FRACTION / interval_h / cl_avi
        free_caz = css_caz * pr.fu_caz
        free_avi = css_avi * pr.fu_avi

        avi_ok = free_avi >= pr.avi_target
        caz_ok = (free_caz[:, None] / mics[None, :]) >= pr.caz_target
        out[reg] = {
            "ekfc_class": cls,
            "caz_pta": 100.0 * caz_ok.mean(axis=0),
            "avi_pta": 100.0 * float(avi_ok.mean()),
            "joint_pta": 100.0 * (caz_ok & avi_ok[:, None]).mean(axis=0),
            "exceedance": 100.0 * float(np.mean(css_caz > pr.tox_threshold)),
            "daily_g": dose_g * 24.0 / interval_h,
        }
    return out


# ------------------------------------------------------------ MIC distributions

def load_mic_distributions(path=MIC_INPUT):
    """MIC weights, aligned to the simulated grid.

    NOTE: this file is RECOVERED, not an original project input. Two distributions
    are documented in Supplementary Table S3b; the other two were recovered by
    inversion. See `code/recovered_inputs/` and REPRODUCTION_CHECK.md.
    """
    mics = np.asarray(P.MIC_GRID, float)
    dists = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            d = dists.setdefault(r["distribution_id"],
                                 {"population": r["population"], "role": r["role"],
                                  "w": np.zeros(len(mics))})
            k = int(np.argmin(np.abs(mics - float(r["mic_mg_l"]))))
            d["w"][k] += float(r["frequency"])
    for did, d in dists.items():
        total = d["w"].sum()
        if abs(total - 1.0) > 1e-9:
            d["w"] /= total
    return dists


def cfr(result, weights):
    """MIC-weighted cumulative fraction of response for one regimen."""
    return float(result["joint_pta"] @ weights)


# --------------------------------------------------------------- verification

# Column names differ between the two layouts' frozen tables (see module docstring
# and _locate_primary_model). Both map to the same four quantities.
_FROZEN_COLUMNS = {
    "lower": dict(key_reg="regimen", key_mic="mic_mg_l", caz_pta="caz_pta_pct",
                  avi_pta="avi_attainment_pct", joint_pta="joint_pta_pct",
                  exceedance="toxicity_pct"),
    "title": dict(key_reg="Regimen", key_mic="MIC (mg/L)", caz_pta="CAZ PTA (%)",
                  avi_pta="AVI attainment (%)", joint_pta="Joint PTA (%)",
                  exceedance="CAZ toxicity (%)"),
}


def verify_against_frozen(tol=None, verbose=True):
    """With zero uncertainty, reproduce the frozen primary PTA table.

    `tol` defaults to `_FROZEN_TOLERANCE`, resolved once at import time by
    `_locate_primary_model`: exactly 0.0 in both known layouts -- see the module
    docstring for why the GitHub layout reads `data/processed/`, not
    `data/reference/`.
    """
    tol = _FROZEN_TOLERANCE if tol is None else tol
    cols = _FROZEN_COLUMNS[_FROZEN_SCHEMA]
    pop = draw_population(P.N_PER_CLASS, P.PRIMARY_SEED)
    res = evaluate(pop, BASE)
    mics = np.asarray(P.MIC_GRID, float)

    frozen = {}
    raw = open(FROZEN_PTA, "rb").read()
    for enc in ("utf-8", "cp1252"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    import io
    for r in csv.DictReader(io.StringIO(text)):
        frozen[(r[cols["key_reg"]], float(r[cols["key_mic"]]))] = r

    worst = {"caz_pta": 0.0, "avi_pta": 0.0, "joint_pta": 0.0, "exceedance": 0.0}
    n = 0
    for reg, d in res.items():
        for k, mic in enumerate(mics):
            f = frozen[(reg, float(mic))]
            n += 1
            worst["caz_pta"] = max(worst["caz_pta"],
                                   abs(d["caz_pta"][k] - float(f[cols["caz_pta"]])))
            worst["avi_pta"] = max(worst["avi_pta"],
                                   abs(d["avi_pta"] - float(f[cols["avi_pta"]])))
            worst["joint_pta"] = max(worst["joint_pta"],
                                     abs(d["joint_pta"][k] - float(f[cols["joint_pta"]])))
            worst["exceedance"] = max(worst["exceedance"],
                                      abs(d["exceedance"] - float(f[cols["exceedance"]])))
    ok = all(v <= tol for v in worst.values())
    if verbose:
        print(f"  REGRESSION TEST against the frozen primary PTA table (layout: "
              f"{'local development' if _FROZEN_SCHEMA == 'lower' else 'GitHub package'}, "
              f"tolerance {tol:g} pp)")
        print(f"    source: {FROZEN_PTA}")
        print(f"    rows compared {n}")
        for k, v in worst.items():
            print(f"    max |delta| {k:12} {v:.3e} pp")
        print(f"    {'PASS' if ok else 'FAIL'}")
    return ok, worst


if __name__ == "__main__":
    print("=" * 74)
    print("MODEL 2 ENGINE — verification")
    print("=" * 74)
    ok, _ = verify_against_frozen()
    d = load_mic_distributions()
    print(f"\n  MIC distributions loaded: {len(d)}")
    for did, x in sorted(d.items()):
        print(f"    {did:28} weights sum {x['w'].sum():.10f}  role {x['role']}")
    raise SystemExit(0 if ok else 1)
