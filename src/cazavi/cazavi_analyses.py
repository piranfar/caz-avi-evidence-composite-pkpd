"""Full analysis suite for the CAZ-AVI evidence-composite model.

Builds on the verified core in `reproduce_primary_run.py`: every parameter is a
published or scenario-defined value, nothing is fitted to the RC1 outputs. Each
analysis can be checked against the frozen RC1 reference tables with `--verify`.

Analyses
    primary      100,000-subject Monte Carlo run over the regimen x MIC grid
    cfr          MIC-weighted cumulative fraction of response, per distribution
    convergence  joint PTA stability from 1,000 to 100,000 subjects
    multiseed    five-seed reproducibility
    gsa          one-at-a-time sensitivity, paired on common random numbers
    psa          Latin-hypercube probabilistic sensitivity analysis

Common random numbers
    The population is drawn once as standard normal deviates `z` plus renal
    function `ekfc`. Scenario variants rebuild the random effects from the same
    `z`, so a scenario difference reflects the parameter change and not Monte
    Carlo noise. This is what makes the paired OAT comparisons meaningful.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass, replace

import numpy as np

from reproduce_primary_run import (
    AVI_CT, CAZ_FRACTION, CAZ_TARGET, CL0_AVI, CL0_CAZ, EKFC_CLASSES, EKFC_REF,
    EXP_AVI, EXP_CAZ, FU_AVI, FU_CAZ, MIC_GRID, N_PER_CLASS, OMEGA_AVI,
    OMEGA_CAZ, PRIMARY_SEED, REGIMENS, RHO, TOX_THRESHOLD, AVI_FRACTION,
)

HERE = os.path.dirname(os.path.abspath(__file__))


def _data_root():
    """Find the data directory, whichever layout this file is sitting in.

    The same module runs from the working folder and from `src/cazavi/` inside
    the repository, so the inputs and frozen reference tables are located by
    looking for them rather than by assuming a fixed relative path.
    """
    candidates = [
        os.path.join(os.path.dirname(os.path.dirname(HERE)), "data"),   # src/cazavi/ -> repo/data
        os.path.join(os.path.dirname(HERE), "data"),                    # src/ -> repo/data
        os.path.join(os.path.dirname(HERE), "CAZ_AVI_Local_First_Reconstruction_v1", "data"),
    ]
    for c in candidates:
        if os.path.isdir(os.path.join(c, "reference")) and os.path.isdir(os.path.join(c, "inputs")):
            return c
    raise FileNotFoundError(
        "Could not locate a data directory containing 'reference' and 'inputs'. Looked in:\n  "
        + "\n  ".join(candidates))


RECON = _data_root()
REFERENCE, INPUTS = os.path.join(RECON, "reference"), os.path.join(RECON, "inputs")

# When the pipeline sits in a repository the results belong in the repository's
# own folders, not beside the source; when it runs from a working directory they
# belong in ./outputs.
_IN_REPO = os.path.basename(HERE) == "cazavi" and \
    os.path.basename(os.path.dirname(HERE)) == "src"
_ROOT = os.path.dirname(os.path.dirname(HERE)) if _IN_REPO else HERE
DEFAULT_OUT = os.path.join(_ROOT, "data", "processed") if _IN_REPO \
    else os.path.join(HERE, "outputs")


# Best permissible regimen for each renal-function class: the decision set that
# carries the reported conclusions, and the scope of the GSA and PSA.
SELECTED_REGIMENS = ("R1", "R8", "R10", "R12", "R13")


@dataclass(frozen=True)
class Scenario:
    """A one-at-a-time or sampled perturbation of the baseline model."""
    name: str = "BASE"
    cl_caz_mult: float = 1.0
    cl_avi_mult: float = 1.0
    omega_caz_mult: float = 1.0
    omega_avi_mult: float = 1.0
    rho: float = RHO
    fu_caz: float = FU_CAZ
    fu_avi: float = FU_AVI
    caz_target: float = CAZ_TARGET
    avi_ct: float = AVI_CT
    tox_threshold: float = TOX_THRESHOLD


# The OAT scenario set reported in the manuscript (Methods 2.13).
GSA_SCENARIOS = [
    Scenario("BASE"),
    Scenario("CAZ_CL_LOW", cl_caz_mult=0.8), Scenario("CAZ_CL_HIGH", cl_caz_mult=1.2),
    Scenario("AVI_CL_LOW", cl_avi_mult=0.8), Scenario("AVI_CL_HIGH", cl_avi_mult=1.2),
    Scenario("IIV_LOW", omega_caz_mult=0.8, omega_avi_mult=0.8),
    Scenario("IIV_HIGH", omega_caz_mult=1.2, omega_avi_mult=1.2),
    Scenario("RHO_050", rho=0.50), Scenario("RHO_000", rho=0.0),
    Scenario("FU_CAZ_LOW", fu_caz=0.80), Scenario("FU_CAZ_HIGH", fu_caz=0.90),
    Scenario("FU_AVI_LOW", fu_avi=0.87), Scenario("FU_AVI_HIGH", fu_avi=0.97),
    Scenario("CAZ_TARGET_2", caz_target=2.0), Scenario("CAZ_TARGET_6", caz_target=6.0),
    Scenario("AVI_CT_2", avi_ct=2.0), Scenario("AVI_CT_6", avi_ct=6.0),
    Scenario("TOX_80", tox_threshold=80.0), Scenario("TOX_128", tox_threshold=128.0),
]


def _cholesky(omega_caz, omega_avi, rho):
    rho = float(np.clip(rho, -0.999, 0.999))
    return np.array([[omega_caz, 0.0],
                     [rho * omega_avi, omega_avi * np.sqrt(1.0 - rho**2)]])


def draw_population(n_per_class, seed):
    """Draw renal function and correlated random effects once, then whiten them.

    The baseline effects are drawn with `multivariate_normal` in the same order
    the original pipeline used (uniform renal function first, then the effects,
    class by class); this reproduces the frozen RC1 tables to within about
    0.06 percentage points. Whitening those draws back to standard normal `z`
    lets every scenario rebuild its own random effects from the identical
    underlying deviates, which is what makes the paired OAT comparisons isolate
    the parameter change rather than Monte Carlo noise. For the baseline
    scenario the round trip is exact.
    """
    rng = np.random.default_rng(seed)
    base_l = _cholesky(OMEGA_CAZ, OMEGA_AVI, RHO)
    cov = np.array([[OMEGA_CAZ**2, RHO * OMEGA_CAZ * OMEGA_AVI],
                    [RHO * OMEGA_CAZ * OMEGA_AVI, OMEGA_AVI**2]])
    population = {}
    for cls, (lo, hi) in EKFC_CLASSES.items():
        ekfc = rng.uniform(lo, hi, n_per_class)
        eta = rng.multivariate_normal(np.zeros(2), cov, size=n_per_class)
        z = np.linalg.solve(base_l, eta.T).T
        population[cls] = (ekfc, z)
    return population


def clearances(ekfc, z, sc: Scenario):
    """Apply a scenario to the shared draws and return individual clearances."""
    l = _cholesky(OMEGA_CAZ * sc.omega_caz_mult, OMEGA_AVI * sc.omega_avi_mult, sc.rho)
    eta = z @ l.T
    cl_caz = CL0_CAZ * (ekfc / EKFC_REF) ** EXP_CAZ * sc.cl_caz_mult * np.exp(eta[:, 0])
    cl_avi = CL0_AVI * (ekfc / EKFC_REF) ** EXP_AVI * sc.cl_avi_mult * np.exp(eta[:, 1])
    return cl_caz, cl_avi


def evaluate(population, sc: Scenario = Scenario()):
    """Return one row per regimen x MIC under the given scenario."""
    rows = []
    for regimen, (cls, dose_g, interval_h) in REGIMENS.items():
        ekfc, z = population[cls]
        cl_caz, cl_avi = clearances(ekfc, z, sc)
        css_caz = dose_g * 1000.0 * CAZ_FRACTION / interval_h / cl_caz
        css_avi = dose_g * 1000.0 * AVI_FRACTION / interval_h / cl_avi
        free_caz, free_avi = css_caz * sc.fu_caz, css_avi * sc.fu_avi

        avi_ok = free_avi >= sc.avi_ct
        toxicity = 100.0 * np.mean(css_caz > sc.tox_threshold)
        avi_pct = 100.0 * avi_ok.mean()

        for mic in MIC_GRID:
            caz_ok = free_caz / mic >= sc.caz_target
            rows.append({
                "scenario": sc.name, "regimen": regimen, "ekfc_class": cls,
                "dose_g": dose_g, "interval_h": interval_h, "mic_mg_l": mic,
                "caz_pta_pct": 100.0 * caz_ok.mean(),
                "avi_attainment_pct": avi_pct,
                "joint_pta_pct": 100.0 * (caz_ok & avi_ok).mean(),
                "toxicity_pct": toxicity,
            })
    return rows


def summarize_exposure(population, sc: Scenario = Scenario()):
    """Steady-state concentration distributions, per component and regimen.

    Section 3.3 of the manuscript describes exposure but tabulates nothing; these
    are the numbers that section needs.
    """
    rows = []
    for regimen, (cls, dose_g, interval_h) in REGIMENS.items():
        ekfc, z = population[cls]
        cl_caz, cl_avi = clearances(ekfc, z, sc)
        css_caz = dose_g * 1000.0 * CAZ_FRACTION / interval_h / cl_caz
        css_avi = dose_g * 1000.0 * AVI_FRACTION / interval_h / cl_avi
        row = {"regimen": regimen, "ekfc_class": cls, "dose_g": dose_g,
               "interval_h": interval_h,
               "daily_dose_g": round(dose_g * 24.0 / interval_h, 4)}
        for tag, css in (("caz", css_caz), ("avi", css_avi)):
            row[f"median_{tag}_css"] = float(np.median(css))
            row[f"{tag}_css_p5"] = float(np.quantile(css, 0.05))
            row[f"{tag}_css_p95"] = float(np.quantile(css, 0.95))
        row["avi_attainment_pct"] = 100.0 * float(np.mean(css_avi * sc.fu_avi >= sc.avi_ct))
        row["toxicity_pct"] = 100.0 * float(np.mean(css_caz > sc.tox_threshold))
        rows.append(row)
    return rows


# --- MIC-weighted cumulative fraction of response ----------------------------

def load_mic_distributions():
    dists = {}
    with open(os.path.join(INPUTS, "mic_distributions.csv")) as fh:
        for r in csv.DictReader(fh):
            d = dists.setdefault(r["distribution_id"], {"population": r["population"],
                                                        "role": r["role"], "weights": {}})
            d["weights"][float(r["mic_mg_l"])] = float(r["frequency"])
    for did, d in dists.items():
        total = sum(d["weights"].values())
        if abs(total - 1.0) > 1e-6:            # acceptance criterion 5
            d["weights"] = {k: v / total for k, v in d["weights"].items()}
    return dists


def compute_cfr(rows, dists):
    """CFR = sum_k w_k * PTA(MIC_k), evaluated per regimen and distribution."""
    by_regimen = {}
    for r in rows:
        by_regimen.setdefault(r["regimen"], {})[r["mic_mg_l"]] = r

    out = []
    for regimen, grid in by_regimen.items():
        any_row = next(iter(grid.values()))
        for did, d in dists.items():
            caz = joint = 0.0
            for mic, w in d["weights"].items():
                if mic not in grid:
                    raise KeyError(f"MIC {mic} of {did} is outside the simulated grid")
                caz += w * grid[mic]["caz_pta_pct"]
                joint += w * grid[mic]["joint_pta_pct"]
            out.append({
                "distribution_id": did, "population": d["population"], "role": d["role"],
                "regimen": regimen, "ekfc_class": any_row["ekfc_class"],
                "dose_g": any_row["dose_g"], "interval_h": any_row["interval_h"],
                "caz_cfr_pct": caz,
                "avi_weighted_pct": any_row["avi_attainment_pct"],  # MIC-independent
                "joint_cfr_pct": joint, "toxicity_pct": any_row["toxicity_pct"],
            })
    return out


# --- Convergence and multi-seed reproducibility ------------------------------

def run_convergence(sizes=(1000, 5000, 10000, 50000, 100000), seed=PRIMARY_SEED):
    """Nested sampling: each smaller run is a prefix of the reference population.

    Redrawing at every sample size would confound convergence with between-run
    noise, so the smaller runs reuse the leading subjects of the 100,000-subject
    reference, as the original pipeline documented.
    """
    full = draw_population(N_PER_CLASS, seed)
    ref = {(r["regimen"], r["mic_mg_l"]): r["joint_pta_pct"] for r in evaluate(full)}
    ref_tox = {r["regimen"]: r["toxicity_pct"] for r in evaluate(full)}

    out = []
    for n_total in sizes:
        per = n_total // len(EKFC_CLASSES)
        subset = {cls: (ekfc[:per], z[:per]) for cls, (ekfc, z) in full.items()}
        rows = evaluate(subset)
        delta = np.array([abs(r["joint_pta_pct"] - ref[(r["regimen"], r["mic_mg_l"])])
                          for r in rows])
        d_tox = np.array([abs(r["toxicity_pct"] - ref_tox[r["regimen"]]) for r in rows])
        out.append({
            "virtual_subjects": n_total, "subjects_per_class": per,
            "mean_abs_delta_joint_pta_pp": float(delta.mean()),
            "p95_abs_delta_pp": float(np.percentile(delta, 95)),
            "max_abs_delta_pp": float(delta.max()),
            "mean_abs_delta_toxicity_pp": float(d_tox.mean()),
        })
    return out


def run_multiseed(seeds=(20260707, 20260708, 20260709, 20260710, 20260711)):
    rows = []
    for s in seeds:
        for r in evaluate(draw_population(N_PER_CLASS, s)):
            rows.append({"seed": s, **r})
    return rows


def summarize_multiseed(rows, mics=(4, 8)):
    keyed = {}
    for r in rows:
        if r["mic_mg_l"] not in mics:
            continue
        keyed.setdefault((r["ekfc_class"], r["regimen"], r["mic_mg_l"]), []).append(r)
    out = []
    for (cls, reg, mic), group in keyed.items():
        for metric, field in (("CAZ PTA", "caz_pta_pct"), ("AVI attainment", "avi_attainment_pct"),
                              ("Joint PTA", "joint_pta_pct"), ("Toxicity", "toxicity_pct")):
            v = np.array([g[field] for g in group])
            out.append({"ekfc_class": cls, "regimen": reg, "mic_mg_l": mic, "metric": metric,
                        "mean": float(v.mean()), "sd": float(v.std(ddof=1)),
                        "min": float(v.min()), "max": float(v.max()),
                        "range_pp": float(v.max() - v.min())})
    return out


# --- One-at-a-time global sensitivity ----------------------------------------

def run_gsa(dists, seed=PRIMARY_SEED, regimens=SELECTED_REGIMENS):
    """Paired OAT on common random numbers; ranked by impact on joint CFR.

    Scope is the selected regimen for each renal-function class crossed with the
    four MIC distributions -- the 20 regimen-distribution pairs that carry the
    reported decisions. Ranking over the full 11-regimen grid inflates every
    maximum, because regimens already discarded on efficacy or toxicity grounds
    swing furthest under perturbation.
    """
    population = draw_population(N_PER_CLASS, seed)
    keep = set(regimens)
    base_rows = [r for r in evaluate(population, GSA_SCENARIOS[0]) if r["regimen"] in keep]
    base_cfr = {(c["regimen"], c["distribution_id"]): c for c in compute_cfr(base_rows, dists)}
    base_pta = {(r["regimen"], r["mic_mg_l"]): r for r in base_rows}

    detail, ranking = [], []
    for sc in GSA_SCENARIOS:
        rows = [r for r in evaluate(population, sc) if r["regimen"] in keep]
        cfr = compute_cfr(rows, dists)
        d_cfr, d_pta4, d_tox = [], [], []
        for c in cfr:
            b = base_cfr[(c["regimen"], c["distribution_id"])]
            d_cfr.append(c["joint_cfr_pct"] - b["joint_cfr_pct"])
            d_tox.append(c["toxicity_pct"] - b["toxicity_pct"])
            detail.append({"scenario": sc.name, "regimen": c["regimen"],
                           "distribution_id": c["distribution_id"],
                           "base_joint_cfr_pct": b["joint_cfr_pct"],
                           "scenario_joint_cfr_pct": c["joint_cfr_pct"],
                           "delta_joint_cfr_pp": c["joint_cfr_pct"] - b["joint_cfr_pct"],
                           "base_toxicity_pct": b["toxicity_pct"],
                           "scenario_toxicity_pct": c["toxicity_pct"]})
        for r in rows:
            if r["mic_mg_l"] == 4:
                d_pta4.append(r["joint_pta_pct"] - base_pta[(r["regimen"], 4)]["joint_pta_pct"])
        a_cfr = np.abs(d_cfr)
        ranking.append({
            "scenario": sc.name,
            "max_abs_delta_joint_cfr_pp": float(a_cfr.max()),
            "mean_abs_delta_joint_cfr_pp": float(a_cfr.mean()),
            "max_abs_delta_joint_pta_mic4_pp": float(np.abs(d_pta4).max()),
            "max_abs_delta_toxicity_pp": float(np.abs(d_tox).max()),
            "direction": ("BASELINE" if sc.name == "BASE" else
                          "LOWER JOINT ATTAINMENT" if np.mean(d_cfr) < 0 else
                          "HIGHER JOINT ATTAINMENT"),
        })
    ranking.sort(key=lambda r: -r["max_abs_delta_joint_cfr_pp"])
    for i, r in enumerate(ranking, 1):
        r["rank"] = i
    return ranking, detail


def toxicity_gate_crossings(detail, gate=15.0):
    """Scenarios that push a regimen across the permissibility ceiling."""
    seen, out = set(), []
    for d in detail:
        k = (d["scenario"], d["regimen"])
        if k in seen or d["scenario"] == "BASE":
            continue
        if d["base_toxicity_pct"] <= gate < d["scenario_toxicity_pct"]:
            seen.add(k)
            out.append({"scenario": d["scenario"], "regimen": d["regimen"],
                        "base_toxicity_pct": d["base_toxicity_pct"],
                        "scenario_toxicity_pct": d["scenario_toxicity_pct"],
                        "crosses_gate": "YES"})
    return out


# --- Probabilistic sensitivity analysis --------------------------------------

PSA_RANGES = {                     # uniform unless noted; scenario bounds, not CIs
    "cl_caz_mult": (0.8, 1.2), "cl_avi_mult": (0.8, 1.2),
    "omega_caz_mult": (0.8, 1.2), "omega_avi_mult": (0.8, 1.2),
    "rho": (0.5, 0.94), "fu_caz": (0.8, 0.9), "fu_avi": (0.87, 0.97),
}
PSA_TOX = (80.0, 128.0, 104.0)     # triangular(low, high, mode)


def load_frozen_draws(n_draws):
    """Reuse the frozen 300-draw design so the PSA is comparable to RC1.

    Regenerating the Latin hypercube would resample the parameter space and move
    the lower-ranked correlations around, which makes a like-for-like comparison
    with the recorded ranking impossible.
    """
    path = os.path.join(INPUTS, "psa_draws_frozen.csv")
    if not os.path.exists(path):
        return None
    fields = list(PSA_RANGES) + ["tox_threshold"]
    with open(path) as fh:
        draws = [{k: float(r[k]) for k in fields} for r in csv.DictReader(fh)]
    return draws[:n_draws] if len(draws) >= n_draws else None


def latin_hypercube(n_draws, seed):
    rng = np.random.default_rng(seed)
    names = list(PSA_RANGES) + ["tox_threshold"]
    cols = {}
    for name in names:
        strata = (rng.permutation(n_draws) + rng.random(n_draws)) / n_draws
        if name == "tox_threshold":
            lo, hi, mode = PSA_TOX
            cols[name] = np.array([
                lo + np.sqrt(u * (hi - lo) * (mode - lo)) if u < (mode - lo) / (hi - lo)
                else hi - np.sqrt((1 - u) * (hi - lo) * (hi - mode)) for u in strata])
        else:
            lo, hi = PSA_RANGES[name]
            cols[name] = lo + strata * (hi - lo)
    return [{k: float(cols[k][i]) for k in names} for i in range(n_draws)]


def run_psa(dists, n_draws=300, inner_per_class=10000, seed=PRIMARY_SEED,
            cfr_target=80.0, tox_gate=15.0, progress=False,
            regimens=SELECTED_REGIMENS):
    keep = set(regimens)
    population = draw_population(inner_per_class, seed)
    draws = load_frozen_draws(n_draws) or latin_hypercube(n_draws, seed)

    acc, ranks = {}, {k: [] for k in list(PSA_RANGES) + ["tox_threshold"]}
    for i, d in enumerate(draws):
        if progress and i % 50 == 0:
            print(f"  PSA draw {i}/{n_draws}", flush=True)
        rows = [r for r in evaluate(population, Scenario(name=f"PSA{i}", **d))
                if r["regimen"] in keep]
        cfr = compute_cfr(rows, dists)
        for c in cfr:
            k = (c["ekfc_class"], c["regimen"], c["distribution_id"])
            acc.setdefault(k, {"joint": [], "tox": [], "params": []})
            acc[k]["joint"].append(c["joint_cfr_pct"])
            acc[k]["tox"].append(c["toxicity_pct"])
            acc[k]["params"].append(d)

    from scipy.stats import rankdata
    probabilities, per_param = [], {k: [] for k in ranks}
    for k, v in acc.items():
        joint, tox = np.array(v["joint"]), np.array(v["tox"])
        p_eff = 100.0 * np.mean(joint >= cfr_target)
        p_safe = 100.0 * np.mean(tox <= tox_gate)
        p_both = 100.0 * np.mean((joint >= cfr_target) & (tox <= tox_gate))
        probabilities.append({
            "ekfc_class": k[0], "regimen": k[1], "distribution_id": k[2],
            "p_joint_cfr_ge_90_pct": 100.0 * np.mean(joint >= 90.0),
            "p_joint_cfr_ge_80_pct": p_eff, "p_toxicity_le_15_pct": p_safe,
            "p_both_pct": p_both,
            "robustness_class": ("Robust" if p_both >= 90 else
                                 "Conditionally robust" if p_both >= 70 else "Fragile"),
        })
        rj = rankdata(joint)
        for name in per_param:
            x = np.array([p[name] for p in v["params"]])
            if x.std() > 0:
                per_param[name].append(abs(np.corrcoef(rankdata(x), rj)[0, 1]))

    ranking = sorted(
        ({"parameter": n, "mean_abs_rank_correlation": float(np.mean(c)),
          "max_abs_rank_correlation": float(np.max(c))} for n, c in per_param.items() if c),
        key=lambda r: -r["mean_abs_rank_correlation"])
    for i, r in enumerate(ranking, 1):
        r["rank"] = i
    return ranking, probabilities


# --- I/O ---------------------------------------------------------------------

def write_csv(rows, path):
    if not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {len(rows):5} rows  {os.path.relpath(path, _ROOT)}")


def compare(label, produced, reference_rows, key_fn, field_pairs):
    """Report agreement between produced rows and a frozen reference table."""
    ref = {key_fn(r): r for r in reference_rows}
    stats = {mine: [] for mine, _ in field_pairs}
    for row in produced:
        k = key_fn(row)
        if k not in ref:
            continue
        for mine, theirs in field_pairs:
            try:
                stats[mine].append(abs(float(row[mine]) - float(ref[k][theirs])))
            except (KeyError, ValueError, TypeError):
                pass
    print(f"\n{label}")
    for mine, vals in stats.items():
        if vals:
            v = np.array(vals)
            print(f"    {mine:28} n={len(v):4}  mean |Δ| {v.mean():7.3f}  max {v.max():7.3f}")


def read_reference(name):
    with open(os.path.join(REFERENCE, name)) as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("analysis", choices=["primary", "cfr", "convergence", "multiseed",
                                         "gsa", "psa", "all"])
    ap.add_argument("--verify", action="store_true", help="compare against frozen RC1 tables")
    ap.add_argument("--psa-draws", type=int, default=300)
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()

    want = ["primary", "cfr", "convergence", "multiseed", "gsa", "psa"] \
        if args.analysis == "all" else [args.analysis]
    dists = load_mic_distributions()

    if "primary" in want:
        rows = evaluate(draw_population(N_PER_CLASS, PRIMARY_SEED))
        write_csv(rows, os.path.join(args.out, "primary_pta_results.csv"))
        if args.verify:
            compare("primary PTA vs RC1", rows, read_reference("primary_pta_results.csv"),
                    lambda r: (r.get("regimen", r.get("Regimen")),
                               float(r.get("mic_mg_l", r.get("MIC (mg/L)", 0)))),
                    [("caz_pta_pct", "CAZ PTA (%)"), ("avi_attainment_pct", "AVI attainment (%)"),
                     ("joint_pta_pct", "Joint PTA (%)"), ("toxicity_pct", "CAZ toxicity (%)")])

    if "cfr" in want:
        rows = compute_cfr(evaluate(draw_population(N_PER_CLASS, PRIMARY_SEED)), dists)
        write_csv(rows, os.path.join(args.out, "cfr_all_distributions.csv"))
        if args.verify:
            compare("CFR vs RC1", rows, read_reference("cfr_summary_primary_three.csv"),
                    lambda r: (r.get("regimen", r.get("Regimen")),
                               r.get("distribution_id", r.get("Distribution ID"))),
                    [("caz_cfr_pct", "CAZ CFR %"), ("joint_cfr_pct", "Joint CFR %")])

    if "convergence" in want:
        rows = run_convergence()
        write_csv(rows, os.path.join(args.out, "convergence_summary.csv"))
        for r in rows:
            print(f"    N={r['virtual_subjects']:>6}  mean |Δ| {r['mean_abs_delta_joint_pta_pp']:.3f} pp"
                  f"  max {r['max_abs_delta_pp']:.3f} pp")

    if "multiseed" in want:
        detail = run_multiseed()
        summary = summarize_multiseed(detail)
        write_csv(detail, os.path.join(args.out, "multiseed_detail.csv"))
        write_csv(summary, os.path.join(args.out, "multiseed_summary.csv"))

    if "gsa" in want:
        ranking, detail = run_gsa(dists)
        write_csv(ranking, os.path.join(args.out, "gsa_sensitivity_ranking.csv"))
        write_csv(detail, os.path.join(args.out, "gsa_detail.csv"))
        write_csv(toxicity_gate_crossings(detail),
                  os.path.join(args.out, "gsa_toxicity_gate_crossings.csv"))
        print("\n    top OAT drivers (max |Δ joint CFR|, pp):")
        for r in ranking[:5]:
            print(f"      {r['rank']}. {r['scenario']:16} {r['max_abs_delta_joint_cfr_pp']:6.2f}")

    if "psa" in want:
        ranking, probs = run_psa(dists, n_draws=args.psa_draws, progress=True)
        write_csv(ranking, os.path.join(args.out, "psa_parameter_ranking.csv"))
        write_csv(probs, os.path.join(args.out, "psa_probabilities.csv"))
        counts = {}
        for p in probs:
            counts[p["robustness_class"]] = counts.get(p["robustness_class"], 0) + 1
        print(f"\n    robustness: {counts}")
        if args.verify:
            compare("PSA parameter ranking vs RC1", ranking,
                    read_reference("psa_parameter_ranking.csv"),
                    lambda r: r.get("parameter"),
                    [("mean_abs_rank_correlation", "mean_abs_rank_correlation")])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
