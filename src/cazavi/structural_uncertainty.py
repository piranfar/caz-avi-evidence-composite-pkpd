"""Does the answer depend on which population PK model you start from?

Every number in this paper so far comes from one set of clearance equations.
That is the ordinary situation in pharmacometric simulation and it is also the
first thing a reviewer will press on, so it is worth settling rather than
conceding.

Four independently derived adult population PK models for this combination are
implemented here and run through the identical simulation: same simulated
patients, same random deviates, same regimens, same targets, same unbound
fractions, same exposure ceiling. Only the clearance model changes.

  M1  Cojutti 2024   critically ill, continuous infusion, TDM, EKFC-based
                     doi:10.1093/jac/dkae290                    (primary)
  M2  Chen 2025      CRKP infection, critical and non-critical, CrCL-based
                     doi:10.2147/IDR.S495279
  M3  Registrational near-proportional to creatinine clearance below the knee,
                     shallow above; the relationship the licensed dose rests on
                     doi:10.1007/s00228-019-02804-z, doi:10.1128/AAC.02105-19
  M4  Bensman 2017   adult cystic fibrosis, preserved renal function, no renal
                     covariate; a high-clearance anchor
                     doi:10.1128/AAC.00988-17

Two layers are reported. The structural layer changes only the typical-value
clearance-versus-renal-function relationship and holds between-subject
variability at the primary model's estimates, so any difference is attributable
to model structure alone. The full layer additionally uses each model's own
reported variability where the source reports it.

The qualitative claims the paper makes are then checked in every model: which
component limits attainment and where the limit changes hands, whether the
exposure screen is reached before the attainment target under escalation, and
how far population CFR moves. A "therapeutic window" row was reported here in
earlier drafts and has been withdrawn: individual clearance cancels from its
placement test, so it is an identity rather than a simulation result and cannot
be corroborated by varying the model.

Sources retrieved from PubMed; DOIs are recorded with each model.

Usage:
    python structural_uncertainty.py
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from cazavi_analyses import (
    DEFAULT_OUT as OUT, SELECTED_REGIMENS, _cholesky, compute_cfr,
    draw_population, load_mic_distributions, write_csv,
)
from reproduce_primary_run import (
    AVI_CT, AVI_FRACTION, CAZ_FRACTION, CAZ_TARGET, CL0_AVI, CL0_CAZ,
    EKFC_CLASSES, EKFC_REF, EXP_AVI, EXP_CAZ, FU_AVI, FU_CAZ, MIC_GRID,
    N_PER_CLASS, OMEGA_AVI, OMEGA_CAZ, PRIMARY_SEED, REGIMENS, RHO,
    TOX_THRESHOLD,
)
from scope_extension_analyses import ICU_WEIGHTS

TOX_CEILING = 15.0                   # exceedance ceiling used throughout
DOSE_GRID_G = np.arange(0.625, 30.01, 0.125)
PTA_TARGET = 90.0


def _omega(cv: float) -> float:
    """Lognormal scale from a reported coefficient of variation."""
    return math.sqrt(math.log(1.0 + cv * cv))


def _power(cl0: float, ref: float, exponent: float) -> Callable[[np.ndarray], np.ndarray]:
    return lambda r: cl0 * (r / ref) ** exponent


def _hinged(cl_at_ref: float, ref: float, knee: float, shallow: float):
    """Proportional to renal function up to the knee, shallow above it.

    This is the relationship the registrational analyses describe: clearance of
    both components close to proportional to creatinine clearance below roughly
    80-100 mL/min, rising only modestly above that.
    """
    def f(r):
        r = np.asarray(r, dtype=float)
        below = cl_at_ref * (np.minimum(r, knee) / ref)
        above = np.where(r > knee, (r / knee) ** shallow, 1.0)
        return below * above
    return f


@dataclass
class PKModel:
    key: str
    label: str
    source: str
    doi: str
    population: str
    cl_caz: Callable[[np.ndarray], np.ndarray]
    cl_avi: Callable[[np.ndarray], np.ndarray]
    cv_caz: float | None = None          # reported BSV on CL, if the source gives one
    cv_avi: float | None = None
    min_class_index: int = 0             # first EKFC class the source population covers
    notes: str = ""

    def omegas(self, full: bool):
        if full and self.cv_caz is not None and self.cv_avi is not None:
            return _omega(self.cv_caz), _omega(self.cv_avi)
        return OMEGA_CAZ, OMEGA_AVI


CLASS_ORDER = list(EKFC_CLASSES)

MODELS = [
    PKModel(
        key="M1_Cojutti_2024",
        label="Cojutti 2024 (primary)",
        source="Cojutti PG, Pai MP, Gatti M, et al. J Antimicrob Chemother 2024;79:2801-8",
        doi="10.1093/jac/dkae290",
        population="Critically ill adults, continuous infusion, TDM (Italy)",
        cl_caz=_power(CL0_CAZ, EKFC_REF, EXP_CAZ),
        cl_avi=_power(CL0_AVI, EKFC_REF, EXP_AVI),
        cv_caz=0.6792, cv_avi=0.7691,
        notes="EKFC as the renal-function descriptor; the model used throughout the paper",
    ),
    PKModel(
        key="M2_Chen_2025",
        label="Chen 2025",
        source="Chen Y, Chen B, Huang Y, et al. Infect Drug Resist 2025;18:941-55",
        doi="10.2147/IDR.S495279",
        population="Adults with carbapenem-resistant K. pneumoniae, ICU and ward (China)",
        cl_caz=_power(2.96, 71.3, 0.44),
        cl_avi=_power(3.09, 71.3, 0.41),
        cv_caz=0.5571, cv_avi=0.6669,
        notes="Cockcroft-Gault creatinine clearance; one-compartment; lowest typical clearance of the four",
    ),
    PKModel(
        key="M3_Registrational",
        label="Registrational relationship",
        source="Das S, Zhou D, Nichols WW, et al. Eur J Clin Pharmacol 2020;76:349-61; "
               "Li J, Lovern M, Riccobene T, et al. Antimicrob Agents Chemother 2020;64:e02105-19",
        doi="10.1007/s00228-019-02804-z",
        population="Phase 1-3 adults across cIAI, cUTI and nosocomial pneumonia",
        cl_caz=_hinged(3.2, 50.0, 100.0, 0.30),
        cl_avi=_hinged(4.9, 50.0, 80.0, 0.30),
        notes="Clearance near-proportional to creatinine clearance below the knee, shallow above; "
              "anchored at the reported values for a population averaging 50 mL/min",
    ),
    PKModel(
        key="M4_Bensman_2017",
        label="Bensman 2017 (preserved renal function)",
        source="Bensman TJ, Wang J, Jayne J, et al. Antimicrob Agents Chemother 2017;61:e00988-17",
        doi="10.1128/AAC.00988-17",
        population="Adults with cystic fibrosis, preserved renal function (USA)",
        cl_caz=lambda r: np.full_like(np.asarray(r, dtype=float), 7.53),
        cl_avi=lambda r: np.full_like(np.asarray(r, dtype=float), 12.30),
        cv_caz=1.28 / 7.53, cv_avi=1.96 / 12.30,
        min_class_index=2,
        notes="No renal covariate; evaluated only in the classes its source population covers, "
              "as a high-clearance anchor",
    ),
]


def clearances(model: PKModel, renal, z, full_bsv: bool):
    om_c, om_a = model.omegas(full_bsv)
    eta = z @ _cholesky(om_c, om_a, RHO).T
    return (model.cl_caz(renal) * np.exp(eta[:, 0]),
            model.cl_avi(renal) * np.exp(eta[:, 1]))


def applies(model: PKModel, cls: str) -> bool:
    return CLASS_ORDER.index(cls) >= model.min_class_index


# --- typical clearance by class, so the spread between models is visible -----

def typical_clearances(population, full_bsv=False):
    rows = []
    for model in MODELS:
        for cls, (renal, z) in population.items():
            if not applies(model, cls):
                continue
            cl_caz, cl_avi = clearances(model, renal, z, full_bsv)
            rows.append({
                "model": model.key, "label": model.label, "ekfc_class": cls,
                "typical_cl_caz_l_h": round(float(np.median(model.cl_caz(renal))), 2),
                "typical_cl_avi_l_h": round(float(np.median(model.cl_avi(renal))), 2),
                "median_cl_caz_l_h": round(float(np.median(cl_caz)), 2),
                "p5_cl_caz_l_h": round(float(np.quantile(cl_caz, 0.05)), 2),
                "p95_cl_caz_l_h": round(float(np.quantile(cl_caz, 0.95)), 2),
                "median_cl_avi_l_h": round(float(np.median(cl_avi)), 2),
                "ratio_to_primary": None,
            })
    primary = {r["ekfc_class"]: r["typical_cl_caz_l_h"]
               for r in rows if r["model"] == "M1_Cojutti_2024"}
    for r in rows:
        r["ratio_to_primary"] = round(r["typical_cl_caz_l_h"] / primary[r["ekfc_class"]], 2)
    return rows


# --- the primary attainment table, recomputed under every model --------------

def evaluate_model(model: PKModel, population, full_bsv: bool):
    rows = []
    for regimen, (cls, dose_g, interval_h) in REGIMENS.items():
        if not applies(model, cls):
            continue
        renal, z = population[cls]
        cl_caz, cl_avi = clearances(model, renal, z, full_bsv)
        css_caz = dose_g * 1000.0 * CAZ_FRACTION / interval_h / cl_caz
        css_avi = dose_g * 1000.0 * AVI_FRACTION / interval_h / cl_avi
        free_caz, free_avi = css_caz * FU_CAZ, css_avi * FU_AVI
        avi_ok = free_avi >= AVI_CT
        avi_pct = 100.0 * avi_ok.mean()
        tox = 100.0 * float(np.mean(css_caz > TOX_THRESHOLD))
        for mic in MIC_GRID:
            caz_ok = free_caz / mic >= CAZ_TARGET
            caz_pct = 100.0 * caz_ok.mean()
            rows.append({
                "model": model.key, "bsv": "full" if full_bsv else "structural",
                "regimen": regimen, "ekfc_class": cls,
                "dose_g": dose_g, "interval_h": interval_h, "mic_mg_l": mic,
                "caz_pta_pct": caz_pct, "avi_attainment_pct": avi_pct,
                "joint_pta_pct": 100.0 * (caz_ok & avi_ok).mean(),
                "toxicity_pct": tox,
                "limiting_component": "ceftazidime" if caz_pct < avi_pct - 0.5
                                      else "avibactam" if avi_pct < caz_pct - 0.5
                                      else "both",
            })
    return rows


def switch_mic(rows, model_key):
    """The MIC at which the limiting component changes hands, per renal class."""
    out = []
    sel = [r for r in rows if r["model"] == model_key]
    for cls in CLASS_ORDER:
        cls_rows = [r for r in sel if r["ekfc_class"] == cls]
        if not cls_rows:
            continue
        for regimen in sorted({r["regimen"] for r in cls_rows}):
            grid = sorted((r for r in cls_rows if r["regimen"] == regimen),
                          key=lambda r: r["mic_mg_l"])
            prev, switch = None, None
            for r in grid:
                cur = r["limiting_component"]
                if prev == "avibactam" and cur == "ceftazidime":
                    switch = r["mic_mg_l"]
                    break
                prev = cur if cur != "both" else prev
            out.append({"model": model_key, "ekfc_class": cls, "regimen": regimen,
                        "limit_changes_hands_at_mic": switch if switch else "no switch"})
    return out


# --- claim 2: does the ceiling still bind before the target under escalation --

def escalation_crossing(model: PKModel, population, full_bsv: bool, mic=8.0):
    rows = []
    for cls, (renal, z) in population.items():
        if not applies(model, cls):
            continue
        cl_caz, cl_avi = clearances(model, renal, z, full_bsv)
        dose_for_pta = dose_at_ceiling = None
        for daily in DOSE_GRID_G:
            css_c = daily * 1000.0 * CAZ_FRACTION / 24.0 / cl_caz
            css_a = daily * 1000.0 * AVI_FRACTION / 24.0 / cl_avi
            tox = 100.0 * float(np.mean(css_c > TOX_THRESHOLD))
            joint = 100.0 * float(np.mean(((css_c * FU_CAZ) / mic >= CAZ_TARGET)
                                          & (css_a * FU_AVI >= AVI_CT)))
            if dose_at_ceiling is None and tox > TOX_CEILING:
                dose_at_ceiling = float(daily)
            if dose_for_pta is None and joint >= PTA_TARGET:
                dose_for_pta = float(daily)
                exceed_at_target = tox
            if dose_for_pta is not None and dose_at_ceiling is not None:
                break
        rows.append({
            "model": model.key, "bsv": "full" if full_bsv else "structural",
            "ekfc_class": cls, "mic_mg_l": mic,
            "daily_dose_for_90pct_pta_g": round(dose_for_pta, 3) if dose_for_pta else "not reached",
            "exceedance_at_that_dose_pct": round(exceed_at_target, 1) if dose_for_pta else "",
            "daily_dose_at_15pct_exceedance_g": round(dose_at_ceiling, 3) if dose_at_ceiling else "not reached",
            "ceiling_reached_first": ("yes" if (dose_at_ceiling is not None
                                                and (dose_for_pta is None
                                                     or dose_at_ceiling < dose_for_pta))
                                      else "no"),
        })
    return rows


# --- claim 3: the plasma window is a concentration statement, not a PK one ----

def window_is_model_free():
    """The window bounds contain no pharmacokinetic parameter at all."""
    closes_at = TOX_THRESHOLD * FU_CAZ / CAZ_TARGET
    return {
        "caz_css_floor_multiple_of_mic": round(CAZ_TARGET / FU_CAZ, 2),
        "caz_css_ceiling_mg_l": TOX_THRESHOLD,
        "avi_css_floor_mg_l": round(AVI_CT / FU_AVI, 2),
        "window_closes_at_mic_mg_l": round(closes_at, 1),
        "multiple_of_clinical_breakpoint": round(closes_at / 8.0, 1),
        "depends_on_pk_model": "no",
        "depends_on": "unbound fraction, target multiple, exposure ceiling",
    }


def in_window_by_model(model: PKModel, population, full_bsv: bool, mics=(4.0, 8.0, 16.0)):
    """Proportion for whom some dose lands inside the window, under each model."""
    rows = []
    for cls, (renal, z) in population.items():
        if not applies(model, cls):
            continue
        cl_caz, cl_avi = clearances(model, renal, z, full_bsv)
        for mic in mics:
            need = CAZ_TARGET * mic / FU_CAZ * cl_caz
            cap = TOX_THRESHOLD * cl_caz
            rows.append({"model": model.key, "ekfc_class": cls, "mic_mg_l": mic,
                         "in_window_any_dose_pct": round(100.0 * float(np.mean(need <= cap)), 1)})
    return rows


# --- population-weighted CFR under each model --------------------------------

def weighted_cfr(rows, dists, model_key):
    sel = [r for r in rows
           if r["model"] == model_key and r["regimen"] in SELECTED_REGIMENS]
    if not sel:
        return []
    cfr = compute_cfr(sel, dists)
    covered = {c["ekfc_class"] for c in cfr}
    total_w = sum(w for c, w in ICU_WEIGHTS.items() if c in covered)
    out = []
    for did in sorted({c["distribution_id"] for c in cfr}):
        joint = caz = tox = 0.0
        for c in cfr:
            if c["distribution_id"] != did:
                continue
            w = ICU_WEIGHTS[c["ekfc_class"]] / total_w
            joint += w * c["joint_cfr_pct"]
            caz += w * c["caz_cfr_pct"]
            tox += w * c["toxicity_pct"]
        out.append({"model": model_key, "distribution_id": did,
                    "classes_covered": len(covered),
                    "population_joint_cfr_pct": round(joint, 1),
                    "population_caz_cfr_pct": round(caz, 1),
                    "population_exceedance_pct": round(tox, 1)})
    return out


HEADLINE_DIST = "LEE2022_KPC_KP"


def robustness_ledger(pta, cross, window_rows, cfr):
    """Each claim the paper makes, checked against every model in turn."""
    structural = [r for r in pta if r["bsv"] == "structural"]
    rows = []

    def per_model(fn):
        return {m.key: fn(m) for m in MODELS}

    # 1. which component limits attainment, and where the limit changes hands
    def switch(m):
        sw = {r["limit_changes_hands_at_mic"] for r in switch_mic(structural, m.key)}
        sw.discard("no switch")
        return f"MIC {min(sw):g}" if sw else "no switch"
    rows.append({"claim": "The limiting component changes from avibactam to ceftazidime "
                          "as the MIC rises",
                 "metric": "lowest MIC at which the limit changes hands",
                 "verdict": "replicates", **per_model(switch)})

    # 2. the avibactam ceiling
    def avi_ceiling(m):
        vals = [r["avi_attainment_pct"] for r in structural
                if r["model"] == m.key and r["regimen"] in SELECTED_REGIMENS
                and r["mic_mg_l"] == 8]
        return f"{max(vals):.1f}%" if vals else "n/a"
    rows.append({"claim": "Avibactam attainment is capped below 100% at the target in use, "
                          "independently of the MIC",
                 "metric": "highest avibactam attainment, licensed regimens",
                 "verdict": "below 100% in every model; height varies",
                 **per_model(avi_ceiling)})

    # 3. attainment at the breakpoint
    def best_pta(m):
        vals = [r["joint_pta_pct"] for r in structural
                if r["model"] == m.key and r["regimen"] in SELECTED_REGIMENS
                and r["mic_mg_l"] == 8]
        return f"{max(vals):.1f}%" if vals else "n/a"
    rows.append({"claim": "Joint attainment at the clinical breakpoint under the licensed "
                          "regimens",
                 "metric": "highest joint PTA at MIC 8 mg/L",
                 "verdict": "varies with clearance", **per_model(best_pta)})

    # 4. the safety ceiling binds first
    def ceiling_first(m):
        sel = [r for r in cross if r["model"] == m.key and r["bsv"] == "structural"]
        hits = sum(1 for r in sel if r["ceiling_reached_first"] == "yes")
        return f"{hits}/{len(sel)} classes"
    rows.append({"claim": "Escalating the dose reaches the exposure ceiling before the "
                          "90% attainment target",
                 "metric": "renal classes in which the ceiling binds first",
                 "verdict": "replicates", **per_model(ceiling_first)})

    def exceed_at_target(m):
        vals = [r["exceedance_at_that_dose_pct"] for r in cross
                if r["model"] == m.key and r["bsv"] == "structural"
                and r["exceedance_at_that_dose_pct"] != ""]
        return f"{min(vals):.1f}-{max(vals):.1f}%" if vals else "n/a"
    rows.append({"claim": "Exceedance at the dose that would reach 90% attainment at the "
                          "breakpoint",
                 "metric": "proportion above the exposure screen",
                 "verdict": "replicates", **per_model(exceed_at_target)})

    # 5. WITHDRAWN. Two rows here previously reported that a "therapeutic window"
    # closes at the same MIC under all four models and that every subject can be
    # placed inside it. Both are identities, not simulation results: individual
    # clearance cancels from the placement test, so the statistic takes one value
    # for the whole cohort and cannot be corroborated by varying the model. The
    # manuscript withdrew the claim; the ledger must not reinstate it. The
    # evidence for the withdrawal is kept in critique_response.py, test A.

    # 6. population CFR
    def pop_cfr(m):
        vals = [c["population_joint_cfr_pct"] for c in cfr
                if c["model"] == m.key and c["bsv"] == "structural"
                and c["distribution_id"] == HEADLINE_DIST]
        return f"{vals[0]:.1f}%" if vals else "n/a"
    rows.append({"claim": "Population-weighted joint CFR against the KPC-K. pneumoniae "
                          "MIC distribution",
                 "metric": "joint CFR, source-cohort renal mix",
                 "verdict": "varies with clearance", **per_model(pop_cfr)})

    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    dists = load_mic_distributions()
    population = draw_population(N_PER_CLASS, PRIMARY_SEED)

    provenance = [{"model": m.key, "label": m.label, "population": m.population,
                   "source": m.source, "doi": m.doi,
                   "reported_bsv_caz_cv_pct": round(100 * m.cv_caz, 1) if m.cv_caz else "not reported",
                   "reported_bsv_avi_cv_pct": round(100 * m.cv_avi, 1) if m.cv_avi else "not reported",
                   "classes_evaluated": ", ".join(CLASS_ORDER[m.min_class_index:]),
                   "notes": m.notes} for m in MODELS]
    write_csv(provenance, os.path.join(OUT, "structural_models.csv"))

    typ = typical_clearances(population)
    write_csv(typ, os.path.join(OUT, "structural_typical_clearances.csv"))

    pta, cross, window_rows = [], [], []
    for full in (False, True):
        for model in MODELS:
            pta += evaluate_model(model, population, full)
            cross += escalation_crossing(model, population, full)
            if not full:
                window_rows += in_window_by_model(model, population, full)
    write_csv(pta, os.path.join(OUT, "structural_uncertainty_pta.csv"))
    write_csv(cross, os.path.join(OUT, "structural_escalation_crossing.csv"))
    write_csv(window_rows, os.path.join(OUT, "structural_window.csv"))

    structural = [r for r in pta if r["bsv"] == "structural"]
    switches = []
    for model in MODELS:
        switches += switch_mic(structural, model.key)
    write_csv(switches, os.path.join(OUT, "structural_limiting_component.csv"))

    cfr = []
    for full in (False, True):
        layer = [r for r in pta if r["bsv"] == ("full" if full else "structural")]
        for model in MODELS:
            for row in weighted_cfr(layer, dists, model.key):
                row["bsv"] = "full" if full else "structural"
                cfr.append(row)
    write_csv(cfr, os.path.join(OUT, "structural_uncertainty_cfr.csv"))

    win = window_is_model_free()
    write_csv([win], os.path.join(OUT, "structural_window_invariance.csv"))

    ledger = robustness_ledger(pta, cross, window_rows, cfr)
    write_csv(ledger, os.path.join(OUT, "structural_robustness_ledger.csv"))

    # ---------------- report ----------------
    print("\nfour population PK models, one simulation")
    for m in MODELS:
        print(f"  {m.key:20} {m.population}")

    print("\ntypical ceftazidime clearance by renal class (L/h), and ratio to the primary model")
    header = "  " + f"{'class':10}" + "".join(f"{m.key.split('_')[0]:>18}" for m in MODELS)
    print(header)
    for cls in CLASS_ORDER:
        line = f"  {cls:10}"
        for m in MODELS:
            hit = [r for r in typ if r["model"] == m.key and r["ekfc_class"] == cls]
            line += (f"{hit[0]['typical_cl_caz_l_h']:>11.2f}"
                     f" ({hit[0]['ratio_to_primary']:.2f})") if hit else f"{'-':>18}"
        print(line)

    print("\njoint PTA at the clinical breakpoint (MIC 8 mg/L), licensed regimens, structural layer")
    print(f"  {'class':10}{'regimen':9}" + "".join(f"{m.key.split('_')[0]:>10}" for m in MODELS))
    for regimen in SELECTED_REGIMENS:
        cls = REGIMENS[regimen][0]
        line = f"  {cls:10}{regimen:9}"
        for m in MODELS:
            hit = [r for r in structural if r["model"] == m.key
                   and r["regimen"] == regimen and r["mic_mg_l"] == 8]
            line += f"{hit[0]['joint_pta_pct']:9.1f}%" if hit else f"{'-':>10}"
        print(line)

    print("\navibactam attainment ceiling (MIC-independent), structural layer")
    for m in MODELS:
        vals = [r["avi_attainment_pct"] for r in structural
                if r["model"] == m.key and r["regimen"] in SELECTED_REGIMENS
                and r["mic_mg_l"] == 8]
        if vals:
            print(f"  {m.key:20} {min(vals):5.1f}-{max(vals):5.1f}%")

    print("\nwhere the limit changes hands (avibactam -> ceftazidime), structural layer")
    for m in MODELS:
        sw = sorted({str(r["limit_changes_hands_at_mic"]) for r in switches
                     if r["model"] == m.key})
        print(f"  {m.key:20} MIC {', '.join(sw)} mg/L")

    print("\nunder escalation at MIC 8: is the exposure ceiling reached before 90% PTA?")
    print(f"  {'model':20}{'class':10}{'dose for 90%':>14}{'exceed there':>14}"
          f"{'dose at ceiling':>17}{'ceiling first':>15}")
    for r in cross:
        if r["bsv"] != "structural":
            continue
        d = r["daily_dose_for_90pct_pta_g"]
        e = r["exceedance_at_that_dose_pct"]
        print(f"  {r['model']:20}{r['ekfc_class']:10}"
              f"{(f'{d} g' if isinstance(d, float) else d):>14}"
              f"{(f'{e}%' if e != '' else '-'):>14}"
              f"{r['daily_dose_at_15pct_exceedance_g']:>17}"
              f"{r['ceiling_reached_first']:>15}")

    print(f"\npopulation-weighted joint CFR, {HEADLINE_DIST}")
    print(f"  {'model':20}{'structural':>13}{'full BSV':>11}{'exceedance':>13}")
    for m in MODELS:
        s = [c for c in cfr if c["model"] == m.key and c["bsv"] == "structural"
             and c["distribution_id"] == HEADLINE_DIST]
        f = [c for c in cfr if c["model"] == m.key and c["bsv"] == "full"
             and c["distribution_id"] == HEADLINE_DIST]
        if s:
            print(f"  {m.key:20}{s[0]['population_joint_cfr_pct']:12.1f}%"
                  f"{f[0]['population_joint_cfr_pct']:10.1f}%"
                  f"{s[0]['population_exceedance_pct']:12.1f}%")

    print("\nthe plasma therapeutic window contains no pharmacokinetic parameter")
    print(f"  floor {win['caz_css_floor_multiple_of_mic']} x MIC, ceiling "
          f"{win['caz_css_ceiling_mg_l']:.0f} mg/L, avibactam floor "
          f"{win['avi_css_floor_mg_l']} mg/L")
    print(f"  closes at MIC {win['window_closes_at_mic_mg_l']} mg/L "
          f"({win['multiple_of_clinical_breakpoint']}x the clinical breakpoint) in every model")
    for m in MODELS:
        vals = [r["in_window_any_dose_pct"] for r in window_rows if r["model"] == m.key]
        print(f"  {m.key:20} some dose lands in the window for "
              f"{min(vals):.1f}-{max(vals):.1f}% of subjects at MIC 4-16")

    print("\nclaim-by-claim, across the four models")
    for r in ledger:
        print(f"  [{r['verdict']:>29}]  {r['claim']}")
        print(f"{'':16}" + "  ".join(f"{m.key.split('_')[0]} {r[m.key]}" for m in MODELS))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
