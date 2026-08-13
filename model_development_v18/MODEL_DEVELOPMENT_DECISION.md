# MODEL_DEVELOPMENT_DECISION.md — Phase 3

> ## ⚠ CORRECTION — superseded numbers
>
> This document quotes the ceftazidime-avibactam clearance correlation as **0.560-0.598**, from
> two-stage analyses. Those estimates are **attenuated by estimation error in the individual
> parameters** (regression dilution). The mixed-effects model in `MODEL1_REPORT.md` estimates the
> comparable quantity at **0.703 (95% CI 0.381 to 0.873)**, which still excludes the assumed 0.94
> but by a smaller margin. **`MODEL1_REPORT.md` is authoritative wherever the two disagree.**
>
> Consequently the "roughly threefold" increase in false reassurance quoted below is
> **2.6-fold** at the fitted value: 3.6% at rho = 0.94 against 9.2% at rho = 0.703.

**Date:** 11 August 2026
**Inputs:** `PROJECT_AUDIT_REPORT.md`, `REPRODUCTION_CHECK.md`, `RESULT_PROVENANCE_MATRIX.csv`,
`PHASE2_DATA_AVAILABILITY_REPORT.md`, `DATA_AVAILABILITY_MATRIX.csv` (8 studies × 35 fields)

> ## SUPERSEDED IN PART — see §6, added after author review
>
> The author has since decided to **abandon the current IJAA submission** rather than revise it, and
> to build **both** models. Removing the submission deadline changed the feasibility of Candidate A,
> and a subsequent identifiability test changed its verdict. **§2 Candidate A is superseded by §6.**
> Everything else in this document stands. See also `NOVELTY_STRATEGY.md`.

---

## 1. The evidence, stated before the recommendation

### 1.1 What exists

| Resource | What it gives | Class |
|---|---|---|
| **Dryad `10.5061/dryad.fxpnvx16s`** — 21 ICU patients on CRRT | 118 paired ceftazidime + avibactam concentrations, 5-7 timepoints each, **CC0 public domain** | **C** |
| **Gatti 2023** (already reference [10]) | 8 patients, 17 TDM occasions, paired clearances, MICs, outcomes, **CC BY-NC-ND** | **C** |
| **Benítez-Cano 2026** (already reference [25]) | Independent non-RRT ICU continuous-infusion cohort; **aggregate** exposures only | **D / E** |
| **Li 2019 CTS** (not currently cited) | Full OMEGA covariance matrices, bootstrap CIs, **executable NONMEM control streams**, n ≈ 2,000 | **D** |
| Cojutti 2024 (the source model) | All parameters confirmed; Ω fully specified | **F / E** |
| Dimelow 2018 | Origin of the 0.85/0.92 unbound fractions; Ω off-diagonals not published | **D / E** |

### 1.2 What does not exist

**No dataset — public, purchasable, or obtainable within a realistic timeframe — contains
individual concentration-time data from critically ill adults receiving continuous-infusion
ceftazidime-avibactam without renal replacement therapy.** That is the primary scenario of this
manuscript, and it is empty. The search covered Dryad, Zenodo, figshare, OSF, Harvard Dataverse,
Mendeley Data, BioModels, DDMoRe, GitHub/GitLab, Vivli, CSDR, YODA, the FDA and EMA document
portals, and the supplementary material of every relevant publication.

The registrational trial data (REPRISE, RECAPTURE, REPROVE, RECLAIM) are reachable **only** through
Vivli, which requires a proposal, a statistical analysis plan, a data use agreement, and
in-platform-only analysis with no raw export. Realistic time to data: **3-9 months.**

### 1.3 The finding that decides the question

Two independent ICU cohorts permit the ceftazidime-avibactam clearance correlation to be estimated
from real paired measurements rather than assumed:

| Source | n | Estimate (log scale) | 95% CI | vs ρ = 0.94 |
|---|---|---|---|---|
| Gatti 2023 (CVVHDF, Italy) | 17 occasions / 8 patients | **r = 0.598** | 0.165 to 0.838 | p ≈ 9 × 10⁻⁵ |
| Dryad / Li (CRRT, China) | 21 patients | **r = 0.560** | 0.169 to 0.799 | p ≈ 3 × 10⁻⁶ |

Two unrelated cohorts, two countries, two renal-replacement modalities, two independent analyses —
**essentially the same answer, and both far below the assumed 0.94.**

The argument is stronger than the numbers alone, because of the direction of the population bias.
**Both cohorts are on renal replacement therapy, where a shared extracorporeal circuit eliminates
both analytes together and should therefore *inflate* the correlation between their clearances.**
That the shared-pathway populations give ≈ 0.57-0.60 makes 0.94 harder, not easier, to justify in a
population with mixed renal and non-renal elimination. Note also that the source's own estimate
carries a relative standard error of **23.8%** — substantial for a parameter bounded at 1.

**Why this matters.** ρ is the parameter behind the manuscript's most clinically actionable
conclusion — that ceftazidime-based prediction of avibactam attainment achieves PPV 95.8% and
NPV 83.6% with only 5.9% misclassification. Re-running the manuscript's own classifier:

| ρ | Specificity | NPV | Wrongly reported as attaining |
|---|---|---|---|
| **0.94** (assumed) | 77.0% | 83.6% | **3.6%** |
| **0.60** (both empirical estimates) | 29.5% | 65.8% | **11.0%** |
| 0.17 (lower CI bound) | 10.9% | 55.5% | 13.9% |

At the empirically supported value the false-reassurance rate **triples** and specificity collapses.
The paper already runs this sensitivity; what it has never had is an empirical anchor for where on
that curve reality sits.

### 1.4 A second, independent check

The primary model was compared against Benítez-Cano's independent non-RRT ICU cohort (6 g/1.5 g per
day by continuous infusion, median eGFR 63 mL/min/1.73 m²). Under continuous infusion at steady
state Css = R/CL, and clearance is lognormal about its typical value, so the median prediction is
exactly R/CL_typ — no simulation assumption enters.

| Analyte | Predicted median Css | Observed median Css | Prediction error |
|---|---|---|---|
| Ceftazidime | 53.8 mg/L | 81.0 mg/L | **−33.5%** |
| Avibactam | 11.6 mg/L | 10.7 mg/L | **+8.7%** |

**The model predicts avibactam well and underpredicts ceftazidime exposure by about a third.** The
direction is favourable to the paper: the finding that ceftazidime becomes limiting at the EUCAST
breakpoint is **conservative**, and the paper's central claim concerns avibactam — the component the
model predicts most accurately.

---

## 2. The four candidate models, assessed

### CANDIDATE A — a new or updated population PK model → **REJECTED**

**Not feasible, and would be indefensible if attempted.**

The only individual concentration-time data available are 21 CRRT patients receiving *intermittent*
8-hourly infusion, plus 8 CVVHDF patients with derived clearances only. Neither matches the primary
scenario on population *or* on administration mode. Fitting a joint two-analyte model with correlated
clearances, renal covariates, and a plasma-ELF link to 29 patients from two incompatible cohorts
would be textbook overparameterisation. The brief's own instruction — *"do not fit an
overparameterised model to a small dataset"* — settles it.

A more specific objection: a model fitted to CRRT patients could not be applied to the manuscript's
non-RRT scenario without an extrapolation larger than anything the model itself estimates.

**Revisit only if** the Vivli application succeeds, or an author releases a non-RRT
continuous-infusion cohort. Both are Phase 8 questions, not Phase 4.

### CANDIDATE B — external validation of the current model → **PARTIALLY ADOPTED, RENAMED**

**Adopt as an external *aggregate-level predictive check*. Do not call it external validation.**

True external validation requires individual observations, and none exist for the primary scenario.
What is available is §1.4: an independent cohort's published aggregate exposures, against which the
model's predicted medians can be compared. That is real external evidence — the cohort had no part
in building the model — but it supports prediction error on a median, not a visual predictive check,
not individual predictions, not calibration of attainment classifications.

**The following must not be attempted or claimed:** prediction-corrected VPCs (no individual data),
individual predictions (no individual data), or ELF observed-versus-predicted (ELF concentrations
are published only as medians). Digitising Benítez-Cano's per-patient figure panels was considered
and **rejected**: logarithmic axes with sparse decade gridlines and heavy overplotting mean
digitisation could not recover exact values, times, or ELF sampling bin. Under the project's own
rule that digitisation must be scientifically defensible, it is not defensible here.

### CANDIDATE C — TU-JAM, target-uncertainty-integrated joint attainment → **ADOPTED, WITH A HARD CONSTRAINT**

**Adopt the decision-analytic machinery. Reject any single fitted distribution over the target.**

The evidence base assembled for this project (Supplementary Table S14, extended by the Phase 2
survey) is not a set of noisy measurements of one underlying quantity. It is a set of **different
constructs**:

| Value | What it actually is |
|---|---|
| 0.15-0.5 mg/L | hollow-fibre **regrowth threshold** (Coleman 2014) — a measured exposure-response quantity |
| 1.0 mg/L | the **regulatory target**, animal-derived (Berkhout 2015/2016, murine thigh and lung, doi:10.1128/AAC.01269-15), expressed as **%fT > C_T**, not a steady-state concentration — and never tested under continuous infusion; see R3 in `NOVELTY_ROUTES.md` |
| 2.5 mg/L | **aztreonam**-avibactam specific — must be **excluded**, it does not belong to a ceftazidime model |
| 4.0 mg/L | the **EUCAST fixed susceptibility-testing concentration**, repurposed as a clinical target |

Three facts from the Phase 2 survey sharpen this further. The pharmacodynamic index is
**time above a threshold**, not AUC or Cmax, established by dose fractionation. The threshold
**tracks the partner β-lactam, not the organism or the MIC** — there is no correlation between C_T
and ceftazidime-avibactam MIC. And the FDA's own Phase 2 clinical exposure-response analysis was
**negative**: attainment values clustered near 100%, formal modelling was infeasible, and the
reduction in microbiological failures was not statistically significant.

**Pooling these into a single probability distribution would be inventing one** — precisely what the
brief forbids. The defensible construction is a **prespecified scenario set of explicitly labelled
distributions**, each carrying its own interpretation, with the results reported across all of them
and sensitivity to the choice reported as a primary output rather than an appendix.

What is genuinely novel here is not the distribution. It is the **decision-analytic layer**: the
probability that each component is limiting once the target is uncertain; the probability of
selecting a different regimen than one would under a single fixed threshold; expected regret; and
robust regimen ranking. Those quantities are new to this literature, they are computable from
machinery this project already has, and they answer the question the manuscript currently only poses.

**Constraints, non-negotiable:**
- Every distribution is labelled by kind — evidence-weighted, triangular, discrete, uniform,
  lognormal, expert-specified — and never presented as a posterior.
- The resulting distribution is **never** described as an established clinical exposure-response
  relationship.
- The 2.5 mg/L aztreonam value is excluded, with the reason stated.
- Results are reported for **all** prespecified distributions, not a preferred one.

### CANDIDATE D — Bayesian model updating → **REJECTED for this revision**

Tempting, and now technically possible: Li 2019 publishes full OMEGA covariance matrices with
bootstrap confidence intervals and executable NONMEM control streams, which would make a genuine
prior specifiable for the first time.

But there is nothing defensible to update **on**. The only individual data are from populations the
primary scenario excludes, and partial pooling across cohorts that differ in renal replacement
status *and* administration mode would produce a posterior belonging to no real population. The
brief's warning — *"do not use Bayesian terminology merely for presentation"* — applies exactly.

**Retain Li 2019 for Candidate C's parameter-uncertainty layer**, where the published bootstrap
intervals give a defensible, non-Bayesian representation of parameter uncertainty. And **cite it**:
it is the richest public parameter package for either drug and the manuscript currently omits it.

---

## 3. Recommendation

> **Adopt a three-part revision. Do not fit a new population pharmacokinetic model.**
>
> **(1) Empirical examination of the clearance correlation.** Report ρ estimated from two
> independent, openly licensed patient-level ICU datasets (r = 0.560 and r = 0.598), state the
> renal-replacement limitation plainly, note that the shared-circuit population should if anything
> inflate the correlation, and re-report the assay operating characteristics across the empirically
> supported range. This converts the paper's assay analysis from a hypothetical sensitivity into an
> evidence-anchored conclusion.
>
> **(2) External aggregate-level predictive check.** Report predicted-versus-observed median
> steady-state exposure against the independent Benítez-Cano cohort: ceftazidime −33.5%, avibactam
> +8.7%. Named precisely, limitations stated, and explicitly **not** called external validation.
>
> **(3) TU-JAM as a prespecified scenario set.** Implement the decision-analytic layer — limiting-
> component probability, probability of regimen misselection under a fixed threshold, expected
> regret, robust ranking — across multiple explicitly labelled target distributions, with no claim
> of a fitted posterior.

**Why this and not more.** It uses every piece of genuine patient data that legally and
scientifically can be used; it adds a real external confrontation the manuscript does not have; it
supplies the novel methodological contribution the brief asks for; and it does none of it by
overclaiming. Each part strengthens a conclusion the paper already draws rather than replacing it.

**Why this and not less.** Parts 1 and 2 together answer the objection that the study is a closed
loop — a simulation calibrated against the same source that supplied its parameters. After this
revision the model will have been confronted with observations from three cohorts that had no part
in building it. That is a substantive change in the study's evidentiary standing, and it is honest.

### Expected contribution

- **First empirical estimate** of the ceftazidime-avibactam clearance correlation from patient-level
  data, from two independent cohorts.
- **First external predictive check** of the Cojutti-derived relationships against an independent
  continuous-infusion ICU cohort.
- A **decision-analytic framework** for target attainment under target uncertainty, with regret and
  misselection probability — quantities not previously reported for this drug combination.
- A **repaired, genuinely reproducible release package** (Phase 1 defects 1-10).

### Assumptions and identifiability risks, stated up front

| Risk | Mitigation |
|---|---|
| Both correlation datasets are RRT populations | State it every time; frame as a bound, never as a replacement value for ρ; retain 0.94 as the primary analysis and report 0.56-0.60 as a sensitivity |
| Trapezoidal AUC from sparse sampling biases clearance high | Bias is common to both analytes and largely cancels in the ratio and correlation; absolute clearances reported as approximate |
| eGFR equation differs between cohorts (CKD-EPI vs EKFC) | State it; report the check as approximate |
| n = 15 medians in the external check | Report as a median-level comparison; do not construct confidence statements the data cannot support |
| TU-JAM distributions could be over-interpreted | Multiple prespecified distributions, explicit labelling, no preferred distribution |

### Validation strategy

Internal: unit tests for dose conversion, infusion rate, free-versus-total, clearance transformation,
correlated random effects, renal classes, MIC weighting, PTA, CFR, limiting-component classification,
exposure constraint, ELF transformation. Deterministic seeds; frozen machine-readable outputs with a
recomputed checksum manifest; Monte Carlo convergence demonstrated. **All internal work is labelled
internal.** Only §1.4 is external, and only at aggregate level.

### Implications for the manuscript

The title, study design and ethics statement **do not change**. The study remains a pharmacometric
simulation; the two patient datasets are already-published, de-identified, openly licensed data used
for secondary analysis, requiring no new ethics approval — but the ethics and data-availability
statements must be rewritten to say exactly that, because the current text states that no individual
patient data were used, and after this revision that will no longer be true.

The abstract gains one sentence on the empirical correlation finding. The discussion gains the
external check. **No existing numerical result changes**, because nothing in the primary analysis is
altered — subject to the Phase 1 corrections (67.1% → 67.2%, and the restated PSA magnitudes once
a frozen design is regenerated).

---

## 4. Decision log

| # | Decision | Rationale |
|---|---|---|
| D1 | Do not fit a new population PK model | No individual data exist for the primary scenario; 29 patients from two mismatched populations cannot support it |
| D2 | Do not digitise Benítez-Cano figure panels | Logarithmic axes, sparse gridlines, heavy overplotting; not scientifically defensible under the project's own rule |
| D3 | Use the Dryad CRRT dataset for one narrow purpose only | CC0 public domain, consented for public release; but CRRT and intermittent infusion make it class C |
| D4 | Use Gatti 2023 as an independent replicate of the same estimate | CC BY-NC-ND; transcription verified against the authors' own reported medians and IQRs |
| D5 | Call §1.4 an external aggregate-level predictive check | It is external, but aggregate; naming it validation would breach the brief's rule 2 |
| D6 | TU-JAM as a scenario set, not a fitted posterior | The evidence comprises different constructs, not noisy measures of one quantity |
| D7 | Exclude the 2.5 mg/L threshold | Aztreonam-avibactam specific; does not belong to a ceftazidime model |
| D8 | Defer Bayesian updating | No defensible data to update on; would be terminology without substance |
| D9 | Cite Li 2019 CTS and Gatti 2025 AAC | Li 2019 is the richest public parameter package; Gatti 2025 (n = 218) is the only positive clinical exposure-response evidence for the joint target, and the manuscript omits both |
| D10 | Retain ρ = 0.94 as the primary analysis | Changing an established result requires rerunning and re-verifying everything; the empirical estimate enters as a sensitivity analysis, not a replacement |

---

## 5. What is needed before Phase 4 begins

**Two questions only the author can answer:**

1. **The manuscript is under active submission** (IJAA-S-26-02124, revision 3). Should this revision
   be prepared as a new version of that submission, or held until the editor responds?
2. **Should the data request to the Benítez-Cano group be sent?** It is drafted and ready. If sent
   and granted, Candidate A and true Candidate B become feasible and the scope changes substantially.
   If not sent, the plan above stands unchanged.

**Everything else can proceed without further input.**

---

# 6. Addendum — Candidate A revisited after an identifiability test

**Added 11 August 2026, after author review. This supersedes the Candidate A verdict in §2.**

## 6.1 What changed

Two things.

**The submission deadline is gone.** The author has abandoned IJAA-S-26-02124. A three-month wait
for author-provided data, previously impossible, is now acceptable — which changes what is worth
attempting.

**An identifiability test was run, and Candidate A passed it.** §2 rejected a new population PK model
on the grounds that 29 patients from two mismatched cohorts could not support one. That reasoning
conflated two different questions: *can we model the primary scenario?* (no — no data exist) and
*can we model the CRRT cohort?* (not previously tested).

Individual one-compartment fits to the 21-patient Dryad cohort:

| | Result |
|---|---|
| Convergence | **21 of 21 subjects, both analytes** |
| Median residual error | **8.7% ceftazidime, 8.4% avibactam** |
| CL and V separately identifiable | yes, from 5-7 samples per analyte |
| Cross-drug clearance correlation (marginal) | **r = 0.563** (95% CI 0.173-0.800) |
| Same, conditional on CRRT covariates | **r = 0.476** (95% CI 0.012-0.772) |

238 observations across 21 subjects with clean individual fits will support a **reduced** joint
mixed-effects model — four fixed effects, four variance terms, and the cross-drug covariance. Not a
full 4 × 4 Ω, and not a covariate model, but enough for the one parameter that matters.

## 6.2 Why the conditional estimate is the one to quote

Cojutti's 0.94 is a **conditional** correlation — between random effects, after the renal-function
covariate. A marginal correlation computed from raw clearances is a different quantity, and when a
shared covariate drives both analytes the marginal is normally the **larger** of the two, because it
is a variance-weighted blend of the near-unity covariate correlation and the residual correlation.

So the like-for-like comparison is **0.476 against 0.94**, and adjusting for covariates moved the
estimate *down*, exactly as the theory predicts. Three independent routes give closely agreeing
marginal values: Dose/AUC 0.560, model-based individual fits 0.563, Gatti cohort 0.598.

## 6.3 Revised verdict

**CANDIDATE A — ADOPTED, with its scope fixed in the title of the model.**

Build a joint two-analyte population PK model of the **CRRT cohort**, estimating the cross-drug
clearance covariance rather than assuming it. Implemented in `code/joint_popk_nlme.py`.

**Its claim to novelty is narrow and real:** the cross-drug clearance correlation has been quantified
exactly once in the literature. No regulatory document contains it — the registrational analyses
resampled η-pairs empirically without ever estimating the covariance. This is the second estimate in
existence and the first from openly available data.

**Its limitations are structural and must appear wherever it appears:**

- The population is on renal replacement therapy; the primary scenario excludes it.
- Administration is intermittent 8-hourly infusion, not continuous infusion.
- 18 of 21 patients have acute pancreatitis; the cohort is not a general ICU population.
- Sparse sampling misses the peak in most patients, so absolute clearances are approximate.
- **Its clearance and volume estimates do not transfer to the primary scenario.** Only the
  correlation structure is carried forward, and only as a sensitivity bound.

**What it does not become:** a replacement for the Cojutti model, a validation of it, or a source of
ρ for the non-RRT population. The primary analysis retains ρ = 0.94; the estimate enters as a
sensitivity analysis.

## 6.4 Revised architecture

```
Model 1  joint population PK, CRRT cohort          -> estimated rho, with uncertainty
Gatti 2023 independent replicate                   -> second estimate of the same quantity
Benitez-Cano external aggregate check              -> prediction error on the primary model
                     |
                     v
Model 2  hierarchical uncertainty-propagating joint attainment model
         parameters as distributions (published RSEs)
         between-study heterogeneity across published models
         rho informed by Model 1, not fixed
         avibactam target as a prespecified scenario set
                     |
                     v
         value of information, limiting-component probability,
         regret, robust ranking, falsifiable prediction
```

## 6.5 What remains impossible, and stays impossible

- A population PK model of the **primary scenario** — non-RRT adults on continuous infusion. No
  individual data exist anywhere for that population.
- True external validation of the primary model.
- Any claim that the study is no longer a simulation built predominantly on published parameters.
- Any claim about bioRxiv eligibility.

**The one action that would remove the first two constraints** is obtaining the Bologna dataset:
112 non-RRT critically ill adults, 185 paired steady-state concentrations under continuous infusion.
That is the primary-scenario population. See `NOVELTY_STRATEGY.md` §6.
