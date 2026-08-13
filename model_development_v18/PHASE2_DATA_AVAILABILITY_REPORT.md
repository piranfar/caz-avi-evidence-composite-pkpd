# PHASE 2 — Real patient-data availability review

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

**Companion file:** `DATA_AVAILABILITY_MATRIX.csv` (6 studies × 35 evidence fields)
**Date:** 11 August 2026
**Method:** every claim below was verified by reading the source — publisher HTML, PMC/Europe PMC
full text, or the downloaded Version of Record PDF — and by retrieving and unpacking every
supplementary file that could be obtained. Fields that could not be verified are marked
NOT ASSESSED rather than inferred.

---

## 1. What was found

| Study | Individual patient-level numeric data | Class |
|---|---|---|
| **Gatti 2023** (CVVHDF case series) | **PUBLISHED IN FULL, open licence** | **C** |
| Cojutti 2024 (the source of every PK parameter) | Not published; paywalled; no data availability statement | F / E |
| Benítez-Cano 2026 (ICU ELF trial) | Not published; figure points only; on-request statement | E / D |
| Dimelow 2018 (healthy volunteers) | Not published; controlled access via Pfizer | D / E |
| O'Jeanson 2025 | None exists — it is a simulation study | F |
| Curtiaud 2024 (ECMO) | Out of scope for the non-RRT, non-ECMO primary scenario | not assessed |

**Classification key:** A suitable for model development · B suitable for external validation ·
C suitable only for sensitivity analysis · D aggregate evidence only · E potentially available by
author request · F not suitable.

**No study qualifies as A or B.** No public repository holds an individual-level ceftazidime or
avibactam pharmacokinetic dataset relevant to this population. That conclusion is unwelcome but it
is the honest one, and it settles the model-development question in Phase 3.

---

## 2. Benítez-Cano 2026 — the study the brief singled out

The citation is correct in every particular: *Critical Care* 2026;**30**:305, published 13 May 2026,
PMID **42129898**, PMCID **PMC13262505**, open access under CC BY-NC-ND.

**The figures in the brief were accurate.** The study randomised 30 patients 1:1; **15 received
ceftazidime-avibactam**. The Results state verbatim that 298 plasma samples were obtained,
comprising **74 each for ceftazidime and avibactam**, and 58 ELF samples comprising **14 each for
ceftazidime and avibactam**. One CAZ/AVI patient died before day 3, so no bronchoalveolar lavage was
obtained from them. Dosing was a 2 g/0.5 g loading dose over 120 minutes followed by 2 g/0.5 g every
8 hours infused over 8 hours — true continuous infusion, 6 g/1.5 g per day. Plasma sampling was at
the end of the loading dose and at 1, 4 and 8 hours; ELF was sampled once per patient on days 3-4 at
one of 1, 4 or 8 hours after a bag change, by standardised bronchoscopy with urea dilution.

### 2.1 The decisive question: are numeric individual concentrations available?

**No.** All three supplementary files were downloaded from Springer's static-content host, unpacked
and read in full:

| Additional file | Size | What it actually contains |
|---|---|---|
| Supplementary Material 1 | 27 kB `.docx` | Bioanalytical methods only — HPLC-UV and UPLC-MS/MS conditions, MRM transitions, gradient programmes, calibration ranges. **No patient data.** |
| Supplementary Material 2 | 16 kB `.docx` | One 4 × 4 table: plasma concentrations on day 1 as **median (IQR) only**, at four nominal times. **No individual rows.** |
| Supplementary Material 3 | 10 MB `.docx` | Population PK narrative, four embedded TIFF figures, Tables S3.1 and S3.2. **No numeric individual data.** |

**There is a trap in Supplementary Table S3.1.** It is titled *"Individual observed and predicted
plasma and intrapulmonary AUC and concentrations at steady state"* — but every cell is a median with
an interquartile range across patients, footnoted as such. "Individual" refers to individual
*model-predicted* empirical-Bayes values that were then summarised. The table has eight rows, not
thirty. It would be easy to cite this as individual data. It is not.

Individual data exist **only as figure points**. Supplementary Figures S3.3 and S3.4 are the only
patient-identifiable graphics: a 5 × 3 grid of panels headed by subject ID, each with a plasma and an
ELF sub-panel. They are digitisable in principle, but the y-axes are logarithmic with sparse decade
gridlines and the day-1 points are tightly overplotted. Digitisation could not recover exact values,
exact times, or which of the 1/4/8-hour ELF bins a patient belonged to. **Under the project's own
rule that digitisation must be scientifically defensible, it is not defensible here**, and no
digitisation was performed.

No repository deposit exists. Keyword scans of the full text and all supplements for *repository*,
*github*, *zenodo*, *figshare*, *code* and *NONMEM* returned nothing. Europe PMC's `hasData: Y` flag
refers to the three `.docx` files above, not to a dataset.

**Data availability statement, verbatim:** *"The datasets used and analyzed during the current study
are available from the corresponding author on reasonable request."*
**Corresponding author:** Dr Luisa Sorlí, `lsorli@hmar.cat`.

A minimal, professionally framed data request has therefore been drafted at
`correspondence/DRAFT_data_request_Benitez-Cano_Sorli.md`. **It has not been sent and will not be
sent without explicit instruction.**

### 2.2 What the published aggregate data already permit

Even without individual data, this study supports something the project does not currently have: an
**external, aggregate-level predictive check against a cohort that had no part in building the
model**. The article publishes median steady-state exposures for both components under a known
continuous-infusion regimen in a cohort of known median renal function. Section 4 below reports that
check.

### 2.3 Three discrepancies worth recording

1. **Three different ELF/plasma ratio pairs are published in the same paper.** Median individual
   AUC(0-8,ss) ratios of **0.41 / 0.44** (the pair this manuscript uses in Table 3), a model-based
   ELF-plasma ratio of **0.51 / 0.60** (Table 2), and observed ratios of **0.46 / 0.40** in
   Supplementary Table S3.1. Using the median individual ratios is defensible and arguably the most
   conservative choice, but the manuscript does not state which pair it used or why, and reports no
   sensitivity to that choice.
2. **A unit inconsistency in the source.** The main text gives calibration ranges in mg/L; the
   bioanalytical supplement gives the identical numbers in ng/mL. The latter is almost certainly a
   typographical error for µg/mL.
3. **An internal inconsistency in the source.** The proportion of simulated subjects exceeding a
   ceftazidime concentration of 78 mg/L on the low-dose regimen is given as 4% in the main text and
   40.1% in the supplement's own narrative — an order-of-magnitude conflict. This does not affect
   the present project, which uses only the penetration ratios, but it argues for caution in citing
   that paper's simulation results.

---

## 3. Cojutti 2024 — the study everything depends on

**Every parameter this project attributes to Cojutti was independently confirmed against the source
publication.** Nothing was refuted:

| Project asserts | Source reports | |
|---|---|---|
| CL ceftazidime 5.0 L/h, exponent 0.70 | 5.0 L/h (RSE 6.36%), 0.70 (RSE 14.0%) | ✓ |
| CL avibactam 5.9 L/h, exponent 0.89 | 5.9 L/h (RSE 7.4%), 0.89 (RSE 12.3%) | ✓ |
| Interindividual variability 67.92% / 76.91% CV | 67.92% (RSE 33.3%) / 76.91% (RSE 30.1%) | ✓ |
| Clearance correlation ρ = 0.94 | 0.94 (RSE 23.8%) | ✓ |
| Unbound fractions 0.85 / 0.92 | *"multiplying by 0.85 and 0.92, respectively"* | ✓ |
| Targets fCss/MIC ≥ 4 and fCss avibactam ≥ 4 mg/L | stated as fCss/MIC ≥ 4 and fCss/CT ≥ 1 with CT = 4 mg/L | ✓ |
| Exposure ceiling: total Css ceftazidime > 104 mg/L | 104 mg/L, neurotoxicity-associated | ✓ |
| 72-row calibration benchmark | **Table 4**: 9 permissible regimens × 8 MIC values = 72 cells | ✓ |

Two facts are worth carrying into the Methods:

- **The source Ω matrix is fully specified.** Interindividual variability is estimated on only two
  parameters, and the single off-diagonal correlation is published. The reimplementation therefore
  inherits no unreported covariance terms — an unusually clean situation that the manuscript should
  state, because it is a genuine strength of the reconstruction.
- **The avibactam target is expressed as a ratio in the source** (fCss/CT ≥ 1 with CT = 4 mg/L),
  not as an absolute concentration. Arithmetically identical, but mirroring the source's own wording
  would be more faithful.

**No individual patient data.** The article is closed access with all rights reserved; Unpaywall
finds no open copy anywhere. Table 1 is aggregate. The supplementary material is paywalled and
consists of model-development tables and diagnostic plots, not patient data. **The article carries no
data availability statement at all** — verified across two independent extractions.

**Consequence:** the model's own source cannot be used to validate the model. Calibration against
Cojutti's Table 4 is an internal reproduction check, exactly as the manuscript already says. That
characterisation must be preserved.

---

## 4. Gatti 2023 — the one genuine individual-patient dataset

This study is **already reference [10] of the manuscript**, and it publishes complete patient-level
data under CC BY-NC-ND:

- **Table 1** — 8 patients, one row each: age, sex, pathogen, **MIC**, infection type, dose and
  adjustment, **average free steady-state concentrations of both components**, target ratios, joint
  target status, treatment and CVVHDF duration, microbiological eradication, 30-day mortality.
- **Table 2** — **17 therapeutic-drug-monitoring occasions**, one row each: weight, full CVVHDF
  circuit settings, dose intensity, residual diuresis, total effluent flow, and **paired ceftazidime
  and avibactam total clearance in L/h**.

The Version of Record PDF is freely downloadable without login from the University of Bologna
repository. The values have been transcribed to
`data_external/Gatti2023_individual_patient_data.csv` with a full provenance record.
**The transcription reproduces the authors' independently reported medians and interquartile ranges
exactly** for both analytes — 2.39 L/h (2.05-2.94) for ceftazidime and 2.56 L/h (2.22-2.96) for
avibactam — which verifies it.

**What it is not.** These are CVVHDF patients; the primary scenario excludes renal replacement
therapy. Total clearance is dominated by the extracorporeal circuit, which compresses between-patient
variability and imposes a shared circuit-driven component in both analytes. It is **class C —
suitable for sensitivity and assumption-testing analysis only**. It cannot validate the primary
model, and nothing fitted to 8 patients could support a population PK model.

**What it is.** The only place where the assumption of near-perfect ceftazidime-avibactam clearance
correlation can be examined against paired measurements in real patients rather than assumed.

---

## 5. Two analyses these findings already make possible

Both were run during this audit. Both are honest about what they are.

### 5.1 External aggregate-level predictive check against Benítez-Cano

Under continuous infusion at steady state, Css = R/CL, and because clearance is lognormal about its
typical value, the **median** predicted concentration is exactly R/CL_typ — no simulation assumption
is needed. Comparing the primary model against the independent Barcelona cohort (6 g/1.5 g per day,
median eGFR 63 mL/min/1.73 m²):

| Analyte | Model CL_typ | Predicted median Css | Observed median Css | Prediction error |
|---|---|---|---|---|
| Ceftazidime | 4.64 L/h | 53.8 mg/L | **81.0 mg/L** | **−33.5%** |
| Avibactam | 5.37 L/h | 11.6 mg/L | **10.7 mg/L** | **+8.7%** |

Equivalently, the model's typical ceftazidime clearance is **1.62-fold higher** than the value
independently estimated in that cohort (4.64 vs 2.86 L/h), while for avibactam the two agree to
within 12% (5.37 vs 6.08 L/h).

**Interpretation.** The model predicts avibactam exposure well in an independent cohort and
underpredicts ceftazidime exposure by about a third. The direction matters: the manuscript's finding
that **ceftazidime becomes the limiting component at the EUCAST breakpoint is therefore
conservative** — if real ceftazidime exposure is higher than modelled, ceftazidime attainment at
MIC 8 mg/L would be better than reported. Meanwhile the manuscript's central claim concerns
avibactam, and avibactam is the component the model predicts most accurately. **The external check
supports the paper's principal conclusion and qualifies its secondary one.**

**Limitations that must travel with this result:** renal function was measured by CKD-EPI in the
comparison cohort and enters the model as EKFC; the comparison is of medians in a 15-patient cohort;
the two cohorts differ in country, case mix and era; and the observed interquartile range for
avibactam (4.7 mg/L) is much narrower than predicted (15.0 mg/L), which 15 patients cannot resolve.
**This is an external aggregate-level predictive check. It is not external clinical validation, and
it must never be described as such.**

### 5.2 Empirical examination of the clearance correlation

From the 17 paired clearances in Gatti Table 2, the correlation of log clearances is
**r = 0.598 (95% CI 0.165 to 0.838)**, significantly below the assumed 0.94 (p ≈ 9 × 10⁻⁵). At the
patient-mean level, r = −0.34 with a confidence interval spanning nearly the entire range (n = 8).

Re-running the manuscript's own second-assay classifier across an extended correlation grid shows
what this would mean if it held in the non-RRT population:

| ρ | Specificity | NPV | Wrongly reported as attaining |
|---|---|---|---|
| **0.94** (assumed) | 77.0% | 83.6% | **3.6%** |
| 0.84 (upper CI) | 58.9% | 75.5% | 6.4% |
| **0.60** (empirical estimate) | 29.5% | 65.8% | **11.0%** |
| 0.17 (lower CI) | 10.9% | 55.5% | 13.9% |

At the empirically estimated correlation, the rate at which a patient would be wrongly reported as
attaining the avibactam target on the basis of ceftazidime alone **triples**, and specificity
collapses from 77% to 30%.

**This does not overturn the manuscript — it sharpens it.** The paper already reports this
sensitivity; what it lacks is any empirical anchor for where on that curve reality sits. The honest
statement is that ρ = 0.94 comes from a single cohort, has never been checked against paired patient
measurements outside it, and that the only paired data available — in a different population, where
circuit clearance is shared between analytes and would if anything *inflate* correlation — give a
substantially lower value. That argues for measuring avibactam rather than inferring it, which is a
clinically actionable conclusion the paper is currently one step short of making.

---

## 6. Answer to the Phase 2 question

> **Can this study be strengthened by incorporating genuine individual patient data?**
>
> **Partly, and not in the way the brief anticipated.** No dataset exists — public or obtainable —
> that could support fitting a new population pharmacokinetic model, and none could serve as an
> external validation set for the primary non-RRT scenario. Every candidate is class C, D, E or F.
>
> But one genuine individual-patient dataset is publicly available under an open licence, in a paper
> this manuscript already cites, and it bears directly on the assumption underpinning the paper's
> most clinically actionable conclusion. And an independent ICU cohort's published aggregate
> exposures permit a real external predictive check that the manuscript does not currently have.
>
> **Neither of these is external clinical validation, and neither may be described as such.** Both
> are worth adding, and together they change the manuscript from an entirely closed-loop simulation
> into one that has been confronted with observations from cohorts that had no part in building it.
