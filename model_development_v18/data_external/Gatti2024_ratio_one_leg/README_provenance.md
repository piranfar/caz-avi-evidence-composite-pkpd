# Provenance record — Gatti 2024: the empirical test of this project's central question

## Source

Gatti M, Viale P, Pea F. *Therapeutic drug monitoring of ceftazidime/avibactam: why one leg is not
enough to run.*
**J Antimicrob Chemother. 2024;79(1):195-199.** doi:10.1093/jac/dkad367 · PMID 38019676

Bologna. Same group as `Cojutti2024_ANCHOR_PopPK/`, `Gatti2023_*` and `Gatti2025_outcome_R5/`.

## Why this is the most directly relevant paper archived so far

**This paper is the empirical, clinical version of exactly the question Model 2's monitoring layer
asks.** Its stated purpose, verbatim from the Background:

> "Some authors hypothesized that the PK/PD target attainment of ceftazidime/avibactam could be
> assessed by means of the TDM of solely ceftazidime, since avibactam concentrations might be
> extrapolated based on the fixed 4:1 ceftazidime-to-avibactam ratio present in the vial. The
> reliability of this hypothesis could be called into question if a wide interindividual variability in
> the ceftazidime-to-avibactam ratio would exist among patients."

That is Model 2's "measure avibactam, or infer it?" decision, posed as a clinical measurement question
and answered with real TDM data.

**Answer, from 188 paired TDM assessments in 107 patients:** the ceftazidime-to-avibactam ratio ranges
from **1.29:1 to 13.46:1** against a vial ratio of 4:1, with 41.0% above 5:1 and 19.1% above 6:1.
Conclusion: **"both ceftazidime and avibactam concentrations should be measured."**

## A published controversy the project did not know about

This paper is one side of a direct exchange in *JAC*, and **the project has already archived the other
side without knowing it was contested**:

1. **Fresan et al. 2023** (`data_external/Fresan2023_CI_TDM/`) measured ceftazidime only, stating the
   assumption explicitly: *"only ceftazidime concentrations were measured, on the assumption that the
   concentration of avibactam would always be sufficient to exert its action."* Cited here as ref 3.
2. **Fresán et al., authors' response, 2023** — *J Antimicrob Chemother* 78:2385-6,
   doi:10.1093/jac/dkad217 — argued the fixed 4:1 vial ratio permits extrapolation. Cited here as
   ref 4. **Not yet archived; worth retrieving to see the argument in their own words.**
3. **This paper (2024)** answers with data: the ratio is not fixed, so extrapolation is unreliable.

The Fresan folder's README should be read alongside this one. Its note that Fresan "cannot speak to the
clearance-correlation question" remains correct, but the fuller point is that Fresan's design choice
became the subject of a published rebuttal — which is itself evidence for how live this project's
central question is in the clinical literature.

## Provenance classification

**DIRECTLY REPORTED — published cohort statistics, n=107 patients / 188 paired assessments.**
`Gatti2024_CAZ_AVI_ratio.csv` transcribes the cohort description, both drugs' median Css and total
clearance, the full ratio distribution, the renal-function split, all four ROC analyses, and —
valuably — **median total clearance for BOTH drugs across seven incremental CrCL strata**, which is
paired renal-function-resolved clearance data of a kind little else in this directory provides.

## An arithmetic check worth recording — the ratio spread does NOT contradict ρ = 0.94

The headline range "1.29:1 to 13.46:1" is easy to read as evidence that the two clearances are only
loosely coupled, which would sit awkwardly against this same group's published ρ = 0.94. **It does
not, and the arithmetic is worth writing down** (all inputs are published; the calculation is
elementary and approximate).

- The **min–max over 188 measurements is not a measure of typical spread.** The **IQR, 3.93:1 to
  5.70:1**, is the informative statistic: on the log scale that is SD ≈ 0.276.
- If the *entire* IQR spread were between-subject variability in clearance, then with Cojutti's
  ω_CAZ = 0.6159 and ω_AVI = 0.6817 it would imply **ρ ≈ 0.92** — remarkably close to the published
  0.94, from a completely independent statistic in a largely overlapping cohort.
- Since part of that spread must be measurement/residual error rather than true between-patient
  differences, the **true ρ implied is higher still**. So this paper's data are *consistent with, and
  weakly supportive of, the anchor value* — not in tension with it.

**A puzzle this raises, flagged rather than resolved.** Cojutti 2024's proportional residual terms are
b1 = 0.31 (ceftazidime) and b2 = 0.33 (avibactam). If the two analytes' residual errors were
*independent*, the ratio would inherit variance 0.31² + 0.33² = 0.205 (SD ≈ 0.45) **from residual
error alone** — which is already almost three times the total observed variance of 0.076. The observed
ratio is therefore **far more stable than independent residual errors would permit.** The natural
explanation is that the two errors are strongly positively correlated — same patient, same blood
sample, same LC-MS/MS run — so they largely cancel in the ratio. That is mechanistically very
plausible and, if true, matters for Model 2: **inferring one component from the other is more reliable
than an independent-assay-error model assumes.**

**Treat this as an observation, not a result.** The confounds are real: the two cohorts overlap but are
not identical; Cojutti's b terms are model residuals that absorb structural misspecification rather
than pure assay imprecision; the 188 assessments are repeated measures on 107 patients, not
independent draws; and the observed median is 4.67:1, not 4:1. It is worth a proper look if the
monitoring layer is ever revisited — `model2_monitoring.py` currently treats assay error as
independent across the two analytes.

## What this dataset is, and is not

**It is** direct clinical evidence that the 4:1 vial ratio does not hold in plasma, with a quantified
renal-function gradient (patients with CrCL > 80 are far more likely to have high ratios — 59.3% vs
23.8% above 5:1), a mechanistic explanation (avibactam undergoes active tubular secretion of ~27% via
OAT1/OAT3 in addition to filtration, ceftazidime filtration only), and independent replication of the
gradient Cojutti 2024 found across eGFR classes (3.47 → 5.47).

**It is not** patient-level data, and not an analysis of attainment or outcome. It reports a
distribution of ratios, not whether patients hit targets. Total clearances are back-calculated as
infusion rate / Css. Concentrations are **total, not free**. CRRT patients were excluded. Single-centre,
and the authors themselves call the findings "preliminary", noting that renal function changed during
treatment in about half the ICU patients and that the ratio spread under *intermittent* infusion
remains unconfirmed.

**Its conclusion is directionally consistent with this project's R7 breaking-point result** (measuring
beats inferring at every correlation from 0.30 to 0.99) — but note that R2's triage analysis reaches a
more refined answer than "always measure both": at ρ = 0.94, ~12.5% of patients measured captures 90%
of the benefit, and selective measurement can beat measuring everyone. This paper argues the binary
case; R2 answers the "which patients" question the paper leaves open. Its ROC cut-offs (CrCL > 75–78
mL/min/1.73 m², serum urea ≤ 45–51 mg/dL) are a clinically-derived triage rule that could be compared
directly against R2's model-derived concentration windows — **that comparison has not been done and
would be a genuinely novel check.**

## Reuse conditions to honour

1. Cite Gatti et al. 2024 wherever these values are used.
2. **Do not quote "1.29:1 to 13.46:1" as a measure of typical variability** — it is a min–max over 188
   single measurements. Quote the IQR (3.93:1–5.70:1) alongside it.
3. State that concentrations are total, not free, and that clearances are back-calculated.
4. Published under OUP's Standard Journals Publication Model, "All rights reserved" — not openly
   licensed. Cite and link; no PDF is archived here and none should be redistributed.
