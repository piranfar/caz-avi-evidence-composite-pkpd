# Provenance record — Wu 2025 population PK / AKI risk model

## Source

Wu T, Ding Q, Huang S, Zhang S, Yang R, Qin Y, Liu J, Pei Q. *Model-informed individualized
administration of ceftazidime-avibactam in critically ill patients: population pharmacokinetics
studies and a parametric time-to-event analysis.*
**J Antimicrob Chemother. 2025;80(10):2693-2704.** doi:10.1093/jac/dkaf275 · PMID 40795198

## Legal basis for use — subscription access, same footing as Tian2025_CVVH/

Retrieved 12 August 2026 via the user's own New York University institutional library subscription
(Oxford Academic / British Society for Antimicrobial Chemotherapy journal), through the user's
personal, already-authenticated browser session. No credentials were entered or handled by Claude.

**This PDF must be excluded from every GitHub push**, exactly like
`Tian2025_CVVH/Tian2025_EJCMID_..._SUBSCRIPTION-ACCESS.pdf` — see the top-level `README.md`'s
exclusion list. The extracted numeric data below (facts, not the copyrighted text) may be shared and
cited; the PDF itself may not be redistributed.

## Provenance classification

**DIRECTLY REPORTED — a real fitted population PK model, not a summary-statistics table.**
This is qualitatively different from every other file in `data_external/`: those report NCA summary
statistics (Cmax, Cmin, CL, AUC as mean/median with a range or CI) from small case series. This paper
reports an actual **structural + covariate population PK model**, fitted in NONMEM with FOCE-I,
externally validated on an independent 32-patient cohort, with bootstrap 95% CIs on every parameter —
methodologically the same class of evidence as the primary model this project audits and extends, just
from a different population.

`Wu2025_PopPK_parameters.csv` — every value read from the article's own tables and equations:

| `record_type` | Source | Content |
|---|---|---|
| `structural` | Table 2, Table 3, and the printed model equations | CL, V typical values, CrCL and CRP covariate exponents, interindividual variance (ω²), residual variance (σ²), all with RSE% and bootstrap 95% CI |
| `exposure` | Table 4 | Monte Carlo-simulated steady-state AUCss,24h, Css,max, Css,min for both drugs, median (range), across the 164 patients' own actual regimens |
| `aki_tte` | Results text and the printed AKI hazard equation | the sigmoidal-Emax EC50 (7450 μg·h/mL) linking ceftazidime AUCss,24h to AKI hazard, and the Gompertz baseline-hazard time-to-event structure; **both RSEs are genuinely high (45.6%, 52.7%) and the authors themselves caution this is a risk-alert tool, not a decision rule** — carry that forward with any use |
| `dosing_recommendation` | Table 5 | the paper's own renal-function- and ventilation-stratified dosing table |
| `demographics` | Table 1, Results text | cohort sizes and composition for the PopPK derivation (n=31), validation (n=32), and AKI analysis (n=159) sub-cohorts |

**Not transcribed**: full goodness-of-fit / VPC figure data (Figures 1-3) and the supplementary tables
S1-S8 (bioanalytical detail, base-model parameters before covariate selection, AKI covariate screening
steps) — available in the archived PDF if needed later.

## What this model is, and is not

**It is** an independently-fitted structural population PK model for ceftazidime and avibactam in a
different population (critically ill Chinese adults, n=31 development + 32 validation, median CrCL
47.8 mL/min — a broad renal-function spectrum, NOT restricted to any single RRT modality) from a
different research group, using different software conventions (NONMEM, one-compartment structural
model for both drugs) than this project's own primary model. Its clearance estimates are markedly
lower than Western PopPK models (Lodise 2022, Li 2019) the paper itself cites and compares against,
and the authors attribute this to their cohort's renal impairment — a genuine, citable data point about
between-population PK variability, independent of anything already in this project's evidence base.

**It is not** a validation set for the primary model, and should not be described as one. It is a
different model fitted to different data for a different (if related) clinical question — dosing
against an AKI-risk endpoint, not against a joint attainment target. Sample size (31 development
patients, single centre, single country) is the authors' own stated limitation.

## Reuse conditions to honour

1. Cite Wu et al. 2025 wherever these values are used.
2. Do not describe the AKI EC50 as precise — the authors explicitly flag high parameter uncertainty
   (RSE 45.6-52.7%) and recommend it as a risk-alert tool only.
3. **Do not push the PDF to GitHub or any public location.** The extracted CSV and this README may be
   shared; the PDF may not.
