# Provenance record — Chen 2025 population PK, critically vs non-critically ill

## Source

Chen Y, Chen B, Huang Y, Li X, Wu J, Lin R, Chen M, Liu M, Qiu H, Cheng Y. *Population
Pharmacokinetics-Based Evaluation of Ceftazidime-Avibactam Dosing Regimens in Critically and
Non-Critically Ill Patients With Carbapenem-Resistant Klebsiella pneumoniae.*
**Infect Drug Resist. 2025;18:941-955.** doi:10.2147/IDR.S495279 · PMID 39990787 · PMC11846486

## Legal basis for use — open access via PMC, exact licence not confirmed

Retrieved 12 August 2026 as **full text from PubMed Central** (PMC11846486) via the NCBI E-utilities
API — PMC only serves full text for articles the publisher has deposited for public access, and Dove
Medical Press (the publisher of *Infection and Drug Resistance*) is a fully open-access publisher whose
articles normally carry a Creative Commons licence (typically CC BY-NC).

**However, the specific licence could not be confirmed programmatically:** PubMed's copyright endpoint
returned only "© 2025 Chen et al." with a null licence type and `is_open_access: false`, which is a
metadata gap rather than a positive statement that the article is closed. **Verify the licence at
https://doi.org/10.2147/IDR.S495279 before redistributing anything beyond the extracted facts below.**

No PDF is archived in this folder — only the extracted numeric parameters, which are facts and not
subject to copyright regardless of how the licence question resolves. This folder therefore raises no
redistribution question of the kind attached to the Tian/Wu/Lanini PDFs.

## Provenance classification

**DIRECTLY REPORTED — a fitted population PK model.**
`Chen2025_PopPK_parameters.csv` transcribes the structural and covariate model parameters, the cohort
description, the assumed protein-binding fractions, the two PK/PD target definitions, and the
healthy-volunteer comparator values the authors cite — all read from the article's Results and
Discussion text and its printed model equations.

**Not transcribed:** the bootstrap confidence intervals and the full parameter tables (the PMC text
extraction rendered table contents as figure/table *references* without their numeric bodies, so the
per-parameter RSEs and bootstrap CIs are not available from the source used). The article states all
fixed-effect RSEs were <30% and that all final estimates fell within bootstrap 95% CIs with <15%
deviation from bootstrap medians, but the individual values would need the PDF or the publisher HTML
tables. **This is a real gap versus the Wu 2025 folder, which does have per-parameter RSEs and CIs.**

## What this dataset is, and is not

**It is** a second independently-fitted population PK model for both ceftazidime and avibactam (after
Wu 2025), from a different Chinese centre, on a different population, with one property none of the
other datasets in this directory has: **it deliberately spans and distinguishes critically ill
(APACHE II > 15) and non-critically ill (APACHE II ≤ 15) patients**, and assigns them different PK/PD
targets by design. Its structural findings are directly comparable to Wu 2025's — both chose a
one-compartment model for both drugs, both found CrCL the only significant covariate on clearance —
which makes the pair a genuine (if small) between-study replication:

| | Chen 2025 | Wu 2025 |
|---|---|---|
| CAZ CL (L/h) | 2.96 | 1.60 |
| AVI CL (L/h) | 3.09 | 2.74 |
| CAZ V (L) | 17.76 | 15.90 |
| AVI V (L) | 18.25 | 28.00 |
| CrCL exponent, CAZ | 0.44 | 0.705 |
| CrCL exponent, AVI | 0.41 | 0.792 |
| Cohort median CrCL (mL/min) | 71.3 | 47.8 |

The two disagree substantially on the CrCL exponent and on avibactam volume — a useful, citable
illustration of between-study parameter heterogeneity in exactly the quantities this project's Layer 2
heterogeneity analysis is about.

**A notable negative finding worth carrying forward:** 15 of 45 patients (33.3%) were on CRRT, and CRRT
was **not** a significant covariate on clearance in either drug's model. The authors attribute this to
their CRRT patients retaining substantial residual renal function (median CrCL 69.3 mL/min), and note
that other work suggests CRRT affects clearance mainly below CrCL ~10 mL/min. This cuts against the
assumption — implicit in treating RRT cohorts as a separate population throughout this project — that
RRT status is always the dominant determinant.

**It is not** a validation set for the primary model. Single centre, n=45, 91 concentrations, sparse
sampling (1-3 per patient), 64.4% of the cohort with renal insufficiency (the authors' own stated
generalisability limitation), no MIC data for the isolated pathogens (so PTA could not be linked to
clinical outcome), and free fractions **assumed** (90% CAZ, 92% AVI) rather than measured.

## Reuse conditions to honour

1. Cite Chen et al. 2025 wherever these values are used.
2. Confirm the licence before redistributing anything beyond the extracted numeric facts.
3. State that free fractions are assumed literature values, not measured.
4. If the per-parameter RSEs or bootstrap CIs are needed, they must be sourced from the PDF or the
   publisher's HTML tables — they are not in this CSV, and no estimate of them should be invented.
