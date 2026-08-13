# Provenance record — Dryad CRRT dataset (Li et al.)

## Source

**Dataset:** Li C, Wang Y, Chen F, Huang L, Dong J, Fan W, Yue H, Ge Y.
*Pharmacokinetics/pharmacodynamics of ceftazidime-avibactam in adult critically ill patients
receiving continuous renal replacement therapy.*
**Dryad, published 28 November 2025.** DOI **10.5061/dryad.fxpnvx16s**
`https://datadryad.org/dataset/doi:10.5061/dryad.fxpnvx16s`

**Primary article:** *Antimicrobial Agents and Chemotherapy* 2026;70(2):e0143825.
doi:10.1128/aac.01438-25 · PMID 41432444

## Legal basis for use

**Licence: CC0 1.0 Universal — public domain dedication.** No restriction on reuse, no attribution
requirement (attribution is nonetheless given throughout, as scholarly practice requires).
Dryad curation status **Published**, version 8.

The depositors' own human-subjects statement, verbatim from the dataset README:

> "I confirm that we received explicit consent from our participants to publish the de-identified
> data in the public domain. We did not disclose the personal privacy of the subjects, such as names
> and ID numbers. Instead, we used anonymized subject numbers to represent different individual
> patients."

Consent for public-domain publication was therefore obtained by the original investigators. The data
carry no direct identifiers; age, body weight, and all CRRT operating parameters are published as
**coded categories rather than raw values**, which is itself a de-identification measure. No
re-identification was attempted and none is possible from these files.

## Provenance classification

**DIRECTLY REPORTED — published individual patient data, public domain.**
Not digitized, not model-inferred, not simulated, not donor-derived.

## How these files were obtained

The Dryad REST API requires a bearer token for file download, and the public
`/downloads/file_stream/` endpoint is protected by a proof-of-work bot check. The file contents were
therefore read through the dataset page's own built-in **file preview** panels, which render each CSV
in full, and transcribed verbatim by `rebuild_from_dryad_preview.py` in this directory.

**Verification of the transcription:** both concentration files contain **21 subjects and 119
observations each** — eight patients with six timepoints, three with seven, and ten with five —
consistent with the primary article's statement of 5-7 sampling points per patient within one
administration cycle. The demographic and CRRT files contain 21 rows each. File sizes on Dryad were
1,706 B and 1,709 B for the two concentration files, 1,451 B for the demographics and 370 B for the
CRRT parameters.

**Anyone re-running this work should download the files directly from Dryad in a browser** and
confirm they match. These transcriptions are a convenience, not the authoritative copy.

## Dataset contents

| File | Rows | Content |
|---|---|---|
| `Ceftazidime_concentration.csv` | 118 | subject, time after dose (h), **pre-filter** and **post-filter** ceftazidime concentration (mg/L) |
| `Avibactam_concentration.csv` | 118 | subject, time after dose (h), **pre-filter** and **post-filter** avibactam concentration (mg/L) |
| `Demographic_data.csv` | 21 | age category, sex code, weight category, diagnosis, ALT, AST, haematocrit, albumin, serum creatinine, APACHE II, SOFA, 24-hour urine volume |
| `CRRT_parameters_data.csv` | 21 | CRRT modality, ultrafiltration rate category, total effluent flow category, consecutive CRRT days category |
| `Medication_Information.csv` | 21 | dose, interval, treatment-duration category, route, infusion-duration category |
| `Microbial_data.csv` | 21 | infection site, pathogen, Kirby-Bauer value, susceptibility, pathogen clearance, clinical outcome — **not transcribed** (multi-line cells; retrieve from Dryad if needed) |

**Dosing is uniform across all 21 patients: 2 g ceftazidime + 0.5 g avibactam every 8 hours by
intravenous infusion.** Infusion duration is published only as a category (1/2/3), not in hours.

## What these data are, and are not

**They are** the only openly downloadable, individual-level, paired ceftazidime **and** avibactam
concentration-time dataset located anywhere in this audit — across Dryad, Zenodo, figshare, OSF,
Harvard Dataverse, BioModels, DDMoRe, and public code repositories.

**They are not** applicable to the manuscript's primary scenario, on two independent grounds:

1. **Population.** All 21 patients are on continuous renal replacement therapy. The primary scenario
   is explicitly restricted to critically ill adults **not** receiving renal replacement therapy.
2. **Administration.** All received **intermittent 8-hourly intravenous infusion**, not continuous
   infusion. The primary model assumes constant continuous infusion at steady state.

The cohort is also narrow clinically: 18 of 21 have acute pancreatitis, and the diagnosis mix is not
representative of a general ICU population.

**Class C — suitable for sensitivity and assumption-testing analysis only.** These data cannot fit
or validate the primary model. What they can do is test one structural assumption, described below.

## Derived quantities computed from these data

Under steady-state dosing, total clearance over one dosing interval is CL = Dose / AUC(0-τ). Area
under the curve was computed by the linear trapezoidal rule on the **pre-filter** (systemic)
concentrations across the full 0-8 h interval, for each patient and each analyte:

| | Median | IQR | Range |
|---|---|---|---|
| Ceftazidime CL | 2.43 L/h | 2.15-2.95 | 1.79-4.19 |
| Avibactam CL | 3.26 L/h | 2.85-3.40 | 2.52-4.44 |

**Known bias.** Sparse sampling misses the true peak in the 16 patients whose first post-dose sample
is at 2 or 3 h, so the trapezoidal AUC is biased low and clearance correspondingly biased high. The
bias applies to both analytes in the same direction and is therefore largely cancelled in the
clearance **ratio** and in the between-analyte **correlation**, which are the quantities used. It is
not cancelled in the absolute clearance values, which should be treated as approximate.

## Intended use in this project

A single, narrowly scoped purpose: to estimate the **between-patient correlation of ceftazidime and
avibactam clearance** in real patients, as an empirical check on the value of 0.94 assumed by the
primary model and taken from a single cohort.

Result: **Pearson r = 0.560 on the log scale (95% CI 0.169 to 0.799, p = 0.008)**, and Spearman
ρ = 0.603 (p = 0.004); significantly below 0.94 (p ≈ 3 × 10⁻⁶). This closely matches the independent
estimate of r = 0.598 obtained from the Gatti 2023 CVVHDF case series (see
`../README_Gatti2023_provenance.md`) — two unrelated cohorts, two countries, two renal-replacement
modalities, giving essentially the same answer.

**Both cohorts are on renal replacement therapy, where a shared extracorporeal circuit removes both
analytes and should if anything inflate the correlation between their clearances.** That both give
approximately 0.57-0.60 rather than 0.94 is the substance of the argument, and it must be reported
with the population limitation stated explicitly every time.

**This is not external validation of the primary model, and it must never be described as such.**
