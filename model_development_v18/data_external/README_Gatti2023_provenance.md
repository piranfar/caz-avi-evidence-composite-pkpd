# Provenance record — Gatti 2023 individual patient data

## Source

Gatti M, Rinaldi M, Gaibani P, Siniscalchi A, Tonetti T, Giannella M, Viale P, Pea F.
*A descriptive pharmacokinetic/pharmacodynamic analysis of continuous infusion ceftazidime-avibactam
for treating DTR gram-negative infections in a case series of critically ill patients undergoing
continuous veno-venous haemodiafiltration (CVVHDF).*
**J Crit Care. 2023 Aug;76:154301.** doi:10.1016/j.jcrc.2023.154301 · PMID 37059003

## Legal basis for use

**Open access, CC BY-NC-ND 4.0.** Verbatim from the published article:

> "0883-9441/© 2023 The Authors. Published by Elsevier Inc. This is an open access article under the
> CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/)."

The Version of Record PDF was obtained without login or subscription from the University of Bologna
institutional repository (record `https://cris.unibo.it/handle/11585/929418`, file
`1-s2.0-S0883944123000503-main.pdf`, 520 kB) and is archived here as
`Gatti2023_JCritCare_154301_CC-BY-NC-ND.pdf`. Reuse for non-commercial research with attribution is
permitted by the licence; **the data must be cited, not redistributed as if newly generated**.

Ethical approval for the original study: local ethics committee, IRCCS Azienda Ospedaliero-Universitaria
di Bologna, **No. EM 232-2022_308/2021/Oss/AOUBo, 16 March 2022.** No new human-subjects
approval is required for secondary analysis of these already-published, de-identified,
non-identifiable aggregate-per-patient values.

## Provenance classification

**DIRECTLY REPORTED — published individual patient data.**
Not digitized, not model-inferred, not simulated, not donor-derived. Every value in
`Gatti2023_individual_patient_data.csv` was read from the printed numeric tables of the article:

| Rows | Source | Content |
|---|---|---|
| `record_type = tdm` (17 rows) | **Table 2** | one row per TDM occasion: patient ID, weight, CVVHDF operating settings, dose intensity, residual diuresis, total effluent flow, **ceftazidime CL (L/h)**, **avibactam CL (L/h)**, reduced-dose flag |
| `record_type = patient` (8 rows) | **Table 1** | one row per patient: age, sex, pathogen, **MIC**, infection type, dose and adjustment, **average free Css for both components**, target ratios, joint target status, treatment and CVVHDF duration, 30-day mortality |

**No figure digitization was performed and none was needed** — Figure 1 of the source plots exactly
the 17 (dose intensity, CL) pairs already printed in Table 2.

## Transcription verification

The transcription was verified against summary statistics the authors report independently in their
Results text. Recomputing from the transcribed table reproduces both exactly:

| Statistic | Published in the article | Recomputed from this file |
|---|---|---|
| Ceftazidime CL, median (IQR) | 2.39 L/h (2.05–2.94) | **2.39 L/h (2.05–2.94)** ✓ |
| Avibactam CL, median (IQR) | 2.56 L/h (2.22–2.96) | **2.56 L/h (2.22–2.96)** ✓ |
| n occasions / n patients | 17 / 8 | 17 / 8 ✓ |

An exact match on both medians and both interquartile ranges is strong evidence that no
transcription error was introduced.

## What these data are, and are not

**They are** the only publicly available, openly licensed, patient-level ceftazidime **and**
avibactam clearance pairs in critically ill adults receiving continuous infusion that this audit
could locate. They are the only dataset in which the ceftazidime–avibactam clearance correlation —
the parameter that drives the manuscript's "incremental value of an avibactam assay" analysis — can
be examined empirically rather than assumed.

**They are not** a validation set for the manuscript's primary scenario. The primary scenario is
explicitly restricted to critically ill adults **not** receiving renal replacement therapy; every
patient here is on CVVHDF. In this population, total clearance is dominated by the extracorporeal
circuit, which compresses between-patient variability (apparent CV of patient-mean log clearance
≈ 10% for ceftazidime and ≈ 13% for avibactam, against 67.9% and 76.9% in the non-RRT source
cohort) and induces a shared circuit-driven component in both analytes. Any inference drawn from
these data must be reported with that limitation stated plainly.

**Sample size is 8 patients and 17 occasions.** Nothing fitted to this dataset can support a
population PK model, and no analysis of it may be described as external clinical validation of the
primary model.

## Reuse conditions to honour

1. Cite Gatti et al. 2023 wherever these values are used.
2. Do not describe the values as author-provided or as raw data supplied by the investigators — they
   are values published in the article's own tables.
3. Non-commercial use only; no derivative redistribution of the PDF.
4. Report the CVVHDF population restriction alongside every result derived from them.
