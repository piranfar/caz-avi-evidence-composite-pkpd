# Provenance record — Fresan 2023 continuous-infusion TDM individual patient data

## Source

Fresan D, Luque S, Benítez-Cano A, Sorlí L, Milagro Montero M, De-Antonio M, Prim N, Vega V,
Horcajada JP, Grau S. *Pharmacokinetics/pharmacodynamics and therapeutic drug monitoring of
ceftazidime/avibactam administered by continuous infusion in patients with MDR Gram-negative bacterial
infections.*
**J Antimicrob Chemother. 2023;78(3):678-683.** doi:10.1093/jac/dkac439 · PMID 36626402

**This is the Benítez-Cano / Sorlí / Hospital del Mar (Barcelona) group** — the same team the project's
unsent second data request is addressed to
(`correspondence/DRAFT_data_request_Benitez-Cano_Sorli.md`). Read that file before any contact: the
standing rule is that no request goes out without explicit instruction.

## Legal basis for use — free-to-read, but NOT openly licensed

The article displays a green **"FREE"** badge on Oxford Academic and is flagged **Editor's Choice**; the
full text and both tables were readable without institutional sign-in (the page still offered "Sign in
through New York University Libraries", i.e. no subscription was consumed to read it). Retrieved 12
August 2026.

**Free to read is not the same as open access.** The article's own footer reads: "© The Author(s) 2023.
Published by Oxford University Press on behalf of British Society for Antimicrobial Chemotherapy. All
rights reserved. For permissions, please e-mail: journals.permissions@oup.com" — published under OUP's
Standard Journals Publication Model, which carries no reuse licence. The extracted numeric data below
(facts, not the copyrighted text) may be cited and shared; the article PDF may not be redistributed.

The Version of Record PDF is archived here as `Fresan2023_JAC_dkac439_FREE-TO-READ.pdf` (303,846 bytes;
its own footer records "Downloaded from academic.oup.com/jac/article/78/3/678/6982766 by guest on 12
August 2026" — *by guest*, confirming no subscription was used). Automated download failed twice; the
user downloaded it manually and supplied it. Because the article is free to read but **not** openly
licensed, the same rule as the Tian/Wu/Lanini folders applies: **exclude this PDF from any GitHub
push.** Anyone re-running this work can obtain it freely from https://doi.org/10.1093/jac/dkac439.

## Provenance classification

**DIRECTLY REPORTED — published individual patient data, n=31.**
Not digitized, not model-inferred, not simulated. `Fresan2023_individual_patient_data.csv` transcribes
**Table 2 in full** — one row per patient, every column: infection type, pathogen, MIC (with a flag for
the 13 patients where 8 mg/L was *assumed* rather than measured), eGFR, total daily dose, free
ceftazidime steady-state concentration, the free-Css/MIC ratio, the exposure band, the TDM-guided dose
recommendation, whether the treating physician accepted it, and 30-day all-cause mortality.

Table 1 (cohort-level summary) is not separately transcribed; its contents are recoverable from the
per-patient rows plus the article text.

## Transcription verification

Recomputed from the transcribed CSV against statistics the authors state independently in their Results
text:

| Statistic | Published in the article | Recomputed from this file |
|---|---|---|
| free-Css/MIC ratio, all 31 patients | (printed per patient) | **all 31 reproduce to <2%** ✓ |
| PK/PD target attained (≥4× MIC) | 26/31 (83.9%) | **26/31 (83.9%)** ✓ |
| Overexposed (>10× MIC) | 15/31 (48.4%) | **15/31 (48.4%)** ✓ |
| Dose reduction / maintenance recommended | 16 / 12 (of 28) | **16 / 12** ✓ |
| 30-day all-cause mortality | 6 (19.4%) | **7 (22.6%)** ✗ — see below |

**One discrepancy — traced to the source article itself, not to this transcription.** Every ratio and
every exposure-band count reproduces exactly. The mortality column of Table 2 carries "Y" for **seven**
patients (nos. 2, 5, 9, 10, 11, 18, 29), while the Results text states "Thirty-day all-cause mortality
was 19.4% (six patients)" — and 19.4% is 6/31, not 7/31 (22.6%), so the text is internally consistent
with six.

This was initially read from the publisher's HTML, so the obvious suspect was a rendering artifact.
**It is not.** The archived typeset PDF was checked row by row and prints the identical seven "Y"
values. **The inconsistency is in the published article: its Table 2 and its Results text disagree with
each other about how many patients died.** The CSV records the table as printed.

Consequence: do not use the mortality column without noting this, and do not quote "six deaths" and the
per-patient table together as if they agreed. The PK columns — which is what this project actually
needs — are verified sound and unaffected.

## What this dataset is, and is not — read this before using it

**The single most important limitation, stated by the authors themselves:** *"only ceftazidime
concentrations were measured, on the assumption that the concentration of avibactam would always be
sufficient to exert its action."* **There are no avibactam concentrations in this study at all.**

That means this dataset **cannot contribute to the clearance-correlation question** — the parameter at
the centre of this project (Model 1, and the manuscript's "incremental value of an avibactam assay"
analysis) requires paired ceftazidime *and* avibactam measurements in the same patient. This dataset has
only one of the two. It is, ironically, a direct real-world illustration of the very gap the manuscript
argues about: a well-run TDM programme that measured only the beta-lactam and assumed the inhibitor away.

**What it is genuinely useful for:**
- Individual patient-level free ceftazidime Css under **true continuous infusion** (12-hourly doses each
  infused over 12 h), a dosing mode matching the manuscript's primary scenario far better than the
  q8h intermittent infusion of most other datasets here.
- A wide, real eGFR spread (9-160 mL/min) with per-patient dose, letting the exposure/renal-function
  relationship be examined at patient level rather than through a fitted model's covariate term.
- Quantified evidence on **overexposure**: 48.4% of patients exceeded 10× MIC at doses often *below* the
  licensed 6 g/1.5 g daily, which is directly relevant to this project's exposure-ceiling/toxicity-screen
  layer.

**It is not** a validation set: retrospective, single centre, n=31, one TDM sample per patient at day 2,
free concentrations *derived* from an assumed 15% protein binding rather than measured, and 8 of 31
patients (25.8%) were on renal replacement therapy — so it is a mixed RRT/non-RRT cohort, not a clean
non-RRT one like Lanini 2024.

## Reuse conditions to honour

1. Cite Fresan et al. 2023 wherever these values are used.
2. **Always state that avibactam was never measured in this study** — omitting that would badly
   misrepresent what the data can support.
3. State that free concentrations are derived using an assumed 15% protein binding, not measured.
4. Do not use the mortality column without disclosing the article's own internal 6-vs-7 inconsistency
   (documented above, confirmed against the typeset PDF).
5. The article is free to read but not openly licensed — cite and link it; **do not push the archived
   PDF to GitHub or redistribute it.**
