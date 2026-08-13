# Provenance record — Lanini 2024 non-RRT CrCl correlation

## Source

Lanini S, Giuliano S, Angelini J, Ferin S, Martini L, Baraldo M, Cossettini S, Roberts J, Tascini C.
*Renal function and its impact on the concentration of ceftazidime-avibactam: A cross-sectional study.*
**Int J Antimicrob Agents. 2024;64(6):107351.** doi:10.1016/j.ijantimicag.2024.107351 · PMID 39362612

A reply to a subsequent Letter to the Editor exists (Angelini et al. 2025, doi:10.1016/j.ijantimicag.2025.107480)
but was not retrieved — it responds to a methodological critique, not additional data.

## Legal basis for use — subscription access, same footing as Tian2025_CVVH/ and Wu2025_PopPK_AKI/

**Confirmed NOT open access** via Elsevier's own Article Retrieval API metadata
(`"openaccess": "0"`, `"openaccessArticle": false`, `"openaccessType": null`) — unlike O'Jeanson 2024
elsewhere in this directory, this article is subscription-only even through Elsevier's official channel.

Retrieved 12 August 2026 via the user's own New York University institutional library subscription. The
first attempt hit a Cloudflare bot-detection CAPTCHA on ScienceDirect; Claude does not solve CAPTCHAs
under any circumstance, and stopped there. **The user solved the CAPTCHA themselves, in their own
browser, and supplied the resulting PDF and its Elsevier supplementary file (`mmc1.docx`) directly.**
No credentials were entered or handled by Claude at any point in this retrieval.

**This PDF must be excluded from every GitHub push**, on the same footing as `Tian2025_CVVH/` and
`Wu2025_PopPK_AKI/` — see the top-level `README.md`'s exclusion list. The extracted numeric data below
(facts, not the copyrighted text) may be shared and cited; the PDF may not be redistributed.

## Provenance classification

**DIRECTLY REPORTED — published cross-sectional correlation statistics, n=52.**
Not digitized, not model-inferred, not simulated. Every value in `Lanini2024_CrCl_correlation.csv` was
read from the article's own text and Table 1 — three simple linear regressions (free-CAZ vs. CrCl,
free-AVI vs. CrCl, CAZ:AVI ratio vs. CrCl) with Pearson r, R², back-transformed slope and 95% CI, plus
a 2×2 risk table (sub-target exposure below vs. above the augmented-renal-clearance threshold of 130
mL/min) for both drugs.

**The Elsevier supplementary file (`mmc1.docx`, Figure 3 at higher resolution) was checked and adds
no data beyond the main text** — its caption is verbatim identical to Figure 3's in the main PDF, and
the embedded image is the same scatter plot already fully summarised by the regression statistics
above. Not separately archived.

## What this dataset is, and is not — the key difference from every other file in this directory

**This is the only non-RRT dataset in `data_external/`.** The study explicitly states: "We exclude
form [sic] that analyses all patients who received renally adjusted dosing and those who undergo
dialysis" — every other file here (Gatti 2023, Li 2025, O'Jeanson 2024, Tian 2025) is RRT-specific and
therefore excluded from the manuscript's primary scenario population by definition. This one is not:
n=52 critically ill patients on **standard-dose continuous-infusion CAZ-AVI**, explicitly excluding
both RRT and renally-adjusted dosing — the closest population match in this project's external evidence
to the manuscript's own primary scenario.

**It is still not a validation set** for the primary model, for several independent reasons stated
plainly by the authors themselves: (1) cross-sectional design — one concentration measurement per
patient at day 3, not a full concentration-time profile; (2) free-fraction values are *derived*, not
measured — total plasma concentration corrected by a single literature protein-binding value (10% CAZ,
7% AVI) applied to every patient, not patient-specific; (3) single-centre design, explicitly flagged by
the authors as limiting generalisability; (4) the correlation quantifies *exposure* vs. renal function,
not *attainment* against a joint PK/PD target the way this project's own analyses do.

**What it is directly useful for:** a genuine, independent, non-RRT confirmation that the qualitative
direction and rough magnitude of the renal-function/exposure relationship this project's own primary
model assumes (higher CrCl -> lower concentration, avibactam more sensitive to CrCl than ceftazidime)
holds in a real non-RRT critically-ill cohort — and a caution, independent of anything already in this
project, that augmented renal clearance (CrCl >130 mL/min, present in 30.8% of this cohort) carries a
substantial, quantified risk of sub-target exposure to both components (68.75% for free-CAZ, 56.25% for
free-AVI) even at the standard dose.

## Reuse conditions to honour

1. Cite Lanini et al. 2024 wherever these values are used.
2. State explicitly that free concentrations are derived (not measured) using a single literature
   protein-binding value, not patient-specific measurements.
3. Do not describe this as validation of the primary model — it is an independent exposure/renal-function
   correlation in a different, if closely related, cohort.
4. **Do not push the PDF (or the `mmc1.docx` original, wherever it ends up) to GitHub or any public
   location.** The extracted CSV and this README may be shared; the PDF may not.
