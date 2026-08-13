# Provenance record — Tian 2025 CVVH pharmacokinetics

## Source

Tian S, Chen Y, Qiu M, Wu W, Dou L, Wang J, Xu L, Zhou Z, Wu M, Li J, Wu X, Ren J.
*PK/PD study of ceftazidime/avibactam in patients with severe intra-abdominal infections treated by
continuous veno-venous hemofiltration.*
**Eur J Clin Microbiol Infect Dis. 2026;45:711-721** (published online 14 November 2025).
doi:10.1007/s10096-025-05343-x

## Legal basis for use — DIFFERENT FROM EVERY OTHER FILE IN THIS DIRECTORY

**This is NOT an open-access article.** It was retrieved 12 August 2026 via the user's own New York
University institutional library subscription, accessed through the user's personal, already
logged-in browser session (Claude only navigated pages the user was already authenticated on; no
credentials of any kind were entered or handled). The page confirmed "Access provided by NEW YORK
UNIVERSITY LIBRARIES" before the PDF was downloaded.

This is legitimate personal access for the user's own scholarly research, exactly as reading the paper
in a library would be — **but it is not a licence to redistribute the PDF.** Consequently:

- The PDF (`Tian2025_EJCMID_s10096-025-05343-x_SUBSCRIPTION-ACCESS.pdf`) is kept in this **local
  development directory only**, for the user's own reference and for building the extracted-data CSV
  below.
- **This PDF, and this PDF alone in `data_external/`, must be excluded from every GitHub push**, on
  the same footing as `audit/extracted/` and `correspondence/` (see the top-level `README.md`). It is
  not open access and Springer's copyright applies in full.
- The **extracted numeric data** (below) is factual data reported in the article — facts are not
  copyrightable — but should still be cited to the source, not presented as if newly generated.

## Provenance classification

**DIRECTLY REPORTED — published summary PK data (n=7 patients, 6 sampling occasions).**
Not digitized, not model-inferred, not simulated. Every value in `Tian2025_PK_parameters.csv` was read
from **Table 2** of the article: steady-state NCA PK parameters (Cmax, Cmin, t1/2, AUC0-8, CL, Vss) for
ceftazidime and avibactam, each reported during vs. outside CVVH periods, mean ± SD.

**Not transcribed into CSV** (available in the PDF if needed later): Table 1 (patient demographics,
n=7 — mean age 56.4±15.5 y, median eGFR 28.81 mL/min/1.73m², E. coli/K. pneumoniae/P. aeruginosa
isolates) and Table 3 (Monte Carlo-simulated CFR by dosing regimen and pathogen — the paper's own
prescriptive analysis, not raw PK data, and out of scope for a "what was measured" extraction).

## What this dataset is, and is not

**It is** a third independent CVVH cohort (after Gatti 2023's CVVHDF, Li 2025's CRRT, and O'Jeanson
2024's CVVHDF, all elsewhere in this directory) — this one CVVH specifically, from Jinling Hospital,
Nanjing, in patients with severe intra-abdominal infection rather than the mixed indications in the
other cohorts. A further class-C consistency check.

**It is not** validation of the primary (non-RRT) scenario — every patient here is on CVVH, same
caveat as every other RRT dataset in this directory. **n=7 patients, 6 sampling occasions, single
centre, single short period** — the authors' own stated limitations, quoted directly: "the study
cohort comprised only seven patients... this limitation would have affected the reliability of our
findings" and results are explicitly labelled "hypothesis-generating," not validated. The article's own
Data Availability statement: *"No datasets were generated or analysed during the current study"* —
unlike O'Jeanson 2024 or the Gatti 2023 case series, there is no path to obtaining patient-level data
for this cohort even on request; the numbers in Table 2 are the ceiling of what is available.

## Reuse conditions to honour

1. Cite Tian et al. 2025 wherever these values are used.
2. State the CVVH population restriction alongside every result derived from them.
3. **Do not push `Tian2025_EJCMID_s10096-025-05343-x_SUBSCRIPTION-ACCESS.pdf` to GitHub or any public
   location.** The extracted CSV and this README may be shared; the PDF may not.
