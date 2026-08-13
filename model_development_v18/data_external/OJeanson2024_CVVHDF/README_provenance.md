# Provenance record — O'Jeanson 2024 CVVHDF pharmacokinetics

## Source

O'Jeanson A, Ioannidis K, Nielsen EI, Galani L, Ginosyan A, Paskalis H, Loryan I, Giamarellou H,
Friberg LE, Karaiskos I. *Ceftazidime-avibactam (CAZ-AVI) pharmacokinetics in critically ill patients
undergoing continuous venovenous hemodiafiltration (CVVHDF).*
**Int J Antimicrob Agents. 2025;65(1):107394.** doi:10.1016/j.ijantimicag.2024.107394 · PMID 39581557 ·
PMC no PMCID assigned at time of retrieval; Scopus EID 2-s2.0-85211985518.

## Legal basis for use

**Open access, CC BY 4.0.** Verbatim from the article's own metadata (Elsevier Article Retrieval API,
`openaccessType: "Full"`, `openaccessUserLicense: http://creativecommons.org/licenses/by/4.0/`):

> "© 2024 The Author(s). Published by Elsevier Ltd. This is an open access article under the CC BY
> license (http://creativecommons.org/licenses/by/4.0/)"

Retrieved 12 August 2026 via Elsevier's official Article Retrieval API
(`api.elsevier.com/content/article/pii/S0924857924003108`), using an API key the user holds and
authorized for this research. No login, no paywall bypass, no scraping — this is the publisher's own
documented API, and the article is genuinely open access regardless of API-key use. Archived here as
`OJeanson2024_IJAA_107394_CC-BY.pdf` (604,027 bytes, matches the size Elsevier's own metadata reports
for the file). CC BY permits redistribution, including commercial, with attribution — no non-commercial
or no-derivatives restriction, unlike the Gatti 2023 file elsewhere in this directory.

## Provenance classification

**DIRECTLY REPORTED — published individual-level (n=4) and summary PK data.**
Not digitized, not model-inferred, not simulated. Every value in `OJeanson2024_PK_parameters.csv` was
read from the article's own printed tables:

| Rows | Source | Content |
|---|---|---|
| `table = Table1` | **Table 1** | demographics/labs, n=4, mean (SD), median, range |
| `table = Table2` | **Table 2** | steady-state NCA PK parameters (Cmax, Cmin, t1/2, CLss, Vss, AUCtau) for CAZ and AVI, at both dose groups (2000/500 mg and 1000/250 mg q8h), median [90% CI], for both the observed cohort and a n=10,000 simulated phase-III (non-RRT) comparator population |

**Figure 1 (concentration-time comparison plot) was visually inspected but not digitized** — Table 2
already reports every NCA-derived parameter Figure 1 illustrates, to full numeric precision. Digitizing
the figure would add estimation error, not information, given the table exists.

## What this dataset is, and is not

**It is** a second, independent (non-Bologna, non-China) CVVHDF cohort, from a European centre
(Athens), with a transparent within-study comparison against a simulated non-RRT phase-III population
using the same Li et al. 2019 PopPK model this project's own primary comparisons are built against —
useful as a further class-C consistency check, not as a validation set.

**It is not** validation of the primary (non-RRT) scenario — every patient here is on CVVHDF, same
caveat as the Gatti 2023 and Li 2025 datasets in this directory. **n=4 patients, 6 dosing occasions.**
Nothing fitted to this dataset can support a population PK model or be described as external clinical
validation of the primary model. The paper's own Discussion notes an implausible CAZ unbound fraction
(median 1.07, capped at 1.0 in analysis) and higher-than-typical Vss values, both flagged by the
authors themselves as needing further mechanistic explanation — carry that caveat forward with any use
of these numbers.

## Reuse conditions to honour

1. Cite O'Jeanson et al. 2024 wherever these values are used.
2. State the CVVHDF population restriction alongside every result derived from them (same as every
   other RRT dataset in this project).
3. CC BY permits redistribution with attribution; still, don't present the transcription as
   author-supplied raw data — it is values read from the published tables.
