# Provenance record — Gatti 2025: the clinical outcome paper behind novelty route R5

## Source

Gatti M, Rinaldi M, Cojutti PG, Bonazzetti C, Siniscalchi A, Tonetti T, Ambretti S, Tedeschi S,
Giannella M, Viale P, Pea F. *A pre-post quasi-experimental study of antimicrobial stewardship
exploring the impact of a multidisciplinary approach aimed at attaining an aggressive joint
pharmacokinetic/pharmacodynamic target with ceftazidime/avibactam on treatment outcome of
KPC-producing Klebsiella pneumoniae infections and on ceftazidime/avibactam resistance development.*
**Antimicrob Agents Chemother. 2025;69(7):e0048825.** doi:10.1128/aac.00488-25 · PMID 40476843 ·
PMC12217479

Bologna group again — same authors as `Cojutti2024_ANCHOR_PopPK/` and `Gatti2023_*`, and the addressee
of the data request sent 11 August 2026.

## Legal basis for use — genuinely open access, confirmed

**CC BY 4.0**, confirmed programmatically against PMC (`license.type: "CC BY 4.0"`,
`is_open_access: true`, https://creativecommons.org/licenses/by/4.0/), "Copyright © 2025 Gatti et al."
Retrieved 12 August 2026 as full text from PubMed Central. CC BY permits redistribution with
attribution and carries **no** non-commercial or no-derivatives restriction — so unlike the Tian, Wu,
Lanini, Fresan and Cojutti folders, **nothing here is excluded from a GitHub push on copyright
grounds.**

## Why this paper: it is exactly what novelty route R5 is about

`NOVELTY_ROUTES.md` R5 describes "Gatti 2025 (*AAC* 69(7):e00488-25, n = 218) ... the only clinical
exposure-response signal for the joint target: microbiological failure odds ratio 0.03, resistance odds
ratio 0.07 for patients attaining it." **Every element of that description is confirmed correct
against the source**: n = 218 (116 pre + 102 post), OR 0.03 (95% CI 0.005–0.20) for microbiological
failure, OR 0.07 (95% CI 0.01–0.69) for 90-day resistance development. The route's premise was
accurately recorded.

R5 remains **on hold** — it needs individual attainment/outcome data, which this paper does not
publish, and obtaining it would mean a *second* request to Bologna while the first is outstanding. The
standing rule (no follow-up, no second request without explicit instruction) is unchanged by archiving
this paper. What this folder does is make the published half of R5 available now, so that "what is
possible without the data" — R5's own stated fallback — can be worked on whenever wanted.

## Provenance classification

**DIRECTLY REPORTED — published cohort outcome statistics and multivariate odds ratios, n = 218.**
Not digitized, not model-inferred. `Gatti2025_outcome_associations.csv` transcribes the cohort
description (Table 1), both microbiological-failure models (Tables 2 and 3), both resistance-development
models (Tables 4 and 5), the clinical-failure models reported in the text, the TDM/target-attainment
figures, and the explicit definition of the aggressive joint PK/PD target.

## What this dataset is, and is not

**It is** the only clinical exposure-response evidence in this project's entire evidence base linking
**joint** (both-component) PK/PD target attainment to hard outcomes. Its target definition —
fCssCAZ/MIC > 4 **and** fCssAVI/CT > 1 with CT = 4 mg/L — is a real-world operationalisation of exactly
the joint-attainment construct Model 2 formalises. Effect sizes are large and consistent across three
independent outcomes (microbiological eradication, clinical cure, 90-day resistance).

**Two findings that bear directly on this project's own analyses:**
- **CRRT was an independent predictor of microbiological failure** in the pre-intervention phase
  (OR 5.20; 1.21–22.34). Set against Chen 2025, which found CRRT *not* a significant covariate on
  clearance, this is a useful reminder that RRT can matter clinically even where it does not move a PK
  parameter — the two findings are not in conflict, they are about different things.
- **Augmented renal clearance predicted clinical failure** (OR 8.34) and was a stated reason for
  target non-attainment — converging with Lanini 2024's finding that CrCl ≥130 mL/min carries a high
  risk of sub-target exposure. Two independent cohorts, same direction.

**It is not** a validation set, and must not be read as causal evidence for the target. This is a
**pre-post quasi-experimental design**: the intervention bundled several changes at once — continuous
infusion rose from 31.9% to 96.1%, combination therapy fell from 67.2% to 15.7%, treatment duration
shortened from 14 to 10 days, all alongside TDM-guided dosing. The outcome improvements cannot be
attributed to joint-target attainment alone, and the authors do not claim otherwise. Secular trend over
2018–2024 is uncontrolled. Single centre, retrospective, ~15% of patients lacked follow-up cultures
(a selection bias the authors flag). Several confidence intervals are very wide (ARC 1.00–69.27;
HAP/VAP 1.55–87.06), reflecting small event counts.

**Free fractions were calculated, not measured** (×0.90 CAZ, ×0.93 AVI). Note these **differ from the
multipliers used by the same group in Cojutti 2024** (×0.85 and ×0.92) — see that folder's README. The
literature, including within one research group, is not internally consistent about avibactam protein
binding, and this project should not present any single pair as settled.

## Reuse conditions to honour

1. Cite Gatti et al. 2025 wherever these values are used.
2. **Never quote the odds ratios as evidence that attaining the target causes better outcomes** — the
   pre-post design bundles multiple co-interventions. Describe them as associations within the
   post-intervention phase.
3. State the free-fraction multipliers used, and that they differ from Cojutti 2024's.
4. CC BY 4.0: redistribution is permitted with attribution. This folder has no push restriction.
