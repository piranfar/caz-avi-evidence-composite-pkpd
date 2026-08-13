# Provenance record — Gatti 2023 (AAC): where the "aggressive joint PK/PD target" comes from

## Source

Gatti M, Rinaldi M, Bonazzetti C, Gaibani P, Giannella M, Viale P, Pea F. *Could an optimized joint
pharmacokinetic/pharmacodynamic target attainment of continuous infusion ceftazidime-avibactam be a way
to avoid the need for combo therapy in the targeted treatment of deep-seated DTR Gram-negative
infections?*
**Antimicrob Agents Chemother. 2023;67(11):e0096923.** doi:10.1128/aac.00969-23 · PMID 37843260 ·
PMC10648963

Bologna, same group as the other four Bologna papers in this directory.

**Note this is a different paper from the Gatti 2023 already archived** as
`Gatti2023_JCritCare_154301_CC-BY-NC-ND.pdf` (J Crit Care 76:154301, the CVVHDF case series). The
Bologna group published several papers in 2023; this one is the AAC study.

## Why this folder exists — it is the definitional root

**This is the paper that defined the "aggressive joint PK/PD target" that the rest of the chain uses.**
Both `Cojutti2024_ANCHOR_PopPK/` (its ref 9) and `Gatti2025_outcome_R5/` build directly on this
definition, and `Gatti2024_ratio_one_leg/` cites it as its ref 8. Until now the definition was in the
project through its downstream users; this records it at source.

The definition, verbatim:

> "The joint PK/PD target was considered optimal when both the *f*C_ss/MIC ratio for ceftazidime ≥4
> (equivalent to 100%T) and the *f*C_ss/C_T ratio for avibactam >1 (equivalent to 100%T > C_T of
> 4.0 mg/L) were simultaneously achieved (quasi-optimal if only one of the two and suboptimal if
> neither of the two was achieved)."

This is the structure Model 2 formalises: a **joint** requirement over two components, where failing
either one fails the pair. The three-level optimal / quasi-optimal / suboptimal grading is this
paper's, and it is what "attainment" means in Gatti 2025's outcome analysis.

## The headline result, and why it needs care

> "In the multivariate analysis, the suboptimal/quasi-optimal joint PK/PD target emerged as the only
> independent predictor of microbiological failure (odds ratio [OR] 11.11; 95% confidence interval [CI]
> 1.31–93.98; P = 0.023), whereas monotherapy was not (P = 0.99)."

**The confidence interval spans a factor of 72 (1.31 to 93.98) and rests on five failure events in 58
patients.** The point estimate is not a reliable magnitude. What the result supports is a *direction*
and a *ranking* — target attainment mattered, and monotherapy-versus-combination did not — not the
number 11.11. Any use of this in the manuscript should quote the CI, and should not present 11.11 as
an effect size.

The secondary claim — that optimised joint target attainment "could render unnecessary combo therapy"
— is an inference from a non-randomised comparison in which combination therapy was preferentially
given to ICU patients and to pneumonia (both P ≤ 0.023), i.e. to sicker patients with a
harder-to-treat site. Confounding by indication is severe and unadjusted for in that claim.

## Legal basis and retrieval limitation

"Copyright © 2023 American Society for Microbiology." **Not openly licensed.** PMC holds the record
but returns an **empty full text** for PMC10648963 — only the abstract is retrievable, exactly as with
Das 2019. **Everything in the CSV therefore comes from the abstract.** The multivariate table, the
per-patient data and the covariates entered into the model have **not** been seen.

**This is a real limitation of this folder**, not a formality: the CSV records the definition and the
headline OR, and nothing about how the multivariate model was specified. If this paper is to carry
weight in the manuscript beyond supplying the target definition, the full text should be obtained from
journals.asm.org (which served Das 2019 successfully when PMC would not).

## Provenance classification

**DIRECTLY REPORTED — abstract only.** Target definition, cohort composition, failure count, and the
single multivariate result. No table-level data.

## What this dataset is, and is not

**It is** the definitional origin of the joint-target construct this project models, and the first
evidence linking that construct to a clinical outcome.

**It is not** strong evidence on its own. n = 58, five events, retrospective, single centre, and the
key comparison (mono vs combo) is confounded by indication. Its role in the chain is definitional; the
outcome evidence was later extended in Gatti 2025 (n = 218), which is the paper R5 should lean on and
which has its own pre-post design limitations recorded in that folder.

## Reuse conditions to honour

1. Cite Gatti et al. 2023 (AAC) for the joint-target definition.
2. **Always quote the CI (1.31–93.98) with the OR** and never present 11.11 as an effect size.
3. Do not repeat the "renders combo therapy unnecessary" inference without noting confounding by
   indication.
4. Flag that only the abstract has been read — the multivariate specification is unverified.
5. Not openly licensed; cite and link, do not redistribute.
