# A four-round published controversy on this project's exact question

## What this folder is

Not a dataset. **A record of a complete, four-round exchange in the *Journal of Antimicrobial
Chemotherapy* (2023–2024) arguing about whether you must measure avibactam or can infer it from
ceftazidime** — which is, precisely, the question Model 2's monitoring layer answers.

The project had archived **round 1 only** (`Fresan2023_CI_TDM/`) and **round 4 only**
(`Gatti2024_ratio_one_leg/`), without knowing they were two ends of a dispute. Rounds 2 and 3 were
found on 12 August 2026 while systematically sweeping JAC with institutional access.

| Round | Citation | Position |
|---|---|---|
| 1 | Fresan et al. 2023, JAC 78:678-83, [dkac439](https://doi.org/10.1093/jac/dkac439) | measure ceftazidime only |
| 2 | Gatti & Pea 2023, JAC 78:1556-7, [dkad108](https://doi.org/10.1093/jac/dkad108) | **measure both** |
| 3 | Fresan et al. 2023, JAC 78:2385-6, [dkad217](https://doi.org/10.1093/jac/dkad217) | measure ceftazidime only |
| 4 | Gatti, Viale & Pea 2024, JAC 79:195-9, [dkad367](https://doi.org/10.1093/jac/dkad367) | **measure both** |

All four are OUP "all rights reserved"; only the extracted arguments are recorded here, no PDFs.

## Why this matters more than any single paper here

**1. It establishes that this project's central question is genuinely contested in the clinical
literature, not invented by the project.** A manuscript arguing about the value of an avibactam assay
can now cite a live four-round dispute rather than asserting the question is open.

**2. The single-analyte position rests on an assumed correlation that is never quantified.** Round 3
states it outright: *"a correlation between ceftazidime target achievement and avibactam target
achievement was assumed in our study."* Nowhere in the exchange does anyone put a number on it. **That
is exactly the parameter Model 1 estimates and Model 2 propagates** — so the project supplies the
missing quantity in a published argument, rather than answering a question nobody asked.

**3. Round 2 raises a mechanism the project does not currently model — MICi.** Gatti & Pea argue that
avibactam concentration *changes the ceftazidime MIC* ("effective MIC with an inhibitor"), citing
Tam et al. 2022 (JAC 77:3130-7) for a 76.1% fT>MICi regrowth-suppression threshold. If true, the
avibactam target is not a fixed 4 mg/L threshold but a modifier of the ceftazidime target — i.e. the
two components' targets are **coupled**, which Model 2 currently treats as independent (joint but
separable). **This is a real structural limitation of Model 2 and is not currently acknowledged
anywhere in the project.** It should be, and Tam 2022 should be read before Phase 6.

## The objection this project can answer, and neither side did

Round 3's strongest argument is **actionability**, verbatim:

> "the measurement of avibactam concentrations for optimizing ceftazidime/avibactam dosing does not
> seem useful because ceftazidime/avibactam is formulated in a fixed combination of 4:1. Consequently,
> in real clinical practice it is not possible to increase only the dose of avibactam."

**Round 4 never rebuts this.** Gatti et al. prove the ratio varies; they do not say what a clinician
should *do* about it given a fixed-ratio product. The exchange ends with the empirical question settled
and the decision question open.

**Model 2 answers the decision question, and this is probably the strongest positioning available for
the manuscript.** The reply the project can make:

- The decision variable was never "titrate avibactam" — it is **which regimen**, and every regimen
  carries the fixed 4:1 ratio. So the fixed formulation is not an argument against measuring; it is
  what makes the measurement *decision-relevant*.
- `limiting_probability()` reports which component is limiting as a **probability**, not a label. If
  avibactam is limiting, escalating the product helps; if ceftazidime is limiting, the same escalation
  is largely wasted.
- Escalating for avibactam drags ceftazidime toward the 104 mg/L neurotoxicity screen. `utility()`
  prices that trade-off explicitly through the λ exchange rate — which is exactly the calculation the
  fixed 4:1 coupling forces on you, and exactly what round 3 assumes cannot be done.
- EVPI puts a number on what resolving the avibactam uncertainty is worth, in attainment percentage
  points, so "is the assay worth it?" stops being a matter of opinion.

**In short: round 3 says the fixed ratio makes the information useless; the project's answer is that
the fixed ratio is precisely why the information has value.** That is a substantive contribution to a
live published debate, and it requires no new data.

## Honest caveats before using this framing

- **The project's answer is model-based.** It shows the information *would* change the regimen choice
  under the model's own utility function, not that it improves patient outcomes. Round 3's objection is
  practical and clinical; the reply is decision-analytic. Say so.
- **Round 3's second objection stands and is not answerable by modelling**: avibactam assays are not
  routinely available in most hospitals. No amount of value-of-information analysis creates an assay.
  R2's triage result (a small measured fraction captures most of the benefit) is the appropriate
  response to that, and it should be offered as such.
- **The MICi mechanism (round 2) is a genuine gap in Model 2**, not a point in its favour. If avibactam
  concentration modifies the ceftazidime MIC, the "joint but separable" target structure is an
  approximation. This should be stated as a limitation in the manuscript rather than discovered by a
  reviewer.

## Reuse conditions to honour

1. Cite all four rounds when describing the controversy; citing only one side misrepresents it.
2. Do not claim the project "resolves" the debate — it supplies a missing quantity (the correlation)
   and a decision framework for an objection neither side addressed.
3. Attribute the MICi concept to Gatti & Pea 2023 and the underlying data to Tam et al. 2022.
4. None of the four is openly licensed; cite and link, no PDFs archived or redistributed.
