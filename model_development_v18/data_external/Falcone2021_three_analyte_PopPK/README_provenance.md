# A second study that had everything needed to compute ρ, and did not

## Source

Falcone M, Menichetti F, Cattaneo D, et al. **Pragmatic options for dose optimization of
ceftazidime/avibactam with aztreonam in complex patients.**
*J Antimicrob Chemother* 2021;76(4):1025-31.
[doi:10.1093/jac/dkaa549](https://doi.org/10.1093/jac/dkaa549) — PMID 33378458.
OUP Standard Journals Publication Model, all rights reserved. Retrieved 12 August 2026 through
institutional access during the systematic JAC sweep. No PDF archived; extracted numbers only.

## What it is

Prospective observational PK study, **n = 41**, Pisa, July 2019 – February 2020. Elderly and complex:
median age 75 (IQR 63–79), 48.8% ICU, 19.5% burns, 36.6% chronic kidney disease, 6 on CRRT, 29.3%
septic shock, 30-day mortality 26.8%. Five samples around the fourth dose, **assayed for all three
drugs — ceftazidime, avibactam and aztreonam — in the same patients**. Population PK by Monolix
(SAEM) then NPAG in Pmetrics, with Monte Carlo simulation to a dosing nomogram.

## Why it is here — the R1 finding, second instance

Falcone 2021 estimated **individual clearance for both ceftazidime and avibactam in the same 41
patients**. It reports them in separate columns of the same table. It **never reports the correlation
between them.**

This is the second independent instance of the exact gap route R1 rests on, after Li 2019. The two are
not identical and the difference matters:

- **Li 2019** built a model that *merged both drugs' random effects per patient* — the joint structure
  existed and the correlation was one line of output away. It was not computed.
- **Falcone 2021** fitted the two drugs *separately*, so no joint structure existed. The correlation is
  still recoverable post hoc from the individual empirical Bayes estimates, which they have for both
  drugs in every patient. It was not computed.

So the claim R1 can now make is stronger and more precisely stated: **the cross-drug clearance
correlation has been computable in at least three published datasets (Li 2019, Falcone 2021,
Cojutti 2024) and reported in exactly one (Cojutti 2024, ρ = 0.94).** That is a documented gap in the
literature, not a question this project invented.

## Paired trough concentrations — the 4:1 ratio again fails

Median (range) C_min after the fifth dose, all three analytes in the same patients:

| analyte | median C_min (mg/L) | range |
|---|---|---|
| aztreonam | 53.6 | 1.5–176 |
| **ceftazidime** | **57.9** | **1.2–175** |
| **avibactam** | **10.9** | **0.4–46.3** |

Ratio of medians = **5.31 : 1**, not the 4:1 of the vial. Independent of Gatti 2024 (n=107, median
4.67:1, range 1.29–13.46:1), in a different centre and population, and pointing the same way. Two
independent cohorts now show the plasma ratio departs from the formulation ratio.

**Ceftazidime C_min reached 175 mg/L** in this cohort — well above the **104 mg/L** neurotoxicity screen
the project applies, in an elderly population at standard or reduced labelled doses. That is
observational support that the screen is set in a clinically reachable range and not a theoretical
concern.

Avibactam C_min ranged 0.4–46.3 mg/L against the project's C_T of 4 mg/L, i.e. spanning from well below
to more than ten times the threshold.

## Population PK parameters, as printed

Median (percentage interindividual variability). CL = θ1 × (CKD-EPI/50)^β + θ2, eGFR in mL/min
transformed with body surface area.

| | aztreonam | ceftazidime | avibactam |
|---|---|---|---|
| **Base** V_d (L) | 45.9 (52.9%) | 26.3 (114%) | 40.2 (98.7%) |
| **Base** CL (L/h) | 3.3 (82.6%) | 3.2 (82.0%) | 4.9 (76.4%) |
| **Final** V_d | 32.0 (61.0%) | 28.5 (67.0%) | 41.6 (88.5%) |
| **Final** θ1 | 1.8 (62.4%) | 2.5 (39.0%) | 3.5 (45.7%) |
| **Final** θ2 | 0.9 (91.2%) | NA | 0.7 (113%) |
| **Final** β | 1.2 (37.7%) | 1.2 (25.0%) | 1.3 (34.7%) |

Between-subject variability on ceftazidime V_d reaches **114%** in the base model — far above anything
the project currently carries. Not directly transferable (different structural model, 1-compartment,
elderly cohort), but a marker that the project's variability assumptions are not conservative in this
population.

## Renal-function equation — relevant to a project choice

Falcone tested Cockcroft–Gault with four weight descriptors, MDRD, and CKD-EPI, and found **CKD-EPI had
the lowest AIC and was selected automatically for all three drugs**. Weight, burns, CRRT and septic
shock were all supplanted once estimated kidney function entered. The project uses **EKFC**; this is a
third equation in the mix (Cojutti 2024 also examines kidney-function equations). No action implied —
EKFC and CKD-EPI are on a similar scale — but it is a live methodological question in this literature
and worth one sentence in the manuscript rather than silence.

## Unbound fractions — a third published pair

Protein binding 10% ceftazidime, 5% avibactam ⇒ **fu = 0.90 / 0.95**. The project uses 0.85 / 0.92
(Dimelow). Both fall inside the project's existing sensitivity ranges (0.80–0.90 and 0.87–0.97), so
this is corroboration rather than a gap. See also
[`Cojutti2026_FN_avibactam_limiting/`](../Cojutti2026_FN_avibactam_limiting/) for the third pair.

## The avibactam threshold was treated as uncertain here too

Falcone simulated avibactam PTA at **C_T = 1, 2 and 4 mg/L** rather than fixing one value — the same
recognition of threshold uncertainty that motivates the project's target-distribution scenarios
(T1–T7). A useful precedent to cite: treating C_T as uncertain is established practice, not a
device invented for this manuscript.

## Honest limits

- **n = 41, single centre, elderly, 85.4% NDM-producing Enterobacterales** treated with the
  ceftazidime/avibactam-plus-aztreonam combination. Aztreonam co-administration is a real difference
  from the project's population.
- Only 6 patients received ceftazidime/avibactam alone.
- **No individual-level data are published** — the paper reports model parameters and summary
  statistics only. Like most sources in this folder, it is class C: it can test assumptions, not
  supply data.
- The paper itself concedes "a correlation of PK/PD with clinical outcomes is not reliable with this
  sample size".

## Reuse conditions

OUP all-rights-reserved. Cite and link; no PDF redistributed. Extracted numbers in
`paired_trough_and_popPK.csv`. If the R1 claim about the uncomputed correlation is used in the
manuscript, cite Li 2019 and Falcone 2021 together and state the difference between them — the
merged-random-effects case and the separate-models case are not the same evidence.
