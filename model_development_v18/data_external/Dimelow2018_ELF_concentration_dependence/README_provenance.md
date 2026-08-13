# The ELF penetration ratios are read off a non-linear curve at the wrong concentration

## Source

Dimelow R, Wright JG, MacPherson M, Newell P, Das S. **Population Pharmacokinetic Modelling of
Ceftazidime and Avibactam in the Plasma and Epithelial Lining Fluid of Healthy Volunteers.**
*Drugs R D* 2018;18(3):221-30. [doi:10.1007/s40268-018-0241-0](https://doi.org/10.1007/s40268-018-0241-0)
— PMID 30054895, PMC6131119, **CC BY-NC 4.0**. Full text retrieved via PubMed Central.

Underlying data: Nicolau DP, Siew L, Armstrong J, et al. *J Antimicrob Chemother* 2015;70:2862-9,
[doi:10.1093/jac/dkv170](https://doi.org/10.1093/jac/dkv170) (NCT01395420, 43 healthy male
volunteers, 2 h infusion q8h). Found on 12 August 2026 during the systematic JAC sweep; it is the
phase 1 study Dimelow 2018 re-analyses, and it was not previously tracked as a project source.

## What the package does

`revision_support/add_icu_elf_scenario.py` and the lung therapeutic-window analysis apply **fixed**
ELF/plasma penetration ratios:

| scenario | ceftazidime | avibactam |
|---|---|---|
| "healthy-volunteer estimate" / "central estimate" | **0.52** | **0.42** |
| "conservative" | 0.30 | 0.30 |
| ICU trial (Benítez-Cano 2026) | 0.41 | 0.44 |

The 0.52/0.42 pair is cited to Dimelow 2018 and **is quoted correctly** — those numbers are in the
abstract. The manuscript also already states, correctly, that this compartmental analysis found
greater penetration than the non-compartmental calculation.

## The problem

**Dimelow's ratios are not constants.** They are the value of a non-linear plasma-ELF relationship
evaluated at **one plasma concentration each**, and the paper says so explicitly:

- Ceftazidime: plasma-ELF link is a **saturable Michaelis-Menten** function.
  E_max = 45.4 mg/L ELF, half-maximal at a plasma concentration of 71.7 mg/L.
  **52% is the ratio at a plasma concentration of 15.3 mg/L** (where ELF reaches the 8 mg/L target).
- Avibactam: plasma-ELF link is a **power** function with exponent < 1.
  **42% is the ratio at a plasma concentration of 2.4 mg/L** (where ELF reaches the 1 mg/L target).

In both cases **the penetration ratio falls as plasma concentration rises**. That is the entire point
of the paper — it contrasts its concentration-specific ratios with the non-compartmental AUC ratios
(≈31–35%) precisely because AUC ratios "average across the studied concentration range".

The project's continuous-infusion regimens operate far above 15.3 and 2.4 mg/L. The ceftazidime
neurotoxicity screen used throughout the package sits at **104 mg/L** — roughly seven times the
concentration at which the 0.52 figure is valid.

## Quantified

`code/elf_penetration_concentration_check.py` rebuilds both published functions from the paper's own
parameters and validates them against **all ten** numeric checkpoints the paper states (ceftazidime
ELF 22.5 mg/L and ratio 32.1% at plasma 70 mg/L; ratio 63.3% as plasma → 0; ELF 22.7 mg/L at plasma
71.7 mg/L; avibactam 47% at 1 mg/L, 42% at 2.4 mg/L, ELF 4.0 mg/L and ratio 33.2% at 12 mg/L; etc.).
All ten reproduce. The avibactam exponent, which the paper does not print, is recovered as 0.8715
from two checkpoints and then confirmed against a third.

Only then does it evaluate the ratio where the project actually operates:

| plasma ceftazidime | Dimelow's own model | applied | applied ratio overstates by |
|---|---|---|---|
| 15.3 mg/L (where 0.52 is valid) | 52.2% | 0.52 | 0% |
| 40 mg/L | 40.6% | 0.52 | 28% |
| 70 mg/L | 32.0% | 0.52 | 62% |
| **104 mg/L (neurotoxicity screen)** | **25.8%** | 0.52 | **101%** |

| plasma avibactam | Dimelow's own model | applied | overstates by |
|---|---|---|---|
| 2.4 mg/L (where 0.42 is valid) | 42.0% | 0.42 | 0% |
| 10 mg/L | 35.0% | 0.42 | 20% |
| 20 mg/L | 32.0% | 0.42 | 31% |

Full grid in `elf_penetration_vs_concentration.csv`.

## Why it matters — a stated conclusion turns on this

`revision_support/outputs/lung_therapeutic_window.csv` reports whether the lung therapeutic window is
open at the clinical breakpoint:

| scenario | window closes at MIC | multiple of breakpoint | open at breakpoint? |
|---|---|---|---|
| plasma (1.0 / 1.0) | 22.1 mg/L | 2.76× | yes |
| **central estimate (0.52 / 0.42)** | 11.5 mg/L | 1.44× | **yes** |
| **conservative (0.30 / 0.30)** | 6.6 mg/L | 0.83× | **no** |

**The qualitative answer flips between the two ELF rows.** And at continuous-infusion concentrations
Dimelow's own model gives ceftazidime ≈26–32% and avibactam ≈32–35% — i.e. the row labelled
**"conservative" is the concentration-appropriate one**, and it is the row where the window is
**closed** at the breakpoint.

So the scenario presented as the central estimate is the one whose input is valid only at
concentrations the regimens do not operate at, and the scenario presented as a pessimistic
sensitivity is the one that matches the cited source at the relevant concentrations.

## What was NOT done, and why

**No file outside `model_development_v18/` was modified.** The ELF scenario code lives in
`revision_support/`, which is frozen under this project's standing constraint. This folder records a
diagnostic, not a fix. Three defensible responses exist and the choice is the author's:

1. **Relabel.** Keep all three scenarios, stop calling 0.52/0.42 the central estimate, and state that
   the conservative row is the one consistent with the cited model at continuous-infusion exposures.
   Cheapest, honest, and changes no numbers.
2. **Re-run with concentration-dependent penetration.** Replace the fixed multiplier with Dimelow's
   actual functions, which are fully specified above and validated here. Most correct; changes
   results; requires touching frozen code.
3. **Report as a limitation.** State that fixed-ratio ELF scaling is an approximation that is
   optimistic at the high steady-state concentrations continuous infusion produces.

Option 1 or 3 can be done in Phase 6 text alone. Option 2 cannot.

## Honest caveats

- This does **not** show the package misquotes its source. It quotes Dimelow correctly; the issue is
  applying a concentration-specific value as a concentration-independent constant.
- Dimelow's subjects are **healthy male volunteers on 2 h intermittent infusion**, not critically ill
  patients on continuous infusion. Extrapolating the *shape* of the curve to CI concentrations is
  itself an assumption — it is, however, the same assumption the package already makes by using the
  number at all, and applying the curve is strictly less of a stretch than applying one point on it.
- The ICU-trial scenario (Benítez-Cano 2026, 0.41/0.44) is measured under continuous infusion and so
  does not carry this defect. Note the separate, already-recorded issue that Benítez-Cano publishes
  **three** different ratio pairs and the sensitivity to that choice is not reported.
- Ratios are on **total** drug in both matrices; Dimelow notes ELF protein binding is lower than
  plasma, so total-based ratios understate free-drug penetration. That pushes the other way and is
  not quantified here because the data to quantify it do not exist in either source.

## Reuse conditions

Dimelow 2018 is CC BY-NC 4.0 — reusable with attribution for non-commercial purposes. Nicolau 2015 is
OUP "all rights reserved"; only extracted numbers are recorded, no PDF archived. Cite both when using
this: Dimelow is the model, Nicolau is the data underneath it.
