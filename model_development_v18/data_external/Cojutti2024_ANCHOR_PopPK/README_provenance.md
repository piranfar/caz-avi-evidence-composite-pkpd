# Provenance record — Cojutti 2024: THE ANCHOR PAPER

## Why this folder matters more than the others

**This is the source of the ρ = 0.94 assumption that this entire project exists to test**, and of every
one of the six published relative standard errors that Model 2's parameter-uncertainty layer propagates.
Until now it was cited throughout the project (as "Cojutti 2024 Table 2") but **never archived or
independently verified**. This folder closes that gap.

## Source

Cojutti PG, Pai MP, Gatti M, Rinaldi M, Ambretti S, Viale P, Pea F. *An innovative population
pharmacokinetic/pharmacodynamic strategy for attaining aggressive joint PK/PD target of continuous
infusion ceftazidime/avibactam against KPC- and OXA-48-producing Enterobacterales and preventing
resistance development in critically ill patients.*
**J Antimicrob Chemother. 2024;79(11):2801-2808.** doi:10.1093/jac/dkae290 · PMID 39159014

Bologna (IRCCS Azienda Ospedaliero-Universitaria) + University of Michigan. **Same group as Gatti 2023
and Gatti 2025 elsewhere in this directory, and the addressee of the data request sent 11 August 2026.**

## Legal basis for use

Retrieved 12 August 2026 as page text from Oxford Academic. Published under OUP's Standard Journals
Publication Model — "© The Author(s) 2024 ... All rights reserved" — so **not openly licensed**. Only
the extracted numeric facts are recorded here; no PDF is archived in this folder. Cite and link; do not
redistribute the article text.

## Verification of the project's own parameters — all seven match exactly

`Cojutti2024_PopPK_parameters.csv` transcribes Table 2 in full. Cross-checking every value this project
already uses, against `outputs/model2_scenario_register.csv` and `code/model2_hujam.py`:

| Project parameter | Project's stated value | Cojutti 2024 Table 2 | Match |
|---|---|---|---|
| `cl0_caz` RSE | 6.36% | CL_CAZ 5.0 L/h, RSE 6.36% | ✓ |
| `exp_caz` RSE | 14.0% | β_CLCAZ 0.70, RSE 14.0% | ✓ |
| `cl0_avi` RSE | 7.4% | CL_AVI 5.9 L/h, RSE 7.4% | ✓ |
| `exp_avi` RSE | 12.3% | β_CLAVI 0.89, RSE 12.3% | ✓ |
| `omega_caz` RSE | 33.3% | Ω CL_CAZ 67.92, RSE 33.3% | ✓ |
| `omega_avi` RSE | 30.1% | Ω CL_AVI 76.91, RSE 30.1% | ✓ |
| `C1_cojutti` | **0.94, RSE 23.8%** | **Correlation CL_CAZ–CL_AVI 0.94, RSE 23.8%** | ✓ |

The primary model's own constants in `revision_support/reproduce_primary_run.py` were checked against
the same table and are equally faithful: `CL0_CAZ 5.0 / EXP_CAZ 0.70`, `CL0_AVI 5.9 / EXP_AVI 0.89`,
`EKFC_REF 70.0`, `CV_CAZ 0.6792 / CV_AVI 0.7691`, `RHO 0.94`, `FU_CAZ 0.85 / FU_AVI 0.92`, and
`TOX_THRESHOLD 104.0` — every one matches Table 2 or the Methods. The file's own header comment
("CL equations, IIV, rho — Cojutti et al. 2024, Table 2 — REPORTED_DIRECT") is accurate.

**Every parameter this project attributes to Cojutti 2024 is faithfully sourced.** The central
assumption — ρ = 0.94 with 23.8% RSE — is exactly as the project has always stated it. Nothing was
misquoted, rounded, or invented anywhere in the chain.

Two further project claims are also confirmed:
- **"non-RRT ICU on continuous infusion"** (the `C1_cojutti` scenario note) — the Methods state
  plainly: *"Patients with renal replacement therapies were excluded from this study."* Correct.
- **The 67.9% / 76.9% between-subject CVs** quoted in `README_Gatti2023_provenance.md` as belonging to
  "the non-RRT source cohort" are Ω CL_CAZ 67.92% and Ω CL_AVI 76.91% from this table. Correct.

## One ambiguity now partly resolved, and one caveat newly surfaced

**The ω-versus-ω² ambiguity — resolved, and the project already handles the main part correctly.**
Table 2 groups both variability terms under the heading **"CV (%) of the random effects"**, and the
Discussion confirms it in words ("The coefficients of variation between subjects ... declined to 68%
and 77%"). So the reported quantity is a **coefficient of variation in per cent** — neither ω nor ω².

**The structural conversion is already right.** `revision_support/reproduce_primary_run.py` reads the
values as CVs and converts them properly:
```python
CV_CAZ, CV_AVI = 0.6792, 0.7691        # reported interindividual variability
OMEGA_CAZ = sqrt(log(1 + CV_CAZ**2))   # 0.6159
OMEGA_AVI = sqrt(log(1 + CV_AVI**2))   # 0.6817
```
So the primary model uses ω = 0.6159 / 0.6817, correctly derived. No error here.

**What remains is a small, quantified, and conservative approximation in Model 2's uncertainty
sampler.** `model2_engine.draw_parameters` perturbs ω multiplicatively using the published RSEs
(33.3%, 30.1%) applied *directly on the ω scale* — but those RSEs are reported on the **CV** scale.
The two are not interchangeable: differentiating ω = √(ln(1+CV²)) gives an elasticity of
(∂ω/∂CV)(CV/ω) = CV²/[(1+CV²)ω²], so

| | published RSE (on CV) | elasticity | equivalent RSE on ω | project applies |
|---|---|---|---|---|
| ceftazidime | 33.3% | 0.832 | **27.7%** | 33.3% (**+20% wider**) |
| avibactam | 30.1% | 0.800 | **24.1%** | 30.1% (**+25% wider**) |

The project therefore propagates *more* variability-parameter uncertainty than the source strictly
implies, by roughly a fifth to a quarter. **The direction is conservative** — it widens the
uncertainty band rather than narrowing it — and `model2_engine.py` already labels the treatment
"SCENARIO ASSUMPTION". **Nothing has been changed here.** This note records the exact magnitude and
sign so the choice is deliberate and documented rather than incidental; whether to refine the sampler
is a modelling decision for the author, and the register's wording ("Whether the reported RSE applies
to omega or to omega-squared is not stated in the source") can now be replaced with the sharper fact
that it is reported on the CV.

**A caveat that should be carried into Model 1's comparison, newly surfaced by reading the source:**
volumes of distribution in this model were **fixed, not estimated** — V1 = 18.0 L and V2 = 18.1 L,
taken from Li 2019. The authors state why: *"As with any study of CI antibiotics, we could not estimate
V and so fixed values to previous population PK estimates and modelled variability as a function of CL
alone."* This matters because **Model 1 estimates volume freely** (and found a shared-volume structure
necessary). The two models are therefore not structurally comparable on volume, and ρ = 0.94 here is a
correlation between clearances in a model where all between-subject variability was forced into
clearance. That is a legitimate difference in modelling choice, not an error — but a Model 1 vs.
Cojutti comparison that ignores it would be misleading.

## What this dataset is, and is not

**It is** the anchor: n = 112 critically ill patients, 185 paired steady-state concentrations with
**both compounds measured** — the property most of this directory's other datasets lack, and the one
the clearance-correlation question requires. Continuous infusion. RRT excluded. Wide renal function
(eGFR 8–215 mL/min). It also supplies, at source, the 104 mg/L neurotoxicity threshold Model 2 uses as
its exposure ceiling, and documents that 21/112 patients exceeded it on 24 occasions.

**It is not** patient-level data — Table 2 is a model summary, not a per-patient listing. The
correlation ρ = 0.94 is a single fitted population parameter with a 23.8% RSE, from one centre, in one
retrospective cohort, with V fixed. That is precisely why this project fitted Model 1 independently
rather than adopting 0.94 on faith, and precisely why the two estimates (0.94 here, 0.703 in Model 1's
CRRT cohort) are reported as **separate scenarios and never pooled** — different populations, different
infusion modes, different structural assumptions.

**A further limitation the authors state, relevant to every free-concentration claim in this project:**
*"only total plasma concentrations were determined for both compounds"* — the free fractions are
calculated by multiplying by 0.85 (CAZ) and 0.92 (AVI), not measured. Note this differs from the
factors used by the same group's Gatti 2025 (0.90 and 0.93) and by Lanini 2024 (0.90 and 0.93) — the
literature is not internally consistent about avibactam protein binding, and the project should not
present any single pair of factors as settled.

## Reuse conditions to honour

1. Cite Cojutti et al. 2024 wherever these values are used — which, given ρ = 0.94, is most of Model 2.
2. State that V was fixed, not estimated, in any comparison against Model 1.
3. State that free fractions are calculated, not measured, and that the multipliers vary across sources.
4. Not openly licensed — no PDF is archived here and none should be redistributed.
