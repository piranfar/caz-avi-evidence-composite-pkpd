# MODEL 1 — joint ceftazidime/avibactam population PK model

**FINAL. Definitive record of the fitted model.** Where this document and any other document in
`model_development_v18/` disagree, this one is correct.

**Code:** `code/joint_popk_nlme.py`, `code/model1_finalise.py`, `code/test_model1.py`
**Outputs:** `outputs/model1_final_parameters.csv`, `model1_profile_likelihood.csv`,
`model1_diagnostics.csv`, `model1_vpc.csv`, `model1_sensitivity.csv`,
`model1_individual_parameters.csv`
**Figures:** `figures/model1_gof.png`, `figures/model1_vpc.png`
**Logs:** `audit/log_model1_finalise.txt`, `audit/log_test_model1.txt`
**Data:** Dryad `10.5061/dryad.fxpnvx16s` (CC0), 21 patients, 238 observations.

---

## ⚠ Correction to earlier numbers in this project

Earlier documents here reported the ceftazidime–avibactam clearance correlation as **0.560–0.598**,
from two-stage analyses (Dose/AUC, and individual nonlinear least-squares fits).

**Those estimates are biased low. The correct value is 0.703.**

Two-stage estimation computes each patient's clearance first, then correlates the estimates. Every
individual clearance carries estimation error, and correlating error-contaminated estimates
attenuates the correlation towards zero — the classic regression-dilution effect. The mixed-effects
model separates residual error (10.0% and 9.3% here) from genuine between-subject variability and
estimates the correlation of the underlying random effects directly. That is the quantity Cojutti
reports as 0.94, and therefore the only one comparable to it.

| Method | Estimate | Comparable to 0.94? |
|---|---|---|
| Dose/AUC two-stage (Dryad) | 0.560 | No — attenuated |
| Individual least-squares two-stage (Dryad) | 0.563 | No — attenuated |
| Two-stage adjusted for CRRT covariates | 0.476 | No — attenuated |
| Gatti 2023, published clearances | 0.598 | No — attenuated |
| **Mixed-effects, this model** | **0.703** | **Yes** |

The correction reduces the magnitude of the finding but not its direction. Recorded here because the
earlier numbers circulated inside this project and should not be quoted.

---

## 1. Model

One compartment per analyte, steady state under repeated 8-hourly intravenous infusion, solved in
closed form and verified against explicit superposition of 400 doses (maximum relative difference
below 10⁻⁸).

Random effects are three standard-normal deviates:

```
z = [z_CL_caz, z_CL_avi, z_V_shared]      corr(z_CL_caz, z_CL_avi) = rho
CL_caz = θ₁ exp(ω₁ z_CL_caz)     CL_avi = θ₂ exp(ω₂ z_CL_avi)
V_caz  = θ₃ exp(ω₃ z_V_shared)   V_avi  = θ₄ exp(ω₄ z_V_shared)
```

**Why the volume deviate is shared.** A first fit with four independent deviates estimated the
ceftazidime–avibactam *volume* correlation at exactly **1.000** — a boundary value, and the standard
signature of an overparameterised variance model. The reduced structure imposes that boundary as a
constraint instead of estimating it. It is also more defensible mechanistically: both agents
distribute into extracellular water, so a patient with an expanded volume has it for both. The
reduction costs nothing — the objective is identical to four decimal places (−722.3558 with four
deviates, −722.3558 with three) with one parameter fewer.

Estimation is first-order conditional with the Laplace approximation; the individual objective is a
penalised nonlinear least-squares problem solved by Levenberg–Marquardt, and curvature at the mode
uses the Gauss–Newton form H = 2(J′J/σ² + Ω⁻¹). **Convergence** alternates a Nelder–Mead pass with an
L-BFGS-B polish and repeats until successive rounds improve the objective by less than 10⁻⁴. Both
fits **stopped on that tolerance**, not on an evaluation limit. Estimation is fully deterministic;
the only stochastic step in the whole analysis is the visual predictive check, seed 20260811.

---

## 2. Results

### Fixed effects

| Parameter | Estimate |
|---|---|
| CL ceftazidime | 2.572 L/h |
| CL avibactam | 3.223 L/h |
| V ceftazidime | 19.98 L |
| V avibactam | 27.03 L |

### Between-subject variability

| Parameter | ω | CV |
|---|---|---|
| CL ceftazidime | 0.2022 | 20.4% |
| CL avibactam | 0.1411 | 14.2% |
| V ceftazidime | 0.2673 | 27.2% |
| V avibactam | 0.1994 | 20.1% |

Variability is far below the 67.9% and 76.9% of the non-RRT source cohort, as expected when a shared
extracorporeal circuit dominates elimination.

### The quantity of interest

> **corr(η_CL ceftazidime, η_CL avibactam) = 0.703**
> **95% profile-likelihood interval 0.380 to 0.874**
> **The interval excludes the assumed value of 0.94.**

The interval is from the **profile likelihood**, not from the Hessian. The two agree almost exactly
(Hessian-based: 0.381 to 0.873), which is itself reassuring about the curvature approximation.

The profile is smooth and unimodal, with its minimum at the reported estimate (ΔOFV = 1 × 10⁻⁶ at
ρ = 0.703, confirming convergence):

| ρ | 0.25 | 0.35 | 0.45 | 0.55 | 0.65 | **0.703** | 0.75 | 0.85 | 0.88 | **0.94** |
|---|---|---|---|---|---|---|---|---|---|---|
| ΔOFV | 6.41 | 4.37 | 2.58 | 1.12 | 0.17 | **0** | 0.17 | 2.46 | 4.18 | **11.78** |

**At ρ = 0.94 the objective is 11.78 units worse — the data reject that value at p = 0.0006.**

Against a model with no cross-drug correlation: **ΔOFV = 12.41 on 1 df, p = 0.00043**. The
correlation is strongly supported; it is its *magnitude* that differs from the assumption.

### Residual error and shrinkage

Proportional residual error 10.0% (ceftazidime) and 9.3% (avibactam). Shrinkage of the clearance
deviates is **−0.4% and 0.7%** — effectively zero, so the individual clearances are determined by
the data, not by the prior. Volume deviate shrinkage 8.7%. Low shrinkage is what makes the
correlation estimate credible.

---

## 3. Model evaluation

### Goodness of fit (`figures/model1_gof.png`)

| Analyte | CWRES mean | CWRES SD | \|CWRES\| > 2 | log-IPRED vs log-DV | log-PRED vs log-DV |
|---|---|---|---|---|---|
| Ceftazidime | −0.067 | 1.003 | 5.9% | r = 0.975 | r = 0.787 |
| Avibactam | +0.024 | 1.002 | 6.7% | r = 0.968 | r = 0.840 |

Conditional weighted residuals are centred on zero with unit standard deviation for both analytes —
close to textbook — and show no trend against time. Individual predictions track observations along
the identity line; population predictions scatter as they should.

### Visual predictive check (`figures/model1_vpc.png`)

1,000 replicates of the observed design. Observed medians track simulated medians closely at every
sampling time for both analytes. Coverage of the simulated 5th–95th interval is **82.7% against a
nominal 90%** (ceftazidime 85.4%, avibactam 79.9%).

**This is a mild under-prediction of variability and it is reported as such.** Two contributors: the
1-hour bin contains only 3 observations, so its coverage is unstable; and with no covariate model,
between-subject differences driven by CRRT settings are absorbed into a single lognormal term that
is slightly too narrow in the tails. It does not affect the correlation estimate, which depends on
the *joint* structure of the deviates rather than on their marginal spread.

### Sensitivity analyses

| Analysis | Variant | ρ |
|---|---|---|
| Reference | infusion duration read as hours | **0.703** |
| Infusion duration | all 1 h | 0.781 |
| Infusion duration | all 2 h | 0.744 |
| Infusion duration | all 3 h | 0.765 |
| Residual error | shared across analytes | 0.704 (ΔOFV +0.23, 1 df — no improvement) |
| Leave one subject out | range over 21 refits | **0.651 to 0.756**, median 0.700 |

**The estimate is robust.** The infusion-duration ambiguity in the source dataset — the column is
named `Infusion duration_h` but the README calls it a category — moves ρ only between 0.70 and 0.78,
and **every variant remains below 0.94.** Leave-one-subject-out refits span 0.651 to 0.756, so no
single patient drives the result; the most influential is subject 7, and removing it still gives
0.756.

### Test suite

`code/test_model1.py` — **119 checks, all passing.** Covers the brief's required list: dose
conversion, infusion rate, free versus total concentration, clearance transformation, correlated
random effects, renal-function classes, MIC weighting, PTA, CFR, limiting-component classification,
the exposure constraint and the ELF transformation — plus structural-model identities (steady state
C(0) = C(τ); AUC/τ = Dose/(CL·τ); dose linearity; peak at the end of infusion), variance-structure
guarantees (Ω symmetric, positive definite, unit diagonal, correlation round-trip), and numerical
guards (positive concentrations, sampling times inside one interval, a singular Ω rejected rather
than crashing). The primary model is imported **read-only** from the v16 package; nothing there is
modified.

---

## 4. What this changes for the manuscript

Re-running the manuscript's own second-assay classifier at the fitted correlation and its
profile-likelihood bounds:

| ρ | Specificity | NPV | Wrongly reported as attaining |
|---|---|---|---|
| **0.94** — assumed in the manuscript | 77.0% | 83.6% | **3.6%** |
| 0.873 — upper bound | 64.6% | 77.8% | 5.5% |
| **0.703 — this model** | **40.7%** | **69.3%** | **9.2%** |
| 0.381 — lower bound | 15.1% | 59.9% | 13.3% |

At the fitted value, specificity roughly **halves** and the rate at which a patient is wrongly
reported as attaining the avibactam target on the basis of ceftazidime alone rises **2.6-fold**.
**Even at the upper bound the false-reassurance rate is 1.5 times the assumed value.**

The clinical reading: inferring avibactam attainment from ceftazidime is less reliable than the
assumed correlation implies, across the whole plausible range.

---

## 5. Limitations — all of which must travel with any use of this model

**Population and administration.** CRRT patients on intermittent 8-hourly infusion. The manuscript's
primary scenario is non-RRT adults on continuous infusion. **The clearance and volume estimates do
not transfer.** Only the correlation is carried forward, and only as a sensitivity bound — never as
a replacement value for ρ in the primary analysis, which retains 0.94.

**Cohort composition.** 18 of 21 patients have acute pancreatitis. Not a general ICU population.

**Sample size.** 21 subjects for 11 parameters. Adequate for this reduced structure and shrinkage is
near zero, but it will not support covariate modelling or a richer variance structure.

**Covariates coded, not measured.** Age, weight and all CRRT operating parameters are published as
categories, so no allometric or covariate model is possible. Serum creatinine, albumin, APACHE II,
SOFA and urine output are continuous and remain available for future work.

**Sparse sampling.** Most patients have no sample before 2 hours, so the absorption-phase peak is
poorly resolved and absolute clearances should be treated as approximate. The correlation is far less
affected, since the bias acts on both analytes in the same direction.

**Variability slightly under-predicted**, as the visual predictive check shows.

**One compartment only.** A two-compartment structure was not tested. With 5–7 samples per analyte
starting at the trough, a distribution phase is unlikely to be identifiable, but this remains an
untested assumption and should be stated as one.

**This is not external validation** of the primary model, and must never be described as such.

---

## 6. Estimator validation (simulation-based calibration)

**Complete: 75/75 replicates, three scenarios, no data-integrity issues after a race-condition fix
(see below).** `code/model1_sbc.py`; results in `outputs/model1_sbc_replicates.csv` and
`model1_sbc_summary.csv`; log in `audit/log_model1_sbc.txt`.

**What this tests, and what it does not.** Each scenario simulates 25 replicate datasets from a
known correlation, reusing the observed design exactly (same 21 patients, same sampling times, same
infusion durations), then refits every replicate with the actual Model 1 estimator. Recovering the
value that was put in is a check on the **software**, not evidence about patients — it cannot and
does not validate the model against reality, and is not described as doing so anywhere in this
project.

| Scenario | True ρ | Mean estimate | Bias | Relative bias | SD | RMSE | 95% range |
|---|---|---|---|---|---|---|---|
| Correctly specified | 0.703 | 0.721 | +0.018 | +2.6% | 0.120 | 0.119 | 0.457 to 0.881 |
| Two-compartment truth (misspecified) | 0.703 | 0.724 | +0.021 | +3.0% | 0.125 | 0.124 | 0.450 to 0.890 |
| At the published value | 0.940 | 0.950 | +0.010 | +1.1% | 0.029 | 0.030 | 0.906 to 1.000 |

**The estimator is essentially unbiased at both true values checked.** A 2.6-3.0% upward bias at
ρ = 0.703 and 1.1% at ρ = 0.94 are both small relative to sampling noise, and in the same direction
at every true value tested — consistent with the ordinary small-sample upward bias of a bounded
correlation estimator rather than with a flaw specific to this implementation.

**Misspecification barely moves the answer.** Simulating from a two-compartment model — which Model 1
does not assume — and fitting the one-compartment model shifts the bias by 0.003, well inside
replicate-to-replicate noise. The one-compartment assumption, if wrong, does not meaningfully distort
the correlation estimate; the untested two-compartment alternative remains an item for future work
(§9), but this result lowers the priority of that item.

**Precision differs sharply by true value, and that difference explains the design analysis.** At
ρ = 0.94 the estimator is far tighter (SD 0.029) than at ρ = 0.703 (SD 0.120), because a correlation
near the boundary is intrinsically easier to pin down than one mid-range. This is exactly why the
design analysis (`model1_design_analysis.py`) found high power to reject 0.94 in favour of ~0.70
without needing high absolute precision on the exact value — precision near the boundary is cheap,
precision mid-range is not.

**A boundary artefact, disclosed rather than hidden.** At the true value 0.94, **2 of 25 replicates
(8%) returned an estimate of exactly 1.0000** — the correlation matrix hit the edge of its valid
range and the optimiser was stopped there. This inflates the mean estimate slightly (0.950 against
the true 0.940) and is the main contributor to the +1.1% bias at that scenario; it is a known
property of estimating a variance-covariance parameter near its boundary and is not a sign of a
coding defect. It has no analogous effect at ρ = 0.703, which is far from the boundary.

### A data-integrity issue during this run, and how it was resolved

Mid-run, two identical copies of `model1_sbc.py` were found running simultaneously against the same
checkpoint file, both writing results for the same replicate indices. This produced duplicate rows —
7 replicates in the two-compartment scenario and 5 in the ρ = 0.94 scenario were each written twice.
Both processes were stopped, the checkpoint file was deduplicated (69 rows to 57, keeping the first
occurrence of each scenario-replicate pair), and a single process completed the remaining replicates
cleanly. The final file was re-verified to contain exactly 25 rows per scenario with no duplicate
replicate indices before the results above were computed.

---

## 7. Novelty claim, stated precisely

The ceftazidime–avibactam cross-drug clearance correlation has been quantified **once** in the
published literature. No regulatory document contains it: the registrational analyses resampled
η-pairs empirically without ever estimating the covariance.

**This is the second estimate in existence, and the first from openly available individual patient
data.** The claim is narrow and checkable, and it is narrow deliberately, because the population
restriction prevents any broader one.

---

## 8. Remaining optional work

Not blocking — the model is complete and defensible as it stands, and §8 has already answered two of
the three items originally listed here.

1. ~~Estimator validation~~ — **done, §8.** Unbiased to within 3% at both true values checked; a
   two-compartment misspecification shifts nothing meaningfully.
2. **Two-compartment comparison** as an actual fitted alternative (rather than only as a
   misspecification-robustness check), now lower priority given §8's result.
3. **Covariate exploration** on serum creatinine and CRRT dose intensity, reported as exploratory
   given the sample size.
