# MODEL 2 — specification

**Hierarchical Uncertainty-integrated Joint Attainment Model (HU-JAM)**

The project brief named this TU-JAM, for *target-uncertainty-integrated*. The scope has grown
beyond the target: the audit found that **four** assumed constants in the primary model are in fact
uncertain quantities with published uncertainty that the model never uses. The name should reflect
that, but it is the author's call.

**Status: specification. Not yet implemented.** Written while Model 1 finalisation runs.

---

## 1. What this model is for

The primary model treats every input as a known constant and then perturbs them one at a time. That
tells you how sensitive the model is to its own inputs. It does not tell you **what to do when you
do not know the inputs** — which is the situation a clinician is actually in.

Model 2 propagates the uncertainty that is already published, and converts it into decisions:
which regimen, how confident, and what it would be worth to find out more.

**It is not a pharmacokinetic model and makes no new pharmacokinetic claim.** It is a
decision-analytic layer over the existing simulation.

---

## 2. The four layers

### Layer 1 — parameter uncertainty (replaces point estimates)

Every parameter becomes a random variable, using the relative standard errors **already published**
in the source and never used by the current model:

| Parameter | Value | Published RSE | Distribution |
|---|---|---|---|
| CL₀ ceftazidime | 5.0 L/h | 6.36% | lognormal |
| Renal exponent, ceftazidime | 0.70 | 14.0% | normal |
| CL₀ avibactam | 5.9 L/h | 7.4% | lognormal |
| Renal exponent, avibactam | 0.89 | 12.3% | normal |
| ω ceftazidime (CV 67.92%) | 0.6159 | 33.3% | lognormal on ω |
| ω avibactam (CV 76.91%) | 0.6817 | 30.1% | lognormal on ω |
| **ρ, clearance correlation** | **0.94** | **23.8%** | **see Layer 3** |
| Unbound fractions | 0.85 / 0.92 | not reported | scenario bounds, as now |

This is the single most defensible upgrade in the whole plan, because **nothing is invented** — the
uncertainty is taken verbatim from the source table. The current model discards it.

### Layer 2 — between-study heterogeneity (replaces "run four models separately")

The manuscript currently runs four published population PK models one at a time and reports the
range of the results (70.2% to 94.5%). A range across four models is not an uncertainty statement.

Instead, treat the typical clearance as drawn from a distribution across published models:

```
log CL_typ(study s) = mu + tau_s ,     tau_s ~ N(0, tau^2)
```

with τ estimated from the four models' typical clearances by random-effects meta-analysis. Joint
attainment is then integrated over τ rather than reported four times. **Prediction intervals, not
ranges.**

Model sources: Cojutti 2024 (primary), Chen 2025, the registrational relationship, Bensman 2017 —
already assembled in `structural_models.csv`. **Li 2019 (Clin Transl Sci) should be added**: it is
the richest public parameter package for either drug, with full Ω covariance matrices, bootstrap
confidence intervals and executable NONMEM control streams, and the manuscript does not cite it.

### Layer 3 — the clearance correlation, informed rather than assumed

Two estimates now exist:

| Source | Population | ρ | Interval |
|---|---|---|---|
| Cojutti 2024 | non-RRT ICU, continuous infusion | 0.94 | RSE 23.8% |
| **Model 1** | CRRT, intermittent infusion | **0.703** | profile likelihood, see `MODEL1_REPORT.md` |

They describe **different populations**, so they must not be pooled into one posterior. Represent
them as a prespecified scenario set:

- **C1** ρ = 0.94 with its published uncertainty — the primary analysis, unchanged.
- **C2** ρ from Model 1 with its profile-likelihood interval — the sensitivity analysis.
- **C3** ρ uniform across the union of both intervals — the agnostic case.

Report all three. **Do not present a pooled posterior**, and do not describe C2 as the value for the
non-RRT population.

### Layer 4 — the avibactam target

The evidence assembled in Phase 2 is not a set of noisy measurements of one quantity. It comprises
**different constructs**, and pooling them would be inventing a distribution:

| Value | What it actually is |
|---|---|
| 0.15-0.5 mg/L | hollow-fibre **regrowth threshold** (Coleman 2014) — a measured exposure-response quantity |
| 1.0 mg/L | the **regulatory target**, animal-derived (Berkhout 2015/2016, murine thigh/lung, doi:10.1128/AAC.01269-15), expressed as **%fT > C_T**, not a steady-state concentration — and never tested under continuous infusion; see R3 in `NOVELTY_ROUTES.md` |
| 4.0 mg/L | the **EUCAST fixed susceptibility-testing concentration**, repurposed as a clinical target |
| 2.5 mg/L | **aztreonam**-avibactam specific — **excluded**, with the reason stated |

Three further facts constrain any distribution: the pharmacodynamic index is **time above a
threshold**, not a steady-state concentration; the threshold **tracks the partner β-lactam, not the
organism or the MIC**; and the FDA's own Phase 2 clinical exposure-response analysis was **negative**
(values clustered near 100%, formal modelling infeasible).

So: a **prespecified scenario set of explicitly labelled distributions**, every one reported, none
preferred.

| ID | Kind | Definition | Interpretation |
|---|---|---|---|
| T1 | point mass | 1 mg/L | the regulatory target |
| T2 | point mass | 4 mg/L | the EUCAST testing convention (current primary) |
| T3 | discrete, evidence-weighted | {0.5, 1, 4} with stated weights | weights are analyst-specified and must be shown |
| T4 | uniform | U(0.5, 4) | agnostic across the defensible range |
| T5 | triangular | Tri(0.25, 4; mode 1) | mass concentrated at the regulatory value |
| T6 | lognormal | median 1, 95th centile 4 | smooth, right-skewed |

**Non-negotiable:** no distribution may be described as an established clinical exposure-response
relationship, sensitivity to the choice is a primary result rather than an appendix, and the
weights in T3 are labelled analyst-specified wherever they appear.

---

## 3. Outputs — where the novelty is

Let *r* index regimens, *θ* the uncertain parameters, and *U(r, θ)* the utility: joint attainment
probability, subject to the exposure constraint.

### 3.1 Uncertainty-integrated attainment
Joint PTA and CFR integrated over θ, reported as a median with a 95% **prediction** interval rather
than a point estimate.

### 3.2 Limiting-component probability
Replace the label *"avibactam is limiting at MIC 2"* with **P(avibactam is limiting | MIC)**,
integrating over parameter and target uncertainty. A quantity, not a verdict. New for this drug
combination.

### 3.3 Probability that each regimen is optimal
P(r = argmax U). A regimen chosen on a point estimate may be optimal in only a minority of the
uncertainty space, which is exactly what a decision-maker needs to know.

### 3.4 Probability of misselection under a fixed threshold
The probability that using a single fixed avibactam target selects a **different regimen** from the
one that maximises expected utility under the full uncertainty. This directly quantifies the harm
done by the practice the manuscript criticises.

### 3.5 Expected regret
Regret(r, θ) = max_{r'} U(r', θ) − U(r, θ). Report the expected regret of the fixed-threshold
decision rule, in percentage points of joint attainment forgone.

### 3.6 Value of information ★ the strongest single contribution

- **EVPI** = E_θ[max_r U(r, θ)] − max_r E_θ[U(r, θ)] — what perfect knowledge of everything is worth.
- **EVPPI** for each of ρ, the avibactam target, and avibactam clearance — what resolving *that one
  parameter* is worth, computed by the Strong–Oakley nonparametric regression method, which needs
  only the single-loop sample already generated.

This answers two questions nobody in this literature has posed quantitatively: **is the second assay
worth buying**, and **is the exposure-response study worth running**. Value of information is
standard in health-technology assessment and essentially absent from beta-lactam target attainment.

### 3.7 A falsifiable prediction
State, before any such cohort exists, the predicted negative predictive value of ceftazidime-based
prediction of avibactam attainment in a non-RRT continuous-infusion population, with an interval. A
modelling paper that can be refuted is worth more than one that cannot.

---

## 4. Implementation

**Two-loop Monte Carlo.** Outer loop: M draws of θ from Layers 1-4. Inner loop: N virtual patients
per draw. Common random numbers across regimens within a draw, so regimen comparisons reflect the
parameter change and not simulation noise — the same device the existing code already uses correctly
for its one-at-a-time analyses.

Target sizes: M = 2,000 outer, N = 10,000 inner, subject to a demonstrated convergence check on
EVPI and on the probability of misselection, which converge more slowly than mean attainment.

**EVPPI** uses the Strong–Oakley regression estimator on the existing single-loop sample; nested
loops are not required and would be prohibitive.

**Reuses unchanged:** the population draw, the scenario engine, the CFR weighting and the exposure
screen from `cazavi_analyses.py`. Model 2 is a layer over that machinery, not a replacement.

**Reproducibility:** deterministic seeds recorded in the output; frozen machine-readable results
with a recomputed checksum manifest; unit tests extending `test_model1.py` to cover the utility
function, the regret identity (regret is zero for the optimal regimen and non-negative everywhere),
the EVPI identity (EVPI ≥ EVPPI ≥ 0), and the scenario-set definitions.

---

## 5. What Model 2 must not claim

- It is **not** a pharmacokinetic model and adds no pharmacokinetic knowledge.
- The target distributions are **not** established exposure-response relationships.
- The ρ scenario set is **not** a pooled posterior across populations.
- No output is a **bedside dosing recommendation**.
- Value-of-information results are in units of **attainment probability**, not money or lives, and
  must not be translated into either.
- **Nothing here converts the study into a clinical study**, and none of it speaks to bioRxiv
  eligibility.

---

## 6. Order of work

1. Layer 1 — parameter uncertainty over the existing engine. Self-contained, immediately checkable
   against the current point-estimate results as a degenerate case (zero uncertainty must reproduce
   the published numbers exactly — a strong regression test).
2. Layer 4 — the target scenario set. Independent of Layers 2 and 3.
3. Outputs 3.1-3.3, which follow directly from Layers 1 and 4.
4. Layer 3 — the ρ scenario set, once Model 1 is final.
5. Outputs 3.4-3.5 — misselection and regret.
6. Output 3.6 — value of information. The highest-value output and the one to protect time for.
7. Layer 2 — between-study heterogeneity. Genuinely useful but the most arguable, so last.
8. Output 3.7 — the falsifiable prediction, written only after everything above is frozen.

**Step 1 carries a decisive regression test:** with all uncertainty set to zero, Model 2 must
reproduce the frozen v16 outputs exactly. If it does not, the wrapper is wrong, not the model.
