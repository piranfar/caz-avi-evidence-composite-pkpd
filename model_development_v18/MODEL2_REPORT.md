# MODEL 2 — HU-JAM results

**Hierarchical Uncertainty-integrated Joint Attainment Model.** Decision-analytic layer over the
primary simulation. **Not a pharmacokinetic model**; it adds no pharmacokinetic knowledge.

**Code:** `code/model2_engine.py`, `model2_hujam.py`, `model2_monitoring.py`,
`model2_heterogeneity.py`, `model2_layer2_compare.py`, `model2_breaking_point.py`, `model2_triage.py`
**Outputs:** `outputs/model2_*.csv` · **Logs:** `audit/log_model2_*.txt`
**Specification:** `MODEL2_SPECIFICATION.md`
**Run:** 2,000 outer draws × 20,000 virtual subjects, seed 20260811; monitoring analysis
1,000 draws.

---

> ## ⚠ CORRECTION — the C1 correlation sampler was biased; affected results were rerun
>
> Writing `code/test_model2.py` during the R6 packaging effort (see `SOFTWARE.md`) found that the
> **C1_cojutti** correlation scenario — the manuscript's own assumed value, 0.94 — sampled
> `Normal(0.94, 0.238×0.94)` directly on the correlation scale and clipped to a valid range. This
> clipped **39.5% of draws** to the 0.999 boundary, which pulled the sampled **mean** down to 0.877
> while simultaneously creating an artificial cluster of draws exactly at the boundary — a
> distribution that resembled neither the true value nor a clean approximation of it. **Fixed** by
> sampling on the Fisher-z scale instead, the same approach the C2 scenario already used
> (`model2_hujam.sample_rho`), which removes the clipping entirely and recovers a median of 0.940.
>
> **Everything that used this sampler was rerun.** Two patterns emerged, and both matter:
>
> - **§3.1–3.3 and §3.6 (Decision A: which regimen) barely moved** — EVPPI(target) 0.638 → 0.630 pp,
>   limiting-component probability at MIC 2 mg/L 52–68% → 51.7–67.5%, Layer 2 prediction interval
>   25.47 → 25.40 pp. This is a *confirmation*, not a correction: it is exactly what §1's finding
>   ("regimen choice is almost insensitive to ρ") predicts, since regimen choice is driven by the
>   avibactam target, not the correlation, regardless of how the correlation is sampled.
> - **§2 (Decision B: measure or infer) moved substantially**, because this is the one place ρ
>   genuinely drives the answer. At the manuscript's 4 mg/L target with no assay error, the gain
>   from measuring rose from **4.3 to 5.6 percentage points**; with 20% assay error and a uniform
>   target, it rose from **0.44 to 0.61 pp**. The largest shift is at the uniform-target, no-assay-
>   error combination: **0.67 → 1.88 pp, nearly triple.** Every number in §2 below is the corrected
>   one. **The qualitative conclusion is unchanged and was re-verified under the fix: EVPPI(ρ) is
>   still exactly 0.0000 pp in every scenario** — measuring still beats inferring across the entire
>   plausible correlation range, and that conclusion still does not depend on resolving ρ.
>
> Pre-correction files are kept as `outputs/*.csv.pre_fix_backup` for anyone who wants to audit the
> change directly. Full account: `SOFTWARE.md`, "A correctness fix this packaging effort found."

---

## 0. The engine is verified against v16

With every uncertainty set to zero, ρ at 0.94 and the avibactam target at 4 mg/L, the engine
reproduces the frozen v16 `primary_pta_results.csv` **exactly**:

```
rows compared 121
max |delta| caz_pta      0.000e+00 pp
max |delta| avi_pta      0.000e+00 pp
max |delta| joint_pta    0.000e+00 pp
max |delta| exceedance   0.000e+00 pp
```

Everything below is therefore the same simulation with uncertainty switched on, not a different one.

---

## 0.1 What this report is asking — and why that is the point

*Added 12 August 2026. The results below are unchanged; what follows is the frame they belong in,
which `NOVELTY_STRATEGY.md` Part I §4.1 identified as "the strongest single candidate for the novelty
peak" and which the work has been carrying implicitly ever since.*

This project cannot answer the question a target-attainment study normally asks. It has no
individual patient data, one anchor cohort, and a correlation reported by exactly one group. Any
claim of the form *"regimen X is best"* would rest on parameters it cannot verify.

**So it asks a different question, and the shortage of data is the reason the question is worth
asking rather than an apology for it:**

> We do not claim to know which regimen is better. We quantify **what it would be worth to find
> out** — which uncertainty, if resolved, would actually change what a clinician or a formulary
> does — and we report that in percentage points of patients correctly classified, not in currency.

This is standard practice in health-technology assessment (Claxton and Sculpher's value-of-information
framework, which now informs research prioritisation at NICE) and **essentially absent from the
beta-lactam target-attainment literature**, which reports how attainment moves when an input moves and
stops there. Sensitivity analysis says *the answer is fragile*. Value of information says *here is
what to go and measure, and here is what not to bother measuring.* The second is a decision; the first
is an observation.

### The two results that frame everything below

The analysis produces one nonzero answer and one zero, and **both are findings**:

| decision | question | EVPPI | what to do about it |
|---|---|---|---|
| **A** — which regimen | Is it worth resolving the **avibactam target**? | **0.630 pp**, the largest of any input and 29× that of ρ *on this same decision* | **Yes.** Among everything unknown, this is what to spend an experiment on. (§3.1) |
| **B** — measure or infer | Is it worth resolving the **clearance correlation ρ** first? | **0.0000 pp**, in every scenario | **No.** Measure avibactam regardless; the decision does not turn on ρ. (§2) |

**The two rows are two different decisions and the ρ figures in them are not the same number.** ρ has
EVPPI 0.021 pp on decision A — small but nonzero, which is what the "29×" compares against — and
0.0000 pp on decision B. Quoting either in place of the other would misstate the result.

The zero is the more useful of the two and the easier to misread. **It is not a failed analysis.** It
says the monitoring decision is *robust to a parameter the field is currently arguing about* — see the
four-round JAC exchange in `data_external/JAC_exchange_measure_one_or_both/`, where both sides dispute
a correlation neither ever quantifies. A study can be worth doing precisely because it shows a
contested quantity does not need to be settled.

Read that way, the two results are a **research-prioritisation statement**: resolve the target, not
the correlation. That is a claim the current literature does not make, and it costs no data anyone
lacks.

### What this frame does not license

- **It does not make the model right.** VOI inherits every assumption underneath it. A well-framed
  wrong model is still wrong, and §3.7's breaking-point analysis is the honest counterweight: regimen
  choice is fragile to the *choice of population PK model*, and no amount of framing fixes that.
- **The units are attainment percentage points, not lives or money**, and must not be converted into
  either (§6).
- **"Worth finding out" is worth-to-this-decision**, under this model's utility function. A different
  utility, or a different decision, gives a different ranking.
- EVPI is a statement about *expected* value under current uncertainty. It is not a promise that an
  experiment resolving the target would change practice.

---

## 1. Two decision problems, and only one of them depends on ρ

The first run exposed something worth stating plainly: **regimen choice is almost insensitive to the
clearance correlation.** Misselection and value of information were nearly identical at ρ = 0.94 and
ρ = 0.703. That is not a bug. ρ shifts the *level* of joint attainment but moves every regimen within
a renal class in the same direction, so it barely changes the *ranking*.

Where ρ is decisive is the **monitoring** decision. So the model addresses two problems:

| | Decision | Driven by |
|---|---|---|
| **A** | Which regimen for this renal class? | the avibactam target |
| **B** | Measure avibactam, or infer it from ceftazidime? | the clearance correlation |

Problem B is the clinically actionable one, and it is where Model 1 feeds in.

---

## 2. Decision problem B — is a separate avibactam assay worth running?

Two strategies. **INFER**: measure ceftazidime only, then predict avibactam attainment from the
conditional distribution of the avibactam random effect. **MEASURE**: assay avibactam directly.
The output is the **accuracy gain** from measuring, in percentage points of correctly classified
patients — a break-even quantity, so no monetary cost has to be invented.

At the manuscript's 4 mg/L target:

| ρ scenario | Assay CV | Infer | Measure | **Gain** | 95% interval | P(gain > 0) |
|---|---|---|---|---|---|---|
| 0.94 (published) | 0% | 94.4% | 100% | **5.6 pp** | 1.6 to 13.4 | 100% |
| 0.94 (published) | 20% | 91.9% | 95.2% | 3.2 pp | −1.2 to 12.1 | 89.1% |
| **0.703 (Model 1)** | 0% | 88.3% | 100% | **11.7 pp** | 3.7 to 21.8 | 100% |
| **0.703 (Model 1)** | 20% | 87.4% | 95.3% | **7.7 pp** | 0.7 to 19.1 | 98.6% |
| Agnostic U(0.38, 0.98) | 0% | 88.7% | 100% | 11.3 pp | 3.6 to 22.4 | 100% |

At the uniform-target scenario, the correction is larger still:

| ρ scenario | Assay CV | Infer | Measure | **Gain** |
|---|---|---|---|---|
| 0.94 (published) | 0% | 98.1% | 100% | **1.9 pp** (was 0.67) |
| 0.94 (published) | 20% | 97.5% | 98.3% | 0.6 pp (was 0.44) |

**At the correlation Model 1 estimates, measuring avibactam is still worth roughly twice as much as
at the assumed correlation** — 11.7 against 5.6 percentage points of correct classification at the
manuscript's fixed target, or 7.7 against 3.2 once realistic assay error is included. The ratio is
smaller than first reported (previously appeared as a factor of ~3, using the biased sampler's
understated denominator), but the direction and the substance are unchanged: measuring is worth
several-fold more if the true correlation is closer to Model 1's estimate than to the assumed value.

### The decision-theoretic result

> **Expected value of partial perfect information on ρ = 0.0000 pp in every scenario.**

That is not a null result, it is the answer: **the monitoring decision does not depend on resolving
ρ.** Measuring beats inferring across the entire plausible range, from 0.38 to 0.98. Resolving the
correlation would change *how much* is gained, not *what to do*.

The practical statement: **avibactam should be measured rather than inferred, and that conclusion
does not require agreement about the correlation.**

---

## 3. Decision problem A — which regimen, when the target is uncertain

### 3.1 The value of information ranking — the strongest single result

Partial EVPI for each uncertain input, averaged over renal classes (correlation scenario C1, target
scenario T4):

| Parameter | EVPPI | Relative to ρ |
|---|---|---|
| **Avibactam target** | **0.630 pp** | **29×** |
| Ceftazidime clearance | 0.088 pp | 4.1× |
| Avibactam variability (ω) | 0.054 pp | 2.5× |
| Avibactam clearance | 0.044 pp | 2.1× |
| Avibactam unbound fraction | 0.038 pp | 1.8× |
| Clearance correlation ρ | 0.021 pp | 1× |

(Recomputed under the corrected C1 sampler — see the correction notice at the top of this document.
Change from the first computation is under 6% in every row; the ranking and every conclusion below
are unchanged.)

**Resolving the avibactam target is worth about seven times more than resolving ceftazidime
clearance and about thirty times more than resolving the correlation.**

This is the manuscript's central claim restated in decision-theoretic terms, and it is a stronger
form of it. The current paper shows that *changing* the target changes the answer — which a reviewer
may fairly call arithmetic. This shows that, among everything not currently known, the target is the
thing most worth finding out. It converts a sensitivity observation into a research-prioritisation
statement.

### 3.2 Misselection and regret

The fixed-threshold rule uses point estimates at a 4 mg/L target — the manuscript's current practice.

| Target scenario | Misselection (5 renal classes) | Expected regret |
|---|---|---|
| T2, fixed at 4 mg/L | 0.0 to 14.2% | 0.00 to 0.19 pp |
| T4, uniform 0.5-4 mg/L | 0.0 to 75.2% | 0.00 to 2.91 pp |
| **T7, Coleman 2014 evidence-derived** | **0.0 to 98.6%** | **0.00 to 5.11 pp** — see §3.8 |

(Recomputed under the corrected C1 sampler; changes from the first computation are under 0.5
percentage points throughout — this section barely moved, consistent with §1.)

**One renal class (121-150, regimen R13) shows 0% misselection and 0 regret in every scenario** —
the highest dose step is optimal regardless of the target, so the fixed rule can never get it wrong
there. Excluding that trivial class, misselection ranges 1.8 to 14.2% at the fixed target and 51.3
to 75.2% under target uncertainty.

Under target uncertainty the fixed rule picks a different regimen from the one that maximises
expected net benefit in **roughly half to three-quarters of draws** across the four non-trivial
classes, and the modal optimum shifts *down* one dose step in most of them (a lower avibactam target
needs less drug, so the exposure penalty favours the smaller regimen).

The regret is nonetheless modest — under 3 percentage points of net benefit even in the worst class —
because the competing regimens are close. **Reporting both matters:** misselection alone would
overstate the problem, and regret alone would hide it.

### 3.3 The limiting component is a probability, not a label

The manuscript states categorically that avibactam is limiting at MIC 2 mg/L and ceftazidime at
MIC 8 mg/L. Integrating over target and parameter uncertainty:

| MIC (mg/L) | P(avibactam limiting) | Verdict |
|---|---|---|
| 1 | 80.5 to 84.9% | uncertain |
| **2** | **51.7 to 67.5%** | **uncertain — close to a coin flip** |
| 4 | 3.5 to 31.5% | uncertain |
| **8** | **0.0 to 0.6%** | **ceftazidime — robust** |

(Recomputed under the corrected C1 sampler; each bound moved by under 1 percentage point.)

**The MIC 8 claim survives; the MIC 2 claim does not.** "Avibactam was the limiting component in all
11 regimens" at MIC 2 holds at a fixed 4 mg/L target but becomes barely better than even once the
target is allowed to be uncertain. This is a finding the manuscript should report, and it is exactly
the kind of thing a fixed-threshold analysis cannot see.

### 3.4 The exchange rate the manuscript never states

Utility is net benefit: joint CFR minus λ × exposure-screen exceedance. λ is the rate at which a
percentage point of exceedance is traded against a percentage point of attainment, and it is
**analyst-specified**.

The manuscript makes the same trade implicitly, through a 15% ceiling that its own source describes
as arbitrary. Making it explicit shows how much rides on it:

| λ | Classes where the explicit rule and the 15% ceiling pick the same regimen |
|---|---|
| 0 | 4 of 5 |
| 0.5 | 4 of 5 |
| 1.0 | 4 of 5 |
| 2.0 | 3 of 5 |
| 4.0 | 2 of 5 |

**No value of λ reproduces the manuscript's selected set exactly.** At λ ≤ 1 the lowest renal class
selects R2, which the 15% rule excludes at 16.2% exceedance; at λ ≥ 2 three classes drop a dose step.
The selected regimens therefore depend materially on an exchange rate that is nowhere stated.

**This is not a correction of the manuscript.** It is a different decision rule, and it is reported
as such. The finding is that the choice of rule matters and should be made explicit.

### 3.5 An earlier utility was discarded, and why

The first implementation used the 15% ceiling as a hard feasibility constraint — utility zero above
it. Expected regret then came out at 20 to 33 percentage points with a 95th centile near 95, because
regret was measuring draws that crossed a cliff where 14.9% exceedance scores 85 and 15.1% scores 0.
That measures the cliff, not the decision. The linear penalty replaced it. The constrained variant is
retained in the code and reported as a secondary analysis.

---

## 3.6 Layer 2 — between-study heterogeneity

The manuscript runs four published population PK models one at a time and reports the range of the
results: population joint CFR **70.2% to 94.5%**, a spread of 24 points. A range across four models
is a description, not an uncertainty statement — it has no coverage interpretation, and it widens
simply by adding another model.

Layer 2 instead treats the typical clearance as drawn from a distribution across studies,
log CL(study, class, analyte) = μ(class, analyte) + u_study with u ~ N(0, τ²), and integrates over u.
τ is estimated by the ANOVA method of moments from a balanced two-way layout: four studies × three
renal classes (those all four cover) × two analytes.

**τ = 0.424 on the log scale — a between-study coefficient of variation of 44.4%.** The
study-by-cell residual is much smaller (SD 0.139, CV 13.9%), which confirms that heterogeneity is a
property of the *study*, shared across analytes, rather than something that varies cell by cell.

| Study | Mean log deviation from the grand mean | Multiplier |
|---|---|---|
| Bensman 2017 (cystic fibrosis) | +0.394 | ×1.48 |
| Cojutti 2024 (the primary model) | +0.126 | ×1.13 |
| Registrational relationship | +0.089 | ×1.09 |
| **Chen 2025** | **−0.609** | **×0.54** |

### What including it does

| Scenario | τ | Median prediction interval | Mean EVPI |
|---|---|---|---|
| **Off — primary analysis** | 0 | **25.4 pp** | 0.96 pp |
| **Estimated τ** | 0.424 | **48.0 pp** | 2.95 pp |
| Leave-one-out low | 0.158 | 29.3 pp | 1.39 pp |
| Leave-one-out high | 0.515 | 56.6 pp | 3.49 pp |

(Recomputed under the corrected C1 sampler; every figure moved by under 0.5 pp — this layer uses C1
as its base rho scenario but, like the rest of Decision A, is driven by the target and by tau, not
by rho.)

**Including between-study heterogeneity roughly doubles the prediction interval, from about 25 to
about 48 percentage points, and triples the value of perfect information.** The manuscript's
four-model range materially understates how uncertain the attainment estimate is once the choice of
source model is itself treated as unknown.

### Why this is a scenario and not a result

**τ is driven almost entirely by one study.** Leave-one-out refits give τ from 0.158 to 0.515 —
dropping Chen 2025 alone collapses it by two-thirds, because Chen is the single low outlier at ×0.54.
The honest statement is therefore a range: **the prediction interval is 30 to 57 percentage points
wide depending on which studies one is willing to treat as exchangeable.**

Three further constraints, all of which must travel with any number from this layer:

1. **Four studies.** Standard meta-analysis guidance is that a between-study variance estimated from
   fewer than about five studies is unreliable. τ here has roughly three degrees of freedom.
2. **The studies are not exchangeable.** Critically ill adults on continuous infusion (Italy), adults
   with carbapenem-resistant *K. pneumoniae* (China), pooled registrational phase 1–3 subjects, and
   adults with cystic fibrosis and preserved renal function. The cystic fibrosis cohort is arguably
   not from the same population at all.
3. **Different renal descriptors.** EKFC, Cockcroft–Gault and creatinine clearance are compared as if
   interchangeable at the same class boundary. They are not exactly interchangeable.

Every result is therefore reported with and without the layer. **The primary analysis remains
Layer 2 off**; the layer is a sensitivity analysis that quantifies how much the choice of source
model matters, which is a question the manuscript currently answers with a range.

---

## 3.7 How wrong would each input have to be? — the breaking-point analysis

A sensitivity analysis asks how much the output moves when an input moves, and leaves the reader to
judge whether that matters. This asks the inverse: **for each input, what is the smallest error that
would change the decision, and is an error that large plausible?**

Each input is moved on its own across a range from 0.30× to 4× its point estimate, and the first
value that changes the recommended regimen is recorded. Two comparators make the answer
interpretable: the required error expressed in **published standard errors**, and — for the
clearances — expressed against the **between-study standard deviation** of 0.424 from Layer 2.

### Decision A — which regimen

| Input | Smallest error that changes any choice | In published SEs | Against between-study SD |
|---|---|---|---|
| **Clearance correlation ρ** | **none in the whole range** | — | — |
| **Ceftazidime unbound fraction** | **none in the whole range** | — | — |
| **Ceftazidime variability ω** | **none in the whole range** | — | — |
| Avibactam variability ω | −45% | 2.0 | — |
| **Ceftazidime clearance** | **−22%** | **3.9** | **0.59** |
| **Avibactam clearance** | **−22%** | 3.4 | **0.59** |
| **Avibactam target** | **−22%** (4 → 3.1 mg/L) | — | — |
| Exposure screen threshold | −22% | — | — |
| Renal exponents | +54% to +81% | 3.5 to 4.2 | — |

**The important line is the last column.** A 22% error in ceftazidime clearance is **3.9 published
standard errors** — implausible as estimation error. But it is only **0.59 of the between-study
standard deviation**, meaning **the difference between two published population PK models is more
than enough to change which regimen is recommended.**

That is the sharpest statement this project can make about its own fragility, and it is honest in
both directions: the regimen choice is robust to *parameter uncertainty within a model* and fragile
to *the choice of model*. It is also an independent argument for Layer 2, arrived at from a
completely different direction.

Three inputs never change the decision at all across a range from 0.30× to 4×: the clearance
correlation, the ceftazidime unbound fraction, and the ceftazidime variability. For ρ this restates
§1 in a stronger form — regimen choice does not depend on it under *any* value it could take.

### Decision B — measure avibactam, or infer it

| Assay CV | Correlation at which inferring catches up | Gain from measuring |
|---|---|---|
| 0% | **never** (scanned 0.30 to 0.99) | 15.0 pp at ρ = 0.30 → **3.5 pp at ρ = 0.99** |
| 10% | **never** | 12.5 → 1.8 pp |
| 20% | **never** | 10.1 → 1.4 pp |
| 30% | **never** | 7.7 → 1.4 pp |

**There is no breaking point.** Measuring avibactam beats inferring it at every correlation from 0.30
to 0.99 and at every assay imprecision from 0% to 30%. **Even at a correlation of 0.99 — higher than
anyone has claimed — measuring still gains 3.5 percentage points of correct classification**, because
inference propagates through a model whose other parameters remain uncertain.

This is the strongest robustness statement available for the paper's most actionable recommendation,
and it is much stronger than the value-of-information result in §2. That said the decision does not
depend on resolving ρ; this says the decision does not depend on ρ's value *at all*.

### Decision B, restated — the decision boundary of a published dispute

The table above is a sensitivity analysis. It is also, once the JAC exchange is on the table
(`data_external/JAC_exchange_measure_one_or_both/`), something considerably more useful: **the
decision boundary of a live published disagreement.**

Fresan et al. monitor ceftazidime alone, and their stated justification is a correlation they say
they *assumed*: "a correlation between ceftazidime target achievement and avibactam target achievement
was assumed in our study." Neither side ever puts a number on it. So the question the exchange leaves
open is precisely: **what would that correlation have to be for the single-analyte position to be
right?**

**One thing has to be fixed before that question can be answered, and it is not cosmetic.** Fresan's
assumed correlation is between **target achievements** — two binary outcomes. The ρ scanned above is
between **clearances**. These are linked but they are not the same number, and treating them as
interchangeable would be exactly the slippage this project exists to avoid.
`code/model2_dispute_boundary.py` computes the mapping.

| clearance ρ | induced attainment φ | ceiling on φ | % of ceiling | inference catches up? |
|---|---|---|---|---|
| 0.30 | 0.217 | 0.522 | 41% | no |
| 0.703 | 0.407 | 0.520 | 80% | no |
| 0.90 | 0.500 | 0.519 | 97% | no |
| 0.94 | 0.508 | 0.517 | 99% | no |
| **0.99** | **0.519** | **0.520** | **100%** | **no** |

**The result is sharper than the ρ scan alone.** A clearance correlation of 0.99 — higher than anyone
in this literature has claimed — induces an attainment correlation of only **0.52**. The mapping is
strongly compressive, and it does not compress because of anything in this model. It compresses
because **two binary outcomes with different prevalences cannot be strongly correlated**: at the
breakpoint, ceftazidime attains in ~58% of patients and avibactam in ~83%, and the Fréchet–Hoeffding
limit for those margins is **φ ≤ 0.52**. That bound holds for *any* joint distribution, under any
pharmacokinetics, however tightly the clearances co-vary. The simulation runs into it and never
crosses it (verified at every ρ on the grid).

So the single-analyte position does not fail because the correlation happens to be low in this
dataset. **It fails because the correlation it assumes is bounded above by the mismatch in attainment
rates between the two components** — and at that ceiling, inference still misclassifies.

That is the answer to the question the exchange left open, it is derived rather than measured, and it
requires no data neither side had.

**What this does not settle, and must not be claimed to:**

- **Fresan's second objection is untouched.** They also argue that an avibactam assay is not routinely
  available. No bound answers that; the triage rule in §3.9 does, and it should be offered as the
  response to that leg rather than this one.
- The `measure` column in the script's output is 100% only because that table uses a **noiseless
  assay**, deliberately, to isolate the ceiling on inference. The gain figures to quote are the ones
  in the Decision B table above, which carry assay CV up to 30%.
- "Accuracy" is correct classification under this model, not a clinical outcome.
- φ depends on the MIC at which it is evaluated. The breakpoint is a defensible choice because it is
  where the empirical-therapy decision is taken, not the only one.
- The prevalences that set the ceiling are this model's, so the *numerical* bound is model-dependent
  even though the *inequality* is not. What is general is that unequal attainment rates cap the
  correlation; what is specific is that they cap it at 0.50 here.

---

## 3.8 An evidence-derived target scenario (Coleman 2014 hollow-fibre threshold)

T1-T6 (§3.1-3.3) are, by design, either fixed regulatory/testing conventions or analyst-specified
distributions expressing a *range* of plausible targets without evidence weighting one point in that
range over another. **T7 is different: the first target scenario in Model 2 built directly from
measured data rather than specified by the analyst.**

**Correcting an earlier internal error.** `NOVELTY_ROUTES.md`'s original description of this route
said Coleman et al. 2014 reported "strain-specific regrowth thresholds (approximately 0.15 to 0.28
mg/L across eight strains)." That is wrong on the strain count, caught only once the primary source
was actually read for this implementation (rather than relied on from an earlier paraphrase). Coleman
et al. 2014 (*Antimicrob Agents Chemother* 58(6):3366-72, doi:
[10.1128/AAC.00080-14](https://doi.org/10.1128/AAC.00080-14), PMC4068505) used **eight strains** for
an unrelated single-dose killing experiment elsewhere in the same paper. The regrowth-threshold
estimate — the avibactam concentration below which bacterial regrowth resumed during continuous-
infusion ceftazidime with a single avibactam bolus — was estimated in a separate experiment using
only **three strains**, reported in the paper's Table 2:

| Strain | Enzyme | Critical concentration (CT) |
|---|---|---|
| *E. cloacae* 293HT96 | serine β-lactamase | ≤0.15 mg/L (at 12 h) |
| *K. pneumoniae* 283CF5 | serine β-lactamase | ≤0.22 mg/L |
| *K. pneumoniae* Tunisie K4 | serine β-lactamase | ≤0.28 mg/L |
| *E. cloacae* 293HT96 (re-estimate) | same strain, 18-20 h | ~0.2 mg/L |

**T7 uses the first three rows only** — the fourth is a second measurement on the same strain
(293HT96) at a later timepoint, not an independent fourth strain; including it would double-count
that strain's contribution. It is retained here only as an informal consistency check: 0.15 mg/L at
12 h versus ~0.2 mg/L at 18-20 h for the same isolate is the same order of magnitude and the same
direction (rising slightly as the exponential-decline extrapolation is pushed later), which is what a
noisy but non-contradictory re-estimate should look like.

**T7 = {0.15, 0.22, 0.28} mg/L, equal weight.** This is deliberately *not* presented as a
random-effects synthesis with an estimated between-strain variance — three points, one of them not
fully independent of a fourth, cannot support one. It is the plain empirical distribution over the
three strains actually measured: a new strain from the same source population is modelled as equally
likely to resemble any of the three, no smoother and no more precise than the data justify. All three
values are themselves upper bounds (extrapolated from exponential-decline curves at the last
pre-regrowth timepoint sampled, not exact measurements), so treating them as point values is itself a
conservative simplification of an already small dataset. **This remains a distribution over an
in-vitro regrowth threshold in a hollow-fibre model, not over the clinical target** — the gap between
the two is exactly the problem T1-T6 already document, and T7 does not close it; it only anchors one
scenario to the specific published numbers instead of an analyst-chosen approximation of them (the
existing T3 discrete scenario had approximated this construct with a placeholder value of 0.5 mg/L,
noticeably higher than any of the three actual measurements).

**Result.** T7's target (0.15-0.28 mg/L) sits far below the 4 mg/L the fixed-threshold decision rule
is calibrated against (§3.2). Consequently, across the four renal classes where the decision is not
already trivial, the fixed rule selects a **more aggressive regimen than optimal in 93.0 to 98.6% of
draws** — because so little avibactam exposure is needed to clear a target this low that the extra
dose only adds exposure-screen risk without adding attainment. The expected cost of that
over-selection is nonetheless modest: **3.66 to 5.11 percentage points of expected regret**, because
the fixed choice and the T7-optimal choice are adjacent dose steps, not far apart ones. EVPI is
correspondingly low (0.004 to 0.735 pp) for the same reason regret is low: resolving the target
would rarely change the choice by more than one step. The top renal class (121-150, regimen R13)
is unaffected, as in every other target scenario — the highest dose step remains optimal regardless
of target.

**Read together with §3.1-3.3, not instead of them.** T7 does not settle which target is correct; it
adds one more internally consistent scenario, now traceable to specific measured numbers rather than
an analyst's placeholder, to a picture that already shows the choice of avibactam target dominates
every other source of uncertainty in this decision.

---

## 3.9 The triage rule against a clinically-derived comparator (R2 vs Gatti 2024)

R2's triage rule is model-derived: measure avibactam when the conditional probability of attainment
given the observed ceftazidime concentration is near one half. **A completely independent, clinically
derived rule for the same decision now exists**, and comparing them is the first external check any
part of Model 2 has had.

Gatti, Viale & Pea 2024 (*J Antimicrob Chemother* 79:195-9,
doi:[10.1093/jac/dkad367](https://doi.org/10.1093/jac/dkad367); archived at
`data_external/Gatti2024_ratio_one_leg/`) measured both components in 107 patients and found the
ceftazidime-to-avibactam ratio ranges 1.29:1 to 13.46:1 against the 4:1 vial ratio. Their ROC analysis
gives a rule based on renal function alone: **CrCL > 75 mL/min/1.73 m²** identifies ratios > 5:1
(AUC 0.694), **> 78** identifies ratios > 6:1.

`code/model2_triage_vs_gatti.py` puts both rules on the same population, the same model, the same 20%
assay imprecision, and the same budget axis. Three results.

### 3.9.1 R2's rule dominates at every assay budget

| % of patients measured | R2 rule | Gatti rule | difference |
|---|---|---|---|
| 5% | 93.32% | 92.35% | **+0.97 pp** |
| 12.5% | 94.35% | 92.62% | **+1.73 pp** |
| 25% | 94.88% | 93.01% | **+1.87 pp** |
| 50% | 95.00% | 93.71% | **+1.29 pp** |

(ρ = 0.94, renally-adjusted dosing; both converge at 0% and 100% by construction. At ρ = 0.703 the
margin is wider, peaking at **+3.38 pp** at a 40% budget.)

### 3.9.2 They select almost entirely different patients

At a 12.5% budget the two rules pick the **same patient only 16.0% of the time** (Jaccard 0.087) at
ρ = 0.94. These are not variants of one rule; they are different rules that happen to answer the same
question.

### 3.9.3 The reconciliation — Gatti's rule encodes a dosing policy

This is the result worth reporting, and it was not the one predicted.

Gatti's rule **does** do what it was built to do: the patients it selects have a median CAZ:AVI ratio
of 5.14:1 against 4.04:1 in the rest — their finding reproduces cleanly in this model. But under **this
project's renally-adjusted dosing grid**, those same patients need no assay: their median conditional
probability of attainment is 0.999, and inference is already correct in **91.6%** of them against
**92.7%** of everyone else. A 1.1-point separation is no discrimination at all.

Re-run with the dose held **fixed at 2.5 g q8h in every renal class** — the policy Gatti's own cohort
was largely on, since 85% started full dose — and the rule comes alive:

| Dosing policy | inference already correct, Gatti-selected | in the rest | separation |
|---|---|---|---|
| Renally-adjusted grid (ρ=0.94) | 91.6% | 92.7% | **1.1 pp** |
| **Fixed 2.5 g q8h (ρ=0.94)** | **88.7%** | **97.9%** | **9.2 pp** |
| **Fixed 2.5 g q8h (ρ=0.703)** | **81.4%** | **97.1%** | **15.7 pp** |

Rule agreement rises correspondingly, from 16.0% to 30.3%.

**Neither rule is wrong.** Gatti's encodes a real mechanism — high renal clearance means low avibactam
means risk — which holds when the dose does not compensate, and largely stops holding once it does.
The project's grid escalates from 1.25 g/day in the lowest EKFC class to 10 g/day in the highest, and
that escalation is what removes renal function's informativeness about who needs a second assay.

**What this is worth saying in the manuscript:** a clinically-derived triage rule validated in one
dosing context does not automatically transfer to another. That is a general point about
model-informed monitoring, it is demonstrated here with two independent rules rather than asserted, and
it is the first time any Model 2 output has been checked against something derived outside the model.

### 3.9.4 What this comparison does not establish

- **Accuracy is against the model's own definition of attainment, not a clinical outcome.** R2 winning
  on this metric means it better predicts what the model calls attainment. It is not evidence about
  patients.
- **The renal-function measures differ.** Gatti used CKD-EPI creatinine clearance; this model carries
  EKFC. Different equations on a similar scale, so the comparison is of **rule shape and ranking**, not
  of a specific mL/min cut-off. The cut-off is not transferable as a number.
- **Gatti's rule is being extended beyond its own claim.** They derived it to show the ratio varies
  with renal function, not to triage assays. Treating it as a triage rule is this project's extension,
  and its underperformance here is not a criticism of their paper.
- **R2's rule needs the ceftazidime result first**; Gatti's needs only renal function and is therefore
  cheaper to apply. The comparison assumes ceftazidime TDM is already being done, which is the premise
  of the whole question.

---

## 4. Convergence

EVPI and misselection converge more slowly than mean attainment, so both were checked across three
independent seeds:

| Draws | EVPI (mean) | EVPI (range across seeds) | Misselection | Range |
|---|---|---|---|---|
| 250 | 0.956 | 0.137 | 49.3% | 3.9 |
| 500 | 0.964 | 0.068 | 48.3% | 4.6 |
| 1,000 | 0.979 | 0.095 | 47.8% | 5.9 |
| **2,000** | **1.004** | **0.080** | **47.5%** | **4.8** |

EVPI is stable to about ±0.07 pp at the reported sample size. **Misselection remains noisy at ±5
percentage points and should be quoted as a rounded range, not a point estimate.**

---

## 5. Scenario register — what is evidence and what is assumption

Every distribution is recorded in `outputs/model2_scenario_register.csv` with its kind and its
provenance. In summary:

**Taken from the source, not invented:** all six pharmacokinetic relative standard errors
(6.36%, 14.0%, 7.4%, 12.3%, 33.3%, 30.1%) are published in Cojutti 2024 Table 2 and are simply
unused by the primary model.

**Analyst-specified, and labelled as such:** the exchange rate λ; the unbound-fraction bounds (no
uncertainty is published); the target distributions T3-T6, including the equal weights in T3, which
express no preference rather than an evidence synthesis; the agnostic correlation range C3; and the
assumption that the reported RSE on the variability terms applies to ω rather than ω².

**Evidence-derived, not analyst-specified:** T7, the three Coleman 2014 Table 2 hollow-fibre
regrowth-threshold strains — see §3.8. Still not a distribution over the clinical target; a
distribution over an in-vitro threshold that happens to be built from the actual measured values
rather than from an approximation of them.

**Excluded, with the reason stated:** the 2.5 mg/L avibactam target, which belongs to
aztreonam-avibactam and not to a ceftazidime model.

**Never done:** the two correlation estimates are never pooled. They describe different populations
and are reported as separate scenarios throughout.

---

## 5.1 A structural limitation found on 12 August 2026 — the targets may not be separable

Model 2 treats the joint target as **joint but separable**: ceftazidime must clear its MIC-based
threshold, avibactam must clear a fixed 4 mg/L threshold, and the two conditions are evaluated
independently. That is the structure Gatti 2023 defined and the whole Bologna chain uses.

**A published argument says the two are coupled, and Model 2 does not represent that.** Gatti & Pea
(*J Antimicrob Chemother* 2023;78:1556-7, [doi:10.1093/jac/dkad108](https://doi.org/10.1093/jac/dkad108);
see `data_external/JAC_exchange_measure_one_or_both/`) argue for the **effective MIC with an inhibitor
(MICi)**: avibactam concentration *changes the ceftazidime MIC*, so avibactam's role is not to clear a
threshold of its own but to move ceftazidime's. They cite Tam et al. 2022 (*JAC* 77:3130-7) for a
76.1% fT>MICi threshold associated with suppression of bacterial regrowth, and their own case series
in which lower avibactam Css produced smaller ceftazidime MIC reductions and a trend to microbiological
failure.

**If MICi is right, the separable structure is an approximation** — the avibactam target would be a
modifier of the ceftazidime target rather than a parallel condition, and "which component is limiting"
would not be a well-posed question in the form §3.3 asks it.

This is recorded as a limitation rather than fixed, for two reasons: implementing MICi needs a
concentration→MIC-shift relationship that this project has not obtained (Tam 2022 has not been read),
and doing so would change the target definition that the entire Bologna evidence chain — including the
outcome evidence in Gatti 2025 — is built on, so the comparison to that literature would break.

**It should be stated in the manuscript as a limitation**, not left for a reviewer to raise. Reading
Tam et al. 2022 is the obvious next step if it is to be addressed rather than acknowledged.

---

## 5.2 A second finding from the same sweep — the ELF penetration ratios are evaluated at the wrong concentration

This one is outside Model 2 (it concerns the lung scenarios in `revision_support/`), but it was found
by the same JAC sweep and it is the more consequential of the two, because **a stated conclusion turns
on it**.

The package applies **fixed** ELF/plasma penetration ratios of **0.52** (ceftazidime) and **0.42**
(avibactam), cited to Dimelow 2018. The citation is correct — those numbers are in the paper. But they
are not constants. Dimelow's ceftazidime plasma-ELF link is a **saturable Michaelis-Menten** function
and avibactam's is a **power** function with exponent < 1; in both cases **the penetration ratio falls
as plasma concentration rises**. The published figures are the ratio at **one plasma concentration
each**: 52% at 15.3 mg/L ceftazidime, 42% at 2.4 mg/L avibactam.

The continuous-infusion regimens simulated here run far above those concentrations — the ceftazidime
neurotoxicity screen alone is at **104 mg/L**, about seven times where 0.52 is valid.

`code/elf_penetration_concentration_check.py` rebuilds both functions from the paper's own parameters,
validates them against **all ten** numeric checkpoints the paper states (all ten reproduce), and then
evaluates them where the project operates:

| plasma ceftazidime | Dimelow's own model | applied | overstatement |
|---|---|---|---|
| 15.3 mg/L (where 0.52 is valid) | 52.2% | 0.52 | 0% |
| 70 mg/L | 32.0% | 0.52 | 62% |
| 104 mg/L (neurotoxicity screen) | **25.8%** | 0.52 | **101%** |

Avibactam is milder but the same direction: 35.0% at 10 mg/L against the 0.42 applied.

**Why it matters.** `revision_support/outputs/lung_therapeutic_window.csv` reports whether the lung
therapeutic window is open at the clinical breakpoint. Under the "central estimate" (0.52/0.42) it is
**open** (1.44× breakpoint); under the "conservative" scenario (0.30/0.30) it is **closed** (0.83×).
The qualitative answer flips between those rows — and at continuous-infusion concentrations Dimelow's
own model gives ≈26–32% (ceftazidime) and ≈32–35% (avibactam), i.e. **the row labelled "conservative"
is the concentration-appropriate one**, and it is the row where the window is closed.

Nothing outside `model_development_v18/` was modified; `revision_support/` is frozen. Three responses
are defensible and the choice is the author's: **relabel** (stop calling 0.52/0.42 the central
estimate), **re-run** with the concentration-dependent functions (fully specified and validated in the
check script), or **report as a limitation**. The first and third are Phase 6 text changes; the second
is not. Full detail, caveats, and the counter-argument on free-drug binding are in
`data_external/Dimelow2018_ELF_concentration_dependence/README_provenance.md`.

---

## 6. What Model 2 does not claim

- It is not a pharmacokinetic model and adds no pharmacokinetic knowledge.
- The target distributions are **not** established exposure-response relationships. They are
  prespecified scenarios, and the sensitivity across them is a primary result.
- The ρ scenarios are **not** a pooled posterior.
- The joint target is treated as **separable** (two independent conditions). If the MICi mechanism in
  §5.1 holds, that is an approximation and the limiting-component analysis in §3.3 is conditional on it.
- Value-of-information results are in units of **attainment probability**, not money or lives, and
  must not be translated into either.
- No output is a bedside dosing recommendation.
- Nothing here converts the study into a clinical study, and none of it speaks to bioRxiv
  eligibility.

---

## 7. Outstanding

1. ~~Layer 2, between-study heterogeneity~~ — **implemented, see §3.6.**
2. ~~Unit tests for the decision layer~~ — **implemented, 44 checks, `code/test_model2.py`,
   documented in `SOFTWARE.md`.** Writing them found and fixed the C1 sampler bug corrected
   throughout this document.
3. ~~The breaking-point analysis~~ — **implemented as R7**, see §3.7 above; not originally listed
   here but delivered as part of the same push.
4. ~~An evidence-derived target scenario~~ — **implemented as R4**, see §3.8 above; corrects an
   earlier internal citation error along the way (Coleman 2014's regrowth threshold rests on three
   strains, not eight).
5. **Figures** for the EVPPI ranking, the monitoring gain against ρ, and the limiting-component
   probability surface.
6. **The falsifiable prediction** (specification §3.7 of `MODEL2_SPECIFICATION.md`), to be written
   only after everything above is frozen.
