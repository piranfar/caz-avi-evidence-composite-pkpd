# Pre-registered prediction — locked before any non-RRT continuous-infusion cohort is held

**Status: PRE-DATA.** As of 12 August 2026 this project holds **no** individual patient data from a
non-RRT continuous-infusion cohort. The Bologna data request was sent 11 August 2026 and is
outstanding; no reply has been received and no follow-up sent. The Barcelona request is drafted and
unsent. That is what makes this document a prediction rather than a description.

**Lock:** this file is committed to git on 12 August 2026. Its commit hash and timestamp are the
evidence of priority. If it is edited after data arrive, the diff will show it. Any revision must be
made as an appended amendment with its own date, never by rewriting what is below.

Produced by `code/prereg_npv_prediction.py`; numbers reproducible from that script.

---

## 1. What was originally proposed, and why it was changed

`NOVELTY_STRATEGY.md` Part I §4.3 asked for:

> "the predicted **negative predictive value** of ceftazidime-based prediction in a non-RRT
> continuous-infusion population, with an interval."

**That endpoint was tested for power and rejected.** It is recorded here rather than quietly swapped,
because the reason it fails is itself the useful finding.

NPV is estimated only on patients **predicted not to attain** the avibactam target. In a cohort of
112 that is about **10 patients**. The two rival correlations imply NPVs of 69.6% (ρ = 0.703) and
84.7% (ρ = 0.94) — a 15 percentage-point signal — but the binomial interval on 10 patients runs from
roughly 39% to 89%. The intervals overlap almost completely.

| endpoint | signal (ρ=.703 vs ρ=.94) | separable at N = 112? | separable at all? |
|---|---|---|---|
| **NPV**, 0% assay CV | 69.6% vs 84.7% | **no** | only at N ≈ 1600 |
| **NPV**, 20% assay CV | 65.4% vs 75.1% | **no** | not even at N = 1600 |
| **Overall accuracy**, 0% CV | 88.0% vs 94.7% | **no** | gap only 6.7 pp |
| **Overall accuracy**, 20% CV | 87.2% vs 91.8% | **no** | gap only 4.6 pp |
| **ρ estimated directly** | 0.703 vs 0.94 | **YES** | **yes, even at N = 50** |

**Pre-registering NPV would have produced an interval no realistic result could contradict — the
appearance of falsifiability without the substance.** That is worse than not pre-registering at all,
because it buys credibility that the design does not earn.

The endpoint that works is the one that uses every paired sample rather than a tenth of them.

---

## 2. The prediction

> **Primary, falsifiable:** In a cohort of ≥50 non-RRT critically ill adults receiving
> ceftazidime/avibactam by continuous infusion, with paired steady-state concentrations of both
> analytes, the between-subject clearance correlation ρ, estimated as the random-effect correlation
> after adjustment for renal function, will be **below 0.90**, and its 95% confidence interval will
> **exclude 0.94**.

**Point prediction: ρ = 0.75, 95% prediction interval 0.55 to 0.87.**

### Why this number

Model 1 estimates ρ = 0.703 (95% CI 0.381–0.873) from the two openly licensed patient-level datasets,
both of which are **on renal replacement therapy**. Cojutti 2024 reports 0.94 in a non-RRT
continuous-infusion cohort — the target population itself.

The prediction is deliberately placed **above** Model 1's estimate and **below** Cojutti's. Two
reasons, stated so they can be judged:

1. RRT imposes a shared extracorporeal clearance pathway on both analytes, which should *inflate*
   their correlation, not deflate it — so 0.703 from RRT cohorts is arguably an upper bound for what
   a non-RRT population would show, pushing the prediction down.
2. Against that, 0.94 comes from the target population itself and from a group whose other reported
   parameters have been verified against their own tables in this project without discrepancy
   (`MODEL2_REPORT.md` §0). It is unverified, not discredited.

**This prediction says Cojutti's 0.94 will not replicate.** That is a real risk, and it is the point.

### What refutes it

| observed | verdict |
|---|---|
| 95% CI for ρ excludes 0.94 **and** point estimate < 0.90 | **prediction holds** |
| 95% CI for ρ includes 0.94, or point estimate ≥ 0.90 | **prediction refuted** |
| 95% CI so wide it includes both 0.75 and 0.94 | **uninformative** — report as such, do not claim support |

The third row matters. A wide interval is not a pass.

### Consequences if refuted

Route **R1** — the argument that the correlation is assumed rather than established — is weakened but
not destroyed: R1's evidential claim is that ρ was *computable in three published datasets and
reported in one*, which stays true whatever value it takes. What would be lost is the stronger
implication that the assumed value is *probably too high*.

**The monitoring recommendation would survive.** §2 shows EVPPI(ρ) = 0.0000 pp on the measure-or-infer
decision, and §3.7 shows measuring beats inferring at every ρ from 0.30 to 0.99. That conclusion was
built not to depend on ρ, and this is the test of whether that was true.

---

## 3. Analysis plan, fixed in advance

1. **Population.** Non-RRT adults on continuous-infusion ceftazidime/avibactam, paired steady-state
   concentrations of both analytes. Patients on any form of RRT are excluded.
2. **Estimand.** The between-subject random-effect correlation of ceftazidime and avibactam clearance
   **after** adjustment for renal function — the same conditional quantity Cojutti reports, not a
   marginal correlation of raw concentrations. This distinction is load-bearing: `NOVELTY_STRATEGY.md`
   Part I §3 shows adjustment moves the estimate *down* (0.563 marginal → 0.476 conditional in the
   Dryad cohort), so comparing a marginal estimate against Cojutti's conditional 0.94 would be
   rigged in the prediction's favour.
3. **Method.** Joint population PK model with a full random-effect covariance block, as in
   `MODEL1_REPORT.md`. If the joint model fails to converge, a two-stage estimate may be reported
   **only** with the regression-dilution caveat already recorded at the top of `NOVELTY_STRATEGY.md`,
   and it does not count as the primary test.
4. **Interval.** 95% confidence interval from the model's covariance estimate; bootstrap if the
   asymptotic interval is unreliable near the boundary.
5. **No subgroup fishing.** The primary test is the whole eligible cohort. Any subgroup result is
   exploratory and must be labelled so.
6. **Reported regardless of direction.** Including if it refutes the prediction.

---

## 4. Secondary quantities — descriptive, explicitly NOT falsifiable

Recorded so they can be compared against, but they are **underpowered by design** and no claim of
support or refutation may rest on them.

| quantity | ρ = 0.703 | ρ = 0.94 | 95% credible interval (parameter uncertainty) |
|---|---|---|---|
| NPV, 0% assay CV | 69.6% | 84.7% | 58.4–79.1 / 72.1–92.4 |
| NPV, 20% assay CV | 65.4% | 75.1% | 54.1–74.6 / 59.1–85.3 |
| Overall accuracy, 0% CV | 88.0% | 94.7% | — |
| Overall accuracy, 20% CV | 87.2% | 91.8% | — |
| Share predicted negative | ~9% | ~14% | — |

If a future cohort reports NPV near 70%, that is **consistent with** ρ = 0.703 and does not establish
it. The intervals above overlap.

---

## 5. Honest limits of this pre-registration

- **A model predicting its own parameter is a weak test.** The prediction is not derived from
  independent theory; it is an extrapolation from two RRT cohorts to a non-RRT one. It is falsifiable,
  which is the claim being made here — not that it is a strong prior.
- **N = 112 is the Bologna cohort as described in the data request.** If the cohort delivered differs
  in size or composition, the power table in §1 must be recomputed *before* unblinding, not after.
- **Fisher's z is an approximation** for a simple correlation; ρ here is a random-effect correlation
  from a hierarchical fit, whose precision depends on sampling density per subject as well as N. The
  approximation is sanity-checked against Model 1's achieved precision (95% CI 0.381–0.873 at n = 21;
  scaling by √(112/21) = 2.31 gives a width near 0.21, matching the Fisher-z row), so it is not
  flattering itself — but it remains an approximation, and a sparsely sampled cohort would do worse.
- **This does not convert the study into a clinical study**, and nothing here speaks to outcomes.
- **The prediction may simply be wrong.** That is what it is for.
