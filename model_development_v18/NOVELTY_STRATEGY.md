# Novelty strategy — how to make this paper distinctive

> ## ⚠ CORRECTION — superseded numbers
>
> This document quotes the ceftazidime-avibactam clearance correlation as **0.560-0.598**, from
> two-stage analyses. Those estimates are **attenuated by estimation error in the individual
> parameters** (regression dilution). The mixed-effects model in `MODEL1_REPORT.md` estimates the
> comparable quantity at **0.703 (95% CI 0.381 to 0.873)**, which still excludes the assumed 0.94
> but by a smaller margin. **`MODEL1_REPORT.md` is authoritative wherever the two disagree.**
>
> Consequently the "roughly threefold" increase in false reassurance quoted below is
> **2.6-fold** at the fitted value: 3.6% at rho = 0.94 against 9.2% at rho = 0.703.

> ## Status, 12 August 2026
>
> **Part I below is the original strategy note and is preserved verbatim.** Several of its
> recommendations have since been carried out; those are annotated in place rather than rewritten, so
> the reasoning that led to them stays readable.
>
> **Part II is new** and addresses a question Part I raises in §7 but does not answer: *given that a
> simulation study on published parameters has a novelty ceiling, what have comparable bodies of work
> actually done about it?*

**Decision taken (author, this session):** abandon the current IJAA submission
(IJAA-S-26-02124) rather than revise it; build both new models; maximise novelty.

---

# Part I — the original strategy note

---

## 1. Diagnosis: where the novelty ceiling actually is

The current central claim is:

> *"The choice of avibactam target had a major effect on the estimated target attainment."*

A fair reviewer will observe that this is close to arithmetic. If you raise a threshold, fewer
patients clear it. The 20.29 pp sensitivity to the avibactam target is a real and well-executed
result, but it demonstrates that the model is sensitive to its own input — not that anything about
ceftazidime-avibactam was previously unknown.

**This is the ceiling, and it has nothing to do with how many components the model has.** The paper
already contains 21 supplementary tables, 6 figures, 6 tables, ELF scenarios, four structural
models, PSA, GSA, multi-seed replication, calibration, decision grids, dose escalation, resistance-
suppression targets and an ARC subgroup. Adding more will read as unfocused, not as thorough.

**Adding model components lowers the score in this field.** At AAC, JAC, CPT:PSP and IJAA the most
common reason a modelling paper is rejected is complexity that the data cannot identify. Every added
compartment or module invites the question *"what identifies this?"* With 21 subjects, each addition
is a liability.

---

## 2. The claim the paper should be built around instead

> The two components of ceftazidime-avibactam are assumed to move together — a random-effect
> clearance correlation of 0.94 — and **that assumption is what makes single-analyte therapeutic drug
> monitoring appear sufficient.** Using the only two openly licensed patient-level datasets in
> existence, the correlation is estimated at approximately 0.48-0.60. At that value, ceftazidime-based
> prediction misclassifies avibactam attainment roughly three times as often as assumed. Avibactam
> should therefore be measured, not inferred.

Why this is a stronger spine:

- It is a **clinical claim with a decision consequence**, not a sensitivity observation.
- It is **anchored in real patient data**, from two independent cohorts.
- It is **falsifiable** — a non-RRT continuous-infusion cohort could refute it tomorrow.
- It **attacks a load-bearing assumption**, which is what makes methodological work matter.

The target-uncertainty analysis becomes the second act rather than the first: having shown that one
assumed constant does not survive contact with data, the paper asks what happens when the other
assumed constant — the avibactam target — is also allowed to be uncertain.

---

## 3. A technical point that makes the correlation argument rigorous

Cojutti's 0.94 is a **conditional** correlation: between random effects, after the renal-function
covariate. A marginal correlation computed from raw clearances is not comparable to it, and with a
shared covariate driving both analytes the marginal is normally the **larger** of the two.

Both quantities were therefore computed on the Dryad cohort:

| Quantity | Estimate |
|---|---|
| Marginal (unadjusted) | 0.563 |
| **Conditional** on CRRT modality, ultrafiltration, effluent flow and log serum creatinine | **0.476** (95% CI 0.012-0.772) |

The like-for-like comparison is therefore **0.476 against 0.94**, and adjustment moves the estimate
*down*, exactly as theory predicts. Three independent routes — Dose/AUC, model-based individual
fits, and the Gatti cohort — give 0.560, 0.563 and 0.598 marginally.

**The limitation stands and must be stated every time:** both cohorts are on renal replacement
therapy, and the primary scenario excludes it. This is an argument that 0.94 is unverified and
probably high, not a measured replacement value for the non-RRT population.

---

## 4. Three additions that raise novelty — none of which is a model component

### 4.1 Value-of-information analysis  ★ highest value

Expected value of perfect information on the two parameters that matter: the clearance correlation,
and the avibactam target. This answers the question a clinician and a hospital formulary actually
have — *is it worth buying the second assay, and is it worth running the trial?* — in the currency of
patients correctly classified.

Value of information is standard in health-technology assessment and essentially absent from the
beta-lactam target-attainment literature. It is computable from machinery this project already has.
This is the strongest single candidate for the novelty peak.

> **✅ Built** (`MODEL2_REPORT.md` §2, EVPI/EVPPI, `code/model2_hujam.py`, `model2_monitoring.py`)
> **and ✅ promoted, 12 Aug 2026** — `MODEL2_REPORT.md` §0.1 now states VOI as the organising frame
> before any result appears, with the shortage of data as the *reason the question is worth asking*
> rather than an apology. This section was right that it was "the strongest single candidate for the
> novelty peak"; what was missing was never the analysis, only where it sat.

### 4.2 Limiting-component probability under uncertainty

Replace the categorical label *"avibactam is limiting at MIC 2"* with
P(avibactam is limiting | MIC), integrating over target and parameter uncertainty. A quantity, not a
verdict. New for this drug combination.

> **✅ Built** (`limiting_probability()`, `MODEL2_REPORT.md` §3.3). Since corroborated externally:
> Cojutti 2026 measured both analytes in 17 patients and found avibactam was the failing component in
> 6 of 7 empirical failures and every targeted failure
> (`data_external/Cojutti2026_FN_avibactam_limiting/`).

### 4.3 A pre-registered, falsifiable prediction

State, before any such cohort is analysed, the predicted negative predictive value of
ceftazidime-based prediction in a non-RRT continuous-infusion population, with an interval. A
modelling paper that exposes itself to refutation is worth more than one that cannot be wrong.

> **✅ Done 12 Aug 2026 — `PREREGISTRATION.md` — but the endpoint proposed here does not work.**
>
> NPV was tested for power before being registered and **rejected**. It conditions on the patients
> predicted negative — about 10 of 112 — so the binomial interval (≈39% to 89%) swamps the 15 pp
> signal between ρ = 0.703 and ρ = 0.94. It separates only at N ≈ 1600, and at a realistic 20% assay
> CV **not even then**. Overall accuracy fails too (gap 6.7 pp). Registering NPV would have bought the
> appearance of falsifiability without the substance, which is worse than not registering.
>
> **ρ estimated directly is the endpoint that works** — it uses every paired sample and separates the
> rivals even at N = 50. The registration is built on that instead. The idea in this section was
> right; the endpoint was not, and the power check is what caught it.

---

## 5. What to remove

Seriously consider cutting or demoting to the supplement: the resistance-suppression target analysis,
the ARC subgroup, the dose-escalation crossing analysis, and two of the four MIC distributions
(including `FR2024_OXA484_SERINE_ONLY`, which has no citation and must go regardless). A shorter paper
with one sharp claim outranks a long paper with six soft ones.

---

## 6. The highest-leverage action available

**Request the Bologna dataset: Cojutti / Gatti / Pea, 112 non-RRT critically ill adults, 185 paired
steady-state concentrations under continuous infusion.**

This is precisely the population of the primary scenario. If obtained:

- The correlation can be estimated **in the right population**, converting "we question this
  assumption" into "we re-estimated this assumption".
- A joint population PK model becomes fittable to the correct population — Candidate A moves from
  rejected to feasible.
- True external validation becomes possible.

That combination is unambiguously a top-tier contribution. Two things make it plausible: the same
group already publishes complete patient-level tables in their 2023 paper, so they are demonstrably
willing to share; and **abandoning the current submission removes the time pressure** that made a
2-3 month wait impossible.

Priority order for data requests:
1. **Bologna** (Gatti `milo.gatti2@unibo.it`, Pea) — the primary-scenario population.
2. **Barcelona** (Sorlí `lsorli@hmar.cat`) — ELF and an independent non-RRT continuous-infusion cohort.
3. Vivli / Pfizer — registrational data, 3-9 months, in-platform analysis only.

**No message has been sent. Both drafts require explicit authorisation.**

> **Partly overtaken.** The Bologna request was sent on 11 August 2026 and is outstanding; no
> follow-up, no second request. Barcelona (Sorlí) is drafted but unsent. R5 remains on hold pending a
> reply. The judgement that this is the highest-leverage action **still stands** — nothing in Part II
> substitutes for the data, and Part II says so.

---

## 7. Honest ceiling

A simulation study built predominantly on published parameters will not reach "10/10 novelty" on the
strength of methodology alone, and no amount of added model structure will change that. What is
achievable is a paper with **one distinctive, data-anchored, falsifiable claim** plus a decision
framework that is genuinely new to the field.

With the Bologna data, the ceiling rises substantially, because the study would then re-estimate a
published parameter in its own population rather than argue about it from adjacent ones.

**This will not convert the study into a clinical study, and it does not guarantee bioRxiv
eligibility.** What has genuinely changed is that the work now incorporates individual patient data
from two independent cohorts rather than published aggregates alone.

---
---

# Part II — how comparable work escaped the "no new data" ceiling

*Added 12 August 2026. Part I §7 states the ceiling honestly and stops there. This part asks what
bodies of work facing the same ceiling actually did, and what transfers.*

## 1. The problem class

Every input is a published summary statistic any reviewer can look up, and the primary output is a
Monte Carlo target-attainment simulation — a genre reviewers rate as competent-but-ordinary. Part I
§1 diagnoses this correctly and concludes that adding model components makes it worse. It is right.
But it leaves open what makes it *better*.

Several bodies of work have escaped this exact position. **None of them compensated for the data they
lacked. All of them changed the object of study.**

| from | to |
|---|---|
| the phenomenon | the literature about the phenomenon |
| one error | a class of error |
| "which is right?" | "what would change the answer?" |
| "what is the answer?" | "what is it worth to find out?" |

## 2. Five precedents

**2.1 Shift the object of study — Ioannidis 2005.** *"Why Most Published Research Findings Are False"*
used no new data: an analytical argument over published base rates and power. The move was to study
the epistemics of a field rather than its subject matter.

*Licenses:* the project's two strongest findings are already of this type and are currently buried —
R1 (the cross-drug correlation was computable in three published datasets and reported in one) and
§5.2 (a concentration-specific penetration ratio applied as a constant). Neither is a simulation
result. Both should lead.

**2.2 Find an error *class*, not an *instance* — Westreich & Greenland 2013.** The "Table 2 fallacy"
paper used no new data; it identified a mistake most of a field makes. **This was checked rather than
assumed.** Fixed ELF/plasma penetration ratios are standard practice across the beta-lactam PTA
literature — Layios 2022 ([10.1128/AAC.02052-21](https://doi.org/10.1128/AAC.02052-21)) holds a ratio
constant "whatever the level of creatinine clearance"; Bader 2019
([10.1128/AAC.00318-19](https://doi.org/10.1128/AAC.00318-19)) uses a single ratio as the Monte Carlo
input; Lepak 2023 ([10.1128/aac.01452-22](https://doi.org/10.1128/aac.01452-22)) uses literature
ratios **to set susceptibility breakpoints**. Meanwhile Dimelow 2018 shows the relationship is
*saturable* for CAZ-AVI, and Kawaguchi 2022 ([10.1002/jcph.1986](https://doi.org/10.1002/jcph.1986))
shows that building a proper intrapulmonary model for cefiderocol changed the answer by 1.4×.
*(Retrieved from PubMed.)*

*Licenses:* a claim that the fixed ELF ratio is a systematically biased input wherever the plasma-ELF
relationship is non-linear.

> **⚠ Checked 12 Aug 2026, and the broad version is dead.** §5 below required verifying that other
> drugs are non-linear before any error-class paper. They are not. **imipenem** was explicitly fitted
> with a *time-independent* penetration coefficient of 0.44 across three populations (van Hasselt
> 2016); lefamulin, ceftaroline and sulbactam/durlobactam are all modelled linearly. Writing "the
> field's fixed ratios are systematically biased" would be **false**, and the imipenem paper alone
> would refute it.
>
> **What survives is narrower and closer to this project's actual subject:** for **BL/BLI
> combinations**, the two components do not penetrate alike, and the joint target is evaluated by
> applying one fixed ratio per component — usually from healthy volunteers, at concentrations unlike
> continuous infusion. Untested in most analyses, and demonstrably wrong in at least two combinations
> (CAZ/AVI saturable; piperacillin/tazobactam 49.3% vs 121.2%, tracking pulmonary permeability rather
> than dose). Sharpened by the fact that **the inhibitor is usually the limiting component**, so the
> component least carefully characterised is the one that decides attainment.
>
> Realistic novelty **7–8/10**, not 9–10. Full evidence both ways:
> `data_external/Dimelow2018_ELF_concentration_dependence/is_nonlinearity_general.md`.

**2.3 Make the dispute decidable — Bell 1964.** No data; an inequality that made a decades-old
disagreement decidable.

*Licensed, and now executed.* See §4 below — this produced a result, not just a reframing.

**2.4 Convert an unanswerable empirical question into an answerable decision question — Claxton &
Sculpher.** Health economics could not answer "which treatment is better?" without trials it lacked,
so it asked "what is it worth to find out?" VOI became a subfield and now informs NICE research
prioritisation.

*Licenses:* exactly what Part I §4.1 already identified and built. The remaining gap is presentational
and is the next action (§3 below).

**2.5 The publishing precedent — reanalysis of a live controversy.** Bayesian reanalyses of contested
trials publish in high-impact journals with no new data, because the controversy already exists and
the paper supplies the missing frame.

*Licenses:* the four-round JAC exchange, where round 3's actionability objection was never answered.

## 3. Ranked options

| option | novelty | feasible now | cost |
|---|---|---|---|
| ~~ELF error-class paper, broad version~~ | ~~9–10/10~~ | **✗ refuted 12 Aug — imipenem is linear** | — |
| **Narrow paper: fixed-ratio ELF scaling in BL/BLI combinations** | **7–8/10** | yes — evidence assembled both ways | new paper |
| ~~Reframe §3.7 as the dispute's decision boundary~~ | 8/10 | **✅ done — and it yielded a bound (§4)** | — |
| ~~Promote VOI from add-on to organising frame~~ | 7/10 | **✅ done — MODEL2_REPORT.md §0.1** | — |
| R1 as "computable in three, reported in one" | 7/10 | yes — evidence assembled | text only |
| ~~Pre-registered falsifiable prediction (Part I §4.3)~~ | 7/10 | **✅ done — `PREREGISTRATION.md`** | — |
| A well-executed PTA simulation | **4/10** | already done | — |

## 4. What the Bell move actually produced

Reframing §3.7 required closing a gap that had not been noticed: Fresan et al. assume a correlation
between **target achievements** (binary), while Model 2 parameterises one between **clearances**.
Computing the mapping (`code/model2_dispute_boundary.py`) showed it is strongly compressive — a
clearance correlation of **0.99 induces an attainment correlation of only 0.52**.

The reason is structural. At the breakpoint ceftazidime attains in ~58% of patients and avibactam in
~85%, and two binary outcomes with unequal prevalences cannot be strongly correlated: the
Fréchet–Hoeffding limit for those margins is **φ ≤ 0.52**, for any joint distribution under any
pharmacokinetics.

**So the single-analyte position does not fail because the correlation is low in this dataset. It
fails because the correlation it assumes is bounded above by the mismatch in attainment rates** — and
at that ceiling, inference still misclassifies. Derived, not measured.

*Caveats that must travel with it:* the **inequality** is general, the **number** is not (0.50 is this
model's prevalences); and Fresan's second objection, that avibactam assays are not routinely
available, is untouched — that is answered by the triage rule in §3.9.

## 5. Honest limits of Part II

- **These are framings, not results** — with one exception (§4). Better framing exposes weak analysis
  faster, not slower.
- **The ELF error-class claim is verified as a *practice*, not as a *quantified field-wide bias*.**
  Fixed ratios are standard and the relationship is saturable for CAZ-AVI; whether other drugs'
  plasma-ELF relationships are non-linear has **not** been checked. Without that, the claim is about
  ceftazidime/avibactam only. Temocillin's constant 0.73 is *suspicious*, not *shown wrong*.
- **Nothing here substitutes for data.** Part I §6 is still right: the Bologna dataset remains the
  highest-leverage action available, and every clinical source archived so far is class C.
- **Precedent is not permission.** Ioannidis and Bell are cited as illustrations of a *move*. The
  manuscript should never invoke them.

## 6. Next actions

1. ~~Reframe §3.7 as the decision boundary of the Fresan–Gatti dispute.~~ **✅ done — yielded §4.**
2. ~~Promote VOI from technical add-on to the organising frame.~~ **✅ done — MODEL2_REPORT.md §0.1.**
   Framing only; no number changed. The two VOI results are now stated as a research-prioritisation
   claim — *resolve the target, not the correlation* — with the zero EVPPI on ρ presented as a
   finding rather than a null. Checking the numbers while writing it caught a real hazard: ρ has
   EVPPI 0.021 pp on decision A and 0.0000 pp on decision B, and the draft table had silently mixed
   the two decisions. Now labelled.
3. **Lead with the meta-research findings (R1, §5.2) rather than simulation output.** *(next)*
4. Decide the ELF finding's home — separate paper, or a limitation paragraph. Not equivalent.
5. ~~Revisit Part I §4.3, the pre-registered prediction.~~ **✅ done — `PREREGISTRATION.md`.** The
   power check inverted the proposal: NPV cannot be falsified by any realistic cohort, so ρ is
   registered instead. Prediction: ρ < 0.90 with a 95% CI excluding 0.94, point 0.75 (0.55–0.87).
   **This says Cojutti's 0.94 will not replicate**, and if refuted, R1 is weakened while the
   monitoring recommendation should survive — which is itself a test of the claim that it never
   depended on ρ.
6. ~~Before any ELF paper, check whether other drugs' plasma-ELF relationships are non-linear.~~
   **✅ done — they are not.** The broad error-class paper should not be written; its central claim
   would be false. The narrow BL/BLI version is defensible at 7–8/10. **`MODEL2_REPORT.md` §5.2 is
   unaffected** — it concerns this package's own misuse of Dimelow's ratios and stands independently.
