# Routes to higher novelty — a live register

Companion to `NOVELTY_STRATEGY.md`, which diagnosed where the ceiling is. This file lists the routes
that could raise it, ranked, with what each would cost and what could kill it. **Kept updated as
routes are opened, tested or abandoned.**

**Status key:** ✅ done · 🔄 in progress · ⏸ available, not started · ❓ blocked on a fact · ❌ rejected

---

## Already banked

| # | Contribution | Status |
|---|---|---|
| B1 | Second published estimate of the ceftazidime–avibactam clearance correlation; first from openly available individual data | ✅ |
| B2 | External aggregate-level predictive check against an independent ICU cohort | ✅ |
| B3 | Decision-analytic layer with value of information — new to beta-lactam target attainment | ✅ |
| B4 | Design analysis: what a future dataset must look like to settle the question | ✅ |
| B5 | Falsifiable prediction, stated before the data exist | ✅ |

These make a good paper. They do not by themselves make a landmark one, because every element is
about a single drug pair.

---

## R1 — Generalise to the whole BL/BLI class ⚠️ *checked — PARTIALLY CONFIRMED, claim must be reworded*

**The idea.** The problem is not specific to ceftazidime-avibactam. Every
beta-lactam / beta-lactamase-inhibitor product has two components, two different pharmacodynamic
targets, and clearances that may or may not move together — and in every case the cross-drug
clearance correlation determines whether measuring one component tells you about the other.

If it turns out that **almost no published model for any BL/BLI combination has ever estimated that
correlation**, then the paper stops being about one drug and becomes: *a parameter that silently
decides therapeutic drug monitoring adequacy across an entire drug class is, with one exception,
unestimated everywhere.*

**Why this is the highest ceiling.** It converts a drug-specific analysis into a methodological
contribution with class-wide reach, and it reframes the finding from "this number looks wrong" to
"this number is missing, and its absence has consequences." Reviewers reward the second far more.

**What it needs.** Verification, currently under way, across ceftazidime/avibactam,
ceftolozane/tazobactam, meropenem/vaborbactam, imipenem/relebactam, aztreonam/avibactam,
piperacillin/tazobactam, and the newer agents. The critical distinction is between a *within*-drug
correlation such as corr(η_CL, η_V) — common, and not the point — and a *between*-drug correlation
linking the two components.

**Result of the check.** The literal claim — "Cojutti 2024 is the only published estimate" — is
**REFUTED**. Two more exist, both verified through independent full-text reads:

- **Piperacillin/tazobactam, ρ(CL) = 0.93 (RSE 11.5%)** — Cojutti/Pai et al., *AAC* 2024;68(4):e0140423,
  PMID 38411995. **The same Bologna/Michigan group**, published months before the anchor paper.
- **Aztreonam/avibactam, ρ(CL) = 0.976, ρ(Vc) = 0.986** — Xie et al. (Pfizer), *AAC* 2025;69(8):e0195024,
  PMID 40530972. An independent sponsor model, higher magnitude than Cojutti's estimate.

A related but weaker case: Wallenburg et al. 2022 (piperacillin/tazobactam, Clin Pharmacokinet,
PMID 35377133) forces near-unity between-drug correlation through a structural CL-link and shared
random effects rather than reporting a rho. And Kong et al. 2024 (pooled piperacillin/tazobactam,
n=415, PMID 39722108) explicitly TESTED for the CL-CL correlation and could not identify it, while
succeeding for the volume terms — informative negative evidence.

**The broader pattern survives, and is reinforced.** Two clean confirmatory non-hits were found:
sulbactam/durlobactam (Cammarata 2024, PMID 39569973 — fully joint 4-compartment model, reports only
within-drug CL-Vc covariance, no cross-drug term anywhere) and imipenem/cilastatin/relebactam
(Patel/Bellanti 2022, Merck, PMID 34704389 — same pattern). Ceftolozane/tazobactam is fitted as a
FULLY SEPARATE model in one recent paper (Gatti 2025, PMID 40980909) despite having TDM data for
both drugs — the cleanest illustration of the underlying problem. Meropenem/vaborbactam and the
minor combinations remain genuinely unverified (paywalled or sparse literature), not confirmed
non-hits.

**Revised, defensible claim:** *"A between-component clearance correlation has been published for
only three of more than ten marketed BL/BLI combinations (ceftazidime/avibactam 0.94,
piperacillin/tazobactam 0.93, aztreonam/avibactam 0.976-0.986), and even where reported it is
typically an incidental by-product of a joint model built for another purpose rather than a routinely
estimated component. For the majority — including at least two combinations where a fully joint model
WAS fitted with a full OMEGA block available — it is not estimated at all."* Three isolated estimates
out of a large model literature is still a striking and citable fact; "the only one" is not.

**A third category, added 12 August 2026 — and it is the strongest single example R1 has.**
Archiving the registrational model (`data_external/Li2019_registrational_PopPK/`) showed that
**Li et al. 2019 acknowledges the ceftazidime–avibactam correlation, preserves it, and never computes
a value.** The two drugs are fitted as separate models with no cross-drug OMEGA block. Full-text
search: **"OMEGA" 0 hits, "omega" 0 hits, "off-diagonal" 0 hits, "rho" 0 hits** — every
covariance/correlation hit is within-drug.

**The supplement has since been obtained, and it makes the point far sharper.** Data S2 gives the
method in full:

> "Each patient's random effects from the ceftazidime and avibactam models were **merged into a single
> data file, in which there was a single record for each patient containing all his/her random effects
> for both compounds**. Patient-level random effect records were bootstrapped ... **thereby preserving
> any underlying correlations.**"

So the sponsor **built the exact data structure a cross-drug correlation is computed from** — one row
per patient, both drugs' random effects side by side — and then resampled intact patients to carry the
correlation implicitly rather than estimating it. The number was not merely left unpublished: **it was
never calculated**, by the group with 1,975 / 2,249 patients across five phase III trials, with the
data already in the right shape, in the analysis that set the licensed dose.

R1's claim can now distinguish three states rather than two: *reported* (3 combinations),
**acknowledged, structurally available, and still never computed (Li 2019, ceftazidime/avibactam)**,
and *not addressed at all* (the rest). "They had it one line of code away and did not take it" is a
much stronger sentence than "not estimated" — and it is now evidenced verbatim from the sponsor's own
supplementary methods, not inferred.

### Second instance found by the JAC sweep, 12 August 2026 — Falcone 2021

Falcone et al. (*JAC* 2021;76:1025-31, [dkaa549](https://doi.org/10.1093/jac/dkaa549)) assayed
**ceftazidime, avibactam and aztreonam in the same 41 patients** and estimated individual clearance
for all three. The two clearances are reported in adjacent columns of the same table. **The
correlation between them is never computed.**

It is not the same evidence as Li 2019, and the difference should be stated rather than blurred:

| | structure | how far from ρ |
|---|---|---|
| **Li 2019** | both drugs' random effects merged per patient | one line of output |
| **Falcone 2021** | drugs modelled separately | post hoc correlation of individual EBEs |

R1's claim is now sharper and countable: **the cross-drug clearance correlation has been computable in
at least three published ceftazidime/avibactam datasets — Li 2019 (n≈2,000), Falcone 2021 (n=41),
Cojutti 2024 (n=142) — and reported in exactly one** (Cojutti 2024, ρ = 0.94). One reported value out
of three opportunities is a documented gap in the literature, not a question this project invented.

Corroborating detail from the same paper: median C_min 57.9 mg/L ceftazidime against 10.9 mg/L
avibactam, a **5.31:1** ratio of medians rather than the 4:1 of the vial — an independent replication
of Gatti 2024's central finding in a different centre and population. Full extraction in
`data_external/Falcone2021_three_analyte_PopPK/`.

**Cost if pursued:** moderate, and lower than expected — the evidence is already assembled above.
**Full audit trail:** `PHASE2_DATA_AVAILABILITY_REPORT.md` companion note, or ask for the full
subagent report to be written out as a standalone file.

---

## R2 — A triage monitoring rule ✅ *done*

**The idea.** Model 2 answered a binary question: measure avibactam, or infer it. The clinically
useful question is narrower — **which patients need the second assay?**

Inferring avibactam from ceftazidime fails mainly near the decision boundary. For a patient whose
ceftazidime concentration is far above or far below the implied threshold, the inference is safe;
for one near it, the conditional probability of attainment is close to one half and the inference is
worthless. So the rule swept a confidence band and measured only patients inside it, with the band
expressed as an actual ceftazidime concentration window per renal class — a form a clinician can act
on directly.

**Done — `code/model2_triage.py`, results in `outputs/model2_triage_curve.csv`.** Two findings:

- At ρ = 0.94 (the manuscript's assumption), **90% of the benefit of measuring everyone is reached
  by measuring only 12.5% of patients** — those with total ceftazidime Css in a window of roughly
  13–29 mg/L depending on renal class (e.g. 12.9–18.8 mg/L in the lowest class, 19.7–28.7 mg/L in
  the highest). **Selective measurement in fact beats measuring everyone outright** — 95.3% accuracy
  at 28.9% of patients measured, against 94.9% when every patient is measured — because assaying a
  patient whose inference was already confident only adds assay noise.
- At ρ = 0.703 (Model 1's estimate), the same 90% coverage needs a much larger **39.8% of patients**,
  with correspondingly wider windows (roughly 6–44 mg/L depending on class). This is a second,
  independent way the corrected correlation matters clinically: not just whether to measure, but how
  selectively.

This result uses fixed ρ values directly rather than the scenario samplers, so it was **not affected
by the C1 sampler bug** fixed during R6 and needed no rerun.

**Not a dosing recommendation** — a monitoring-strategy analysis under the model, evaluated against
the model's own definition of attainment, not a clinical outcome. Stated as such in the code.

**A clinically-derived comparator now exists, and comparing against it is an open, novel check.**
Gatti, Viale & Pea 2024 (`data_external/Gatti2024_ratio_one_leg/`) asks R2's question with real TDM
data in 107 patients and answers the *binary* version: the ceftazidime-to-avibactam ratio ranges
1.29:1 to 13.46:1 against a 4:1 vial ratio, so avibactam cannot be extrapolated from ceftazidime and
**both should be measured**. That is directionally consistent with R7's breaking-point result, but R2
goes further than the paper does — it answers *which* patients, which the paper leaves open.

The paper also offers its own triage rule, derived from ROC analysis rather than from a model:
**CrCL > 75–78 mL/min/1.73 m²**, or **serum urea ≤ 45–51 mg/dL**, identify the patients whose ratios
exceed 5:1 and 6:1 (AUCs 0.685–0.717).

**✅ That comparison has now been run** — `code/model2_triage_vs_gatti.py`, full results in
`MODEL2_REPORT.md` §3.9. It is the first external check any Model 2 output has had, and the
disagreement turned out to be the interesting part:

1. **R2's rule dominates at every assay budget** — +1.73 pp at a 12.5% budget (ρ = 0.94), up to
   +3.38 pp at ρ = 0.703.
2. **The two rules pick almost entirely different patients** — the same patient only 16.0% of the time
   at a 12.5% budget (Jaccard 0.087). Not variants of one rule.
3. **The reconciliation, and the publishable point: Gatti's rule encodes a dosing policy.** It does
   reproduce their own finding in this model (selected patients' median CAZ:AVI 5.14:1 vs 4.04:1), but
   under this project's renally-adjusted grid those patients need no assay — inference is already
   correct in 91.6% of them vs 92.7% of the rest. Hold the dose **fixed** at 2.5 g q8h, which is what
   85% of Gatti's own cohort started on, and the rule comes alive: 88.7% vs 97.9% at ρ = 0.94, and
   81.4% vs 97.1% at ρ = 0.703 — a 9–16 point separation where there had been ~1.

**Neither rule is wrong.** Gatti's encodes a real mechanism that holds when the dose does not
compensate for renal clearance and largely stops holding once it does. **A clinically-derived triage
rule validated in one dosing context does not automatically transfer to another** — demonstrated with
two independent rules rather than asserted. That is a genuine, citable contribution and it cost no new
data.

A prediction written into the script before running it — that Gatti's patients would sit at
conditional attainment probabilities near *zero* — was **wrong in direction** (they sit near one under
adjusted dosing). The docstring records the error rather than hiding it.

**One caution before that comparison is attempted:** the paper's ratio spread is *not* evidence against
the ρ = 0.94 anchor. Its IQR (3.93:1–5.70:1) implies ρ ≈ 0.92 or higher once measurement error is
allowed for — see the arithmetic in that folder's README. The min–max range is a 188-measurement
extreme and should not be read as typical variability.

---

## R3 — The index mismatch nobody has addressed ⚠️ *checked — GAP CONFIRMED, not bridgeable with the published literature*

**The idea.** The regulatory avibactam target is a *time-above-threshold* index — **50% fT > 1 mg/L** —
derived from murine infection models. Continuous-infusion practice uses **Css ≥ 4 mg/L** — a
steady-state concentration.

**Under continuous infusion, %fT > C_T is either 0% or 100%.** The index the target was defined on
degenerates entirely. So the two targets are not two values of one quantity; they are two different
quantities, and no published work states what steady-state concentration is pharmacodynamically
equivalent to a given %fT > 1 mg/L under intermittent dosing.

**⚠️ A correction, then a correction of the correction — both recorded, because the second one matters.**

On 12 August 2026 this entry was edited to say: *"This entry originally quoted the target as a flat
'50% fT > 1 mg/L.' No source anywhere in this project, and none of the papers read to check this route,
actually states that number."* **That edit was wrong, and has been reverted above.** It was based on
reading Berkhout 2015/2016 — the primary murine dose-fractionation study — which indeed reports only
strain- and site-specific stasis values (lung 0–21.4%, thigh 14.1–62.5%) and never a flat 50%. The
error was generalising from the animal study without having read the *registrational* paper that
actually sets the target.

**Das et al. 2019** (*Antimicrob Agents Chemother* 63(4):e02187-18,
doi:[10.1128/AAC.02187-18](https://doi.org/10.1128/AAC.02187-18) — now archived at
`data_external/Das2019_dose_selection/`) states it explicitly: *"50% fT>CT of 1 mg/liter was considered
a robust avibactam target for use in dose selection ... the joint PK/PD target for dosage selection was
defined as ceftazidime 50% fT > 8 mg/liter and avibactam 50% fT > 1 mg/liter."* So the original
wording was correct.

**The distinction that resolves it, and that R3 should actually be making:**

| | What it is | Value |
|---|---|---|
| Berkhout 2015/2016 | *empirical* murine stasis thresholds, strain- and site-specific | lung 0–21.4%; thigh 14.1–62.5% |
| Das 2019 | the *regulatory dose-selection target* — a deliberately conservative round value | **50% fT > 1 mg/L** |

Neither refutes the other. The 50% is a **regulatory choice**, not a measurement: Das rounds up to
cover both organism groups and to align with ceftazidime's own 50% fT>MIC convention. That the
headline target is a rounded administrative convention rather than an observed threshold *strengthens*
R3's argument rather than weakening it — but it must be stated accurately, and the empirical range must
be attributed to Berkhout, never the 50% to Berkhout.

**What this route originally proposed** (a bridging calculation through Aubry 2025 and Kroemer 2023)
**turned out not to be executable, and here is why, checked against the actual primary sources rather
than assumed:**

1. **Aubry et al. 2025** (*Antimicrob Agents Chemother* 69(5):e01797-24,
   doi:[10.1128/aac.01797-24](https://doi.org/10.1128/aac.01797-24), PMC12057351) models
   ceftazidime/avibactam's inoculum effect with avibactam held at a **fixed 4 mg/L background**, and
   states explicitly that "it was assumed that avibactam concentration was sufficient to inhibit
   β-lactamases and had no bactericidal effect on its own." Avibactam is not a variable in this
   model's kill equation at all.
2. **Kroemer et al. 2023** (*Microbiol Spectr* 12(1):e0331823,
   doi:[10.1128/spectrum.03318-23](https://doi.org/10.1128/spectrum.03318-23), PMC10783110) — already
   cited in the manuscript as ref [22] — folds avibactam into ceftazidime's potency (a reduction of
   ceftazidime's EC50 "by >99%") under an explicit assumption of "a permanent inhibition of the
   beta-lactamases by avibactam," stating directly that "a concentration-dependent inhibition was not
   implemented in the model." Same problem: no independent avibactam concentration-response term.
3. Searching further (neither of the above was, in hindsight, ever going to work — both are
   combination-synergy studies built around ceftazidime, not avibactam-focused dose-fractionation
   studies) surfaced the actual paper that defined the %fT>C_T avibactam index in the first place:
   **Berkhout et al. 2015/2016** (*Antimicrob Agents Chemother* 60(1):368-75,
   doi:[10.1128/AAC.01269-15](https://doi.org/10.1128/AAC.01269-15), PMC4704241), murine neutropenic
   thigh and lung infection, dose-fractionation of avibactam (q2h/q8h) against a fixed ceftazidime
   background. This is the right paper — but it has two gaps of its own that block the bridging
   calculation: (a) the exposure-response relationship is reported only as scatter plots with *r*²
   values, with **no fitted closed-form Emax/Hill equation and no published EC50/slope parameters**;
   and (b) **every regimen in the study is intermittent (q2h, q8h, or q12h) — continuous infusion was
   never tested.** There is no data point, published or derivable, at which the %fT>C_T index was
   ever evaluated under the exposure regime the manuscript's own primary scenario uses.

**What Berkhout 2015/2016 does establish, precisely:** stasis was reached at %fT>C_T (1 mg/L) ranging
from 0-21.4% across four strains in the lung model (mean 20.1%, range 16.1-23.5%) and 14.1-62.5%
across six strains in the thigh model — the thigh model consistently the more stringent of the two.
**These are empirical, strain- and site-dependent stasis values; the regulatory 50% is a separate,
later, deliberately conservative rounding by Das 2019 (see the correction note above).** Both belong in
the write-up, clearly distinguished. The T1 scenario in `model2_hujam.py` cites Berkhout as the
index-defining study, which is correct — but its note should also name Das 2019 as the source of the
50% itself.

**The full derivation chain, now traced end to end** (from Das 2019, which names its own sources):
Coleman et al. hollow-fibre → minimum C_T of 0.5 mg/L for Enterobacteriaceae; **Berkhout et al.**
neutropenic murine thigh/lung, 7 ceftazidime-resistant *P. aeruginosa* strains → C_T of 1 mg/L is the
best predictor and the index is %fT>C_T rather than Cmax or AUC; **Das 2019** → rounds to 50% fT>C_T
"taking a conservative approach", pairs it with ceftazidime 50% fT>8 mg/L, and uses the joint target
for dose selection and breakpoint setting. Note one discrepancy worth not smoothing over: Das
summarises Berkhout's efficacy range as "20% to 50%", which matches neither of Berkhout's own published
bounds — **do not attribute "20–50%" to Berkhout.**

**The conclusion this route lands on, checked rather than assumed:** the degeneracy of %fT>C_T under
continuous infusion is real and structural — a time-above-threshold index cannot, by construction,
distinguish among concentrations once they clear the threshold for the whole interval. What is new
here is confirming that **no published avibactam pharmacodynamic model bridges that degeneracy**: not
the two combination-synergy models this route originally proposed (avibactam is not their free
variable), and not even the index-defining dose-fractionation study itself (no continuous-infusion
arm, no fitted equation). The manuscript's comparison between the murine %fT>1 mg/L target and its own
Css ≥ 4 mg/L primary scenario is therefore a comparison between two quantities that no existing
published model connects — not a gap this project failed to close, but a gap that does not yet have a
published closing.

**What this route does NOT deliver, and should not be written up as delivering:** a bridging
calculation, a numeric equivalence, or a defended claim about whether 4 mg/L is more or less
appropriate than the murine-derived range. Building one would require new dose-fractionation data
with a continuous-infusion arm — an experimental study, not a literature synthesis — and is out of
scope here.

**Cost:** moderate, as expected — most of it was verification, not modelling. **Risk:** the risk that
materialised was not "the answer is model-dependent" (the originally anticipated risk) but "no model
answers the question at all," which turned out to be the more defensible and more citable finding.

---

## R4 — Evidence-derived target distribution from hollow-fibre data ✅ *done*

**The idea.** Model 2's target distributions T3 to T6 are analyst-specified, and labelled as such.
One subset of the evidence is genuinely homogeneous in kind: Coleman et al. 2014 reports a
**hollow-fibre regrowth threshold** for avibactam. Repeated measurements of the same quantity in
different organisms are exactly what an evidence-derived (as opposed to analyst-chosen) scenario is
for.

**Correction before implementation.** This entry originally said the threshold rested on
"approximately 0.15 to 0.28 mg/L across eight strains," estimated by random-effects synthesis. Both
parts were wrong, caught only once the primary source (Coleman K et al., *Antimicrob Agents
Chemother* 2014;58(6):3366-72, doi:[10.1128/AAC.00080-14](https://doi.org/10.1128/AAC.00080-14),
PMC4068505) was actually read rather than relied on from an earlier paraphrase. The **eight** strains
in that paper belong to a *different*, single-dose killing experiment. The regrowth-threshold
estimate itself (Table 2) comes from a separate continuous-infusion experiment using only **three**
strains — *E. cloacae* 293HT96 (CT ≤0.15 mg/L), *K. pneumoniae* 283CF5 (CT ≤0.22 mg/L), *K.
pneumoniae* Tunisie K4 (CT ≤0.28 mg/L), plus a fourth Table 2 entry that re-estimates the same *E.
cloacae* strain at a later timepoint (~0.2 mg/L) rather than a fourth independent strain. Three points
— one of them not fully independent of a fourth — cannot support a random-effects variance estimate;
that framing has been dropped.

**Done — `code/model2_hujam.py`, scenario `T7_coleman_evidence`.** Implemented instead as the plain
empirical distribution over the three actually-measured strains ({0.15, 0.22, 0.28} mg/L, equal
weight) — no smoother and no more precise than three real data points justify. Full results,
including the correction above, are in `MODEL2_REPORT.md` §3.8. Headline finding: because T7's
evidence-derived target sits far below the 4 mg/L the manuscript's fixed decision rule assumes, that
rule picks a more aggressive regimen than the T7-optimal one in 93.0-98.6% of draws across the four
non-trivial renal classes, at a modest cost (3.66-5.11 pp expected regret, low EVPI) because the two
choices are adjacent dose steps. `code/test_model2.py` gained 4 checks for the new sampler (44 total,
up from 40).

**What it must not do, and still does not do:** be presented as a distribution over the *clinical*
target. It is a distribution over an in-vitro regrowth threshold, and the gap between that and a
clinical target is the whole problem — T7 anchors one scenario to real measured numbers instead of an
analyst's approximation of them (T3 previously stood in with a placeholder of 0.5 mg/L, noticeably
higher than any of the three real values), but it does not resolve which target is clinically correct.

**Cost:** low, as expected. **Risk:** low — the one risk that materialised (an inaccurate source
citation) was caught before anything was built on it, not after.

---

## R5 — Invert the clinical outcome data ⏸

**The idea.** Gatti 2025 (*AAC* 69(7):e00488-25, n = 218) is the only clinical exposure-response
signal for the joint target: microbiological failure odds ratio 0.03, resistance odds ratio 0.07 for
patients attaining it. In principle the threshold that best separates success from failure could be
estimated rather than assumed.

**Premise verified, 12 August 2026.** The paper is now archived and transcribed at
`data_external/Gatti2025_outcome_R5/` (CC BY 4.0 — fully open, no push restriction). Every figure this
entry cited was checked against the source and is correct: n = 218 (116 pre + 102 post), microbiological
failure OR 0.03 (95% CI 0.005–0.20), 90-day resistance OR 0.07 (95% CI 0.01–0.69). Clinical cure adds a
third: OR 0.08 (95% CI 0.02–0.34).

**One thing this entry should have said and did not.** The study is a **pre-post quasi-experimental
design in which the intervention bundled several changes at once** — continuous infusion rose 31.9% →
96.1%, combination therapy fell 67.2% → 15.7%, treatment duration shortened 14 → 10 days, all alongside
TDM-guided dosing. The odds ratios are therefore associations *within the post-intervention phase*, not
causal estimates of what attaining the target does. Any inversion built on them inherits that
limitation, and R5 must not be written up as though the ORs established causation.

**Still blocked on data** — inverting the threshold needs individual attainment and outcome data, which
this paper does not publish. That would be a *second* request to the same Bologna group while the first
(sent 11 August 2026) is outstanding. **The standing rule is unchanged: no second request without
explicit instruction.**

**What is possible right now, without new data:** the published ORs and the target definition
(fCssCAZ/MIC > 4 AND fCssAVI/CT > 1, CT = 4 mg/L) are archived, so the bound on how much a better
threshold could achieve can be computed from them whenever wanted.

**Cost:** low without data, high with. **Risk:** the request could sour the first one. **Hold.**

---

## R6 — Release a reusable software artifact ✅ *done*

**The idea.** An open, documented, tested package implementing joint two-analyte target attainment
with correlated clearance, target uncertainty and value of information. Nothing like it exists.

CPT:PSP and similar journals treat a usable artifact as a contribution in itself, and it converts
the reproducibility work already done into something other groups can apply to their own drug pair —
which multiplies the reach of R1.

**Cost:** moderate, and mostly packaging rather than science. Much of it is already done: the engine
is verified against the frozen outputs, there are 119 passing tests, and the release package now
runs from a clean checkout.
**Risk:** low.

**Done.** `pyproject.toml` + `code/__init__.py` install the analysis under the import name `hujam`
(`pip install -e .`); verified working both as an installed package from an unrelated directory and
as direct scripts, which is how it was originally developed. `code/interface.py` writes down the
actual two-key contract the decision layer depends on (`joint_pta`, `exceedance`) and provides a
runtime conformance checker — discovered by reading `model2_hujam.py`'s `run()` and `utility()` line
by line, not by guessing what "generic" should mean. `code/test_model2.py` adds 40 checks: engine
conformance, the zero-uncertainty regression, all six target samplers against their documented
support and moments, both correlation samplers, and the decision identities (regret ≥ 0 and zero for
the best-in-hindsight choice; EVPI ≥ EVPPI ≥ 0). Full documentation, including exactly what a
different drug pair would need to write versus reuse unchanged, is in `SOFTWARE.md`.

**Writing the sampler tests found a real bug**, which is exactly what packaging-with-tests is for:
the primary (C1) correlation scenario sampled Normal(0.94, RSE) directly on the correlation scale and
clipped to a valid range, clipping 39.5% of draws to the boundary and pulling the sampled mean down
to 0.877 instead of 0.94. Fixed by sampling on the Fisher-z scale, matching the approach the C2
scenario already used. Every Model 2 result that depended on this sampler was rerun — see
`MODEL2_REPORT.md` for the corrected numbers.

---

## R7 — Ask what the model would have to be wrong about ✅ *done*

**The idea.** Instead of asking how sensitive the conclusion is to each parameter, ask the inverse:
**how large would each error have to be to overturn the conclusion?** For each input, find the value
at which the recommendation flips, and state whether that value is plausible.

"Ceftazidime clearance would have to be wrong by 60% — more than the difference between any two
published models — before the recommendation changes" is a far stronger sentence than a tornado
plot, and it is computable from the machinery already built.

**Cost:** low. **Risk:** none.

**Done — and it returned more than expected.** See `MODEL2_REPORT.md` §3.7. Two results stand out.
The regimen choice is unmovable by the clearance correlation, the ceftazidime unbound fraction or the
ceftazidime variability across a 0.30x to 4x range. And a 22% error in ceftazidime clearance — which
would flip the choice — is 3.9 published standard errors but only 0.59 of the between-study standard
deviation, so **the difference between two published models is more than enough to change the
recommendation while parameter uncertainty within a model is not.** For the monitoring decision there
is no breaking point at all: measuring beats inferring at every correlation from 0.30 to 0.99 and
every assay CV from 0% to 30%.

---

## Recommended order

1. ~~**R7**~~ — ✅ done, and it strengthened the case more than expected.
2. ~~**R2**~~ — ✅ done. Selective measurement beats measuring everyone outright at ρ = 0.94.
3. ~~**R1**~~ — ✅ checked, partially confirmed, claim reworded.
4. ~~**R6**~~ — ✅ done, and found + fixed a real sampler bug along the way.
5. ~~**R4**~~ — ✅ done, and corrected a real citation error (three strains, not eight) along the way.
6. ~~**R3**~~ — ✅ checked. The proposed bridging papers don't support the calculation; confirmed no
   published avibactam model does, including the one that actually defines the index. Written up as
   a literature gap, not a simulated equivalence.
7. **R5** — hold until the Bologna request resolves. Nothing else is queued behind it.

## Honest ceiling, restated

A simulation study built on published parameters will not reach "10/10" on methodology alone. The
two things that genuinely change the ceiling are **new data** (the Bologna request, sent 11 August
2026) and **class-wide generality** (R1). Everything else on this list is craft — it makes the paper
better and harder to dismiss, but it does not change what kind of paper it is.
