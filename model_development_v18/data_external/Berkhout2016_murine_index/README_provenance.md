# Provenance record — Berkhout 2016: the study that defined the avibactam PK/PD index

## Source

Berkhout J, Melchers MJ, van Mil AC, Seyedmousavi S, Lagarde CM, Schuck VJ, Nichols WW, Mouton JW.
*Pharmacodynamics of Ceftazidime and Avibactam in Neutropenic Mice with Thigh or Lung Infection.*
**Antimicrob Agents Chemother. 2016;60(1):368-375.** (published online 2 November 2015)
doi:10.1128/AAC.01269-15 · PMID 26525790 · PMC4704241.
Radboud UMC / Canisius-Wilhelmina Hospital / Erasmus MC + AstraZeneca.

## Why this folder exists

This is **the study that established the avibactam PK/PD index itself** — that avibactam efficacy
tracks *time above a threshold concentration* (%fT>C_T) rather than Cmax or AUC, and that the
threshold is 1 mg/L. Everything downstream depends on it:

- **Model 2's `T1_point_1` scenario** (1 mg/L, "regulatory target") — its inline note in
  `model2_hujam.py` cites this paper.
- **Route R3** in its entirety — the index-degeneracy argument is about *this* index.
- **Das 2019's regulatory target** — Das names this study (its ref 64) as the source of C_T = 1 mg/L.

Like Coleman 2014, it was cited throughout the project but never archived. This folder closes that gap.

## Legal basis for use

American Society for Microbiology, 2016. Read from the PMC page (PMC4704241) on 12 August 2026.
**Not openly licensed.** Only extracted numeric facts are recorded here; no PDF is archived and none
should be redistributed.

## Provenance classification

**DIRECTLY REPORTED — published murine dose-fractionation results.**
`Berkhout2016_fT_over_CT.csv` records the index finding, the C_T value, the per-model stasis ranges,
the q2h-vs-q8h comparison that demonstrates the index, the dosing intervals actually tested, and two
explicit *absences* that route R3 turns on.

## The two absences that make R3 what it is

Both are properties of this paper, and both were verified by reading it rather than assumed:

1. **No fitted equation.** The exposure-response relationship is reported only as scatter plots with
   *r*² values. There is **no closed-form Emax/Hill model, no EC50, no slope parameter** — so the
   relationship cannot be re-simulated or inverted by anyone downstream.
2. **No continuous infusion.** Every regimen in the study is intermittent — q2h, q8h or q12h. The
   index was never evaluated under continuous infusion, which is the exposure mode the manuscript's own
   primary scenario uses.

Together these are why R3 concludes that no published avibactam PD model can bridge %fT>C_T to a
steady-state concentration: the index-defining study itself provides neither the equation nor the
exposure regime that a bridge would need.

## The distinction that took two corrections to get right

**Berkhout reports empirical, strain- and site-specific stasis thresholds. It does not report 50%.**

| Model | Strains | %fT>C_T (1 mg/L) for stasis |
|---|---|---|
| Lung | 4 *P. aeruginosa* | 0–21.4% (mean 20.1%, range 16.1–23.5%) |
| Thigh | 6 *P. aeruginosa* | 14.1–62.5% — consistently more stringent |

The familiar **"50% fT > 1 mg/L" is Das 2019's regulatory rounding**, not a Berkhout result — Das
rounds up to cover both organism groups and to align with ceftazidime's own 50% fT>MIC convention.
See `data_external/Das2019_dose_selection/`.

`NOVELTY_ROUTES.md` R3 recorded this wrongly in both directions on 12 August 2026 before settling: it
first claimed the 50% appeared in no source (wrong — Das states it), after an edit that had been
prompted by reading only this paper. Both the empirical range and the regulatory value are real; they
are different kinds of quantity. **Attribute the range to Berkhout and the 50% to Das, never the
reverse.**

**A numerical discrepancy to preserve, not smooth over:** Das summarises this study's efficacy range as
"20% to 50%", which matches **neither** of Berkhout's published bounds. Das may be quoting a different
endpoint (1-log kill rather than stasis) or a strain subset. Quote Berkhout directly for the empirical
range; quote Das only for the regulatory target.

## What this dataset is, and is not

**It is** the empirical foundation of the avibactam target — a proper dose-fractionation study, the
design that can actually identify which PK/PD index drives effect. Its q2h-vs-q8h result is the direct
evidence: a 2.7- and 10.1-fold difference in daily dose for the same static effect, yet nearly the same
%fT>C_T (21.6% vs 18.5%), which is what identifies time-above-threshold as the driving index.

**It is not** a source for a clinical target, a continuous-infusion target, or a usable
concentration-effect model. It is **neutropenic mice**, *P. aeruginosa* only (the Enterobacteriaceae
half of the regulatory target comes from Coleman 2014), with a fixed ceftazidime background chosen so
that ceftazidime alone produced growth rather than stasis.

## Reuse conditions to honour

1. Cite Berkhout et al. 2016 for the index (%fT>C_T), the threshold (1 mg/L), and the empirical ranges.
2. **Never attribute the 50% to this paper** — that is Das 2019.
3. Never attribute "20–50%" to this paper either — that is Das's summary and does not match.
4. State that no fitted equation and no continuous-infusion arm exist, whenever the paper is invoked
   for anything quantitative.
5. Not openly licensed; cite and link, do not redistribute.
