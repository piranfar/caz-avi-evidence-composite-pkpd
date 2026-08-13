# Provenance record — Das 2019: where the regulatory target actually comes from

## Source

Das S, Li J, Riccobene T, Carrothers TJ, Newell P, Melnick D, Critchley IA, Stone GG, Nichols WW.
*Dose Selection and Validation for Ceftazidime-Avibactam in Adults with Complicated Intra-abdominal
Infections, Complicated Urinary Tract Infections, and Nosocomial Pneumonia.*
**Antimicrob Agents Chemother. 2019;63(4):e02187-18.** doi:10.1128/AAC.02187-18 · PMID 30670413 ·
PMC6437548. Minireview. AstraZeneca / Allergan authors.

## Legal basis for use

"Copyright © 2019 American Society for Microbiology." **Not openly licensed.** Note that PMC holds only
the abstract for this article — the full text was read from the publisher site (journals.asm.org) on
12 August 2026. Only extracted numeric facts are recorded here; no PDF is archived and none should be
redistributed.

## Why this paper matters: it corrects a correction I made earlier

**This paper is the registrational dose-selection document — the one that defined the licensed
regimen. It states the avibactam target explicitly, and derives it.** Reading it shows that a
"correction" recorded in `NOVELTY_ROUTES.md` R3 on 12 August 2026 was itself partly wrong.

R3's correction had asserted: *"This entry originally quoted the target as a flat '50% fT > 1 mg/L.'
No source anywhere in this project, and none of the papers read to check this route, actually states
that number."* That was based on reading Berkhout 2015/2016 (the primary murine dose-fractionation
study), which indeed reports only strain- and site-specific stasis values, never a flat 50%.

**But Das 2019 does state it, verbatim:**

> "Based on these studies, and taking a conservative approach to ensure that the CT was appropriate for
> both *Enterobacteriaceae* and *P. aeruginosa*, **50% fT>CT of 1 mg/liter was considered a robust
> avibactam target for use in dose selection**. In combination with the previously established
> ceftazidime target, **the joint PK/PD target for dosage selection was defined as ceftazidime 50% fT >
> 8 mg/liter and avibactam 50% fT > 1 mg/liter**."

So the original `NOVELTY_ROUTES.md` wording ("The regulatory avibactam target is 50% fT > 1 mg/L") was
**correct**, and my correction introduced an error by generalising from the animal study without having
read the registrational paper that sets the target. R3 has been re-corrected accordingly.

**The distinction that both statements need, and that is the actually useful finding:**

| | What it is | Value |
|---|---|---|
| Berkhout 2015/2016 | *empirical* murine stasis thresholds, strain- and site-specific | lung 0–21.4%; thigh 14.1–62.5% (%fT>1 mg/L) |
| Das 2019 | the *regulatory dose-selection target*, a deliberately conservative round value | **50% fT > 1 mg/L** |

These are not in conflict and neither refutes the other. The 50% is a **regulatory choice**, not an
observed threshold — chosen, in the paper's own words, to be conservative across two organism groups
and to be "consistent with the 50% fT>MIC required for efficacy of ceftazidime alone." That framing —
a real, sourced, deliberately-rounded regulatory target rather than an empirical measurement — is what
R3 should have said all along, and now does.

## The full derivation chain, now traced end to end

Das 2019 names its own sources, which closes a chain this project had only partly assembled:

1. **Coleman et al.** (Das ref 63) — hollow-fibre, ceftazidime-resistant Enterobacteriaceae
   (*K. pneumoniae* SHV-5 / CTX-M-15 / KPC-2; *E. cloacae* derepressed AmpC; *C. freundii* stably
   derepressed AmpC), CAZ-AVI MICs ≤0.125–4 mg/L → **minimum CT of 0.5 mg/L** for Enterobacteriaceae.
   This is the same Coleman 2014 study that route R4 built scenario T7 from.
2. **Berkhout et al.** (Das ref 64) — neutropenic murine thigh and lung, **7** ceftazidime-resistant
   *P. aeruginosa* strains (stably derepressed AmpC and/or TEM-24), dose fractionation → **CT of
   1 mg/L was the best predictor of efficacy**, and the index is %fT>CT rather than Cmax or AUC.
   This independently confirms R3's identification of Berkhout as the index-defining study.
3. **Das 2019** → rounds up to **50% fT > CT = 1 mg/L**, pairs it with ceftazidime 50% fT > 8 mg/L,
   and uses the joint target for dose selection and breakpoint setting.

**One numerical discrepancy worth recording rather than smoothing over:** Das summarises Berkhout's
efficacy range as "20% to 50% across the thigh and lung infection models". Berkhout's own published
stasis ranges are 0–21.4% (lung) and 14.1–62.5% (thigh) — which is neither the same lower bound nor
the same upper bound. Das may be quoting a different endpoint (e.g. 1-log kill) or a subset of strains.
**Do not quote "20–50%" as Berkhout's own figure**; quote Berkhout directly if the empirical range is
what is meant, and quote Das only for the regulatory 50%.

## Provenance classification

**DIRECTLY REPORTED — regulatory target definition, derivation, and approved dosing table.**
`Das2019_target_derivation.csv` transcribes the joint target definition verbatim, the full derivation
chain with attributions, and **Table 3** — the approved renal dosage adjustments together with the
joint PTA achieved at MIC 8 mg/L in each band (94.9% to 99.6%).

Note that Table 3 records both the *original* phase-3 protocol regimens and the *approved modified*
ones, which differ in three of five renal bands (moderate: q12h → q8h; severe upper: 1000+250 q24h →
750+187.5 q12h; severe lower: 500+125 q24h → 750+187.5 q24h). Anyone comparing this project's regimen
grid against "the licensed regimen" should be explicit about which of the two they mean.

## What this dataset is, and is not

**It is** the authoritative source for what the licensed target and regimen actually are, and for why.
It is a **minireview of the sponsor's own dose-selection programme**, not new primary data.

**It is not** an independent validation of the target, and not patient-level data. All authors are
AstraZeneca or Allergan employees; the paper reviews the analyses that supported the sponsor's own
approved label. The 50% figure is a conservative regulatory rounding, not an empirical measurement —
which is precisely the point route R3 is making about index degeneracy under continuous infusion, and
this paper does not address continuous infusion at all (every regimen here is a 2-h intermittent
infusion).

## Reuse conditions to honour

1. Cite Das et al. 2019 for the regulatory target and the approved dosing table.
2. **Always distinguish the regulatory 50% (Das) from the empirical murine ranges (Berkhout).** They
   are different kinds of quantity and conflating them is what produced two successive errors here.
3. Do not attribute "20–50%" to Berkhout — it is Das's summary and does not match Berkhout's published
   ranges.
4. Not openly licensed; cite and link, do not redistribute.
