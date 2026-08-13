# Provenance record — Coleman 2014: the primary source behind Model 2's T7 scenario

## Source

Coleman K, Levasseur P, Girard AM, Borgonovi M, Miossec C, Merdjan H, Drusano G, Shlaes D,
Nichols WW. *Activities of ceftazidime and avibactam against β-lactamase-producing Enterobacteriaceae
in a hollow-fiber pharmacodynamic model.*
**Antimicrob Agents Chemother. 2014;58(6):3366-3372.** doi:10.1128/AAC.00080-14 · PMID 24687507 ·
PMC4068505. Novexel SA (Romainville, France) + Drusano + Shlaes + AstraZeneca.

Cited in the manuscript as reference **[20]**.

## Why this folder exists

**`model2_hujam.py`'s `T7_coleman_evidence` scenario — a live scenario in the running model, used in
all 21 scenario combinations — is built directly from this paper's Table 2.** Until now the paper was
cited in the code and in `NOVELTY_ROUTES.md` R4 but was **not archived in `data_external/`**, unlike
every downstream paper that depends on it. That was a provenance gap in a load-bearing input, and this
folder closes it.

## Legal basis for use

"Not Open Access" — confirmed via the NCBI OA service, which returns
`<error code="idIsNotOpenAccess">identifier 'PMC4068505' is not Open Access</error>`. The article is
readable on PMC but is not openly licensed. Table 2 and the abstract were read from the PMC page on
12 August 2026. Only extracted numeric facts are recorded here; **no PDF is archived and none should
be redistributed.**

## Provenance classification

**DIRECTLY REPORTED — published per-strain in-vitro thresholds.**
`Coleman2014_regrowth_thresholds.csv` records the three independent strain values that T7 uses, the
fourth same-strain re-estimate that T7 deliberately excludes, the pooled abstract figure, the pulsed-
exposure result, and — importantly — the separate eight-strain single-dose experiment that is *not*
the source of the threshold.

## The error this paper corrected, recorded so it cannot recur

R4's original entry described the threshold as "approximately 0.15 to 0.28 mg/L across **eight**
strains" and proposed a random-effects synthesis over them. Reading Table 2 directly showed:

- The **range 0.15–0.28 mg/L is correct.**
- The **strain count is not.** Those values come from a continuous-infusion experiment using **three**
  strains, plus a fourth row that re-measures one of them at a later timepoint.
- The **eight strains belong to a different experiment entirely** — a single-dose killing study in
  which a 1 g/250 mg profile sufficed for 7 of 8 strains and 2 g/500 mg was needed for a high-level
  AmpC producer. That experiment did not estimate a regrowth threshold at all.

Three points, one of them not fully independent, cannot support a between-strain variance estimate, so
the random-effects framing was dropped and T7 was implemented as the plain empirical distribution over
the three measured strains.

**Where the mistake most likely came from — now traceable.** Li 2019 (archived at
`data_external/Li2019_registrational_PopPK/`) summarises this paper as *"CT values of 0.15–0.28 mg/L
were sufficient to restore ceftazidime activity"* — quoting the range correctly but **without stating
a strain count**. The "eight strains" was almost certainly imported from Coleman's *other* experiment
and welded onto the range. This is a good illustration of why the project reads primary sources rather
than paraphrases: the number survived the paraphrase, the provenance did not.

## A downstream discrepancy worth knowing about

The two registrational papers summarise this same study differently, and both are defensible:

| Source | What it says Coleman established | Which experiment |
|---|---|---|
| **Li 2019** | CT values of **0.15–0.28 mg/L** restore ceftazidime activity | the continuous-infusion regrowth experiment (Table 2) |
| **Das 2019** | a **minimum CT of 0.5 mg/L** is appropriate for Enterobacteriaceae | the conservative bound, consistent with the pulsed-exposure result (>0.25 and <0.5 mg/L) |

So the "avibactam threshold" attributed to Coleman is anywhere from 0.15 to 0.5 mg/L depending on
which experiment and which reading is meant — a three-fold spread inside a single paper. **T7 uses the
0.15/0.22/0.28 regrowth values specifically**, and any write-up should say so rather than referring
vaguely to "Coleman's threshold".

Das 2019 also supplies the strain characterisation Coleman's own abstract does not: *K. pneumoniae*
producing SHV-5, CTX-M-15 or KPC-2; an *E. cloacae* with derepressed AmpC; a *C. freundii* with stably
derepressed AmpC; CAZ-AVI MICs ≤0.125 to 4 mg/L.

## What this dataset is, and is not

**It is** the empirical foundation of T7 — the only scenario in Model 2 built from measured data rather
than analyst specification.

**It is not** a clinical target, and T7 does not claim it is. These are **in-vitro hollow-fibre regrowth
thresholds**: the avibactam concentration below which bacterial regrowth resumed while ceftazidime was
infused continuously. The gap between that and a clinical target is the substance of route R3. All
three values are also **upper bounds** ("≤", extrapolated from exponential-decline curves at the last
pre-regrowth sampling point), so using them as point values is itself a conservative simplification.

Sample size is three strains. Nothing fitted to them can support a variance estimate, and T7 does not
attempt one.

## Reuse conditions to honour

1. Cite Coleman et al. 2014 wherever these values are used.
2. **Never describe the 0.15–0.28 range as resting on eight strains.** It rests on three.
3. Say which experiment is meant — the 0.15–0.28 regrowth values, or the ≤0.5 pulsed/conservative bound.
4. State that the values are upper bounds and in-vitro, not clinical.
5. Not openly licensed; cite and link, do not redistribute.
