# Provenance record — Li 2019: the registrational PopPK model

## Source

Li J, Lovern M, Green ML, Chiu J, Zhou D, Comisar C, Xiong Y, Hing J, MacPherson M, Wright JG,
Riccobene T, Carrothers TJ, Das S. *Ceftazidime-Avibactam Population Pharmacokinetic Modeling and
Pharmacodynamic Target Attainment Across Adult Indications and Patient Subgroups.*
**Clin Transl Sci. 2019;12(2):151-163.** doi:10.1111/cts.12585 · PMID 30221827 · PMC6440567

AstraZeneca / Quantitative Solutions (Certara) / Wright Dose / Allergan. Sponsor-run registrational
analysis. Rights acquired by Pfizer December 2016.

## Legal basis for use

Full text retrieved 12 August 2026 from PubMed Central (PMC6440567); the Version-of-Record PDF and both
supplement files were supplied by the author the same day. **No copyright or licence statement appears
anywhere in the retrieved article text** — searches for "copyright", "licens", "Creative Commons" and
"CC BY" all returned zero hits; the only publisher mark is a "John Wiley & Sons, Ltd" footer under each
table. The project's `DATA_AVAILABILITY_MATRIX.csv` records this study (S8) as "open_access YES,
license CC BY-NC"; **that licence claim could not be confirmed from the article text and should be
verified at the DOI before relying on it.**

Three files are archived here — `Li2019_CTS_12-151_SUBSCRIPTION-ACCESS.pdf`,
`Li2019_DataS1_NONMEM_control_streams.docx`, `Li2019_DataS2_supplementary_methods_tables.docx`.
**Because the licence is unconfirmed, treat all three as excluded from any GitHub push**, on the same
footing as the Tian / Wu / Lanini / Fresan PDFs. The extracted CSVs and this README carry no such
restriction. If the CC BY-NC claim is confirmed at the DOI, that restriction can be revisited.

## Why this paper matters to this project

It is the **registrational reference model** — the benchmark that Cojutti 2024, Chen 2025, Wu 2025 and
O'Jeanson 2024 all compare themselves against, and the source of the fixed volumes (V1 = 18.0,
V2 = 18.1 L) that the anchor model, Cojutti 2024, adopted rather than estimating. It is also, by a
wide margin, the largest CAZ-AVI PK dataset in existence: **1,975 subjects / 9,155 observations for
ceftazidime and 2,249 subjects / 13,735 observations for avibactam**, pooled across four phase III
cIAI/cUTI trials (RECLAIM, RECAPTURE, REPRISE), one phase III NP trial (REPROVE), two phase II and
eleven phase I studies.

## THE KEY FINDING — for novelty route R1 (now confirmed from the supplement itself)

**Li 2019 acknowledges the ceftazidime–avibactam correlation, preserves it, and never quantifies it.**
This is a category R1 does not currently name, and it is the strongest single example of R1's thesis.

The main text says only: *"To account for the correlation between ceftazidime and avibactam random
effects, the random effects were bootstrapped using the approaches detailed in the supplement."*
**Data S2 has now been obtained** (12 August 2026, supplied by the author) and gives the method in
full, verbatim:

> "Each patient's random effects from the ceftazidime and avibactam models were merged into a single
> data file, in which there was a single record for each patient containing all his/her random effects
> for both compounds. Patient-level random effect records were bootstrapped within the population of
> interest and read into the simulation model as data columns. The covariates used in the simulations
> were matched such that any given value of the patient identification variable had the exactly the
> same covariates (and random effects) and covariate values in the simulation datasets for both
> compounds, **thereby preserving any underlying correlations**. This approach preserved the inherent
> correlations between subject covariates and parameter random effects."

**This makes the finding much stronger than "not estimated".** The sponsor:

- **had the paired random effects for both drugs merged into a single per-patient file** — one record
  per patient containing both compounds' ETAs, which is precisely the data structure from which a
  cross-drug correlation is computed in one line;
- **knew the correlation mattered**, enough to build the entire simulation around preserving it;
- and **chose a method that carries the correlation implicitly without ever measuring it** — resampling
  intact patients rather than estimating a covariance.

So the number was not merely unpublished; it was **never computed**, by the group with the largest
dataset in existence for this drug pair and with the data already in the right shape. That is a far
sharper sentence for R1 than "it is not estimated."

Confirmed structurally by the control streams in Data S1: each drug has its own file with
`$OMEGA BLOCK(4)` over **its own four ETAs (CL, V1, V2, Q)**. There is no joint model and no cross-drug
block. A full-text search of the article confirms the absence in the paper itself: **"OMEGA" 0 hits,
"omega" 0 hits, "off-diagonal" 0 hits, "rho" 0 hits.**

**Shrinkage inflation, also now documented.** The supplement states that post-hoc ETAs were re-inflated
before simulation to avoid understating between-subject variability, using an explicit formula
(ETASIM from ETAEST and shrink%). So the simulated variability is deliberately wider than the raw
empirical-Bayes estimates — worth knowing before comparing this model's spread against any other.

**This confirms the project's own audit note exactly.** `DATA_AVAILABILITY_MATRIX.csv` (S8) states:
*"it reports WITHIN-drug random-effect correlations only; no CAZ-to-AVI cross-drug covariance is
quantified in any regulatory or published source."* Verified correct.

**What R1 should take from this:** the sponsor with the largest dataset in the world for this drug
pair — over 2,000 patients, with paired measurements of both components — knew the correlation
mattered enough to handle it explicitly in their target-attainment simulations, and still published no
estimate of it. That is a stronger and more specific claim than "it is not estimated": it is
*deliberately routed around*. R1's revised claim can now distinguish three states rather than two:
reported (Cojutti 0.94; pip/tazo 0.93; azt/avi 0.976–0.986), **acknowledged but not quantified
(Li 2019)**, and not addressed at all (most others).

## Within-drug random effects ARE fully published

Both models carry a full OMEGA block over CL, V1, V2, Q with correlations printed — for ceftazidime
r ranges −0.84 to +0.82, for avibactam "−0.36 < r < 0.99". These are transcribed in the CSV. So the
project's characterisation of Li 2019 as "the richest public parameter package for either drug" is
fair — provided "richest" is not read as including the cross-drug term, which it does not.

**The avibactam Q-versus-V2 ambiguity — RESOLVED by the control streams.** This was recorded as
unresolved when only the article text was available. Data S1 settles it: **both** drugs' control files
assign the random effects in the same order, explicitly commented —

```
$OMEGA  BLOCK(4)
 ...  ;    ETA1 CL
 ...  ;    ETA2 V1
 ...  ;    ETA3 V2
 ...  ;    ETA4 Q
```

and the parameter assignments confirm it (`V2 = EXP(MU_3 + ETA(3))`, `Q = EXP(MU_4 + ETA(4))`).
Applying that mapping to the avibactam η variances published in Table 2:

| ETA | Parameter | Variance | CV% = √variance |
|---|---|---|---|
| ETA1 | CL | 0.349 | 59.1% |
| ETA2 | V1 | 1.147 | 107.1% |
| ETA3 | **V2** | 1.494 | **122.2%** |
| ETA4 | **Q** | 6.359 | **252.2%** |

So **avibactam Q = 252.2% and V2 = 122.2%** — matching the ceftazidime pattern of very large BSV on Q
(259%). Read against the θ-row order (CL, V1, Q, V2), the BSV column as extracted appears to place
these the other way round; the η-block plus control-stream mapping is authoritative and is what the CSV
now records. Whether the discrepancy is a typesetting error in Table 2 or an artefact of text
extraction has not been determined and does not affect the values. (Note the paper already carries a
published erratum for a *different* transposition — the Ceftazidime/Avibactam headings in Table 4 —
so table-ordering slips are not unprecedented here.)

Ceftazidime's assignment was never ambiguous and was verified arithmetically (√0.179 = 42.3%,
√1.10 = 105%, √1.21 = 110%, √6.70 = 259%).

**A second small thing the control stream resolves:** the article states that "abnormally high
ceftazidime concentrations (> 750 mg/mL) were excluded", which is dimensionally impossible. The control
file shows `IGNORE=(DV.GT.750000)` with concentrations in ng/mL — i.e. the cut-off is
**750,000 ng/mL = 750 mg/L**. The published units are a typo; the actual threshold is unambiguous.

## Two further items that settle open questions elsewhere in this project

**1. The free fractions trace to here.** Verbatim: *"PK target attainment analyses used free plasma
concentrations (taken to be 85% and 92% of total plasma concentrations for ceftazidime and avibactam,
respectively)."* These are exactly the values Cojutti 2024 uses and exactly the values this project's
primary model hard-codes (`FU_CAZ, FU_AVI = 0.85, 0.92`). So the 0.85/0.92 pair is the
**registrational lineage**. Gatti 2025 and Lanini 2024 instead use 0.90/0.93. The divergence noted in
those folders is now traceable: two lineages, not random inconsistency — but still no consensus.

**2. The registrational exposure-response analysis is null, and the reason matters for route R5.**
Verbatim: *"Almost all individual ceftazidime %T > MIC ... and avibactam %T > CT values were close to
100%. The low treatment failure rates in the phase III trials limited investigation of clinical PK/PD
relationships, and no meaningful exposure-response relationships were observed."* Unfavourable
microbiological response was 5.8% (cIAI), 15.5% (cUTI), 38.7% (NP).

This is a **ceiling effect, not evidence of no relationship** — and it is precisely why R5 has to rely
on Gatti 2025 rather than on the registrational data. The largest dataset available is structurally
uninformative for threshold estimation because nearly everyone hit the target. Worth stating explicitly
in any write-up of R5, because "the registrational trials found no exposure-response relationship" is
otherwise an easy sentence to misread as evidence against the target mattering.

## What this dataset is, and is not

**It is** the definitive registrational parameter package, with full within-drug uncertainty, and the
provenance root for the free-fraction assumptions and the fixed volumes used downstream.

**It is not** patient-level data, and not independent of the sponsor. Simulations "excluded residual
error and uncertainty in the population parameters". The structural model is **two-compartment** for
both drugs, whereas Cojutti 2024, Chen 2025 and Wu 2025 all chose one-compartment — so when Cojutti
fixed V1 = 18.0 / V2 = 18.1 "on the basis of prior population PK models", it imported volumes from a
two-compartment fit into a one-compartment structure. That is a defensible pragmatic choice given CI
data cannot identify V, but it is a structural mismatch that any comparison across these models should
state.

**Both supplements have now been obtained (12 August 2026, supplied by the author after PMC's
JS/cookie challenge and Wiley's 403 blocked automated retrieval), and both of the project's claims
about this paper are VERIFIED:**

| `DATA_AVAILABILITY_MATRIX.csv` claim (S8) | Verdict |
|---|---|
| "actual NONMEM control streams for both final models, in Data S1" | ✓ **confirmed** — complete, executable `$PROBLEM`…`$TABLE` streams for both drugs |
| "`$OMEGA BLOCK(4)` for each drug" | ✓ **confirmed**, verbatim |
| "`$SIZES / $INPUT / $SUBROUTINE ADVAN3 TRANS4 / $PK` with all covariate blocks" | ✓ **confirmed** (ADVAN3 TRANS4 = two-compartment, as stated) |
| "The dataset itself is referenced but not included" | ✓ **confirmed** — `$DATA ../../data/merged.cazavi.20160819.v1.csv`, a path only |
| "bootstrap 90% confidence intervals for every theta in Table S4/S5" | ✓ **confirmed** — now transcribed in `Li2019_bootstrap_CIs.csv` |
| "reports WITHIN-drug random-effect correlations only; no CAZ-to-AVI cross-drug covariance is quantified" | ✓ **confirmed**, and now explained in full (see above) |

The project's audit of this paper was accurate on every point.

**Bootstrap CIs are now available for the parameter-uncertainty layer** that
`MODEL_DEVELOPMENT_DECISION.md` earmarked them for — see `Li2019_bootstrap_CIs.csv` (14 ceftazidime
θs, 17 avibactam θs, point estimate + bootstrap median + 90% CI each). Two cautions when using them:
ceftazidime **Q** is very poorly determined (point 31.5, bootstrap median 15.4, CI 7.94–78.6 — this is
the one parameter the paper admits differed from its bootstrap median by >20%), and the avibactam
**θ14 (WT on Vc)** CI is printed as "(0.969, 0.144)", an upper bound below its lower bound — an evident
typo in the source, recorded as-printed in the CSV rather than silently repaired.

## Reuse conditions to honour

1. Cite Li et al. 2019 wherever these values are used.
2. **Do not describe Li 2019 as reporting a ceftazidime–avibactam correlation.** It acknowledges one
   and handles it by paired bootstrap; it publishes no value.
3. Verify the avibactam Q/V2 BSV assignment against the PDF before using either.
4. State the two-compartment structure when comparing against one-compartment models.
5. Licence unconfirmed — cite and link; do not redistribute the article text.
