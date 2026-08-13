# `data_external/` — index and honest assessment

Sixteen external sources, each in its own folder with a `README_provenance.md` and an extracted CSV.
Assembled 11–12 August 2026.

**Read this before using any of them.** The single most important fact about this collection is stated
first, because it is easy to lose in the volume:

> **None of these is validation data for the manuscript's primary scenario, and assembling them has not
> changed that.** The gap identified in Phase 2 — no individual patient data, with both analytes
> measured, in non-RRT patients on continuous infusion — is exactly where it was. What this collection
> does is make every number in the project traceable, catch several errors, and materially strengthen
> one novelty claim. It adds no new evidence about patients.

**Relation to `DATA_AVAILABILITY_MATRIX.csv`:** that Phase 2 artefact formally classifies 8 studies
(S1–S8) against 35 evidence fields and remains the authoritative classification for those. This index
covers all 16 and is descriptive, not a replacement. Where both exist, the matrix's A–F classification
governs.

---

## What each source is for

### The anchor and its lineage
| Folder | What it is | Why it matters here |
|---|---|---|
| `Cojutti2024_ANCHOR_PopPK/` | n=112, 185 paired Css, RRT excluded | **THE ANCHOR.** Source of ρ = 0.94 and all six RSEs Model 2 propagates. All seven parameters verified exact against source. |
| `Li2019_registrational_PopPK/` | n=1975/2249, registrational | The reference model everyone benchmarks against; Cojutti fixed its volumes from here. **Never computed the cross-drug correlation** — key R1 evidence. Includes NONMEM control streams + bootstrap 90% CIs. |
| `Das2019_dose_selection/` | registrational dose selection | Source of the licensed regimen and of the "50% fT > 1 mg/L" target, with its full derivation chain. |

### Preclinical primary sources (load-bearing for live model scenarios)
| Folder | What it is | Why it matters here |
|---|---|---|
| `Coleman2014_hollowfibre_T7/` | hollow-fibre, 3 strains | **Model 2's live T7 scenario is built from this.** 3 strains, not 8 — corrected here. |
| `Berkhout2016_murine_index/` | murine dose-fractionation | Defined the %fT>C_T index and C_T = 1 mg/L. Underpins T1 and all of route R3. |

### The project's own central question
| Folder | What it is | Why it matters here |
|---|---|---|
| `Gatti2024_ratio_one_leg/` | n=107, 188 paired TDM | **The empirical test of "can you infer avibactam from ceftazidime?"** Answer: no — ratio 1.29:1 to 13.46:1 vs 4:1 vial. Supplies a clinical triage rule to compare against R2's. |
| `Gatti2023_joint_target_origin/` | n=58, abstract only | Defined the aggressive joint PK/PD target the whole Bologna chain uses. |
| `Gatti2025_outcome_R5/` | n=218, CC BY 4.0 | The clinical outcome paper behind route R5. Pre-post design — associations, never causal. |
| `JAC_exchange_measure_one_or_both/` | 4 letters, 2023–24 | **The published four-round dispute over this project's exact question.** The single-analyte side rests on a correlation it *assumes and never quantifies* — the quantity Model 1 estimates. Also contains the actionability objection neither side answered. |

### Found by the JAC sweep of 12 August 2026 — affects results, not just provenance
| Folder | What it is | Why it matters here |
|---|---|---|
| `Dimelow2018_ELF_concentration_dependence/` | CC BY-NC full text + validated re-derivation | **The 0.52/0.42 ELF ratios are correct quotations of a *non-linear* function evaluated at 15.3 and 2.4 mg/L, applied as constants at up to 104 mg/L.** At the neurotoxicity screen the same model gives 25.8%, not 52%. The lung therapeutic-window conclusion flips between the affected scenarios. |
| `Cojutti2026_FN_avibactam_limiting/` | n=17 CZA of 256, both analytes | **Independent clinical evidence that avibactam is the limiting component** — 6/7 empirical and all targeted failures were avibactam, not ceftazidime. Single-analyte TDM would have missed every one. Also an independent EUCAST basis for C_T = 4 mg/L. |
| `Falcone2021_three_analyte_PopPK/` | n=41, CAZ+AVI+ATM assayed | **Second study with both clearances estimated and the correlation never computed** (after Li 2019) — R1 is now one-reported-in-three. Median C_min ratio 5.31:1, replicating Gatti 2024 independently. Ceftazidime C_min reached 175 mg/L, above the 104 mg/L screen. |

### Clinical cohorts (all class C — assumption-testing only)
| Folder | n | Population | Both analytes? |
|---|---|---|---|
| `Gatti2023_*` (top level) | 8 pts / 17 occasions | CVVHDF | yes |
| `dryad_Li2025_CRRT/` | 21 | CRRT | yes — **raw individual data, CC0** |
| `OJeanson2024_CVVHDF/` | 4 | CVVHDF | yes |
| `Tian2025_CVVH/` | 7 | CVVH | yes |
| `Lanini2024_nonRRT_CrCl/` | 52 | **non-RRT** | yes, but derived free concs |
| `Fresan2023_CI_TDM/` | 31 | mixed, true CI | **NO — ceftazidime only** |
| `Wu2025_PopPK_AKI/` | 31+32 | mixed critically ill | model summary |
| `Chen2025_PopPK_CRKP/` | 45 | crit + non-crit | model summary |

---

## Was it worth it? An honest accounting

### What this genuinely bought

**1. The anchor is now verified at source.** ρ = 0.94 and all six RSEs were cited throughout the
project but had never been checked against Cojutti 2024's Table 2. They are now, and every one matches
exactly — as do the primary model's own constants in `reproduce_primary_run.py`. The foundation is
sound. That is worth knowing with certainty rather than by assumption.

**2. Route R1 got materially stronger.** Li 2019's supplement shows the sponsor **merged both drugs'
random effects into one per-patient file and still never computed their correlation** — resampling
intact patients to preserve it implicitly instead. "They had it one line of code away and did not take
it," evidenced verbatim, is a much sharper claim than "it is not estimated."

**3. Four real errors were caught before they reached the manuscript.**
- R4's "eight strains" (Coleman's threshold rests on three; the eight belong to a different experiment)
- R3's "50% appears in no source" — **my own correction, which was itself wrong**; Das 2019 states it
- The free-fraction split: 0.85/0.92 (registrational lineage) vs 0.90/0.93 (Gatti/Lanini) — no consensus
- Fresan 2023's mortality column disagrees with its own text (7 vs 6), confirmed against the typeset PDF

**4. A quantified, documented conservatism.** Model 2 applies the published RSEs on the ω scale though
they are reported on the CV scale, carrying ~20–25% more variability-parameter uncertainty than
strictly implied. Conservative in direction, now documented with exact magnitude rather than left as an
open question in the register.

**5. One genuinely novel check is now possible.** Gatti 2024 supplies a clinically-derived triage rule
(CrCL > 75–78, urea ≤ 45–51). Comparing it against R2's model-derived concentration windows — two
independent routes to "which patients need the second assay" — has not been done and would be
publishable either way it comes out.

**6. A usable parameter-uncertainty package.** Li 2019's bootstrap 90% CIs for 31 parameters, which
`MODEL_DEVELOPMENT_DECISION.md` earmarked for Candidate C and which were previously unobtained.

### What it did NOT buy — stated plainly

**1. No new data about patients.** All sixteen are published aggregate or summary results. Not one is
a new observation.

**2. The validation gap is completely unchanged.** Model 1 still cannot be validated. The closest
population match (Lanini 2024, non-RRT) is cross-sectional with derived rather than measured free
concentrations. The one dataset with true continuous infusion and individual patients (Fresan 2023)
**never measured avibactam at all**. No source here provides paired individual clearance estimates for
both analytes in non-RRT patients — which is exactly what Model 1 would need.

**3. The correlation question is still empirically thin.** There remains exactly **one** published
CAZ–AVI clearance correlation (Cojutti's 0.94). Model 1's 0.703 comes from a CRRT cohort and is not
comparable. Knowing that Li 2019 declined to compute one strengthens the *argument* about the gap; it
does not fill it.

**4. Most of the output is provenance, not science.** These folders make the project defensible under
review and traceable end to end. That is real value for a manuscript facing reviewers, and it should
not be mistaken for new findings.

> **Superseded in part on 12 August 2026.** The sentence originally here — "They do not change a single
> number in any result" — was true when written and is no longer. The JAC sweep found that the ELF
> penetration ratios are concentration-specific values applied as constants
> (`Dimelow2018_ELF_concentration_dependence/`), and the lung therapeutic-window conclusion flips
> between the scenarios this affects. That is a result-level finding, not a provenance one. It is the
> first time this archiving exercise changed an answer rather than a citation.

**5. The manuscript is still untouched.** Phase 6 has not started.

### The honest bottom line

This was **worth doing, but it is insurance rather than discovery.** It converts a project that cited
its sources into one that has verified them, and it caught four errors — including one I introduced
myself. The single most valuable output is probably the R1 finding, which is a genuine strengthening of
the paper's novelty claim. The single most important thing it did *not* do is close the data gap, and
no amount of further literature archiving will: that needs either the Bologna reply, a second request,
or new measurements.

**If the aim is a stronger paper, the highest-value remaining moves are the R2-vs-Gatti-2024 triage
comparison (novel, computable now, no new data needed) and Phase 6.** Further paper-archiving has
reached diminishing returns.

### Revised after the JAC sweep, 12 August 2026

The judgement above — that archiving had hit diminishing returns — was wrong, and the sweep that
tested it is what showed so. A systematic pass over all 49 JAC avibactam PK/PD papers returned three
things no earlier round did:

1. **A four-round published dispute** on exactly this project's question, of which the project had
   archived only the two endpoints without knowing they were opposing sides. It supplies both a
   citation that the question is contested and an unanswered objection the project can answer.
2. **A structural limitation in Model 2** (MICi coupling, §5.1 of `MODEL2_REPORT.md`) — against the
   project's interest, and better self-reported than reviewer-discovered.
3. **A result-level error in the ELF scenarios** (§5.2) that flips a stated conclusion.

The lesson is narrower than "keep archiving": ad-hoc retrieval had indeed exhausted itself, but a
**systematic, journal-complete sweep** is a different instrument and found what ad-hoc retrieval could
not. Doing the same for *AAC*, *CID*, *Clin Pharmacokinet* and *IJAA* is now the obvious next move, and
it should happen before Phase 6 rather than after — §5.2 shows the sweep can change what the text needs
to say.
