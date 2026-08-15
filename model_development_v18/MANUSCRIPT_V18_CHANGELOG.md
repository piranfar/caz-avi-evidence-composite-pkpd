# v17 → v18: what changed before resubmission to IJAA

Applied 2026-08-13 to produce `Piranfar_CAZ-AVI_IJAA_Original_Article_v18.docx`.
`v17.docx` is untouched. Word count 3,833 → 4,221.

The document was edited at the XML level, validated against the OOXML schema
(`validate.py` — all checks passed), and every renumbered citation was verified to still
point at the same reference. Note: a visual render was not possible on this machine
(no LibreOffice), so **open the file in Word once before submitting.**

## Scientific changes

### 1. The clearance-correlation assumption is now disclosed (§2.1, §3.4, §4.3, abstract)

The headline diagnostic figures — PPV 95.8%, NPV 83.6% — are computed at a clearance
correlation of 0.94, taken from the source model. `MODEL1_REPORT.md` estimates that
correlation at **0.703 (95% CI 0.381–0.873)** from individual patient data, and **0.94 lies
outside that interval**. v17 hedged this in a single unquantified sentence and said nothing
in the abstract.

- **§2.1** — the sentence "No individual patient data were used, and we did not fit a new
  population pharmacokinetic model" is replaced. It was both the sentence that matched the
  stated bioRxiv/medRxiv rejection reason and, since the joint NLME fit, factually false.
  The replacement states the fit, the estimate, the confidence interval, and that the
  estimator was validated by simulation-based calibration.
- **§3.4** — quantified: at 0.703 the NPV falls to 69.3%, specificity falls from 77.0% to
  40.7%, and the proportion wrongly predicted to attain rises from 3.6% to 9.2%; at the
  upper bound 0.873 it is 5.5%.
- **§4.3** — carries the honest caveat: the Dryad cohort was on renal replacement therapy
  while the primary scenario excludes it, so the estimate shows the assumed 0.94 is
  **unverified in this population**, not that 0.703 replaces it.
- **§4.2** — one sentence noting the case for direct measurement is strengthened, not
  weakened, by a lower true correlation.
- **Abstract** — Methods names the NLME fit; Results states the 0.94 assumption, the 84.1%
  attainment prevalence against which PPV/NPV must be read, and the 0.703 estimate with its
  interval; Conclusions add one sentence.

### 2. The ELF healthy-volunteer scenario is labelled as an optimistic bound (§3.3)

Dimelow's ceftazidime plasma–ELF link is saturable and avibactam's is a power function with
exponent < 1, so the penetration ratio **falls as plasma concentration rises**. The applied
0.52 is the value at 15.3 mg/L. Re-verified by
`code/elf_penetration_concentration_check.py`, which rebuilds both published functions and
reproduces all ten stated checkpoints:

| plasma | Dimelow's own model | applied |
|---|---|---|
| 15.3 mg/L | 52.0% | 52% |
| 70 mg/L | 32.0% | 52% |
| 104 mg/L (the §2.3 exposure screen) | **25.8%** | 52% |

Continuous-infusion regimens operate in that upper range. Text-only fix applied; the
scenario is not re-run, so the 44.0% figure stands but is now explicitly an optimistic
bound. **If the therapeutic-window result appears in the supplement, its "conservative"
label needs revisiting — that scenario is the concentration-appropriate one.**

### 3. Title (§ title)

> Avibactam Target Selection ~~Drives~~ **Determines Estimated** Joint Target Attainment
> During Continuous-Infusion Ceftazidime-Avibactam**: A Pharmacometric Simulation**

"Drives" reads as causal; "Estimated" locates the claim in the model. The design subtitle
is restored.

## Mechanical corrections

### 4. References renumbered into order of first citation

v17 cited [24] in §2.4 before [19]–[23], and [25]–[26] before [20]–[23]. All 27 references
were renumbered and the list reordered:

| v17 | v18 | | v17 | v18 |
|---|---|---|---|---|
| 19 | 27 | | 24 | 19 |
| 20 | 22 | | 25 | 20 |
| 21 | 24 | | 26 | 21 |
| 22 | 23 | | 27 | 26 |
| 23 | 25 | | 1–18 | unchanged |

Verified: citation order in v18 is exactly 1…27, and each renumbered citation resolves to
the same reference text as before. Five blank spacer paragraphs inside the reference list
were dropped; no content was lost.

### 5. Reference [19] removed from the §2.4 model list

v17 cited "three additional published adult population pharmacokinetic models [10-14,19]".
Barreto (*Crit Care Explor* 2021) is a therapeutic-range paper, not a population PK model.
Now `[10-14]`. Barreto remains cited in §4.3.

### 6. AI-use declaration rewritten

v17 read: "the author used OpenAI ChatGPT 5.6 Sol and Codex … Also, Montecarlo test ran by
Anthropic Clade Opus 5.0."

Two typos ("Clade", "Montecarlo") and one substantive risk: stating that an AI *ran the
analysis* is a different claim from Elsevier's policy scope, and could raise a desk query.
The replacement describes the tools as assisting with language, drafting, figure formatting
and code implementation/checking, and states explicitly that study design, model and target
selection, interpretation and all scientific conclusions are the author's own.

**The author must confirm this description is accurate before submitting.**

### 7. Affiliation and reference style

- "Iran University of Medical **Science**" → "**Sciences**"; stray space in "Farname Inc ,"
  removed; the two affiliations separated with a semicolon.
- Reference [20] (Benítez-Cano) expanded from three authors to six + *et al.*, and volume,
  issue and month added to [20], [21] and [26] to match the style of [1]–[19].

## Verified and left alone

- All 27 DOIs resolve in PubMed. No fabricated reference.
- Reference [20] (Benítez-Cano, *Crit Care* 2026;30(1):305, PMID 42129898) confirmed real,
  with correct title, journal, volume and author order.
- All 27 references are cited; none orphaned.
- `https://github.com/piranfar/caz-avi-evidence-composite-pkpd` returns HTTP 200.
- Internal arithmetic: 5 × 20,000 = 100,000; 1.71 → 8.94 g/day is 5.2-fold ("more than
  five-fold" is correct); the 84.3% plasma joint CFR sits inside the stated 73.5–88.6% range.

## Not done

The restructuring in `NOVELTY_STRATEGY.md` Part II — leading with the Fréchet bound, the
value-of-information framing, and the classification-performance analysis — is a different
paper and was deliberately left out. See `MANUSCRIPT_CORRECTIONS_v17.md` §5.
