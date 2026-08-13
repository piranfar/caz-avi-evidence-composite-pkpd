# PROJECT_AUDIT_REPORT.md — Phase 1

**Manuscript:** Avibactam Target Selection Drives Joint Target Attainment During Continuous-Infusion
Ceftazidime-Avibactam: A Pharmacometric Simulation
**Author:** Vahhab Piranfar · **Target journal:** International Journal of Antimicrobial Agents
**Package audited:** `E:\Github Project\IJAA_submission_package_v16` (156 files, 66.9 MB)
**Audit date:** 11 August 2026
**Companion documents:** `FILE_INVENTORY.csv`, `RESULT_PROVENANCE_MATRIX.csv`, `REPRODUCTION_CHECK.md`

**Nothing in the v16 package was modified.** All reproduction runs were executed on a sandbox copy.
Every new artifact is inside `model_development_v18/`.

---

## 1. Executive summary

The science is sound and the numbers are reproducible. The **release package is not**.

| Dimension | Verdict |
|---|---|
| Are the reported results reproducible? | **Yes** — 60/60 frozen outputs reproduce to ≤ 1.2 × 10⁻¹⁰; PSA reproduces rank order 8/8 and robustness class 20/20 |
| Do the manuscript's numbers trace to frozen outputs? | **Yes** — 25 of 30 traced claims match exactly; 1 has a 0.1 pp typographical discrepancy |
| Are the source parameters correctly transcribed? | **Yes** — all seven Cojutti 2024 parameters independently confirmed against the source publication |
| Can a third party run the shipped code? | **No** — a required input directory is absent; 10 of 12 analysis scripts abort on import |
| Is the published checksum manifest valid? | **No** — 4 of 55 artifacts fail |
| Can all figures be regenerated? | **Partly** — 6/7 matplotlib figures are bit-identical; the Node toolchain is undeclared and one generator crashes |
| Is there version control? | **No** — not a Git repository |

**The single most important finding of this audit is not a defect.** It is that genuine, openly
licensed, individual patient-level ceftazidime and avibactam clearance data **do exist** in the
literature already cited by this manuscript, and they permit an empirical test of the assumption
that most strongly drives one of the paper's clinical conclusions. See §6.

---

## 2. Package structure and authoritative files

The package has no README and no manifest at its root. Authorship of each deliverable had to be
established from timestamps, content and checksums.

### 2.1 Authoritative deliverables

| Deliverable | Authoritative file | Modified | Note |
|---|---|---|---|
| **Manuscript** | `Piranfar_CAZ-AVI_IJAA_Original_Article_v17.docx` | 11 Aug 10:07 | **v17, despite the folder being named v16** |
| Manuscript PDF | `Piranfar_CAZ-AVI_IJAA_Original_Article_v17.pdf` | 11 Aug 10:08 | |
| **Editorial Manager submission** | `IJAA-S-26-02124_3.pdf` | 11 Aug 10:16 | 24 pp, system-generated; carries manuscript number **IJAA-S-26-02124**, revision `_3` |
| Tables | `Tables_IJAA.docx` | 5 Aug 17:26 | Tables 1-6 with data, not legends only |
| Supplementary | `Supplementary_Tables_S1-S21_v13.xlsx` | 4 Aug 14:43 | 22 sheets |
| Figures | `Figures/Figure_1..6_*.png` | 5 Aug | byte-identical to `revision_support/figures/` sources |
| Cover letter | `Cover_Letter_IJAA_v16.docx` / `.pdf` | 5 Aug 17:31 | |
| Highlights | `Highlights_IJAA.docx` / `.pdf` | 5 Aug | |
| Graphical abstract | `Graphical_Abstract_..._v21_300dpi_Print.png` / `.pdf` | 5 Aug 15:03 | |
| Analysis code | `revision_support/*.py`, `*.mjs` | 2-4 Aug | 32 Python + 6 Node scripts |
| Frozen outputs | `revision_support/outputs/*.csv` | 2-4 Aug | 60 CSVs + `VERIFICATION_LOG.txt` |

**The manuscript is under active submission to IJAA** (manuscript number IJAA-S-26-02124, at least
three versions). Any revision must be coordinated with that submission; this is not a
pre-submission draft.

### 2.2 Superseded and duplicate material — retained, not authoritative

- `revision_support/manuscript_package/` — the v3 "Evidence-Composite" generation, including
  `Supplementary_Tables_S1-S9 V5.xlsx` (9 sheets, superseded by the 22-sheet v13 workbook).
- `revision_support/outputs/v3_before_*.docx`, `v3_final.md`, `v3_readable.md` — pre-edit snapshots.
- `revision_support/مقاله_فارسی_v5.md` — a 67 kB Persian-language version, with its builder
  `build_persian_docx.py`. Not part of the IJAA submission.
- `revision_support/apply_*.py`, `fix_v9.py`, `trim_v9.py`, `sync_*.py`, `organize_repo.py` — 15
  one-shot manuscript-editing scripts. **These are document-mutation history, not analysis code**,
  and should not be shipped as "analysis code" in a reproducibility repository.
- `_review_render/`, `revision_support/qa_current_front/`, `revision_support/review_current/`,
  `revision_support/node_modules/` — **all four directories are empty**.
- `revision_support/__pycache__/` — compiled bytecode for Python 3.12 and 3.14; should be excluded.

### 2.3 Two files that are not what their names suggest

- **`revision_support/primary_pta_results.csv`** sits beside the code and looks like the frozen
  reference table. It carries the exact column schema `reproduce_primary_run.py` writes and is
  timestamped *after* the real frozen output in `outputs/`. It is almost certainly the script's own
  output from a prior run. Verifying the model against it is self-referential — which is why that
  check returns a suspiciously perfect 0.000 pp.
- **`Supplementary_Tables_S1-S21_v13.xlsx` sheet S7** lists checksums keyed on bare filenames. Two
  distinct files are named `primary_pta_results.csv`; only the one in `outputs/` matches the
  published checksum.

---

## 3. Analysis architecture

The pipeline is a clean two-layer design, and the layering is genuinely good practice:

```
reproduce_primary_run.py      core model; constants carry inline provenance comments
        │
        └── cazavi_analyses.py    population draw, scenario engine, CFR, convergence,
                │                 multiseed, OAT sensitivity, Latin-hypercube PSA
                │
                ├── reviewer_response_analyses.py   avibactam-threshold sweep, renal boundary
                ├── prescriptive_analyses.py        individualised dose, decision grid
                ├── structural_uncertainty.py       four alternative population PK models
                ├── dose_escalation_analyses.py     escalation, resistance-suppression targets
                ├── scope_extension_analyses.py     ARC subgroup, protein binding, population CFR
                ├── critique_response.py            window identity, RRT exclusion, variance
                ├── critique2_response.py           free-vs-total, second assay, weighting
                ├── v10_analyses.py                 penetration dependence and variability
                ├── make_structural_figure.py       Supplementary Figure S12
                └── make_v9_figures.py
        │
        └── add_icu_elf_scenario.py   ELF scenarios (imports the core only)
```

Two design decisions deserve explicit credit. **Common random numbers**: the population is drawn
once as standard-normal deviates and whitened, so every one-at-a-time scenario rebuilds its random
effects from identical underlying draws — a scenario difference therefore reflects the parameter
change and not Monte Carlo noise. **Nested convergence sampling**: smaller runs reuse the leading
subjects of the 100,000-subject reference rather than redrawing, so convergence is not confounded
with between-run noise. Both are correct choices that many published target-attainment analyses get
wrong.

Output-to-script attribution for all 66 outputs is recorded in `RESULT_PROVENANCE_MATRIX.csv`. Four
outputs (`v3_*.docx/.md`) have no identifiable producing script; all four are superseded manuscript
snapshots, not results.

---

## 4. Reproduction outcome

Full detail in `REPRODUCTION_CHECK.md`. Summary:

- **The package as shipped does not run.** `cazavi_analyses.py` needs `data/inputs/` and
  `data/reference/`; neither exists in the package or anywhere on the drive. Ten downstream scripts
  fail on import.
- **Once the missing MIC-distribution input was recovered, everything reproduced.** The weights were
  recovered exactly — by non-negative least squares inversion of CFR = **P**·**w** against the frozen
  PTA matrix, residual **0.0000 pp**, and cross-validated against the two distributions whose weights
  *are* published in Supplementary Table S3b.
- **60 of 60 frozen outputs reproduce** to ≤ 1.2 × 10⁻¹⁰. The manuscript's headline sensitivity
  result reproduces exactly: the avibactam target at 6 mg/L shifts joint CFR by **20.29 pp** and at
  2 mg/L by **17.75 pp**, against **8.30 pp** for a 20% error in avibactam clearance. The claim that
  the analyst's target choice dominates the pharmacokinetics is correct and reproducible.
- **The PSA is the one gap.** `psa_draws_frozen.csv` is missing and a 300 × 8 Latin-hypercube design
  is not recoverable from 20 summary correlations. Regenerating it preserves the parameter rank
  order (8/8) and the robustness classification (13 fragile / 7 conditionally robust, exactly as
  published) but moves correlation magnitudes by up to 0.066.

---

## 5. Discrepancies requiring disclosure before any result is relied upon

Ordered by severity. Items 1-4 must be fixed before resubmission; 5-10 should be.

**1 — Critical. The reproducibility claim is not currently satisfiable by a reader.**
The manuscript states that code and frozen outputs are available for reproduction. A reader who
obtains this package cannot run ten of the twelve analysis scripts. Either the GitHub repository
named in the Data Availability statement contains the missing `data/` directory — which must be
verified — or the claim is not met.

**2 — Critical. `psa_draws_frozen.csv` is unrecoverable.** The published PSA correlation magnitudes
cannot be exactly reproduced by anyone. The fix is to generate a frozen design, ship it, and restate
the PSA numbers from it. Conclusions do not change.

**3 — Major. The Supplementary Table S7 checksum manifest is stale.** Four of 55 artifacts fail:
`cazavi_analyses.py`, `make_figures.py`, `population_weighted_cfr.csv`,
`structural_uncertainty_cfr.csv`. All four were modified after the manifest was frozen, and both
CSVs reproduce exactly from the shipped code — so this is a documentation failure, not a numerical
one. But a reader who checks the manifest will see four failures with no way to tell they are benign.
**A published checksum manifest that fails is worse than none.**

**4 — Major. `FR2024_OXA484_SERINE_ONLY` has no source.** This MIC distribution is used in the PSA
and the robustness ledger (Supplementary Table S17) and appears in the frozen CFR outputs, but **no
citation for it exists anywhere** in the manuscript, its reference list, or the supplement. It must
be sourced or removed. Its weights are recoverable by inversion but that establishes what was used,
not where it came from.

**5 — Moderate. Figures 2 and 3 are 220 dpi**, below Elsevier's 300 dpi minimum for combination
artwork. Both are matplotlib outputs and Figure 2 regenerates bit-identically, so this is a one-line
`dpi=` change affecting no number.

**6 — Moderate. Figure 5's generator crashes** on a fresh installation. Root cause diagnosed:
`write_csv()` does not pin `encoding="utf-8"`, so on a Windows locale it emits cp1252, while the
Node scripts read UTF-8 — the en-dash in the EKFC class labels breaks the lookup key. One-line fix.
Separately, the Node toolchain has **no `package.json` and no lockfile**; `@napi-rs/canvas` and
`sharp` had to be identified from import statements and installed at guessed versions.

**7 — Moderate. No version control.** `git rev-parse` returns *not a git repository*. There is no
commit history, so the post-freeze modifications inferred in item 3 have no audit trail.

**8 — Minor but must be fixed. A 0.1 pp mismatch between the manuscript and its own frozen output.**
Results §3.1 and Table 1 print joint PTA at MIC 4 as **67.1%** for R1; the frozen table holds
**67.2%**. In a paper whose central methodological claim is full traceability to machine-readable
outputs, every printed digit must match.

**9 — Minor. Three inconsistent author affiliations.** The manuscript states *"Department of
Microbiology, Iran University of Medical Science, Tehran, Iran; Farname Inc, Th, Canada"* — with
"Th" evidently a placeholder — the cover letter states *"Independent Researcher, New York, United
States"*, and the Editorial Manager record states *"Farname Inc, Jersey City, NJ, United States"*.
These must be reconciled.

**10 — Minor. The ELF penetration ratio choice is undocumented.** Table 3 uses 0.41/0.44 from the
ICU pneumonia trial. That source publishes **three** different ratio pairs: median individual
AUC ratios 0.41/0.44 (used), a model-based ELF-plasma ratio 0.51/0.60, and a supplementary observed
ratio 0.46/0.40. Using the median individual ratios is defensible, but the choice is not stated and
no sensitivity to it is reported — and since ELF joint CFR (45.9%) is one of the paper's more
striking numbers, the alternative that would raise it should be shown.

---

## 6. What the audit found beyond the defects

Two findings materially change what this project can become. Both are developed in
`DATA_AVAILABILITY_MATRIX.csv` and `MODEL_DEVELOPMENT_DECISION.md`.

### 6.1 Every source parameter is confirmed, and the source Ω matrix is fully specified

Independent verification against Cojutti 2024 confirmed all seven pharmacokinetic constants —
CL 5.0 L/h with exponent 0.70, CL 5.9 L/h with exponent 0.89, CV 67.92% and 76.91%, ρ = 0.94 — plus
the unbound fractions 0.85/0.92, both PK/PD targets, and the 104 mg/L total-Css exposure ceiling.
Nothing was refuted.

A useful structural fact emerged: because interindividual variability in that model is estimated on
**only two parameters** and the single off-diagonal correlation is published, **the source Ω matrix
is completely specified** — unusual, and worth stating explicitly in the Methods, because it means
the reimplementation inherits no unreported covariance terms.

### 6.2 Genuine individual patient data exist, and they test the assumption that matters most

Gatti et al. 2023 — **already reference [10] of this manuscript** — publishes complete patient-level
data under an open CC BY-NC-ND licence: Table 1 gives 8 patients with MICs, free steady-state
concentrations and outcomes; **Table 2 gives 17 therapeutic-drug-monitoring occasions with paired
ceftazidime and avibactam clearances.** The Version of Record is freely downloadable. These values
have been transcribed to `data_external/Gatti2023_individual_patient_data.csv`; the transcription
reproduces the authors' independently reported medians *and* interquartile ranges exactly for both
analytes, which verifies it.

This matters because **ρ = 0.94 is the parameter behind the manuscript's assay conclusion** — the
finding that ceftazidime-based prediction of avibactam attainment achieves PPV 95.8% and NPV 83.6%
with only 5.9% misclassification. In these real paired data the correlation of log clearances is
**r = 0.598 (95% CI 0.165 to 0.838, n = 17 occasions)**, significantly below 0.94 (p ≈ 9 × 10⁻⁵). At
the patient-mean level, r = −0.34 with a CI spanning nearly the whole range (n = 8).

**This is a comparison to be reported carefully, not a refutation.** These are CVVHDF patients; the
primary scenario excludes renal replacement therapy. Extracorporeal clearance compresses
between-patient variability — apparent CV ≈ 10% and 13% here against 67.9% and 76.9% in the non-RRT
cohort — and imposes a shared circuit-driven component. The correlation structure in this population
is not the correlation structure in the primary scenario. What the comparison does establish is that
**ρ = 0.94 has never been checked against paired patient data outside the cohort that produced it**,
and that the manuscript's own sensitivity analysis already shows the assay conclusion is fragile to
it: at ρ = 0.5, NPV falls from 83.6% to 63.2% and misclassification rises from 5.9% to 12.3%.

### 6.3 Where an independent aggregate-level check is available

Benítez-Cano 2026 — **already reference [25]** — is an independent randomized PK trial in 15
non-RRT ICU patients on continuous-infusion ceftazidime-avibactam, with published median and
interquartile steady-state exposures for both components. Its individual concentrations are **not**
published; they exist only as points in per-patient supplementary figure panels, and the numeric data
are gated behind an on-request statement. But its aggregate exposure summaries permit a genuine
**external, aggregate-level predictive check** of the primary model against a cohort that had no
part in building it — the strongest form of external evidence available to this project without
obtaining data from authors.

---

## 7. What the audit did not find

- **No fabricated or untraceable results.** Every principal number traces to a frozen output.
- **No evidence of tuning to outputs.** The core model's constants are published values, and the
  reimplementation reproduces the frozen tables without any fitted parameter.
- **No overstatement of validation status.** The manuscript, the supplement and
  `VERIFICATION_REPORT.md` all correctly describe the 72-row calibration as an *internal reproduction
  check* and explicitly state it is not external validation. Supplementary Table S8, limitation L9,
  says so directly. This is honest reporting and it should be preserved verbatim through any revision.
- **No patient data of any kind in the current package**, consistent with the manuscript's ethics
  statement.

---

## 8. Recommended immediate actions

**Before any new modelling work:**

1. Initialise Git in `IJAA_submission_package_v16`, commit the current state untouched as the
   baseline, and work only in `model_development_v18/` thereafter.
2. Verify whether `https://github.com/piranfar/caz-avi-evidence-composite-pkpd` contains the missing
   `data/inputs/` and `data/reference/` directories. If it does, this audit's item 1 is a packaging
   problem only. If it does not, the Data Availability statement is not currently satisfiable.
3. Source or remove the `FR2024_OXA484_SERINE_ONLY` MIC distribution.

**Fixes that change no number and can be made immediately:**

4. Correct 67.1% → 67.2% in Results §3.1 and Table 1.
5. Re-emit Figures 2 and 3 at 300-600 dpi.
6. Pin `encoding="utf-8"` in `write_csv()`; add `package.json` and a lockfile.
7. Recompute the Supplementary Table S7 checksum manifest, keyed on path rather than filename.
8. Reconcile the three author affiliations.

**Requires a rerun and re-verification:**

9. Generate and freeze a new Latin-hypercube PSA design; restate the PSA correlation magnitudes.

None of items 4-9 alters a single scientific conclusion.

---

## 9. Phase 1 verdict

> The manuscript's numerical results are **reproducible and correctly sourced**. The reimplementation
> regenerates all 60 frozen outputs to floating-point precision, every source parameter was
> independently confirmed against the original publication, and the paper's characterisation of its
> own calibration as an internal check rather than external validation is accurate and should be
> preserved.
>
> The **release package is incomplete** in ways that defeat the reproducibility claim the paper makes
> for itself: a required input directory is missing, the published checksum manifest fails on four
> artifacts, one figure generator crashes, and one MIC distribution used in the analysis has no cited
> source. These are packaging and documentation defects, all fixable, none of which casts doubt on a
> reported result.
>
> Separately, the audit establishes that **genuine individual patient data relevant to this
> manuscript are publicly and legally available** in a study it already cites, and that they bear
> directly on the assumption underpinning one of its clinical conclusions. That is the most promising
> route to strengthening the work, and it is developed in Phase 2 and Phase 3.
