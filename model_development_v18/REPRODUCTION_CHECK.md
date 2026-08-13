# REPRODUCTION_CHECK.md

**Project:** Avibactam Target Selection Drives Joint Target Attainment During Continuous-Infusion
Ceftazidime-Avibactam: A Pharmacometric Simulation
**Package audited:** `E:\Github Project\IJAA_submission_package_v16`
**Date of check:** 11 August 2026
**Auditor environment:** Windows 11 Pro 26200, Python 3.14.6, NumPy 2.5.0, SciPy 1.18.0,
pandas 3.0.3, Matplotlib 3.11.1 (installed during audit), Node.js v24.16.0,
`@napi-rs/canvas` ^1.0.5, `sharp` ^0.35.3 (installed during audit).
**Rule observed:** nothing in `IJAA_submission_package_v16` was modified. Every reproduction run was
performed on a copy of `revision_support/` in a scratch sandbox. The only additions to the package
are inside `model_development_v18/`.

---

## 1. Headline result

**The numerical results of this study are reproducible.** Of the 60 frozen output CSVs, **60 reproduce
to within floating-point noise (max |Δ| ≤ 1.2 × 10⁻¹⁰)** once one missing input directory is supplied.
The single exception is the probabilistic sensitivity analysis, which reproduces its **rank ordering
exactly (8/8 parameters) and its robustness classification exactly (13 fragile / 7 conditionally
robust)** but not its correlation magnitudes, because the frozen Latin-hypercube design file is
**absent from the package**.

**But the package as shipped does not run.** Ten of the twelve analysis scripts abort on import
because the data directory they depend on is not present. This is a packaging failure, not a
modelling failure — the results are right, the release is incomplete.

---

## 2. Blocking defect: the input data directory is missing

`cazavi_analyses.py` locates its inputs and frozen reference tables through `_data_root()`
(lines 40-61), which searches three candidate paths for a directory containing both `reference/`
and `inputs/`. **None of the three exists anywhere in the package, or anywhere on the E: drive.**

```
FileNotFoundError: Could not locate a data directory containing 'reference' and 'inputs'. Looked in:
  E:\Github Project\data
  E:\Github Project\IJAA_submission_package_v16\data
  E:\Github Project\IJAA_submission_package_v16\CAZ_AVI_Local_First_Reconstruction_v1\data
```

Four files are missing:

| Missing file | Consumed by | Consequence |
|---|---|---|
| `data/inputs/mic_distributions.csv` | `load_mic_distributions()` | every CFR result |
| `data/inputs/psa_draws_frozen.csv` | `load_frozen_draws()` | PSA is not bit-reproducible |
| `data/reference/primary_pta_results.csv` | `--verify` | calibration check |
| `data/reference/cfr_summary_primary_three.csv` | `--verify` | calibration check |

Because `cazavi_analyses.py` is imported by ten downstream scripts
(`critique_response`, `critique2_response`, `dose_escalation_analyses`, `make_structural_figure`,
`make_v9_figures`, `prescriptive_analyses`, `reviewer_response_analyses`, `scope_extension_analyses`,
`structural_uncertainty`, `v10_analyses`), **the entire secondary analysis suite is dead on arrival**
for anyone who downloads this package. Only `reproduce_primary_run.py` and `add_icu_elf_scenario.py`
are self-contained.

### 2.1 How the missing input was recovered

The MIC weights were recovered two ways and cross-checked:

1. **Documented (authoritative).** Supplementary Table S3b publishes isolate counts and denominators
   for `LEE2022_KPC_KP` (n = 379) and `LEE2022_OXA_KP` (n = 236). These were used directly.
2. **Recovered by inversion (forensic).** For `INDIA2022_OXA48_KP` and `FR2024_OXA484_SERINE_ONLY`,
   whose weights are published nowhere, the weight vector **w** was recovered by non-negative least
   squares from the identity CFR = **P**·**w**, where **P** is the 11 × 11 matrix of regimen-by-MIC
   PTA in `primary_pta_results.csv` and CFR is the corresponding column of
   `cfr_all_distributions.csv`, under the constraint Σw = 1.

The inversion is exact: residual **0.0000 pp** on the ceftazidime CFR used to fit, and **0.0000 pp**
on the joint CFR held out from the fit. For the two documented distributions the inversion reproduces
the published S3b weights exactly (e.g. `LEE2022_KPC_KP` at MIC 0.5 → 0.17414, published 66/379 =
0.174142). The three lowest MIC bins (0.0625, 0.125, 0.25) are not separately identifiable because
ceftazidime PTA is 100% at all three; the inversion correctly returns their sum in the lowest bin.

The recovered file is written to `model_development_v18/code/recovered_inputs/mic_distributions.csv`
with a `weight_provenance` column labelling each row `documented_S3b` or `recovered_by_inversion`.
**These are recovered values, not original inputs, and are labelled as such.**

`psa_draws_frozen.csv` **cannot be recovered** — a 300 × 8 Latin-hypercube design is not identifiable
from 20 summary correlations.

---

## 3. Reproduction results

### 3.1 Core model — `reproduce_primary_run.py`

```
python reproduce_primary_run.py --verify primary_pta_results.csv
rows compared            121
joint PTA  mean |Δ|      0.000 pp   max 0.000 pp
CAZ PTA    mean |Δ|      0.000 pp   max 0.000 pp
toxicity   mean |Δ|      0.000 pp   max 0.000 pp
```

**Caveat — this check is self-referential.** The file it verifies against,
`revision_support/primary_pta_results.csv`, carries the exact column schema that
`reproduce_primary_run.py` itself writes (`Regimen`, `EKFC class`, …) and is timestamped
**after** `outputs/primary_pta_results.csv`. It is almost certainly the script's own output from a
prior `python reproduce_primary_run.py` invocation, not the original frozen RC1 table. A perfect
0.000 pp agreement is the signature of a file comparing itself. The genuine reference tables lived in
the missing `data/reference/` directory.

The meaningful check is against the **shipped frozen outputs**, reported below, which is what a
reader of the paper would actually verify.

### 3.2 Full analysis suite — `cazavi_analyses.py all --verify`

Run with the recovered `mic_distributions.csv`, comparing against the shipped
`revision_support/outputs/` tables, keyed on identifier columns:

| Analysis | Rows | Agreement with shipped frozen output |
|---|---|---|
| Primary PTA | 121 | max \|Δ\| **0.000 pp** |
| CFR, all distributions | 44 | max \|Δ\| **1.1 × 10⁻¹⁰ pp** |
| Convergence | 5 | **byte-identical** |
| Multi-seed detail / summary | 605 / 88 | max \|Δ\| **0.000** |
| OAT sensitivity ranking | 19 | max \|Δ\| **2.1 × 10⁻¹¹ pp**; rank order identical |
| OAT detail | 380 | max \|Δ\| **8.8 × 10⁻¹¹ pp** |
| Toxicity-gate crossings | 9 | **byte-identical** |
| PSA parameter ranking | 8 | **rank order 8/8 identical**; mean \|Δ\| 0.023, max 0.066 |
| PSA probabilities | 20 | **robustness class 20/20 identical**; max \|Δ\| 4 pp |

The convergence ladder reproduces the published values exactly:
N = 1,000 → 2.014 pp; 5,000 → 0.806; 10,000 → 0.619; **50,000 → 0.169 pp** (the figure quoted in
`VERIFICATION_REPORT.md`); 100,000 → 0.000.

The OAT driver ranking reproduces exactly, including the paper's central sensitivity claim:

```
1. AVI_CT_6          20.29 pp      <- avibactam target moved to 6 mg/L
2. AVI_CT_2          17.75 pp      <- avibactam target moved to 2 mg/L
3. AVI_CL_HIGH        8.30 pp      <- 20% higher avibactam clearance
4. AVI_CL_LOW         7.44 pp
5. IIV_LOW            5.12 pp
```

**PSA discrepancy explained and bounded.** With `psa_draws_frozen.csv` absent, `run_psa()` falls
through to `latin_hypercube(300, seed=20260707)`, which resamples the design. The consequence is
confined to correlation magnitudes; every conclusion drawn from the PSA survives:

| Parameter | Regenerated | Frozen | Rank |
|---|---|---|---|
| `cl_avi_mult` | 0.830 | 0.814 | 1 = 1 |
| `omega_avi_mult` | 0.507 | 0.499 | 2 = 2 |
| `fu_avi` | 0.184 | 0.251 | 3 = 3 |
| `omega_caz_mult` | 0.123 | 0.096 | 4 = 4 |
| `rho` | 0.054 | 0.090 | 5 = 5 |
| `cl_caz_mult` | 0.045 | 0.047 | 6 = 6 |
| `tox_threshold` | 0.044 | 0.040 | 7 = 7 |
| `fu_caz` | 0.043 | 0.016 | 8 = 8 |

Robustness classification: **13 fragile, 7 conditionally robust, 0 robust** — identical to the
supplement's stated counts.

### 3.3 Secondary analysis scripts

All ten scripts that had been blocked ran to completion once the data directory was supplied.
**49 output CSVs were genuinely rewritten; all 49 reproduce exactly** (max |Δ| ≤ 1.1 × 10⁻¹³, the
largest being a rounding-level difference in `median_css_caz` in the dose-escalation tables).

| Script | Exit | Outputs regenerated | Result |
|---|---|---|---|
| `reviewer_response_analyses.py` | 0 | 5 | all reproduce |
| `scope_extension_analyses.py` | 0 | 4 | all reproduce |
| `prescriptive_analyses.py` | 0 | 7 | all reproduce |
| `dose_escalation_analyses.py` | 0 | 5 | all reproduce |
| `structural_uncertainty.py` | 0 | 8 | all reproduce |
| `critique_response.py` | 0 | 6 | all reproduce |
| `critique2_response.py` | 0 | 4 | all reproduce |
| `v10_analyses.py` | 0 | 3 | all reproduce |
| `add_icu_elf_scenario.py` | 0 | 4 | all reproduce |
| `avibactam_evidence_table.py` | 0 | 1 | all reproduce |

`structural_uncertainty.py` reproduces the manuscript's structural range verbatim:
population-weighted joint CFR **M1 84.3%, M2 94.5%, M3 84.3%, M4 70.2%** — the reported 70.2-94.5%.

### 3.4 Figures

**Matplotlib figures — 6 of 7 bit-identical.**

| Figure | Result |
|---|---|
| `fig_avi_threshold.png` (**manuscript Figure 2**) | **byte-identical** |
| `fig_pta_vs_mic.png` | byte-identical |
| `fig_cfr_distributions.png` (Suppl. S6) | byte-identical |
| `fig_oat_tornado.png` (Suppl. S7) | byte-identical |
| `fig_toxicity_gate.png` (Suppl. S8) | byte-identical |
| `fig_second_assay.png` | byte-identical |
| `fig_dose_response.png` (Suppl. S9) | differs (rendering only) |

**Node/canvas figures — 3 of 4 regenerate, none bit-identical, 1 crashes.**
Figures 1, 4, 5 and 6 of the manuscript are produced by four `.mjs` scripts requiring
`@napi-rs/canvas` and `sharp`. **The package ships no `package.json`, no lockfile, and an empty
`node_modules/`** — the dependency set and versions are undeclared and had to be guessed. With
current versions installed, three scripts run and produce output of the correct pixel dimensions and
300 dpi density, but not byte-identical output (font rasterisation differs between canvas versions).

`make_individualised_dose_figure.mjs` (**manuscript Figure 5**) **crashes**:

```
TypeError: Cannot read properties of undefined (reading 'within')
    at drawBarPanel (make_individualised_dose_figure.mjs:191:27)
```

**Root cause identified.** `write_csv()` in `cazavi_analyses.py` (line 474) opens output files with
`open(path, "w", newline="")` and **does not pin an encoding**, so Python uses the platform locale
encoding — cp1252 on this machine. The `.mjs` scripts read with `fs.readFile(csvPath, "utf8")`. The
en-dash in the EKFC class labels (`0–30`) therefore arrives as a replacement character, the lookup
key never matches, the group map stays empty, and the bar panel dereferences `undefined`.

Converting `prescriptive_decision_grid.csv` to UTF-8 makes the script run and emit a correctly sized
2100 × 1120 image at 300 dpi. **All 60 shipped CSVs are valid UTF-8**, so the original author's
environment defaulted to UTF-8; the bug only surfaces on a fresh Windows Python installation. It is
nevertheless a real portability defect and the one-line fix (`encoding="utf-8"`) belongs in v18.

### 3.5 Figure resolution against journal requirements

| Manuscript figure | Pixels | Density |
|---|---|---|
| Figure 1 | 1080 × 1600 | 300 dpi |
| **Figure 2** | 1870 × 1122 | **220 dpi** |
| **Figure 3** | 1789 × 1139 | **220 dpi** |
| Figure 4 | 2400 × 1220 | 300 dpi |
| Figure 5 | 2100 × 1120 | 300 dpi |
| Figure 6 | 2100 × 1240 | 300 dpi |

**Figures 2 and 3 are below the 300 dpi minimum** Elsevier specifies for combination artwork. Both
are matplotlib outputs and Figure 2 regenerates bit-identically, so both can be re-emitted at 300 or
600 dpi with a one-line `dpi=` change and no change to any number.

---

## 4. Checksum manifest verification

Supplementary Table S7 publishes SHA-256 checksums for 55 artifacts. Recomputed against the shipped
files: **51 match, 4 do not.**

| Artifact | Status | Interpretation |
|---|---|---|
| `cazavi_analyses.py` | **MISMATCH** | file modified after the manifest was frozen (mtime 2 Aug 22:41) |
| `make_figures.py` | **MISMATCH** | file modified after the manifest was frozen (mtime 3 Aug 16:36) |
| `population_weighted_cfr.csv` | **MISMATCH** | regenerated after the manifest was frozen |
| `structural_uncertainty_cfr.csv` | **MISMATCH** | regenerated after the manifest was frozen |

**This is a documentation defect, not a numerical one.** Both mismatched CSVs reproduce *exactly*
from the shipped code in §3.3, and `population_weighted_cfr.csv` still carries the 84.3% that the
manuscript reports. The manifest is simply stale: it was computed, then the code and two outputs were
updated, and S7 was never recomputed. The manifest must be regenerated before any release, because a
reader who checks it will find four failures and has no way to tell that they are benign.

A separate ambiguity: S7 lists `primary_pta_results.csv` once, but two files of that name exist
(`revision_support/` and `revision_support/outputs/`) with different schemas. Only the `outputs/`
one matches the published checksum. Release manifests must use paths, not bare filenames.

---

## 5. Manuscript numbers traced to frozen outputs

Every principal number in the abstract, results and tables was traced to a frozen machine-readable
output and recomputed. **All fifteen check.** Full detail in `RESULT_PROVENANCE_MATRIX.csv`;
representative verifications:

| Manuscript claim | Frozen source | Recomputed |
|---|---|---|
| Joint CFR 93.4-98.6% at AVI 1 mg/L | `avi_threshold_sweep_cfr.csv` | 93.4-98.6 ✓ |
| Joint CFR 73.5-88.6% at 4 mg/L | `avi_threshold_sweep_cfr.csv` | 73.5-88.6 ✓ |
| Joint CFR 44.0-60.7% at 8 mg/L | `avi_threshold_sweep_cfr.csv` | 44.0-60.7 ✓ |
| AVI attainment 99.2-100.0% → 44.2-61.0% | `avi_threshold_sweep_cfr.csv` | ✓ |
| Joint PTA at MIC 8 = 10.1-68.1% | `primary_pta_results.csv`, all 11 regimens | 10.1-68.1 ✓ |
| Joint PTA at MIC 4 = 67.1-88.3% | `primary_pta_results.csv`, 5 selected | 67.2-88.3 ✓ (see below) |
| Daily dose 1.71 → 8.94 g/day at MIC 8 | `individualised_attainment.csv` | ✓ |
| Within 10 g/day cap 99.2% → 57.2% | `individualised_attainment.csv` | ✓ |
| PPV 95.8%, NPV 83.6%, false reassurance 3.6% | `critique2_second_assay_operating.csv` | ✓ |
| 5.9% misclassified | derived as 100 − 94.1 accuracy | ✓ |
| Population-weighted joint CFR 84.3% (plasma) | `lung_penetration_icu_trial_summary.csv` | ✓ |
| ELF scenarios 45.9 / 44.0 / 26.7% | `lung_penetration_icu_trial_summary.csv` | ✓ |
| Structural CFR range 70.2-94.5% | `structural_uncertainty_cfr.csv` | ✓ |
| Free-vs-total: +2.0-3.0 pp; MIC 8 ≤ 0.1 pp | `critique2_free_vs_total.csv` | ✓ |

**One rounding discrepancy, immaterial.** Results §3.1 states joint PTA at MIC 4 ranged from
**67.1%** to 88.3%; the frozen table gives **67.2%** for R1 (Table 1 also prints 67.1%). A 0.1 pp
disagreement between the manuscript and its own frozen output. Harmless numerically, but it should
be corrected so that every printed digit matches the machine-readable source.

---

## 6. Version control

**The project is not a Git repository.** `git rev-parse` returns
`fatal: not a git repository`. There is therefore no commit history, no uncommitted-changes check
was possible, and no audit trail exists for the file modifications inferred from timestamps in §4.
The manuscript's Data and Code Availability statements point to
`https://github.com/piranfar/caz-avi-evidence-composite-pkpd`; whether the contents of that
repository match this package could not be verified locally and must be checked before release.

---

## 7. Defects, in priority order

| # | Severity | Defect | Fix |
|---|---|---|---|
| 1 | **Critical** | `data/inputs/` and `data/reference/` missing; 10 of 12 scripts cannot run | ship the recovered `mic_distributions.csv`; regenerate and ship a frozen PSA design |
| 2 | **Major** | `psa_draws_frozen.csv` unrecoverable; PSA not bit-reproducible | regenerate a frozen design, ship it, restate PSA numbers from it |
| 3 | **Major** | S7 checksum manifest stale for 4 of 55 artifacts | recompute manifest at release; key on path, not filename |
| 4 | **Major** | No `package.json`/lockfile; Node figure toolchain undeclared; Figure 5 generator crashes | add `package.json` + lockfile; pin `encoding="utf-8"` in `write_csv()` |
| 5 | **Moderate** | Figures 2 and 3 at 220 dpi, below journal minimum | re-emit at 300-600 dpi (no numbers change) |
| 6 | **Moderate** | `FR2024_OXA484_SERINE_ONLY` MIC distribution has no cited source anywhere | cite the source or drop the distribution |
| 7 | **Moderate** | Not a Git repository; no provenance for post-freeze edits | initialise Git; commit the frozen state before further work |
| 8 | **Minor** | Manuscript prints 67.1% where the frozen output holds 67.2% | correct to match the frozen output |
| 9 | **Minor** | `revision_support/primary_pta_results.csv` is a regenerated output masquerading as the frozen reference | remove or rename; keep references in `data/reference/` only |
| 10 | **Minor** | Three inconsistent author affiliations across manuscript, cover letter and Editorial Manager record | reconcile before resubmission |

---

## 8. Verdict

> The principal numerical results of this manuscript **can be reproduced**, and were reproduced here
> to floating-point precision across all 60 frozen outputs, with the probabilistic sensitivity
> analysis reproducing its rank ordering and robustness classification exactly but not its
> correlation magnitudes. Every number printed in the abstract, results and tables traces to a frozen
> machine-readable output, with one 0.1 pp rounding discrepancy.
>
> The **release package, however, is not self-contained**: a required input directory is absent, so a
> third party downloading it cannot run ten of the twelve analysis scripts, and the published
> checksum manifest fails on four of fifty-five artifacts. These are packaging and documentation
> defects. They do not cast doubt on the reported results, but they do defeat the reproducibility
> claim the paper makes for itself, and they must be fixed before any resubmission.
