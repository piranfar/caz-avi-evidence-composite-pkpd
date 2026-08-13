# model_development_v18 — index

Development directory for the continuation of the CAZ-AVI IJAA manuscript.
**Nothing in `IJAA_submission_package_v16` outside this directory has been modified** — verified by
checksum against `FILE_INVENTORY.csv` after every change in this directory.

**Status, updated 12 August 2026:**
- Phases 1–3 (audit, data-availability review, model-selection decision): **complete.**
- **Model 1** (joint clearance-correlation estimate): **complete and finalised** — see `MODEL1_REPORT.md`.
- **Model 2** (HU-JAM decision-analytic layer, all four uncertainty layers): **complete** — see
  `MODEL2_REPORT.md`.
- **Software packaging (R6):** **complete** — see `SOFTWARE.md`.
- **Data request to the Bologna group:** **sent** 11 August 2026. Awaiting reply; no follow-up planned.
- **GitHub:** merged to `main` —
  [PR #2](https://github.com/piranfar/caz-avi-evidence-composite-pkpd/pull/2) (`3ad1573`),
  [PR #3, R4](https://github.com/piranfar/caz-avi-evidence-composite-pkpd/pull/3) (`43268bc`), and
  [PR #4, R3](https://github.com/piranfar/caz-avi-evidence-composite-pkpd/pull/4) (`797ed7c`), all
  12 August 2026.
- **NOVELTY_ROUTES.md:** every route except R5 (blocked on the Bologna reply) is checked or done —
  R7, R2, R1, R6, R4, R3.
- **Manuscript text:** **not yet revised.** Phase 6 has not started.

---

## Read in this order

| # | File | What it answers |
|---|---|---|
| 1 | `PROJECT_AUDIT_REPORT.md` | What is in the package, what is authoritative, what is broken |
| 2 | `REPRODUCTION_CHECK.md` | Do the published results reproduce? |
| 3 | `RESULT_PROVENANCE_MATRIX.csv` | Every manuscript number traced to a frozen output |
| 4 | `PHASE2_DATA_AVAILABILITY_REPORT.md` | Do genuine patient data exist, and can they be used? |
| 5 | `DATA_AVAILABILITY_MATRIX.csv` | 8 studies × 35 evidence fields, each classified A–F (Phase 2; still authoritative for those 8) |
| 5b | **`data_external/README.md`** | index of all **16** archived sources, what each is for, and an honest accounting of what the collection did and did not buy |
| 6 | `MODEL_DEVELOPMENT_DECISION.md` | Which model was selected, and why |
| 7 | **`MODEL1_REPORT.md`** | The clearance-correlation estimate: 0.703, excludes the assumed 0.94 |
| 8 | **`MODEL2_REPORT.md`** | The decision layer: value of information, misselection/regret, triage |
| 9 | `SOFTWARE.md` | The installable package, and how to adapt it to another drug pair |
| 10 | `NOVELTY_ROUTES.md` | What else could raise the paper's ceiling, ranked, live-updated |

---

## Directory layout

```
model_development_v18/
├── PROJECT_AUDIT_REPORT.md, REPRODUCTION_CHECK.md, RESULT_PROVENANCE_MATRIX.csv, FILE_INVENTORY.csv
│                                     Phase 1: audit of the v16 package
├── PHASE2_DATA_AVAILABILITY_REPORT.md, DATA_AVAILABILITY_MATRIX.csv
│                                     Phase 2: real patient-data availability
├── MODEL_DEVELOPMENT_DECISION.md    Phase 3: which model, and why
├── MODEL1_REPORT.md                 Model 1: the clearance-correlation estimate (FINAL)
├── MODEL2_REPORT.md                 Model 2: the decision-analytic layer
├── MODEL2_SPECIFICATION.md          Model 2's design, written before implementation
├── SOFTWARE.md                      the installable package (R6)
├── NOVELTY_STRATEGY.md              why the paper's spine changed, and to what
├── NOVELTY_ROUTES.md                live register of further novelty routes, ranked
├── pyproject.toml, LICENSE          packaging metadata (installs as `hujam`)
│
├── audit/
│   └── log_*.txt                    stdout from every reproduction and analysis run
│   (audit/extracted/ — unpublished manuscript text extraction — exists locally only;
│    deliberately excluded from every GitHub push)
│
├── code/
│   ├── __init__.py, interface.py    package entry point + the engine/decision-layer contract
│   ├── model2_engine.py             the CAZ-AVI-specific engine, regression-tested against v16
│   ├── model2_hujam.py              the portable decision layer (VOI, misselection, regret)
│   ├── model2_monitoring.py, model2_triage.py, model2_heterogeneity.py,
│   │   model2_layer2_compare.py, model2_breaking_point.py
│   │                                 Model 2's four uncertainty layers and analyses
│   ├── joint_popk_nlme.py, model1_finalise.py, model1_sbc.py, model1_design_analysis.py
│   │                                 Model 1: fitting, diagnostics, validation, design analysis
│   ├── test_model1.py               119 checks
│   ├── test_model2.py               44 checks: decision-layer identities, samplers, interface
│   ├── clearance_correlation_analysis.py, external_aggregate_check_BenitezCano.py
│   │                                 the two-dataset correlation analysis and the external check
│   ├── docx2txt.py, build_*.py, append_*.py
│   │                                 Phase 1/2 report generators
│   └── recovered_inputs/
│       └── mic_distributions.csv    RECOVERED, not original — see below
│
├── outputs/                         every frozen, machine-readable result — one CSV per analysis
├── figures/                         Model 1 goodness-of-fit and visual predictive check
│
├── data_external/
│   ├── Gatti2023_JCritCare_154301_CC-BY-NC-ND.pdf, Gatti2023_individual_patient_data.csv
│   ├── README_Gatti2023_provenance.md
│   ├── dryad_Li2025_CRRT/           CC0 public-domain dataset, 21 patients
│   │   ├── *.csv, README_provenance.md, rebuild_from_dryad_preview.py
│   ├── OJeanson2024_CVVHDF/         CC BY 4.0, n=4, retrieved via Elsevier's own API
│   │   ├── OJeanson2024_IJAA_107394_CC-BY.pdf, OJeanson2024_PK_parameters.csv, README_provenance.md
│   ├── Tian2025_CVVH/               subscription access (NYU), n=7 — PDF EXCLUDED from GitHub, CSV is not
│   │   ├── Tian2025_EJCMID_..._SUBSCRIPTION-ACCESS.pdf, Tian2025_PK_parameters.csv, README_provenance.md
│   ├── Wu2025_PopPK_AKI/            subscription access (NYU), n=31+32 — a real fitted PopPK+AKI-TTE
│   │   │                            model, not a summary table; PDF EXCLUDED from GitHub, CSV is not
│   │   ├── Wu2025_JAC_dkaf275_SUBSCRIPTION-ACCESS.pdf, Wu2025_PopPK_parameters.csv, README_provenance.md
│   ├── Lanini2024_nonRRT_CrCl/      subscription access (NYU), n=52 — the ONLY non-RRT dataset here;
│   │   │                            PDF EXCLUDED from GitHub, CSV is not
│   │   ├── Lanini2024_IJAA_107351_SUBSCRIPTION-ACCESS.pdf, Lanini2024_CrCl_correlation.csv, README_provenance.md
│   ├── Fresan2023_CI_TDM/           free-to-read but NOT OA-licensed, n=31 INDIVIDUAL patients, true
│   │   │                            continuous infusion. NO avibactam measured — cannot inform the
│   │   │                            correlation question. PDF EXCLUDED from GitHub, CSV is not
│   │   ├── Fresan2023_JAC_dkac439_FREE-TO-READ.pdf, Fresan2023_individual_patient_data.csv, README_provenance.md
│   ├── Chen2025_PopPK_CRKP/         open access via PMC, n=45 — 2nd fitted PopPK model; the only
│   │   │                            dataset separating critically vs non-critically ill. No PDF; CSV only
│   │   ├── Chen2025_PopPK_parameters.csv, README_provenance.md
│   ├── Cojutti2024_ANCHOR_PopPK/    *** THE ANCHOR *** source of rho=0.94 and all six published RSEs
│   │   │                            this project propagates. n=112, both compounds measured, RRT
│   │   │                            excluded. All 7 project parameters verified exact. No PDF; CSV only
│   │   ├── Cojutti2024_PopPK_parameters.csv, README_provenance.md
│   ├── Gatti2025_outcome_R5/        CC BY 4.0 (fully open), n=218 — the clinical outcome paper behind
│   │   │                            novelty route R5. No push restriction. No PDF; CSV only
│   │   ├── Gatti2025_outcome_associations.csv, README_provenance.md
│   ├── Li2019_registrational_PopPK/ the registrational reference model, n=1975/2249 — largest CAZ-AVI
│   │   │                            dataset in existence. Merged both drugs' random effects per patient
│   │   │                            yet NEVER computed their correlation (key R1 evidence). Includes
│   │   │                            Data S1 NONMEM control streams + Data S2 + bootstrap 90% CIs.
│   │   │                            PDF and both .docx EXCLUDED from GitHub (licence unconfirmed)
│   │   ├── Li2019_PopPK_parameters.csv, Li2019_bootstrap_CIs.csv, README_provenance.md
│   │   ├── Li2019_CTS_12-151_SUBSCRIPTION-ACCESS.pdf, Li2019_DataS1_NONMEM_control_streams.docx,
│   │   └── Li2019_DataS2_supplementary_methods_tables.docx
│   ├── Das2019_dose_selection/      the registrational dose-selection paper — source of the licensed
│   │   │                            regimen AND of the "50% fT > 1 mg/L" target. No PDF; CSV only
│   │   ├── Das2019_target_derivation.csv, README_provenance.md
│   ├── Coleman2014_hollowfibre_T7/  the PRIMARY SOURCE Model 2's live T7 scenario is built from
│   │   │                            (3 strains, not 8). Not OA. No PDF; CSV only
│   │   ├── Coleman2014_regrowth_thresholds.csv, README_provenance.md
│   ├── Berkhout2016_murine_index/   the study that DEFINED the %fT>C_T index and C_T = 1 mg/L;
│   │   │                            underpins T1 and all of route R3. Not OA. No PDF; CSV only
│   │   ├── Berkhout2016_fT_over_CT.csv, README_provenance.md
│   ├── Gatti2023_joint_target_origin/  the AAC paper that DEFINED the aggressive joint PK/PD target
│   │   │                            Cojutti 2024 and Gatti 2025 both use. ABSTRACT ONLY (PMC empty)
│   │   ├── Gatti2023_joint_target.csv, README_provenance.md
│   └── Gatti2024_ratio_one_leg/     *** the empirical test of THIS PROJECT'S CENTRAL QUESTION ***
│       │                            n=107/188 paired TDM: CAZ:AVI ratio 1.29:1-13.46:1 vs 4:1 vial.
│       │                            Rebuts Fresan 2023's single-analyte assumption. Not OA; CSV only
│       ├── Gatti2024_CAZ_AVI_ratio.csv, README_provenance.md
│
└── correspondence/
    └── DRAFT_data_request_Bologna_Gatti_Pea.md      SENT 11 Aug 2026 — no follow-up planned
    └── DRAFT_data_request_Benitez-Cano_Sorli.md     DRAFT ONLY — NOT SENT
```

---

## Things to know before using anything here

**1. `code/recovered_inputs/mic_distributions.csv` is recovered, not original.**
The v16 package's `data/inputs/` directory is missing, so most of its analysis scripts cannot run
from a clean checkout. The MIC weights here were recovered — two distributions from the published
Supplementary Table S3b, two by non-negative least-squares inversion of the frozen CFR against the
frozen PTA matrix (residual 0.0000 pp, cross-validated against the two documented distributions).
Every row carries a `weight_provenance` column. **Do not treat this file as an original project
input.**

**2. Every patient dataset used is class C — assumption-testing only, not validation.**
There are now **sixteen** entries in `data_external/`. Fourteen are clinical datasets, none of which is
a validation set; two (**Coleman 2014**, **Berkhout 2016**) are *preclinical* primary sources
— in-vitro hollow-fibre and murine — archived because Model 2's T7 and T1 scenarios and the whole of
route R3 rest on them directly. **Gatti 2024 ("why one leg is not enough to run") is the closest thing
here to a direct empirical test of this project's own central question** — whether measuring one
component tells you about the other — and it answers no: the ceftazidime-to-avibactam ratio ranges
1.29:1 to 13.46:1 against a 4:1 vial ratio. It also turns out to be one side of a published exchange
rebutting Fresan 2023, which this project had already archived without knowing it was contested.
Of the clinical datasets, four
(Gatti 2023, Li 2025, O'Jeanson 2024, Tian 2025) are RRT-only cohorts, which the manuscript's primary
scenario excludes by definition; Model 1 itself is likewise fitted to RRT patients, so its
clearance/volume estimates do not transfer to the primary scenario — only its correlation estimate is
carried forward, and only as a sensitivity bound (`MODEL1_REPORT.md` §5). Of the other eight:
**Lanini 2024** explicitly excludes RRT and renally-adjusted dosing — the closest population match to
the primary scenario — but is cross-sectional with *derived*, not measured, free concentrations.
**Wu 2025** and **Chen 2025** are mixed critically-ill cohorts (Chen alone separates critically from
non-critically ill, and found CRRT *not* to be a significant covariate on clearance — worth knowing
before assuming RRT status always dominates). **Fresan 2023** is individual-patient data under true
continuous infusion, but **measured ceftazidime only — no avibactam at all**, so it cannot speak to the
clearance-correlation question that is central to this project. **Cojutti 2024 is the anchor** — the
published source of ρ = 0.94 and of all six RSEs Model 2 propagates; archiving it verified all seven
against the source exactly (see that folder), but it is a model summary, not patient-level data.
**Gatti 2025** is the clinical-outcome paper behind route R5, and is association evidence from a
pre-post design with bundled co-interventions — never causal. **Li 2019 and Das 2019 are the two
registrational reference papers** — the parameter package and the dose-selection document behind the
licensed label; sponsor-run, not independent, and neither contains patient-level data. Read each
folder's own `README_provenance.md` before using any of them; the limitations differ substantially and
are not interchangeable.

**3. What has actually been sent or pushed, and what has not.** The Bologna data request was sent —
see the correspondence file for the exact text and the standing no-follow-up rule. The Benítez-Cano
request is still a draft. GitHub work has been merged to `main` (PR #2, merge commit `3ad1573`, 12
August 2026), so `main` on `caz-avi-evidence-composite-pkpd` now reflects it — this is what the
manuscript's data-availability statement points readers to. Nothing has
been pushed from the *local* copy of this whole package (which includes the manuscript text and
correspondence) — only a curated subset, built specifically to exclude those, has gone to GitHub.
See `SOFTWARE.md` and the PR description for exactly what that subset is. **Four more PDFs join that
exclusion list as of 12 August 2026** — none is openly licensed, so all four stay local-only:
`data_external/Tian2025_CVVH/Tian2025_EJCMID_s10096-025-05343-x_SUBSCRIPTION-ACCESS.pdf`,
`data_external/Wu2025_PopPK_AKI/Wu2025_JAC_dkaf275_SUBSCRIPTION-ACCESS.pdf`,
`data_external/Lanini2024_nonRRT_CrCl/Lanini2024_IJAA_107351_SUBSCRIPTION-ACCESS.pdf` (retrieved after
a ScienceDirect CAPTCHA the user solved themselves — Claude never attempts that), and
`data_external/Fresan2023_CI_TDM/Fresan2023_JAC_dkac439_FREE-TO-READ.pdf` (**free to read is not the
same as openly licensed** — this one carries OUP's "All rights reserved" standard model, so it is
excluded on copyright grounds even though no subscription was needed to obtain it). See each folder's
`README_provenance.md`. **Three more join it on 12 August 2026** —
`data_external/Li2019_registrational_PopPK/` holds the Version-of-Record PDF plus both supplement
`.docx` files (Data S1 NONMEM control streams, Data S2 methods and bootstrap tables); the project's
audit records this article as CC BY-NC but **that licence could not be confirmed from the article
text**, so all three are excluded pending verification. `OJeanson2024_CVVHDF/` (CC BY 4.0) and
`Gatti2025_outcome_R5/` (CC BY 4.0) have no such restriction, and `Chen2025_PopPK_CRKP/`,
`Cojutti2024_ANCHOR_PopPK/`, `Das2019_dose_selection/`, `Coleman2014_hollowfibre_T7/` and
`Berkhout2016_murine_index/` hold no PDFs at all — only extracted facts.

**4. The manuscript text is untouched.** Phase 6 (manuscript revision) has not started. Every number
above is verified and frozen, but none of it is in the manuscript yet.

---

## Reproducing the analyses in this directory

```bash
pip install -e .                              # installs as `hujam`; see SOFTWARE.md
cd code
python model2_engine.py                       # regression test against the frozen v16 outputs
python test_model1.py                         # 119 checks
python test_model2.py                         # 44 checks
python joint_popk_nlme.py && python model1_finalise.py   # Model 1, full fit + diagnostics
python model2_hujam.py                        # Model 2, full decision-analytic pipeline
```

Environment: Windows 11 Pro 26200, Python 3.14.6, NumPy 2.5.0, SciPy 1.18.0, pandas 3.0.3,
Matplotlib 3.11.1, openpyxl 3.1.5, pypdf, Node.js v24.16.0 with `@napi-rs/canvas` ^1.0.5 and
`sharp` ^0.35.3.
