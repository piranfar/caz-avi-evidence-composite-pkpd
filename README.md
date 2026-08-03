# Continuous-Infusion Ceftazidime–Avibactam: What the Avibactam Target Costs

Analysis code and data for the manuscript:

**What the Avibactam Target Costs: A Testing Convention, and What Individualised
Dosing Recovers in Continuous-Infusion Ceftazidime–Avibactam**

An uncertainty-aware pharmacometric evaluation of continuous-infusion
ceftazidime–avibactam in critically ill adults without renal replacement
therapy, built from published evidence with explicit provenance labelling.

> **Not a clinical tool.** This repository exists for reproducibility and
> methodological transparency. It is not a dosing guideline and must not be used
> for patient-level decisions.

## What the analysis finds

1. The single most influential quantity in the model is not pharmacokinetic. The
   assumed avibactam critical concentration moves joint CFR by more than 40
   percentage points across the range in published use, and the 4 mg/L value
   adopted in continuous-infusion practice is a susceptibility-testing
   convention rather than an exposure–response result.
2. The daily dose reaching target at the clinical breakpoint ranges from 1.71 to
   8.94 g/day across renal function, and the licensed maximum covers 99.2% of
   the lowest renal class but 57.2% at augmented clearance.
3. A single ceftazidime assay confirms adequate avibactam exposure well and
   excludes inadequate exposure poorly. Under the source model's clearance
   correlation of ρ = 0.94 it reaches 97.2% sensitivity but 77.0% specificity,
   wrongly reassuring 3.6% of the cohort; specificity falls to 21.4% at ρ = 0.50.
   Accuracy alone (94.1%) is a misleading summary because attainment is common.

## Reproducing the results

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python src/cazavi/cazavi_analyses.py all --verify
python src/cazavi/reviewer_response_analyses.py
python src/cazavi/dose_escalation_analyses.py
python src/cazavi/avibactam_evidence_table.py
python src/cazavi/scope_extension_analyses.py
python src/cazavi/prescriptive_analyses.py
python src/cazavi/structural_uncertainty.py
python src/cazavi/critique_response.py
python src/cazavi/critique2_response.py
python src/cazavi/make_figures.py
python src/cazavi/make_structural_figure.py
python src/cazavi/make_v9_figures.py
python src/cazavi/v10_analyses.py
```

`--verify` compares every analysis against the frozen RC1 tables in
`data/reference/`. The primary run reproduces those tables to within 0.06
percentage points of joint target attainment across 121 regimen–MIC rows.

## What each script does

| Script | Purpose |
|---|---|
| `reproduce_primary_run.py` | The primary Monte Carlo model, implemented from the published equations with nothing fitted |
| `cazavi_analyses.py` | CFR, convergence, multi-seed, deterministic and probabilistic sensitivity |
| `reviewer_response_analyses.py` | Analyses added at peer review |
| `dose_escalation_analyses.py` | Dose sweep at the breakpoint, suppression proxy, accumulation half-life |
| `avibactam_evidence_table.py` | Where the experimental support for an avibactam threshold actually sits |
| `scope_extension_analyses.py` | ICU renal mix, augmented clearance, variable protein binding, lung penetration |
| `prescriptive_analyses.py` | The dose each subject needs, and the grid of renal class against MIC |
| `structural_uncertainty.py` | The same simulation under four published population PK models |
| `critique_response.py` | Tests of the first eight objections raised in review, including the identity check |
| `critique2_response.py` | Free versus total avibactam target, three definitions of the limiting component, the second assay as a classifier, and the population weighting |
| `v10_analyses.py` | Penetration drawn per subject, and the second-assay operating-characteristic figure |

## A note on one withdrawn result

An earlier draft reported a "therapeutic window" and the proportion of subjects
placeable inside it. Individual clearance cancels from both sides of that
placement test, so the proportion is an identity in the MIC — 100% or 0% for the
whole cohort at once — and carries no simulated information. `critique_response.py`
demonstrates this directly (`data/processed/critique_a_window_identity.csv`,
where the distinct-value count is 1 in every row). The claim, its figure and its
supplementary sheet were withdrawn; the robustness analyses that replaced it are
in `critique_b`–`critique_f`.

## Licence

Code is MIT (`LICENSE`). Data and documentation are CC BY 4.0
(`LICENSE-DATA-DOCS.md`). Source article PDFs are not redistributed.
