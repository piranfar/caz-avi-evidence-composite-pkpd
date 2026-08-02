# Continuous-Infusion Ceftazidime–Avibactam: Dosing Window and Uncertainty

Analysis code and data for the manuscript:

**Continuous-Infusion Ceftazidime–Avibactam in Critically Ill Adults: A Narrow Dosing Window Bounded by Avibactam Exposure and Ceftazidime Toxicity**

An uncertainty-aware pharmacometric evaluation of continuous-infusion
ceftazidime–avibactam in critically ill adults without renal replacement
therapy, built from published evidence with explicit provenance labelling.

> **Not a clinical tool.** This repository exists for reproducibility and
> methodological transparency. It is not a dosing guideline and must not be used
> for patient-level decisions.

## Reproducing the results

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python src/cazavi/cazavi_analyses.py all --verify
python src/cazavi/reviewer_response_analyses.py
python src/cazavi/make_figures.py
```

`--verify` compares every analysis against the frozen RC1 tables in
`data/reference/`. Expected agreement:

| Analysis | Agreement |
|---|---|
| Primary 100,000-subject run | mean abs difference 0.060 pp joint PTA (121 rows) |
| MIC-weighted CFR | 0.041 pp (33 rows) |
| Convergence at N = 50,000 | 0.169 pp against 0.167 recorded |
| Deterministic sensitivity | same rank order; top driver 20.29 vs 19.75 |
| Probabilistic sensitivity | all 8 parameters in exact rank order; robustness 20/20 |

Two further analyses have no frozen counterpart, because they were added after
the reviewed draft: the sweep across avibactam critical concentrations (1-8
mg/L), and the renal-function boundary sensitivity that separates the lowest
class at 15 mL/min/1.73 m2.

## What is and is not reproducible

The original analysis code was not preserved. The model in `src/cazavi/` was
reimplemented from the published equations and parameter values, with nothing
fitted to the frozen outputs, and reproduces them to within Monte Carlo noise.
Agreement is statistical, not bitwise: the original random-number stream cannot
be recovered.

## Layout

```text
src/cazavi/          verified pipeline: model, analyses, figures
src/verify_checksums.py
src/legacy/          superseded figure script, kept for reference
data/inputs/         MIC distributions, published calibration anchors, PSA draws
data/reference/      frozen RC1 tables; the reproduction target
data/processed/      outputs of the current pipeline
data/legacy_rc1/     original RC1 summary CSVs
figures/             manuscript figures 1-6
figures/legacy_rc1/  earlier figures; five are now supplementary figures S1-S5
outputs/rc1/         frozen release-candidate workbooks and manifests
references/          source bibliography; article PDFs are not redistributed
docs/                reproducibility notes
```

## Model summary

Critically ill adults without renal replacement therapy, continuous infusion,
ceftazidime and avibactam simulated as separate components and evaluated
jointly. 100,000 virtual subjects across five EKFC renal-function classes, 11
regimens, 11 MIC values.

| Parameter | Value | Source |
|---|---|---|
| CL ceftazidime | 5.0 x (EKFC/70)^0.70 L/h | Cojutti 2024 |
| CL avibactam | 5.9 x (EKFC/70)^0.89 L/h | Cojutti 2024 |
| Interindividual variability | 67.92% / 76.91% CV | Cojutti 2024 |
| Random-effect correlation | 0.94 | Cojutti 2024 |
| Unbound fractions | 0.85 ceftazidime, 0.92 avibactam | Cojutti 2024 |
| Ceftazidime target | fCss/MIC >= 4 | Cojutti 2024 |
| Avibactam target | fCss >= 4 mg/L (1-8 in sensitivity) | continuous-infusion TDM target; 1 mg/L is the registrational threshold |
| Exposure screen | total Css > 104 mg/L | Cojutti 2024 |

## Provenance policy

Inputs are labelled as directly reported, inherited, donor-derived,
user-specified, model-inferred, or scenario-generated. Missingness is retained
as an evidence state rather than imputed, and donor evidence is restricted to
sensitivity roles.

## Manuscript

The manuscript and supplementary workbook are not included here while they are
in preparation. They will be added, with a persistent DOI, on publication.

## Citation

Cite the manuscript and the source articles listed in
`references/source_articles.md`.

## License

Code under the MIT License; data and documentation under CC BY 4.0 unless
otherwise stated. Source article PDFs are not included.
