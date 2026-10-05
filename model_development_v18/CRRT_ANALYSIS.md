# CRRT analysis

The ceftazidime–avibactam clearance correlation during continuous renal replacement therapy (CRRT), and
what it means for monitoring that measures ceftazidime only. The analysis is built on Model 1
(`MODEL1_REPORT.md`, section 0).

## Questions

1. How strongly do the clearances of the two drugs correlate between patients during CRRT?
2. What joint target attainment do common regimens give in a virtual CRRT population?
3. How often does classifying avibactam attainment from one measured ceftazidime trough wrongly
   reassure, and how many avibactam measurements does a ceftazidime trough replace?

## Data

| Source | Used for | Location |
|---|---|---|
| Li et al, Dryad doi:10.5061/dryad.fxpnvx16s (CC0): 21 critically ill adults on CRRT, 2 g/0.5 g every 8 h, 119 concentrations per drug | Model 1, virtual population, cohort checks | `data_external/dryad_Li2025_CRRT/` |
| Gatti et al, J Crit Care 2023;76:154301, Table 2: published clearances, 17 occasions in 8 patients on CVVHDF | second cohort | `data_external/Gatti2023_individual_patient_data.csv` |
| Cojutti et al, J Antimicrob Chemother 2024;79:2801-8 | comparator without renal replacement therapy (ρ = 0.94) | `data_external/Cojutti2024_ANCHOR_PopPK/`, `src/cazavi/` |

## Definitions

- Targets: free ceftazidime trough at least 4 × MIC; free avibactam trough at least 4 mg/L. Also recorded:
  free avibactam above 1 mg/L for 50% of the interval (registrational target), and total ceftazidime above
  104 mg/L (exposure screen).
- Unbound fractions 0.85 (ceftazidime) and 0.92 (avibactam); 0.73 for avibactam as a sensitivity analysis.
- Regimens, 2-h infusions: 2.5 g every 8 h, 2.5 g every 12 h, 1.25 g every 8 h, 0.94 g every 12 h;
  100,000 virtual patients each.
- Classifier: probability that avibactam attains its target given one ceftazidime trough measured with
  the model's proportional residual error, computed from Model 1; attaining if at least 0.5. "Wrongly
  reassured": classified as attaining with the avibactam trough below target, as a percentage of all
  patients.

## Run order

From `model_development_v18/code/`. Paths are resolved relative to each script.

| Step | Script | Writes (in `outputs/` unless stated) | Run time |
|---|---|---|---|
| 1 | `joint_popk_nlme.py` | `model1_joint_popk_parameters.csv`, `model1_parameter_covariance.csv`, `model1_individual_parameters.csv` | about 4 min |
| 2 | `model1_finalise.py` | `model1_final_parameters.csv`, `model1_profile_likelihood.csv`, `model1_diagnostics.csv`, `model1_vpc.csv`, `model1_sensitivity.csv`; `figures/model1_gof.png`, `model1_vpc.png` | about 23 min |
| 3 | `model1_sbc.py --reps 50` | `model1_sbc_replicates.csv`, `model1_sbc_summary.csv` | hours from scratch; resumes from the stored replicates |
| 4 | `crrt_virtual_tdm.py` | `crrt_attainment.csv`, `crrt_classifier.csv` | about 5 min |
| 5 | `crrt_supplementary_analyses.py` | `crrt_decision_threshold.csv`, `crrt_triage_curve.csv`, `crrt_equal_detection.csv`, `crrt_data_checks.csv`, `crrt_seed_replication.csv` | about 3 min |
| 6 | `crrt_parameter_uncertainty.py` | `crrt_parameter_uncertainty_draws.csv`, `crrt_parameter_uncertainty.csv`, `crrt_unbound_fraction_sensitivity.csv` | about 7 min |
| 7 | `table5_estimated_rho_rows.py` | `table5_estimated_rho_rows.csv` | under 1 min |
| 8 | `make_tdm_figures.py` | `figures/TDM_Figure1_joint_model.pdf`/`.tif`, `figures/TDM_Figure2_virtual_crrt.pdf`/`.tif` | under 1 min |
| 8 | `make_model1_figures_600dpi.py` | `figures/Figure7_model1_goodness_of_fit.pdf`/`.png`, `figures/Figure8_model1_visual_predictive_check.pdf`/`.png` | under 1 min |
| | `test_model1.py` | 140 checks; no output files | under 1 min |

Run times are for a 16-thread Windows 11 machine. `model1_sbc.py` and `crrt_parameter_uncertainty.py` use
parallel processes; the environment variable `SIM_WORKERS` sets how many (default: CPU count minus 2, at
most 14). Results do not depend on it, because every replicate and every parameter draw has its own
random stream.

Run with Python 3.14.7, NumPy 2.5.0, SciPy 1.18.0 and Matplotlib 3.11.1. Rerunning every step in a fresh clone (step 3 from its stored replicates) reproduces every CSV in the table above byte for byte. The figures are identical except that each PDF carries its own creation date. Logs of that rerun are
in `audit/log_<step>.txt`.

## Random streams

| Stream | Seed |
|---|---|
| Visual predictive check; estimator checks (`model1_sbc.py`) | 20260811 |
| Virtual population, regimens 2.5 g every 8 h, 1.25 g every 8 h, 0.94 g every 12 h | SEED = 20261002 |
| Regimen 2.5 g every 12 h | SEED + 1 |
| Patient-cluster bootstrap of the second cohort (20,000 resamples) | SEED + 2 |
| Parameter-uncertainty draws (500 draws × 20,000 patients) | SEED + 3, one stream per draw |
| Avibactam unbound fraction 0.73 | SEED + 4 |
| Monte Carlo replication of the classifier | SEED + 10 to SEED + 14 |

Values are written with 6 decimals and rounded once, when reported.
