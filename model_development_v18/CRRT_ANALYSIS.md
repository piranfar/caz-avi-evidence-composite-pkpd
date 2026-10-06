# CRRT analysis

The ceftazidime–avibactam clearance correlation during continuous renal replacement therapy (CRRT), and
what it means for monitoring that measures ceftazidime only. The analysis is built on Model 1
(`MODEL1_REPORT.md`, section 0).

## Questions

1. How strongly do the clearances of the two drugs correlate between patients during CRRT, and how much does
   the estimate depend on the covariance structure, the estimation method and the coding of the deposited data?
2. How much of each clearance does the circuit provide, and how much of the variability of avibactam clearance
   does ceftazidime clearance leave unexplained?
3. What joint target attainment do common regimens give in a virtual CRRT population?
4. How often does classifying avibactam attainment from one measured ceftazidime trough wrongly
   reassure, and how many avibactam measurements does a ceftazidime trough replace?

## Data

| Source | Used for | Location |
|---|---|---|
| Li et al, Dryad doi:10.5061/dryad.fxpnvx16s (CC0): 21 critically ill adults on CRRT, 2 g/0.5 g every 8 h, 119 prefilter and 119 postfilter concentrations per drug | Model 1 (prefilter), circuit clearance (both sites), virtual population, cohort checks | `data_external/dryad_Li2025_CRRT/` |
| Gatti et al, J Crit Care 2023;76:154301, Table 2: published clearances, 17 occasions in 8 patients on CVVHDF | second cohort | `data_external/Gatti2023_individual_patient_data.csv` |
| Cojutti et al, J Antimicrob Chemother 2024;79:2801-8 | ρ = 0.94 and between-patient variability without renal replacement therapy | `data_external/Cojutti2024_ANCHOR_PopPK/`, `src/cazavi/` |
| O'Jeanson et al, Int J Antimicrob Agents 2025;65:107394 (CC BY) | avibactam unbound fraction 0.73 during CVVHDF (sensitivity analysis) | `data_external/OJeanson2024_CVVHDF/` |

## Definitions

- Targets: free ceftazidime trough at least 4 × MIC; free avibactam trough at least 4 mg/L. Also recorded:
  free avibactam above 1 mg/L for 50% of the interval (registrational target), and total ceftazidime above
  104 mg/L (exposure screen).
- Unbound fractions 0.85 (ceftazidime) and 0.92 (avibactam); 0.73 for avibactam as a sensitivity analysis.
- Regimens, 2-h infusions: 2.5 g every 8 h, 2.5 g every 12 h, 1.25 g every 8 h, 0.94 g every 12 h;
  100,000 virtual patients each.
- Classifier: probability that avibactam attains its target given one ceftazidime trough measured with
  the model's residual error, computed from Model 1; attaining if at least 0.5. "Wrongly
  reassured": classified as attaining with the avibactam trough below target, as a percentage of all
  patients.
- Circuit clearance, with the formulas of Li et al: CL_CRRT = Qp − Qout × AUCpost / AUCpre over the sampled
  interval, Qp = 160 mL/min × (1 − hematocrit), Qout = Qp − Quf (CVVHD) or Qp − Quf − QR (postdilution CVVH),
  QR = 2 L/h. Net ultrafiltration Quf is deposited only as a category: 0 L/h, with 0.1 and 0.2 L/h as
  sensitivity analyses. Clearance outside the circuit = empirical Bayes total − circuit clearance.
- Avibactam clearance variability not explained by ceftazidime clearance: the CV of ω_avi × √(1 − ρ²).
- Independent re-estimation (step 14): nlmixr2, first-order conditional estimation with interaction, with the
  residual correlation written as the conditional likelihood of avibactam given ceftazidime in the same sample;
  every fit starts from the same neutral values and is restarted from its own estimates until the OFV changes by
  less than 0.001. The closed-form concentrations are checked against an ODE solution with steady-state dosing.
- Measurement error in circuit and other clearance (step 18): every prefilter and postfilter concentration gets
  rounding noise of ±0.25 mg/L and a log-normal error with Model 1's residual SD, correlated between the 2 drugs
  in the same tube; 2,000 replicates, and half the residual SD as a sensitivity analysis (the residual SD also
  contains model misfit, so it bounds the measurement error from above).
- Calibrated interval (step 17): the cutoff is the 95th percentile of the OFV change at the true ρ in the 500
  datasets simulated at the estimate (step 15), in place of 3.841.
- Weight category (step 16): the deposited category w (0 to 3) multiplies clearance by exp(b_CL (w − 1.5)) and
  volume by exp(b_V (w − 1.5)) for both drugs.
- Covariance structures: primary (ρ and the volume correlation), Ka (adds the clearance-volume correlation within
  each drug), Kb (all 6 correlations, parameterized by C-vine partial correlations so that every parameter
  vector gives a valid matrix).

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
| 8 | `crrt_circuit_clearance.py` | `crrt_extracorporeal_components.csv`, `crrt_extracorporeal_summary.csv`, `crrt_conditional_avibactam.csv`, `model1_parameter_ci.csv` | under 1 min |
| 9 | `model1_structural_refits.py` | `model1_structural_refits.csv`, `model1_rounding_sensitivity.csv`, `model1_vpc_by_infusion.csv` | about 30 min |
| 10 | `model1_covariance_structure.py` | `model1_covariance_structure.csv`, `model1_covariance_profile.csv` | about 30 min |
| 11 | `crrt_covariance_sensitivity.py` | `crrt_covariance_sensitivity.csv` | about 3 min |
| 12 | `model1_vpc_structures.py` | `model1_vpc_structures.csv` | under 1 min |
| 13 | `model1_bayes_numpyro.py` | `model1_bayes_summary.csv` | about 1 min; needs NumPyro 0.21 and JAX 0.11 (`pip install numpyro==0.21.0 jax==0.11.1`) |
| 14 | `Rscript model1_nlmixr2.R` | `model1_nlmixr2.csv` | about 6 min; needs R 4.6.1 with nlmixr2 7.0.1 (nlmixr2est 7.1.0, rxode2 5.1.7) and Rtools45 (`install.packages("nlmixr2")`) |
| 15 | `model1_sbc_coverage.py --reps 500` | `model1_sbc_coverage_replicates.csv`, `model1_sbc_coverage_summary.csv` | about 7 h with 14 processes; resumes from the stored replicates |
| 16 | `model1_weight_covariate.py` | `model1_weight_covariate.csv` | about 10 min |
| 17 | `model1_profile_checks.py` | `model1_calibrated_interval.csv`, `model1_infusion_profile.csv` | about 10 min |
| 18 | `crrt_circuit_uncertainty.py` | `crrt_circuit_uncertainty.csv`, `crrt_circuit_uncertainty_summary.csv` | under 1 min |
| 19 | `crrt_unexplained_variability.py` | `crrt_unexplained_variability.csv` | under 1 min |
| 20 | `make_tdm_figures.py` | `figures/TDM_Figure1_joint_model.pdf`/`.tif`, `figures/TDM_Figure2_virtual_crrt.pdf`/`.tif` | under 1 min |
| 20 | `make_model1_figures_600dpi.py` | `figures/Figure7_model1_goodness_of_fit.pdf`/`.png`, `figures/Figure8_model1_visual_predictive_check.pdf`/`.png` | under 1 min |
| | `test_model1.py` | 140 checks; no output files | under 1 min |

Steps 11, 12 and 16 read the outputs of step 10, step 17 those of steps 9 and 15, steps 18 and 19 those of steps 1, 6 and 8, and step 20 those of steps 8, 10, 11 and 18. Run times are for a
16-thread Windows 11 machine. Scripts that fit many models use parallel processes; the environment variables
`SIM_WORKERS` (steps 3, 6 and 15) and `REFIT_WORKERS` (steps 9, 10, 16 and 17) set how many. Results do not depend on
them, because every replicate, parameter draw and fit has its own random stream or starting point.

Run with Python 3.14.7, NumPy 2.5.0, SciPy 1.18.0 and Matplotlib 3.11.1. Rerunning steps 1 to 9, 11, 12, 14 and 16 to 20
reproduces their CSVs byte for byte (step 3 from its stored replicates), except the `seconds` columns, which
record run time. The figures are identical except that each PDF carries its own creation date. Logs of the
rerun of steps 1 to 7 are in `audit/log_<step>.txt`.

## Random streams

| Stream | Seed |
|---|---|
| Visual predictive checks (overall, by infusion category, by covariance structure); estimator checks (`model1_sbc.py`, `model1_sbc_coverage.py`, same stream per replicate) | 20260811 |
| Virtual population, regimens 2.5 g every 8 h, 1.25 g every 8 h, 0.94 g every 12 h | SEED = 20261002 |
| Regimen 2.5 g every 12 h | SEED + 1 |
| Patient-cluster bootstrap of the second cohort (20,000 resamples) | SEED + 2 |
| Parameter-uncertainty draws (500 draws × 20,000 patients) | SEED + 3, one stream per draw |
| Avibactam unbound fraction 0.73 | SEED + 4 |
| Avibactam trough given a measured ceftazidime trough | SEED + 5 |
| Monitoring under the full covariance structure | SEED + 6 |
| Measurement error in circuit and other clearance (2,000 + 2,000 replicates per SD) | SEED + 7 |
| Monte Carlo replication of the classifier | SEED + 10 to SEED + 14 |
| Rounding noise (20 datasets); Hamiltonian Monte Carlo | 20261005 |

Values are written with 6 decimals and rounded once, when reported.
