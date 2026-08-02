# Reproduction and verification report

**Scope.** Every analysis reported in the manuscript was reimplemented from the
published model description and checked against the frozen RC1 outputs. No
parameter was fitted, tuned, or back-calculated from those outputs.

## Result

All six analyses reproduce within Monte Carlo noise.

| Analysis | Agreement with RC1 |
|---|---|
| Primary 100,000-subject run | mean \|Δ joint PTA\| **0.060 pp**, max 0.320 (121 rows) |
| MIC-weighted CFR | mean \|Δ joint CFR\| **0.041 pp**, max 0.082 (33 rows) |
| Convergence at N = 50,000 | **0.169 pp** vs 0.167 recorded |
| Multi-seed (5 seeds) | reproduced |
| One-at-a-time sensitivity | same rank order; top driver 20.29 vs 19.75 recorded |
| Probabilistic sensitivity | **all 8 parameters in exact rank order**, mean \|r\| differences ≤ 0.002; robustness classification **20/20** identical |

The manuscript's Methods section is therefore sufficient to regenerate its own
results. That is a reproducibility claim the paper can now make truthfully.

## Parameters, and where each comes from

| Parameter | Value | Source |
|---|---|---|
| CL ceftazidime | 5.0 × (EKFC/70)^0.70 L/h | Cojutti 2024, Table 2 |
| CL avibactam | 5.9 × (EKFC/70)^0.89 L/h | Cojutti 2024, Table 2 |
| IIV | 67.92% / 76.91% CV → ω 0.6159 / 0.6817 | Cojutti 2024 |
| Random-effect correlation | ρ = 0.94 | Cojutti 2024 |
| Unbound fractions | **0.85 ceftazidime, 0.92 avibactam** | Cojutti 2024, verbatim |
| Ceftazidime target | fCss/MIC ≥ 4 | Cojutti 2024 |
| Avibactam target | fCss ≥ 4 mg/L | EUCAST fixed testing concentration |
| Exposure screen | total Css > 104 mg/L | Cojutti 2024 cohort |
| Renal function | uniform within each EKFC class | scenario definition |

The reproduction is itself the evidence for the unbound fractions: at 0.90 and
1.0 — the values currently in the manuscript — the model does not land on the
recorded outputs. At 0.85 and 0.92 it does.

## Two corrections this exercise forces

**Unbound fractions.** Methods §2.7 states ≈0.90 and ≈1.0. Both the input
specification and the source paper state 0.85 and 0.92, and only those values
reproduce the results. The manuscript text is wrong; the run is right.

**Calibration breadth.** The manuscript describes 72 comparison rows "spanning
all six source packages". The populated comparison is 9 regimens × 8 MIC values
against Cojutti Table 4 alone, which the project's own output specification
calls an *internal reproduction check*. It should be described as such.

## Findings the analyses support but the manuscript does not report

**The limiting component switches at the breakpoint.** Avibactam caps attainment
at MIC ≤ 4 mg/L. At the EUCAST clinical breakpoint of 8 mg/L, ceftazidime is
limiting in all 11 regimens and joint PTA falls to 10–68%. The source paper
states this directly; the manuscript never mentions the breakpoint.

**The dominant uncertainty is a choice, not a measurement.** Moving the
avibactam critical concentration to 6 or 2 mg/L shifts joint CFR by 20.3 and
17.8 percentage points — more than twice the effect of a 20% error in avibactam
clearance (8.3 pp). The largest single source of uncertainty in the model is an
analyst-specified threshold. The PSA excludes target thresholds by design; that
exclusion should be stated, because without it the ranking looks complete when
it is not.

**Permissibility is fragile.** A 20% reduction in ceftazidime clearance pushes
R8, R10 and R12 from 9–12% to 17–21% exposure-screen exceedance, across the 15%
ceiling. Raising IIV by 20% does the same to R10 and R12. Four of the five
selected regimens lose permissibility under one or more modest perturbations.
The manuscript reports only that all selected regimens sat below 15%.

**The efficacy–safety trade-off at low renal function.** R2 reaches 86.3% joint
CFR in the 0–30 class but is excluded at 16.2% exceedance, while the retained R1
reaches only 73.4%. In this class, adequate efficacy is only available by
crossing the safety ceiling. This is the most clinically pointed result in the
dataset and it is currently invisible.

## What is not reproduced

Exact historical reproduction is impossible: the original code and its
random-number stream were not preserved. Agreement here is statistical, not
bitwise. The recovered call order — uniform renal function, then correlated
random effects, class by class — reproduces the recorded tables far more closely
than an independent draw would, which suggests it matches the original
structure, but that cannot be confirmed.

The reconstruction package `CAZ_AVI_Local_First_Reconstruction_v1` does not
reproduce RC1 (mean \|Δ joint PTA\| 2.72 pp, max 9.68; several rows exceed the
paper's own 10 pp failure threshold). It replaces the published clearance
equations with per-class constants fitted to the RC1 percentiles and sets
ρ = 0.75. Its scaffolding, provenance documentation and frozen reference CSVs
are worth keeping; its model core should not be published as the analysis code.

## Files

```
reproduce_primary_run.py   core model, self-verifying
cazavi_analyses.py         all six analyses
make_figures.py            figures, from verified outputs only
outputs/                   analysis outputs + VERIFICATION_LOG.txt
figures/                   manuscript figures
```

Reproduce everything:

```
python cazavi_analyses.py all --verify
python make_figures.py
```
