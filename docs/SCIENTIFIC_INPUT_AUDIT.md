# Scientific Input Audit — RC1 Archive

## Purpose

This audit distinguishes archival outputs from executable model inputs and records unresolved discrepancies before the Monte Carlo engine is implemented.

## Package inventory

The uploaded archive contains:

- final manuscript PDF;
- Supplementary Tables S1–S8;
- six RC1 analytical workbooks;
- processed CSV summaries;
- figures 1–8;
- freeze manifest, runtime record, release-gate file, and checksums;
- lightweight scripts for figure regeneration and checksum verification.

The package is therefore a strong archival release candidate, but it is not yet a full computational reproduction package. The current scripts regenerate summary figures from processed outputs; they do not regenerate the virtual population, PTA, CFR, convergence, sensitivity, or PSA analyses from scientific inputs.

## Primary executable specification recovered from RC1 workbook

The `Methods` sheet of `First_Monte_Carlo_Run_100000.xlsx` records:

- population size: 100,000;
- five EKFC classes with 20,000 subjects per class;
- EKFC sampled uniformly within 0–30, 31–60, 61–90, 91–120, and 121–150 mL/min;
- typical CAZ clearance: `5.0 × (EKFC/70)^0.70` L/h;
- typical AVI clearance: `5.9 × (EKFC/70)^0.89` L/h;
- bivariate normal random effects on log clearance;
- omega CAZ: 0.62;
- omega AVI: 0.68;
- correlation: 0.94;
- individual clearance: typical clearance × `exp(eta)`;
- continuous-infusion concentration: dose rate / individual clearance;
- product dose split: 80% CAZ and 20% AVI;
- free fractions: CAZ 0.85 and AVI 0.92;
- CAZ target: `fCss/MIC ≥ 4`;
- AVI target: `fCss/4 mg/L ≥ 1`;
- joint target: both component targets attained;
- toxicity: total CAZ Css >104 mg/L;
- residual error: not applied in the primary latent-exposure PTA run;
- primary seed: 20260707.

## Provenance discrepancy requiring adjudication

Supplementary Table S3 and the RC1 executable workbook do not fully agree on all scientific inputs.

| Component | Supplementary S3 | RC1 executable workbook | Assessment |
|---|---:|---:|---|
| CAZ IIV | 67.92% CV | omega 0.62 | compatible after lognormal parameter mapping, within rounding |
| AVI IIV | 76.91% CV | omega 0.68 | compatible after lognormal parameter mapping, within rounding |
| CAZ free fraction | approximately 0.90 | 0.85 | unresolved material discrepancy |
| AVI free fraction | approximately 1.00 | 0.92 | unresolved material discrepancy |

The IIV pairs are not contradictory if the Supplementary values are CV percentages and the workbook values are log-scale standard deviations. This must be documented as `PARAMETER_MAPPING` rather than treated as two independent reported values.

The free-fraction discrepancy is material because it directly changes component-specific and joint target attainment. Exact RC1 reproduction should use the workbook values unless the manuscript/Supplementary package is revised or an explicit version policy selects different values.

## Freeze-state observations

The manifest identifies the model as `CAZ_AVI_EVIDENCE_COMPOSITE_V1_RC1` and records the freeze status as `TECHNICALLY_FROZEN_PENDING_USER_APPROVAL`.

Outstanding approvals include:

- donor MIC distributions;
- MIC-bin harmonization;
- toxicity-threshold policy;
- PK/PD-target policy;
- release-candidate approval as Model v1.

The manifest also preserves the broad French OXA-48-like 2021 MIC distribution as unresolved, nonblocking missingness.

## Implementation policy

1. Preserve the uploaded workbooks and manifest as immutable reference artifacts.
2. Treat workbook `Methods` values as the executable RC1 specification for initial reproduction.
3. Represent CV-to-omega conversion explicitly in provenance.
4. Do not silently reconcile the free-fraction discrepancy.
5. Build the core engine against `First_Monte_Carlo_Run_100000.xlsx` before implementing CFR, convergence, sensitivity, or PSA.
6. Require reference-output tests for regimen summaries and PTA-by-MIC results.
7. Generate new run metadata and checksums for every reproduction attempt.

## Next implementation unit

The next code change should add:

- typed population and regimen configuration;
- deterministic EKFC-class sampling;
- correlated CAZ/AVI log-clearance simulation;
- continuous-infusion CAZ/AVI exposure calculation;
- component-specific, joint PTA, and toxicity calculations;
- reference-output comparison against the RC1 workbook.
