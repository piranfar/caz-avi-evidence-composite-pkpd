# Final manuscript artifacts

The following two files are the authoritative publication artifacts for this repository:

| Artifact | Canonical repository path | SHA-256 | Status |
|---|---|---|---|
| Manuscript | `manuscript/Draft_V5_Revised.docx` | `ed7f835742454e4e01ab3f36208b30a081eea345b9308e29c686b80eefcd3253` | FINAL |
| Supplementary Tables S1-S8 | `supplement/Supplementary_Tables_S1-S8_V4.xlsx` | `b3a1030d689a5e4c18c597f5601e6e5168ec039a32bf9663277f65c99c4bedeb` | FINAL |

These files supersede earlier manuscript, PDF, and supplementary-table versions when resolving scientific wording, parameter definitions, provenance labels, and release documentation.

## Authority rule

1. The final manuscript and Supplementary Tables S1-S8 V4 define the publication-facing specification.
2. Executable model configurations must reproduce the reported outputs while documenting any parameter-scale conversion explicitly.
3. Earlier RC1 workbooks remain immutable reference outputs and calibration anchors; they do not override the final publication artifacts silently.
4. Any implementation difference from the final artifacts must be recorded as a versioned correction or an explicit reproducibility exception.

## Final publication-facing inputs

- `CL_CAZ = 5.0 * (EKFC / 70)^0.70 L/h`
- `CL_AVI = 5.9 * (EKFC / 70)^0.89 L/h`
- CAZ clearance IIV: `67.92% CV`
- AVI clearance IIV: `76.91% CV`
- CAZ-AVI clearance random-effect correlation: `rho = 0.94`
- CAZ unbound fraction: approximately `0.90`
- AVI unbound fraction: approximately `1.0`
- CAZ target: `fCss / MIC >= 4`
- AVI target: `fCss / 4 mg/L >= 1`
- Toxicity screen: total CAZ `Css > 104 mg/L`
- Primary population: 100,000 virtual subjects, 20,000 per each of five EKFC classes
- Primary seed: `20260707`

## Required parameter mapping

The final supplement reports IIV as coefficient of variation. The simulation engine uses log-normal standard deviations:

`omega = sqrt(log(1 + (CV/100)^2))`

This gives approximately:

- CAZ: `CV 67.92% -> omega 0.620`
- AVI: `CV 76.91% -> omega 0.680`

This is a `PARAMETER_MAPPING`, not a change to the reported evidence.

## File inclusion

The binary final artifacts should be committed at the canonical paths above before the repository is tagged for release. Their hashes must match this document and the release manifest.