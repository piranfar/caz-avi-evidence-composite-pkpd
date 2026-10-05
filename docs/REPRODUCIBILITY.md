# Reproducibility notes

## Execution context

The release-candidate analyses were executed in a managed Linux runtime during manuscript development (x86_64; Python 3.13.5 with NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 and python-docx 1.2.0). The runtime record listed in the release-candidate manifest (`CAZ_AVI_Runtime_Environment_v1_RC1.json`) is not distributed in this repository; the details above are taken from it. Its checksum is still listed in `outputs/rc1/CAZ_AVI_SHA256SUMS_v1_RC1.txt`, so that one line does not verify here.

The release-candidate manifest is stored in:

```text
outputs/rc1/CAZ_AVI_Model_Freeze_Manifest_v1_RC1.json
```

The original file-level checksums are stored in:

```text
outputs/rc1/CAZ_AVI_SHA256SUMS_v1_RC1.txt
```

## Determinism

The simulation is intended to be deterministic given identical inputs, random seed, and software/runtime versions. The primary seed was 20260707; multi-seed checks used 20260707 through 20260711.

## Scope

This repository preserves the release-candidate artifacts and provides lightweight scripts to regenerate the summary figures used in the manuscript. It does not claim to provide a full clinical dosing engine, a new population PK model, or a formal systematic review dataset.

## Validation boundary

Calibration against published anchors was treated as a provisional plausibility assessment, not as regulatory-grade validation or clinical dose validation. Prospective patient-level validation would require independent clinical concentration-time data.
