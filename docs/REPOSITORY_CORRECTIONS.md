# Repository corrections derived from the final manuscript and supplement

This document records repository-level corrections required to make the computational release consistent with the authoritative final artifacts. The final DOCX and XLSX are not silently rewritten here.

## C1. Final free-fraction policy

**Correction:** the final publication specification is `fu_CAZ approximately 0.90` and `fu_AVI approximately 1.0`.

Earlier RC1 executable workbooks used 0.85 and 0.92. Those values may be retained only in an explicitly named legacy-reference configuration for comparison with historical outputs. They must not be used by the final scenario configuration.

**Repository action:** implemented in `config/final_primary_scenario.yaml`.

## C2. IIV scale mapping

The manuscript and Supplementary Table S3 report clearance IIV as CV percentages (67.92% and 76.91%). The simulation requires log-normal standard deviations. The conversion is:

`omega = sqrt(log(1 + (CV/100)^2))`

giving approximately 0.620 and 0.680. This is recorded as `PARAMETER_MAPPING`; neither representation should be mislabeled as a directly reported value on the other scale.

## C3. Main-text Table 4 column labels

The final manuscript Table 4 contains two numerical sensitivity columns, but the rendered heading is ambiguous (`Mean` and `rank correlation`). Based on Supplementary Table S6, the intended labels are:

- `Mean |rho|`
- `Maximum |rho|`

Repository-generated tables and figures must use these explicit labels.

## C4. Supplementary Table S7 is a release plan, not the current repository inventory

Supplementary Table S7 includes entries marked `Generated`, `Planned`, and `Pending`. Several listed files were not present in the newly initialized GitHub repository at the time the final supplement was supplied. Therefore:

- its historical checksums must not be presented as checksums of future regenerated files;
- every generated release artifact must receive a fresh SHA-256 checksum;
- the repository release manifest must distinguish `historical_supplement_record` from `current_repository_artifact`;
- no file may be marked generated or final until it exists and its checksum is verified.

## C5. Binary final artifacts

The final manuscript and supplement are identified by hashes in `docs/FINAL_ARTIFACTS.md`. They should be committed to:

- `manuscript/Draft_V5_Revised.docx`
- `supplement/Supplementary_Tables_S1-S8_V4.xlsx`

before tagging a release. Until those binaries are committed, the repository must not claim that the final publication files are hosted in the repository.

## C6. Primary versus sensitivity-only MIC distributions

The final Option A policy is:

- primary: `LEE2022_KPC_KP`, `LEE2022_OXA_KP`;
- sensitivity-only: `INDIA2022_OXA48_KP`, `FR2024_OXA484_SERINE_ONLY`.

Sensitivity-only donors must not be included in the primary CFR summary or counted as locked primary evidence.

## C7. Toxicity semantics

The threshold `total CAZ Css > 104 mg/L` and the 15% impermissibility cutoff remain provisional/user-specified screening rules. Code, tables, and documentation must not describe them as a validated ceftazidime-specific clinical neurotoxicity model.

## C8. Release gate

A release may be tagged only after:

1. the two final binary artifacts are committed and hashes verified;
2. the executable model regenerates the reported primary outputs under the final parameter policy;
3. any inability to reproduce historical RC1 outputs after the free-fraction correction is quantified;
4. current checksums, environment metadata, provenance registry, and tests are generated from the tagged commit;
5. publication-facing tables use corrected sensitivity-column labels.