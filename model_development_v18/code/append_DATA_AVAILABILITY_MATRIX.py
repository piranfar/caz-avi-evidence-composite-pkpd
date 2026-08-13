import csv, sys
sys.stdout.reconfigure(encoding="utf-8")
P = "DATA_AVAILABILITY_MATRIX.csv"
rows = list(csv.reader(open(P, encoding="utf-8")))
H = rows[0]

def mk(**k):
    return [k.get(c, "NOT ASSESSED") for c in H]

new = [
 mk(study_id="S7_LiDryad2025_CRRT",
    citation="Li C, Wang Y, Chen F, Huang L, Dong J, Fan W, Yue H, Ge Y. Dryad dataset, 2025; primary article Antimicrob Agents Chemother 2026;70(2):e0143825",
    doi="10.5061/dryad.fxpnvx16s (dataset); 10.1128/aac.01438-25 (article)", pmid="41432444 (article)", pmcid="NOT ASSESSED",
    open_access="YES - dataset is public domain", license="CC0 1.0 Universal (public domain dedication)",
    population="Critically ill adults on CRRT; 18 of 21 with acute pancreatitis (China)", n_subjects="21",
    infusion_mode="INTERMITTENT intravenous infusion, 2 g + 0.5 g every 8 h - NOT continuous infusion",
    anonymised_patient_id="YES - subjectID 1 to 21",
    concentration_observations="YES - 118 ceftazidime and 118 avibactam observations, 5 to 7 timepoints per patient, PRE- and POST-FILTER",
    sampling_times="YES - time after dose in hours, per observation",
    dose_history="YES - uniform 2 g + 0.5 g q8h for all patients; infusion duration published only as a category",
    infusion_start_stop="PARTIAL - time after dose is given; infusion duration is categorical",
    plasma_conc="YES - individual, per timepoint, both analytes", elf_conc="NOT APPLICABLE",
    renal_function="PARTIAL - serum creatinine and 24-h urine volume per patient; CRRT settings as categories",
    body_weight="CODED - 4 categories, not raw kilograms",
    age_sex="CODED - age in 3 categories; sex as a code",
    rrt_ecmo_status="YES - all on CRRT; modality (CVVH vs CVVHD) per patient",
    mic_data="PARTIAL - Kirby-Bauer zone diameters and S/R interpretation per patient; no numeric MIC",
    isolate_info="YES - per-patient pathogen, infection site, pathogen clearance, clinical outcome",
    blq_observations="NONE apparent", assay_error_lloq="NOT REPORTED in the dataset",
    clinical_outcomes="YES - per patient: pathogen clearance and clinical outcome",
    model_code="NOT PROVIDED (Phoenix used for non-compartmental analysis)",
    model_files_nonmem_monolix="NOT PROVIDED", covariance_matrix="NOT APPLICABLE - non-compartmental analysis",
    supplementary_individual_profiles="NOT APPLICABLE - the raw data themselves are provided",
    public_repository="YES - Dryad, https://datadryad.org/dataset/doi:10.5061/dryad.fxpnvx16s",
    data_availability_statement="Public domain deposit; README confirms explicit patient consent for public-domain publication of de-identified data",
    corresponding_author_route="Not required - the data are already public",
    classification="C (suitable for sensitivity and assumption-testing analysis only) - NOT A or B, because the population is on CRRT and receives intermittent rather than continuous infusion",
    evidence_and_url="Files read in full through the Dryad preview panels and transcribed to model_development_v18/data_external/dryad_Li2025_CRRT/. Row counts and file sizes match the Dryad manifest. THE ONLY OPENLY DOWNLOADABLE INDIVIDUAL-LEVEL PAIRED CAZ AND AVI CONCENTRATION-TIME DATASET FOUND ANYWHERE IN THIS AUDIT."),

 mk(study_id="S8_Li2019_CTS_registrational",
    citation="Li J, Lovern M, Green ML, et al. Clin Transl Sci 2019;12(2):151-163",
    doi="10.1111/cts.12585", pmid="NOT ASSESSED", pmcid="PMC6440567",
    open_access="YES", license="CC BY-NC",
    population="Phase 1 to 3 adults pooled across RECLAIM, RECAPTURE, REPRISE, REPROVE",
    n_subjects="1,975 (ceftazidime) / 2,249 (avibactam)",
    infusion_mode="2-hour intermittent intravenous infusion q8h",
    anonymised_patient_id="NO", concentration_observations="NO - 9,155 and 13,735 observations were modelled but none are published",
    sampling_times="NO", dose_history="NO", infusion_start_stop="NO", plasma_conc="NO", elf_conc="NOT APPLICABLE",
    renal_function="NO - covariate relationships published, not individual values",
    body_weight="NO", age_sex="NO", rrt_ecmo_status="NOT APPLICABLE", mic_data="NOT APPLICABLE",
    isolate_info="NOT APPLICABLE", blq_observations="NOT ASSESSED",
    assay_error_lloq="Residual error stratified by study phase, published",
    clinical_outcomes="NOT APPLICABLE",
    model_code="YES - actual NONMEM control streams for both final models, in Data S1 (CTS-12-151-s001.docx), retrievable from Europe PMC",
    model_files_nonmem_monolix="YES - $SIZES / $INPUT / $SUBROUTINE ADVAN3 TRANS4 / $PK with all covariate blocks and $OMEGA BLOCK(4) for each drug. The dataset itself is referenced but not included.",
    covariance_matrix="YES - FULL OMEGA covariance and correlation matrices published in Tables 1-2, PLUS bootstrap 90% confidence intervals for every theta in Table S4/S5",
    supplementary_individual_profiles="NO",
    public_repository="Supplementary files served by Europe PMC: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6440567/supplementaryFiles",
    data_availability_statement="NOT ASSESSED - individual data would fall under the Pfizer/Vivli controlled-access route",
    corresponding_author_route="Via Pfizer data sharing / Vivli",
    classification="D (aggregate evidence only) - but the RICHEST PUBLIC PARAMETER AND UNCERTAINTY PACKAGE available for either drug",
    evidence_and_url="Not currently cited by the manuscript. Contains the only public full OMEGA covariance matrices WITH bootstrap confidence intervals AND executable model code. Directly usable for structural-uncertainty analysis and for specifying parameter uncertainty. NOTE: it reports WITHIN-drug random-effect correlations only; no CAZ-to-AVI cross-drug covariance is quantified in any regulatory or published source."),
]

rows.extend(new)
with open(P, "w", newline="", encoding="utf-8") as fh:
    csv.writer(fh, lineterminator="\n").writerows(rows)
print(f"matrix now has {len(rows)-1} studies")
ci = H.index("classification")
for r in rows[1:]:
    print(f"  {r[0]:28} {r[ci][:70]}")
