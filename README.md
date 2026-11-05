# K-CTO Registry 2.0 — sex-difference analysis of CTO PCI: analysis code

R scripts underlying "Sex Differences in Clinical Outcomes and Procedural Intensity After Chronic Total Occlusion PCI: K-CTO Registry" (JACC: Advances, manuscript JADV102226-3033-T).

## Contents
| Script | Purpose |
|---|---|
| `01_clean_raw.R` | Raw registry workbook → cleaned patient-level file (variable dictionary, flags) |
| `00_cohort_prep.R` | Analysis cohort (CTO PCI, index year 2014–2023; N = 6,777), anonymised center index, helper functions |
| `02_descriptive_by_sex.R` | Baseline, lesion, procedural and in-hospital tables by sex |
| `03_derive_endpoints.R` | Outcome derivation (MACCE/MACE/NACE as in the parent report; hard composite of death, MI, stroke) |
| `04_primary_analysis.R` | Cox models (center-clustered SE), staged adjustment, age strata, procedural logistic models |
| `05_qa_sensitivity.R` | Sensitivity battery incl. the procedural intensity score engine (`pis_models()`), center models, body-size analyses |
| `06_intensity_robustness.R` | Score construction sensitivity, leave-one-out, decomposition |
| `09_interaction_tests.R` | Sex × lesion-complexity interaction tests |
| `12_v2_macce_landmark.R`, `13_v2_landmark_remedies.R`, `14_v2_landmark_extras.R`, `15_v2_age_stratified_macce.R`, `16_v2_macce_corroborating.R` | Landmark-partitioned analysis of MACCE (0–1, 1–3, >3 years), RMST, remedies, age strata, corroborating estimators (IPTW, entropy balancing, PSM) |
| `17_phase7_reviewer_analyses.R` (+ `18_…` if present) | Analyses added in response to peer review (treated vs untreated registry comparison, single-center follow-up sensitivity and tipping point, bilateral injection eligibility restriction, standardized absolute risks, component cumulative incidence) |
| `07_figures.R`, `10_central_illustration_jacc.R`, `11_figures_pptx.R`, `_theme.R` | Figures and Central Illustration |

## Running
1. Set the project root: `Sys.setenv(KCTO_PROJECT_ROOT = "/path/to/project")` (scripts expect `data/raw/`, `data/processed/`, `output/tables/`, `output/figures/` under it).
2. Run in numeric order; `00_cohort_prep.R` is sourced by the analysis scripts.
3. R ≥ 4.5; packages: survival, lme4, ordinal, MASS, mice, cmprsk/tidycmprsk, WeightIt/cobalt, ggplot2, patchwork, ragg, flextable, readxl, dplyr, tidyr (see `sessionInfo.txt`).

## Data
Patient-level data are held by the K-CTO registry (Korean Society of Cardiology) and cannot be shared publicly; access requests are handled by the registry steering committee. The scripts contain no patient data. Identifier columns (`id`, `id_code`, `patient_key`, `Center_name`) are never printed; centers appear only as an anonymised index.

## Reproducibility note
Parts of the code were drafted with AI assistance (disclosed in the manuscript). For the revision the full pipeline was re-run end-to-end and all published estimates reproduced (`17_phase7_reviewer_analyses.R` contains the assertions).
