# SCAPeSCLC

SCAPeSCLC is a harmonized analytical framework developed from publicly available datasets GSE261345 and GSE261348. These datasets originate from the CANTABRICO and IMfirst cohorts of patients with extensive-stage small cell lung cancer (ES-SCLC).

This repository contains selected analysis scripts and supporting datasets used to generate Bayesian pathway posterior estimates and perform gene-level and pathway-level survival modeling.

---

## Repository Structure

data/
Core datasets used for analysis.

- D5_scaled_gene_expression.csv  
  Patient-level scaled gene expression values.

- D10_ROI_CTA_Zscores.csv  
  ROI-level CTA pathway enrichment Z-scores.

- D13_patient_BP_posteriors.csv  
  Patient-level Bayesian posterior pathway scores.

scripts/
Analysis scripts.

- 01_gene_level_cox_models.R  
  Performs gene-level Cox proportional hazards modeling.

- 01_gene_level_cox_ph_assumptions.R  
  Performs gene-level Cox proportional hazards assumption testing.

- 02_bayesian_patient_level_pathways.R  
  Generates patient-level Bayesian posterior pathway estimates.

- 03_pathway_posterior_cox_models.R  
  Performs pathway-level Cox proportional hazards modeling.

- 03_pathway_posterior_cox_ph_assumptions.R  
  Performs pathway-level Cox proportional hazards assumption testing.
---

## Requirements

R version ≥ 4.2 recommended.

Required packages:

- dplyr  
- survival  
- broom  
- purrr  
- brms  
- tidyr  
- stringr  

Install packages in R:

```r
install.packages(c(
"dplyr",
"survival",
"broom",
"purrr",
"brms",
"tidyr",
"stringr"
))
