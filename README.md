# SCAPeSCLC

**SCAPeSCLC** is a harmonized multi-level transcriptomic and clinical resource derived from the publicly available GEO datasets **GSE261345** and **GSE261348**, originating from the CANTABRICO and IMfirst cohorts of patients with extensive-stage small cell lung cancer (ES-SCLC).

This repository contains the R scripts and supporting datasets used to generate Bayesian pathway posterior estimates, perform gene- and pathway-level survival analyses, assess Cox proportional hazards model assumptions, and generate comprehensive diagnostic atlases.

The repository accompanies the published SCAPeSCLC dataset and data paper.

---

## Associated Resources

- **Zenodo dataset:** https://doi.org/10.5281/zenodo/21523665
- **Data paper:**  
  Shirvaliloo, M. *SCAPeSCLC: An Integrated Spatial Transcriptomic and Bayesian Pathway Enrichment Dataset for Survival Modeling in Extensive-Stage Small Cell Lung Cancer.* **Data** **2026**, *11*, 152. https://doi.org/10.3390/data11070152

---

# Repository Structure

```
SCAPeSCLC
├── data/
│   ├── D5_scaled_gene_expression.csv
│   ├── D10_ROI_CTA_Zscores.csv
│   └── D13_patient_BP_posteriors.csv
│
└── scripts/
    ├── 01_gene_level_cox_models.R
    ├── 01_gene_level_cox_ph_assumptions.R
    ├── 02_bayesian_patient_level_pathways.R
    ├── 03_pathway_posterior_cox_models.R
    ├── 03_pathway_posterior_cox_ph_assumptions.R
    ├── 04_SCAPeSCLC_diagnostic_atlas_generator_for_genes.R
    └── 04_SCAPeSCLC_diagnostic_atlas_generator_for_genes_confounder_adjusted.R
```

---

## Data Files

| File | Description |
|------|-------------|
| **D5_scaled_gene_expression.csv** | Patient-level standardized gene expression matrix. |
| **D10_ROI_CTA_Zscores.csv** | ROI-level Cancer Transcriptome Atlas pathway enrichment Z-scores. |
| **D13_patient_BP_posteriors.csv** | Patient-level Bayesian posterior pathway activity estimates. |

---

## Analysis Workflow

The repository implements the following analytical workflow:

1. Bayesian estimation of patient-level pathway activities from ROI-level Cancer Transcriptome Atlas pathway enrichment scores.
2. Gene-level Cox proportional hazards regression for progression-free, disease-specific, and overall survival.
3. Pathway-level Cox proportional hazards regression using Bayesian posterior pathway activities.
4. Assessment of proportional hazards assumptions using Schoenfeld residuals.
5. Generation of comprehensive diagnostic atlases including:
   - hazard ratios and confidence intervals,
   - Wald test statistics,
   - proportional hazards test results,
   - Martingale residuals,
   - Schoenfeld residuals,
   - Deviance residuals,
   - DFBETA influence diagnostics.

---

## Requirements

R **4.2** or later is recommended.

Required packages:

```r
install.packages(c(
  "dplyr",
  "survival",
  "broom",
  "purrr",
  "brms",
  "tidyr",
  "stringr",
  "ggplot2",
  "patchwork",
  "cowplot",
  "gtable"
))
```

---

## Citation

If you use SCAPeSCLC in your work, please cite both the Zenodo dataset and the accompanying data paper.

---

## License

This project is distributed under the **MIT License**.  Performs pathway-level Cox proportional hazards assumption testing.
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
