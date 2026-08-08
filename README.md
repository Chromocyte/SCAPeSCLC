# SCAPeSCLC

**SCAPeSCLC** is a harmonized multi-level transcriptomic and clinical resource derived from the publicly available GEO datasets **GSE261345** and **GSE261348**, originating from the CANTABRICO and IMfirst cohorts of patients with extensive-stage small cell lung cancer (ES-SCLC).

![SCAPeSCLC workflow](figures/SCAPeSCLC_pipeline.png)

This repository contains the R scripts and supporting datasets used to generate Bayesian pathway posterior estimates, perform gene- and pathway-level survival analyses, assess Cox proportional hazards model assumptions, and generate comprehensive diagnostic atlases for both gene expression and biological pathway activity.

The repository accompanies **SCAPeSCLC v1.3.5**, the published dataset and associated data paper.

---

## Current Release (v1.3.5)

Major additions include:

- Gene-level diagnostic atlases (unadjusted and confounder-adjusted)
- Biological pathway diagnostic atlases (unadjusted and confounder-adjusted)
- Proportional hazards assumption testing for all gene and pathway Cox models
- Bayesian estimation of patient-level Cancer Transcriptome Atlas pathway activities

---

## Associated Resources

- **Zenodo dataset:** https://doi.org/10.5281/zenodo.19897644
- **Data paper:**  
  Shirvaliloo M. SCAPeSCLC: An Integrated Spatial Transcriptomic and Bayesian Pathway Enrichment Dataset for Survival Modeling in Extensive-Stage Small Cell Lung Cancer. *Data*. 2026;11(7):152. https://doi.org/10.3390/data11070152

![SCAPeSCLC Summary](figures/SCAPeSCLC_summary_expanded.png)

---

# Repository Structure

```text
SCAPeSCLC
├── data/
│   ├── D5_scaled_gene_expression.csv
│   ├── D10_ROI_CTA_Zscores.csv
│   └── D13_patient_BP_posteriors.csv
│
├── figures/
│   └── SCAPeSCLC_pipeline.png
│
└── scripts/
    ├── 01_gene_level_cox_models.R
    ├── 01_gene_level_cox_ph_assumptions.R
    ├── 02_bayesian_patient_level_pathways.R
    ├── 03_pathway_posterior_cox_models.R
    ├── 03_pathway_posterior_cox_ph_assumptions.R
    ├── 04_SCAPeSCLC_diagnostic_atlas_generator_for_genes.R
    ├── 04_SCAPeSCLC_diagnostic_atlas_generator_for_genes_confounder_adjusted.R
    ├── 05_SCAPeSCLC_diagnostic_atlas_generator_for_BPs.R
    └── 05_SCAPeSCLC_diagnostic_atlas_generator_for_BPs_confounder_adjusted.R
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
3. Pathway-level Cox proportional hazards regression using Bayesian posterior pathway activity estimates.
4. Assessment of proportional hazards assumptions using Schoenfeld residuals for both gene- and pathway-level Cox regression models.
5. Generation of comprehensive diagnostic atlases for both genes and biological pathways, including:
   - model summary statistics,
   - hazard ratios and 95% confidence intervals,
   - Wald test statistics,
   - proportional hazards test results,
   - Martingale residuals,
   - Schoenfeld residuals,
   - Deviance residuals,
   - DFBETA influence diagnostics.
6. Generation of both unadjusted and confounder-adjusted diagnostic atlases for all survival endpoints (OS, DSS, and PFS).

---

## Data Provenance and Analytical Workflow

The diagram below summarizes the provenance of the major data products included in SCAPeSCLC and their relationships to the analytical pipelines implemented in this repository.

```mermaid
flowchart TD

    A["GEO Datasets<br/>GSE261345 & GSE261348"]

    A --> B["Clinical and Transcriptomic Harmonization"]

    B --> C["Patient-Level Gene Expression<br/>(D5_scaled_gene_expression.csv)"]
    B --> D["ROI-Level CTA Pathway Z-scores<br/>(D10_ROI_CTA_Zscores.csv)"]

    D --> E["Bayesian Pathway Activity Estimation"]
    E --> F["Patient-Level Bayesian Pathway Posteriors<br/>(D13_patient_BP_posteriors.csv)"]

    C --> G["Gene-Level Cox Proportional Hazards Models"]
    G --> H["Gene PH Assumption Testing"]
    H --> I["Gene Diagnostic Atlases"]

    F --> J["Pathway-Level Cox Proportional Hazards Models"]
    J --> K["Pathway PH Assumption Testing"]
    K --> L["Biological Pathway Diagnostic Atlases"]
```
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

If you use **SCAPeSCLC** in your work, please cite both the Zenodo dataset and the accompanying data paper.

---

## License

This project is distributed under the **MIT License**.
