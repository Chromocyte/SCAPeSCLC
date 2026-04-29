############################################################
# SCAPeSCLC – Bayesian Pathway Cox Proportional Hazards
#
# Description:
# This script performs Cox proportional hazards modeling
# using Bayesian posterior pathway scores as covariates.
#
# Both unadjusted and adjusted models are generated:
#
# Unadjusted model:
#   Surv(time, event) ~ pathway_posterior
#
# Adjusted model:
#   Surv(time, event) ~ pathway_posterior
#                        + treatment
#                        + bone_metastasis
#
# Input Requirements:
# - CSV file derived from D13 dataset
# - Posterior columns located at:
#     22, 25, 28, ..., 337
#   (Every 3rd column starting at 22)
#
# Column Structure:
# Each pathway has:
#   Posterior
#   Lower CI
#   Upper CI
#
# Only Posterior columns are used.
#
# Required Clinical Columns:
# OS, DSS, PFS
# OS_event, DSS_event, PFS_event
# treatment
# bone_metastasis
#
# Output:
# Six CSV files:
#   BP_OS_unadjusted.csv
#   BP_DSS_unadjusted.csv
#   BP_PFS_unadjusted.csv
#   BP_OS_adjusted.csv
#   BP_DSS_adjusted.csv
#   BP_PFS_adjusted.csv
#
# Dependencies:
# dplyr, survival, broom, purrr, stringr
############################################################

library(dplyr)
library(survival)
library(broom)
library(purrr)
library(stringr)

# =========================
# User Input Section
# =========================

file_path <- "D13_patient_BP_posteriors.csv.csv"

# =========================
# Load Data
# =========================

df <- read.csv(
  file = file_path,
  stringsAsFactors = FALSE
)

# =========================
# Sanitize Column Names
# =========================

colnames(df) <- colnames(df) %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("&", "and") %>%
  str_replace_all("-", "_") %>%
  str_replace_all("/", "_") %>%
  str_replace_all("\\.", "_")

# =========================
# Identify Posterior Columns
# =========================

posterior_indices <- seq(
  from = 22,
  to = 337,
  by = 3
)

bp_cols <- colnames(df)[posterior_indices]

# =========================
# Cox Model Function
# =========================

run_cox_models <- function(
    pathway,
    data,
    time_var,
    event_var,
    adjusted = FALSE
) {
  
  if (adjusted) {
    
    vars_needed <- c(
      pathway,
      time_var,
      event_var,
      "treatment",
      "bone_metastasis"
    )
    
    formula_text <- paste0(
      "Surv(",
      time_var,
      ", ",
      event_var,
      ") ~ ",
      pathway,
      " + treatment + bone_metastasis"
    )
    
  } else {
    
    vars_needed <- c(
      pathway,
      time_var,
      event_var
    )
    
    formula_text <- paste0(
      "Surv(",
      time_var,
      ", ",
      event_var,
      ") ~ ",
      pathway
    )
    
  }
  
  data_subset <- data[, vars_needed]
  
  surv_formula <- as.formula(formula_text)
  
  model <- survival::coxph(
    surv_formula,
    data = data_subset,
    ties = "efron"
  )
  
  result <- broom::tidy(
    model,
    exponentiate = TRUE,
    conf.int = TRUE
  )
  
  result <- result %>%
    filter(term == pathway) %>%
    mutate(
      pathway_tested = pathway,
      n = nrow(data_subset)
    )
  
  return(result)
  
}

# =========================
# Wrapper Function
# =========================

run_survival_set <- function(
    time_var,
    event_var,
    adjusted_flag
) {
  
  message(
    "Running ",
    ifelse(adjusted_flag,
           "adjusted",
           "unadjusted"),
    " pathway model for ",
    time_var
  )
  
  results <- map_dfr(
    bp_cols,
    run_cox_models,
    data = df,
    time_var = time_var,
    event_var = event_var,
    adjusted = adjusted_flag
  )
  
  results <- results %>%
    mutate(
      p_adj = p.adjust(
        p.value,
        method = "BH"
      )
    )
  
  return(results)
  
}

# =========================
# Run Unadjusted Models
# =========================

results_OS_unadj <- run_survival_set(
  "OS",
  "OS_event",
  adjusted_flag = FALSE
)

results_DSS_unadj <- run_survival_set(
  "DSS",
  "DSS_event",
  adjusted_flag = FALSE
)

results_PFS_unadj <- run_survival_set(
  "PFS",
  "PFS_event",
  adjusted_flag = FALSE
)

# =========================
# Run Adjusted Models
# =========================

results_OS_adj <- run_survival_set(
  "OS",
  "OS_event",
  adjusted_flag = TRUE
)

results_DSS_adj <- run_survival_set(
  "DSS",
  "DSS_event",
  adjusted_flag = TRUE
)

results_PFS_adj <- run_survival_set(
  "PFS",
  "PFS_event",
  adjusted_flag = TRUE
)

# =========================
# Export Results
# =========================

write.csv(
  results_OS_unadj,
  "BP_OS_unadjusted.csv",
  row.names = FALSE
)

write.csv(
  results_DSS_unadj,
  "BP_DSS_unadjusted.csv",
  row.names = FALSE
)

write.csv(
  results_PFS_unadj,
  "BP_PFS_unadjusted.csv",
  row.names = FALSE
)

write.csv(
  results_OS_adj,
  "BP_OS_adjusted.csv",
  row.names = FALSE
)

write.csv(
  results_DSS_adj,
  "BP_DSS_adjusted.csv",
  row.names = FALSE
)

write.csv(
  results_PFS_adj,
  "BP_PFS_adjusted.csv",
  row.names = FALSE
)

message("All Bayesian pathway survival analyses completed.")