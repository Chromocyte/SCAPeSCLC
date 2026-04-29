############################################################
# SCAPeSCLC – Gene-Level Cox Proportional Hazards Analysis
#
# Description:
# This script performs gene-level Cox proportional hazards
# modeling using scaled gene expression values from file D5.
#
# Both unadjusted and adjusted models are generated:
#
# Unadjusted model:
#   Surv(time, event) ~ gene_expression
#
# Adjusted model:
#   Surv(time, event) ~ gene_expression
#                        + treatment
#                        + bone_met
#
# Input Requirements:
# - CSV file derived from D5 dataset
# - Gene expression columns: 22–1743
# - Binary variables formatted as 1/0
#
# Required Columns:
# OS, DSS, PFS
# OS_event, DSS_event, PFS_event
# treatment
# bone_met
#
# Output:
# Six CSV files:
#   OS_unadjusted.csv
#   DSS_unadjusted.csv
#   PFS_unadjusted.csv
#   OS_adjusted.csv
#   DSS_adjusted.csv
#   PFS_adjusted.csv
#
# Dependencies:
# dplyr, survival, broom, purrr
############################################################

library(dplyr)
library(survival)
library(broom)
library(purrr)

# =========================
# User Input Section
# =========================

file_path <- "D5_scaled_gene_expression.csv.csv"

# =========================
# Load Data
# =========================

df <- read.csv(
  file = file_path,
  stringsAsFactors = FALSE
)

# =========================
# Identify Gene Columns
# =========================

gene_cols <- colnames(df)[22:1743]

# =========================
# Cox Model Function
# =========================

run_cox_models <- function(
    gene,
    data,
    time_var,
    event_var,
    adjusted = FALSE
) {
  
  if (adjusted) {
    
    vars_needed <- c(
      gene,
      time_var,
      event_var,
      "treatment",
      "bone_met"
    )
    
    formula_text <- paste0(
      "Surv(",
      time_var,
      ", ",
      event_var,
      ") ~ ",
      gene,
      " + treatment + bone_met"
    )
    
  } else {
    
    vars_needed <- c(
      gene,
      time_var,
      event_var
    )
    
    formula_text <- paste0(
      "Surv(",
      time_var,
      ", ",
      event_var,
      ") ~ ",
      gene
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
    filter(term == gene) %>%
    mutate(
      gene_tested = gene,
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
    " model for ",
    time_var
  )
  
  results <- map_dfr(
    gene_cols,
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
  "OS_unadjusted.csv",
  row.names = FALSE
)

write.csv(
  results_DSS_unadj,
  "DSS_unadjusted.csv",
  row.names = FALSE
)

write.csv(
  results_PFS_unadj,
  "PFS_unadjusted.csv",
  row.names = FALSE
)

write.csv(
  results_OS_adj,
  "OS_adjusted.csv",
  row.names = FALSE
)

write.csv(
  results_DSS_adj,
  "DSS_adjusted.csv",
  row.names = FALSE
)

write.csv(
  results_PFS_adj,
  "PFS_adjusted.csv",
  row.names = FALSE
)

message("All analyses completed successfully.")