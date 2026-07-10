\
############################################################
# SCAPeSCLC – Bayesian Pathway Cox Proportional Hazards
# with Proportional Hazards Assumption Testing
#
# Description:
# This script performs Cox proportional hazards modeling
# using Bayesian posterior pathway scores as covariates.
#
# Both unadjusted and adjusted models are generated.
# For every fitted model, the proportional hazards
# assumption is evaluated using cox.zph().
#
# Unadjusted model:
#   Surv(time, event) ~ pathway_posterior
#
# Adjusted model:
#   Surv(time, event) ~ pathway_posterior +
#                       treatment +
#                       bone_metastasis
#
# Input:
#   D13_patient_BP_posteriors.csv
#
# Posterior columns:
#   22, 25, 28, ..., 337
#
# Required clinical columns:
#   OS, DSS, PFS
#   OS_event, DSS_event, PFS_event
#   treatment
#   bone_metastasis
#
# Outputs:
#   BP_OS_unadjusted_PH.csv
#   BP_DSS_unadjusted_PH.csv
#   BP_PFS_unadjusted_PH.csv
#   BP_OS_adjusted_PH.csv
#   BP_DSS_adjusted_PH.csv
#   BP_PFS_adjusted_PH.csv
############################################################

library(dplyr)
library(survival)
library(broom)
library(purrr)
library(stringr)

# =========================
# User Input
# =========================

file_path <- "D13_patient_BP_posteriors.csv"

# =========================
# Load Data
# =========================

df <- read.csv(file_path, stringsAsFactors = FALSE)

# =========================
# Sanitize Column Names
# =========================

colnames(df) <- colnames(df) %>%
  str_replace_all("\\\\s+","_") %>%
  str_replace_all("&","and") %>%
  str_replace_all("-","_") %>%
  str_replace_all("/","_") %>%
  str_replace_all("\\\\.","_")

# =========================
# Identify Posterior Columns
# =========================

posterior_indices <- seq(from = 22, to = 337, by = 3)
bp_cols <- colnames(df)[posterior_indices]

# =========================
# Cox + PH Function
# =========================

run_zph_model <- function(pathway,
                          data,
                          time_var,
                          event_var,
                          adjusted = FALSE){

  if(adjusted){
    vars_needed <- c(pathway,time_var,event_var,
                     "treatment","bone_metastasis")
    formula_text <- paste0(
      "Surv(",time_var,", ",event_var,
      ") ~ ",pathway,
      " + treatment + bone_metastasis")
  } else {
    vars_needed <- c(pathway,time_var,event_var)
    formula_text <- paste0(
      "Surv(",time_var,", ",event_var,
      ") ~ ",pathway)
  }

  data_subset <- data[,vars_needed]
  surv_formula <- as.formula(formula_text)

  model <- tryCatch(
    coxph(surv_formula,
          data=data_subset,
          ties="efron"),
    error=function(e) NULL
  )

  if(is.null(model)){
    return(data.frame(
      term=pathway,
      estimate=NA_real_,
      std.error=NA_real_,
      statistic=NA_real_,
      p.value=NA_real_,
      conf.low=NA_real_,
      conf.high=NA_real_,
      pathway_tested=pathway,
      n=nrow(data_subset),
      global_zph_p=NA_real_,
      global_zph_chisq=NA_real_,
      global_zph_rho=NA_real_,
      pathway_zph_p=NA_real_,
      stringsAsFactors=FALSE))
  }

  zph_res <- tryCatch(cox.zph(model),
                      error=function(e) NULL)

  global_zph_p <- NA_real_
  global_zph_chisq <- NA_real_
  global_zph_rho <- NA_real_
  pathway_zph_p <- NA_real_

  if(!is.null(zph_res)){
    if("GLOBAL" %in% rownames(zph_res$table)){
      p_col <- colnames(zph_res$table)[ncol(zph_res$table)]
      global_zph_p <- zph_res$table["GLOBAL",p_col]
      if("chisq" %in% colnames(zph_res$table))
        global_zph_chisq <- zph_res$table["GLOBAL","chisq"]
      if("rho" %in% colnames(zph_res$table))
        global_zph_rho <- zph_res$table["GLOBAL","rho"]
    }

    if(pathway %in% rownames(zph_res$table)){
      p_col <- colnames(zph_res$table)[ncol(zph_res$table)]
      pathway_zph_p <- zph_res$table[pathway,p_col]
    }
  }

  result <- tidy(model,
                 exponentiate=TRUE,
                 conf.int=TRUE) %>%
    filter(term==pathway) %>%
    mutate(
      pathway_tested=pathway,
      n=nrow(data_subset),
      global_zph_p=global_zph_p,
      global_zph_chisq=global_zph_chisq,
      global_zph_rho=global_zph_rho,
      pathway_zph_p=pathway_zph_p
    )

  result
}

# =========================
# Wrapper
# =========================

run_survival_set <- function(time_var,
                             event_var,
                             adjusted_flag){

  message("Running ",
          ifelse(adjusted_flag,"adjusted","unadjusted"),
          " models for ",time_var)

  res <- map_dfr(
    bp_cols,
    run_zph_model,
    data=df,
    time_var=time_var,
    event_var=event_var,
    adjusted=adjusted_flag)

  res %>%
    mutate(
      p_adj=p.adjust(p.value,"BH"),
      CI_width=conf.high-conf.low
    ) %>%
    arrange(CI_width)
}

results_OS_unadj  <- run_survival_set("OS","OS_event",FALSE)
results_DSS_unadj <- run_survival_set("DSS","DSS_event",FALSE)
results_PFS_unadj <- run_survival_set("PFS","PFS_event",FALSE)

results_OS_adj  <- run_survival_set("OS","OS_event",TRUE)
results_DSS_adj <- run_survival_set("DSS","DSS_event",TRUE)
results_PFS_adj <- run_survival_set("PFS","PFS_event",TRUE)

write.csv(results_OS_unadj,"BP_OS_unadjusted_PH.csv",row.names=FALSE)
write.csv(results_DSS_unadj,"BP_DSS_unadjusted_PH.csv",row.names=FALSE)
write.csv(results_PFS_unadj,"BP_PFS_unadjusted_PH.csv",row.names=FALSE)

write.csv(results_OS_adj,"BP_OS_adjusted_PH.csv",row.names=FALSE)
write.csv(results_DSS_adj,"BP_DSS_adjusted_PH.csv",row.names=FALSE)
write.csv(results_PFS_adj,"BP_PFS_adjusted_PH.csv",row.names=FALSE)

message("All Bayesian pathway Cox PH assumption analyses completed.")
