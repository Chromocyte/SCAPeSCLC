############################################################
# SCAPeSCLC – Gene-Level Cox PH Assumption Testing
#
# Description:
# This script performs gene-level Cox proportional hazards
# modeling followed by proportional hazards (PH)
# assumption testing using cox.zph().
#
# Both unadjusted and adjusted models are generated.
#
# Input:
# - D5_scaled_gene_expression.csv
# - Gene expression columns: 22–1743
#
# Required clinical variables:
# OS, DSS, PFS
# OS_event, DSS_event, PFS_event
# treatment
# bone_metastasis
#
# Output:
# OS_unadjusted_PH.csv
# DSS_unadjusted_PH.csv
# PFS_unadjusted_PH.csv
# OS_adjusted_PH.csv
# DSS_adjusted_PH.csv
# PFS_adjusted_PH.csv
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

file_path <- "D5_scaled_gene_expression.csv"

# =========================
# Load Data
# =========================

df <- read.csv(file_path, stringsAsFactors = FALSE)

gene_cols <- colnames(df)[22:1743]

# =========================
# PH Assumption Function
# =========================

run_ph_models <- function(gene, data, time_var, event_var, adjusted = FALSE){

  if(adjusted){
    vars_needed <- c(gene,time_var,event_var,"treatment","bone_metastasis")
    formula_text <- paste0(
      "Surv(",time_var,", ",event_var,") ~ ",
      gene," + treatment + bone_metastasis"
    )
  } else{
    vars_needed <- c(gene,time_var,event_var)
    formula_text <- paste0(
      "Surv(",time_var,", ",event_var,") ~ ",gene
    )
  }

  missing_cols <- setdiff(vars_needed, names(data))
  if(length(missing_cols)>0){
    stop("Missing columns: ", paste(missing_cols, collapse=", "))
  }

  data_subset <- data[,vars_needed]
  surv_formula <- as.formula(formula_text)

  model <- tryCatch(
    coxph(surv_formula, data=data_subset, ties="efron"),
    error=function(e) NULL
  )

  if(is.null(model)){
    return(data.frame(
      term=gene, estimate=NA_real_, std.error=NA_real_,
      statistic=NA_real_, p.value=NA_real_,
      conf.low=NA_real_, conf.high=NA_real_,
      gene_tested=gene, n=nrow(data_subset),
      global_zph_p=NA_real_,
      global_zph_chisq=NA_real_,
      global_zph_rho=NA_real_,
      gene_zph_p=NA_real_,
      stringsAsFactors=FALSE
    ))
  }

  zph_res <- tryCatch(cox.zph(model), error=function(e) NULL)

  global_zph_p <- NA_real_
  global_zph_chisq <- NA_real_
  global_zph_rho <- NA_real_
  gene_zph_p <- NA_real_

  if(!is.null(zph_res)){
    p_col <- colnames(zph_res$table)[ncol(zph_res$table)]

    if("GLOBAL" %in% rownames(zph_res$table)){
      global_zph_p <- zph_res$table["GLOBAL",p_col]
      if("chisq" %in% colnames(zph_res$table))
        global_zph_chisq <- zph_res$table["GLOBAL","chisq"]
      if("rho" %in% colnames(zph_res$table))
        global_zph_rho <- zph_res$table["GLOBAL","rho"]
    }

    if(gene %in% rownames(zph_res$table)){
      gene_zph_p <- zph_res$table[gene,p_col]
    }
  }

  result <- tidy(model, exponentiate=TRUE, conf.int=TRUE)

  if(adjusted){
    result <- result %>%
      filter(term == gene)
  }

  result <- result %>%
    mutate(
      gene_tested=gene,
      n=nrow(data_subset),
      global_zph_p=global_zph_p,
      global_zph_chisq=global_zph_chisq,
      global_zph_rho=global_zph_rho,
      gene_zph_p=gene_zph_p
    )

  result
}

# =========================
# Wrapper Function
# =========================

run_ph_set <- function(time_var,event_var,adjusted_flag){

  message("Running ",
          ifelse(adjusted_flag,"adjusted","unadjusted"),
          " PH analysis for ",time_var)

  results <- map_dfr(
    gene_cols,
    run_ph_models,
    data=df,
    time_var=time_var,
    event_var=event_var,
    adjusted=adjusted_flag
  ) %>%
    mutate(
      p_adj=p.adjust(p.value, method="BH"),
      CI_width=conf.high-conf.low
    ) %>%
    arrange(CI_width)

  message(
    "Genes satisfying PH assumption: ",
    sum(results$gene_zph_p>0.05, na.rm=TRUE),
    "/",
    nrow(results)
  )

  results
}

# =========================
# Run Unadjusted Models
# =========================

results_OS_unadj  <- run_ph_set("OS","OS_event",FALSE)
results_DSS_unadj <- run_ph_set("DSS","DSS_event",FALSE)
results_PFS_unadj <- run_ph_set("PFS","PFS_event",FALSE)

# =========================
# Run Adjusted Models
# =========================

results_OS_adj  <- run_ph_set("OS","OS_event",TRUE)
results_DSS_adj <- run_ph_set("DSS","DSS_event",TRUE)
results_PFS_adj <- run_ph_set("PFS","PFS_event",TRUE)

# =========================
# Export Results
# =========================

write.csv(results_OS_unadj,"OS_unadjusted_PH.csv",row.names=FALSE)
write.csv(results_DSS_unadj,"DSS_unadjusted_PH.csv",row.names=FALSE)
write.csv(results_PFS_unadj,"PFS_unadjusted_PH.csv",row.names=FALSE)

write.csv(results_OS_adj,"OS_adjusted_PH.csv",row.names=FALSE)
write.csv(results_DSS_adj,"DSS_adjusted_PH.csv",row.names=FALSE)
write.csv(results_PFS_adj,"PFS_adjusted_PH.csv",row.names=FALSE)

message("All PH assumption analyses completed successfully.")
