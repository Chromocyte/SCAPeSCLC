############################################################
# SCAPeSCLC – Patient-Level Bayesian Posterior Estimation
#
# Description:
# This script estimates patient-level Bayesian posterior
# means from ROI-level CTA pathway enrichment Z-scores
# using hierarchical Gaussian models implemented in brms.
#
# SCAPeSCLC file D10 features:
# - CSV file containing ROI-level CTA pathway Z-scores
# - "Patient_ID" column is present
# - CTA pathway Z-score columns occupy columns 7–112
#
# File naming:
# The input file may have any name. Update the file_path
# variable below to match your file location.
#
# Output:
# A wide-format CSV file containing posterior means and
# 95% credible intervals for each patient and pathway.
#
# Dependencies:
# dplyr, stringr, brms, tidyr
#
# Recommended:
# Run on a system with ≥4 CPU cores.
############################################################

library(dplyr)
library(stringr)
library(brms)
library(tidyr)

# =========================
# User Input Section
# =========================

file_path <- "D10_ROI_CTA_Zscores.csv"  
# Replace with your file name (e.g., ROI_CTA_Zscores.csv)

output_file <- "Patient_Level_Bayesian_Posteriors.csv"

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
  str_replace_all("\\s+", "_") %>%   # spaces -> underscores
  str_replace_all("&", "and") %>%    # ampersand -> 'and'
  str_replace_all("-", "_") %>%      # hyphens -> underscore
  str_replace_all("\\.", "_")

# =========================
# Identify CTA Pathway Columns
# =========================

# Assumes pathway columns are 7 through 112
bp_cols <- colnames(df)[7:112]

# =========================
# Model Priors
# =========================

prior_settings <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

# =========================
# Bayesian Model Fitting
# =========================

all_patient_scores <- lapply(bp_cols, function(bp) {
  
  message("Processing pathway: ", bp)
  
  formula <- as.formula(
    paste(bp, "~ 1 + (1|Patient_ID)")
  )
  
  fit <- brm(
    formula = formula,
    data = df,
    family = gaussian(),
    prior = prior_settings,
    chains = 4,
    iter = 4000,
    warmup = 800,
    cores = 4,
    seed = 123,
    control = list(adapt_delta = 0.95)
  )
  
  # Extract posterior summaries
  posterior <- posterior_summary(
    fit,
    pars = "^r_Patient_ID"
  )
  
  # Convert to tidy format
  data.frame(
    Patient_ID = gsub(
      "r_Patient_ID\\[|,Intercept\\]",
      "",
      rownames(posterior)
    ),
    BP = bp,
    Posterior = posterior[, "Estimate"],
    Lower_CI = posterior[, "Q2.5"],
    Upper_CI = posterior[, "Q97.5"],
    stringsAsFactors = FALSE
  )
  
})

# =========================
# Combine Results
# =========================

patient_bayesian_bp <- do.call(
  rbind,
  all_patient_scores
)

# =========================
# Convert to Wide Format
# =========================

posterior_wide <- patient_bayesian_bp %>%
  select(Patient_ID, BP, Posterior) %>%
  pivot_wider(
    names_from = BP,
    values_from = Posterior,
    names_prefix = "Post_"
  )

lower_wide <- patient_bayesian_bp %>%
  select(Patient_ID, BP, Lower_CI) %>%
  pivot_wider(
    names_from = BP,
    values_from = Lower_CI,
    names_prefix = "Lower_"
  )

upper_wide <- patient_bayesian_bp %>%
  select(Patient_ID, BP, Upper_CI) %>%
  pivot_wider(
    names_from = BP,
    values_from = Upper_CI,
    names_prefix = "Upper_"
  )

# =========================
# Merge Outputs
# =========================

patient_bayesian_wide <- posterior_wide %>%
  left_join(lower_wide, by = "Patient_ID") %>%
  left_join(upper_wide, by = "Patient_ID")

# =========================
# Export Results
# =========================

write.csv(
  patient_bayesian_wide,
  file = output_file,
  row.names = FALSE
)

message("Analysis complete.")
message("Output saved to: ", output_file)
