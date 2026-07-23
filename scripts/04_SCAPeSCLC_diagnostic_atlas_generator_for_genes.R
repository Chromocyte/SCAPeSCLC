# ============================================================================
# SCAPeSCLC Diagnostic Atlas Generator
# ============================================================================
#
# Description:
#   Generates multi-page diagnostic dashboards for gene-level Cox proportional
#   hazards models using standardized gene expression profiles from SCAPeSCLC.
#   Each dashboard includes:
#     • Martingale residuals
#     • Schoenfeld residuals
#     • Deviance residuals
#     • DFBETA influence diagnostics
#     • Model summary panel
#
# Input:
#   D5_scaled_gene_expression.csv
#
# Output:
#   SCAPeSCLC_Diagnostic_Atlas_<endpoint>.pdf
#
# Endpoints:
#   OS
#   DSS
#   PFS
#
# Author: M. Shirvaliloo
# License: MIT
# Version: 1.0.0
# ============================================================================

# ============================================================================
# 1. Setup
# ============================================================================
# Load Libraries
library(ggplot2)
library(grid)
library(patchwork)
library(survival)

# Read Source Data
df <- read.csv(file = "D5_scaled_gene_expression.csv",
               stringsAsFactors = FALSE)

# Gene Columns
GENE_COLUMNS <- 22:1743
gene_cols <- colnames(df)[GENE_COLUMNS]

# ============================================================================
# 2. Function Definitions
# ============================================================================
# ===========================
# 2.1. SCAPeSCLC Theme Layout
# ===========================
# ------------------------------------------------------------------
# Apply SCAPeSCLC Plot Theme
#
# Returns:
#   A customized ggplot2 theme used throughout the diagnostic atlas.
# ------------------------------------------------------------------
theme_scapesclc <- function() {
  
  theme_classic(base_size = 10.5) +
    
    theme(
      
      plot.title =
        element_text(face = "bold",
                     size = 12,
                     family = "Myriad Pro",
                     color = "grey25"),
      
      plot.subtitle = 
        element_text(size = 8.5,
                     family = "Myriad Pro",
                     color = "grey35"),
      
      axis.title =
        element_text(size = 9,
                     family = "Helvetica",
                     color = "grey25",
                     face = "plain"),
      
      axis.text =
        element_text(size = 8.5,
                     family = "Helvetica",
                     color = "grey45"),
      
      legend.position = "none"
    )
}

# ====================
# 2.2. Subtitle Layout
# ====================
# -----------------------------------------------------
# Create Plot Subtitle
#
# Args:
#   endpoint  : Survival endpoint (OS, DSS or PFS).
#   plot_type : Name of the diagnostic plot.
#   extra     : Optional text appended to the subtitle.
#
# Returns:
#   A formatted subtitle string.
# -----------------------------------------------------
make_subtitle <- function(endpoint, data, event_var,
                          plot_type, extra = NULL) {
  subtitle <- paste0(
    endpoint,
    " • ",
    plot_type
  )
  if (!is.null(extra)){
    subtitle <- paste0(
      subtitle,
      " • ",
      extra
    )
  }
  
  subtitle
  
}

# ==================
# 2.3. Summary Table
# ==================
# ------------------------------------------------------------------
# Extract Model Summary
#
# Args:
#   model     : Fitted Cox proportional hazards model.
#   zph       : cox.zph object for proportional hazards diagnostics.
#   data      : Input data frame.
#   event_var : Name of the event indicator variable.
#
# Returns:
#   A data frame containing key model summary statistics.
# ------------------------------------------------------------------
extract_model_summary <- function(model,
                                  zph,
                                  gene,
                                  data,
                                  event_var){
  tidy_model <- broom::tidy(
    model,
    exponentiate = TRUE,
    conf.int = TRUE
  )
  
  patients <- nrow(data)
  events <- sum(data[[event_var]])
  hr <- sprintf("%.3f", tidy_model$estimate)
  ci <- paste0(
    sprintf("%.3f", tidy_model$conf.low),
    "–",
    sprintf("%.3f", tidy_model$conf.high)
  )
  cox_p <- sprintf("%.3f", tidy_model$p.value)
  global_ph_p <- sprintf("%.3f",
                         zph$table["GLOBAL", "p"])
  
  # -----------------------------------------------------------------
  # Create Model Summary Panel
  #
  # Args:
  #   model_summary : Data frame returned by extract_model_summary().
  #   gene          : Gene symbol.
  #
  # Returns:
  #   A grid graphical object (grob) displaying the model summary.
  # -----------------------------------------------------------------
  
  model_summary <- data.frame(
    
    Metric = c(
      "Patients",
      "Events",
      "Hazard Ratio",
      "95% CI",
      "P-value",
      "Global PH p"
    ),
    
    Value = c(
      patients,
      events,
      hr,
      ci,
      cox_p,
      global_ph_p
    ),
    
    stringsAsFactors = FALSE
  )
  
  model_summary
}

make_summary_grob <- function(model_summary, gene) {
  
  # Canvas
  children <- list()
  
  # Gene Title
  children[[length(children) + 1]] <-
    textGrob(
      label = gene,
      x = unit(0.052, "npc"),
      y = unit(0.98, "npc"),
      just = c("left", "top"),
      gp = gpar(
        fontsize = 15,
        fontface = "bold"
      )
    )
  
  # Horizontal Divider
  children[[length(children) + 1]] <-
    segmentsGrob(
      x0 = unit(0.052, "npc"),
      x1 = unit(1.00, "npc"),
      y0 = unit(0.845, "npc"),
      y1 = unit(0.845, "npc"),
      gp = gpar(lwd = 0.25,
                color = "grey75")
    )
  
  # Table Rows
  y <- c(
    0.75,
    0.65,
    0.50,
    0.40,
    0.30,
    0.15
  )
  
  for(i in seq_len(nrow(model_summary))) {
    
    # Metric
    children[[length(children) + 1]] <-
      textGrob(
        label = model_summary$Metric[i],
        x = unit(0.052, "npc"),
        y = unit(y[i], "npc"),
        just = "left",
        gp = gpar(
          fontsize = 9.75,
          fontfamily = "Avenir",
          fontface = "bold",
          color = "grey25"
        )
      )
    
    # Value
    children[[length(children) + 1]] <-
      textGrob(
        label = model_summary$Value[i],
        x = unit(0.556, "npc"),
        y = unit(y[i], "npc"),
        just = "left",
        gp = gpar(
          fontsize = 10,
          fontfamily = "mono"
        )
      )
  }
  
  grobTree(children = do.call(gList, children))
}

# =====================
# 2.4. Diagnostic Plots
# =====================
# ================================
# 2.4.1. Schoenfeld Residuals Plot
# ================================
# ---------------------------------------------------
# Create Schoenfeld Residual Plot
#
# Args:
#   model    : Fitted Cox proportional hazards model.
#   gene     : Gene symbol.
#   endpoint : Survival endpoint.
#
# Returns:
#   A list containing:
#     plot : ggplot object
#     zph  : cox.zph object
#----------------------------------------------------
make_schoenfeld_plot <- function(model, gene, endpoint) {
  
  zph <- cox.zph(model)
  
  df_zph <- data.frame(
    Time = zph$x,
    Residual = zph$y[, 1]
  )
  
  ylim_zph <- max(abs(zph$time))
  
  gene_ph_p <- sprintf("%.3f",
                       zph$table[gene, "p"])
  
  schoenfeld_plot <- ggplot(df_zph,
                            aes(x = Time,
                                y = Residual)) +
    
    geom_point(size = 2.2,
               alpha = 0.65,
               shape = 16) +
    
    geom_smooth(
      method = "loess",
      span = 0.9,
      se = TRUE,
      formula = y ~ x
    ) +
    
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "grey70",
               linewidth = 0.6) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    coord_cartesian(
      ylim = c(-ylim_zph, ylim_zph)
    ) +
    
    labs(
      subtitle = make_subtitle(
        endpoint = endpoint,
        plot_type = "Schoenfeld Residuals",
        extra = paste0(
          "Gene PH p = ",
          gene_ph_p
        )
      ),
      x = "Transformed Time",
      y = "Residual"
    ) +
    
    theme_scapesclc()
  
  list(
    plot = schoenfeld_plot,
    zph = zph
  )
}

# ================================
# 2.4.2. Martingale Residuals Plot
# ================================
# ---------------------------------------------------
# Create Martingale Residual Plot
#
# Args:
#   model    : Fitted Cox proportional hazards model.
#   data     : Input data frame.
#   gene     : Gene symbol.
#   endpoint : Survival endpoint.
#
# Returns:
#   A list containing:
#     plot : ggplot object
# ---------------------------------------------------
make_martingale_plot <- function(model,
                                 gene,
                                 endpoint) {
  
  mart_df <- data.frame(
    Patient = seq_len(nrow(df)),
    Expression = df[[gene]],
    Martingale = residuals(model, type = "martingale")
  )
  
  ylim_mart <- max(abs(mart_df$Martingale))
  
  martingale_plot <- ggplot(mart_df,
                            aes(x = Expression,
                                y = Martingale)) +
    
    geom_point(
      size = 2.2,
      alpha = 0.65,
      shape = 16
    ) +
    
    geom_smooth(
      method = "loess",
      span = 0.9,
      formula = y ~ x,
      se = TRUE,
      linewidth = 1.0
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.6
    ) +
    
    geom_rug(
      alpha = 0.4,
      sides = "b") +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    
    coord_cartesian(
      ylim = c(-ylim_mart, ylim_mart)
    ) +
    
    labs(
      subtitle = make_subtitle(
        endpoint = endpoint,
        plot_type = "Martingale Residuals"
      ),
      x = paste0(gene,
                 " Expression (standardized)"),
      y = "Residual"
    ) +
    
    theme_scapesclc()
  
  list(plot = martingale_plot)
}

# ==============================
# 2.4.3. Deviance Residuals Plot
# ==============================
# ---------------------------------------------------
# Create Deviance Residual Plot
#
# Args:
#   model    : Fitted Cox proportional hazards model.
#   data     : Input data frame.
#   gene     : Gene symbol.
#   endpoint : Survival endpoint.
#
# Returns:
#   A list containing:
#     plot : ggplot object
# ---------------------------------------------------
make_deviance_plot <- function(model,
                               gene,
                               endpoint) {
  
  df_dev <- data.frame(
    Patient = seq_len(nrow(df)),
    Deviance = residuals(model, type = "deviance")
  )
  
  ylim_dev <- max(abs(df_dev$Deviance))
  
  deviance_plot <- ggplot(df_dev,
                          aes(x = Patient,
                              y = Deviance)) +
    
    geom_point(size = 2.2,
               alpha = 0.65,
               shape = 16) +
    
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "grey70",
               linewidth = 0.6) +
    
    geom_hline(yintercept = c(-2, 2),
               linetype = "dotted",
               color = "grey75",
               linewidth = 0.5) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    
    coord_cartesian(
      ylim = c(-ylim_dev, ylim_dev)
    ) +
    
    labs(
      subtitle = make_subtitle(
        endpoint = endpoint,
        plot_type = "Deviance Residuals"
      ),
      x = "Patient",
      y = "Residual"
    ) +
    
    theme_scapesclc()
  
  list(plot = deviance_plot)
}

# ==============================
# 2.4.4. D.DFBETA Influence Plot
# ==============================
# ---------------------------------------------------
# Create DFBETA Influence Plot
#
# Args:
#   model    : Fitted Cox proportional hazards model.
#   data     : Input data frame.
#   gene     : Gene symbol.
#   endpoint : Survival endpoint.
#
# Returns:
#   A list containing:
#     plot : ggplot object
# ---------------------------------------------------
make_dfbeta_plot <- function(model, gene, endpoint) {
  
  df_dfb <- data.frame(
    Patient = seq_len(nrow(df)),
    DFBETA = residuals(model, type = "dfbeta")
  )
  
  ylim_dfb <- max(abs(df_dfb$DFBETA))
  
  dfbeta_plot <- ggplot(df_dfb,
                        aes(x = Patient,
                            y = DFBETA)) +
    
    geom_point(size = 2.2,
               alpha = 0.65,
               shape = 16) +
    
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "grey70",
               linewidth = 0.6) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    
    coord_cartesian(
      ylim = c(-ylim_dfb, ylim_dfb)
    ) +
    
    labs(
      subtitle = make_subtitle(
        endpoint = endpoint,
        plot_type = "DFBETA"
      ),
      x = "Patient",
      y = "DFBETA"
    ) +
    
    theme_scapesclc()
  
  list(plot = dfbeta_plot)
}

# ===================================
# 2.5. Diagnostic Dashboard Generator
# ===================================
# ------------------------------------------------------------------
# Create Diagnostic Dashboard
#
# Args:
#   gene      : Gene symbol.
#   endpoint  : Survival endpoint (OS, DSS or PFS).
#   time_var  : Name of the survival time variable.
#   event_var : Name of the event indicator variable.
#
# Returns:
#   A patchwork object containing the complete diagnostic dashboard.
# ------------------------------------------------------------------
make_diagnostic_dashboard <- function(gene,
                                      endpoint,
                                      time_var,
                                      event_var) {
  
  formula <- as.formula(
    paste(
      "Surv(", time_var, ", ", event_var, ") ~", gene
    )
  )
  
  model <- survival::coxph(
    formula,
    data = df
  )
  
  schoenfeld_results <- make_schoenfeld_plot(
    model = model,
    gene = gene,
    endpoint = endpoint
  )
  
  schoenfeld_plot <- schoenfeld_results$plot
  zph <- schoenfeld_results$zph
  
  martingale_results <- make_martingale_plot(
    model = model,
    gene = gene,
    endpoint = endpoint
  )
  
  martingale_plot <- martingale_results$plot
  
  deviance_results <- make_deviance_plot(
    model = model,
    gene = gene,
    endpoint = endpoint
  )
  
  deviance_plot <- deviance_results$plot
  
  dfbeta_results <- make_dfbeta_plot(
    model = model,
    gene = gene,
    endpoint = endpoint
  )
  
  dfbeta_plot <- dfbeta_results$plot
  
  model_summary <- extract_model_summary(
    model = model,
    zph = zph,
    gene = gene,
    data = df,
    event_var = event_var
  )
  
  summary_grob <-
    make_summary_grob(
      model_summary,
      gene
    )
  
  summary_panel <-
    patchwork::wrap_elements(
      full = summary_grob
    )
  
  dashboard <- summary_panel /
    (martingale_plot | schoenfeld_plot) /
    (deviance_plot | dfbeta_plot)
  
  dashboard
}

# ===============================
# 2.6. Diagnostic Atlas Generator
# ===============================
# ---------------------------------------------------------
# Generate Diagnostic Atlas
#
# Args:
#   endpoint : Survival endpoint (OS, DSS or PFS).
#
# Returns:
#   Generates a multi-page PDF atlas and writes it to disk.
# ---------------------------------------------------------
make_atlas <- function(endpoint) {
  
  start_time <- Sys.time()
  
  time_var <- endpoint
  event_var <- paste0(endpoint, "_event")
  
  cat("=====================================
        Generating Atlas
=====================================\n",
      sep = ""
  )
  
  cat(
    "Endpoint: ", endpoint, "\n",
    "Genes: ", length(gene_cols), "\n",
    "Output: SCAPeSCLC Diagnostic Atlas ", endpoint,
    ".pdf\n\n",
    "Please wait...\n\n",
    sep = ""
  )
  
  cairo_pdf(
    file = paste0("SCAPeSCLC Diagnostic Atlas - ",
                  endpoint,
                  ".pdf"
    ),
    width = 8,
    height = 6.5
  )
  
  for (i in seq_along(gene_cols)) {
    
    gene <- gene_cols[i]
    
    message(
      sprintf("[%d/%d] Building dashboard: %s",
              i,
              length(gene_cols),
              gene)
    )
    
    dashboard <- make_diagnostic_dashboard(
      gene = gene,
      time_var = time_var,
      event_var = event_var,
      endpoint = endpoint
    )
    
    print(dashboard)
    
    rm(dashboard)
  }
  
  dev.off()
  
  elapsed <- Sys.time() - start_time
  
  cat("=====================================
    Atlas completed successfully.
=====================================\n\n",
      "Elapsed: ", round(as.numeric(elapsed,
                                    units = "mins"),
                         2),
      " minutes\n\n",
      sep = ""
  )
  
}

# ================
# 3. Main Program
# ================
cat("
=====================================
 SCAPeSCLC Diagnostic Atlas Generator
=====================================

1. Overall Survival (OS)
2. Disease-Specific Survival (DSS)
3. Progression-Free Survival (PFS)
4. Generate all atlases?\n
")

choice <- readline(prompt = "Selection (1-4): ")

if (choice == "1") {
  
  make_atlas("OS")
} else if (choice == "2") {
  
  make_atlas("DSS")
} else if (choice == "3") {
  
  make_atlas("PFS")
} else if (choice == "4") {
  
  cat(
    "\nGenerating all atlases...\n\n",
    " • Overall Survival\n",
    " • Disease-Specific Survival\n",
    " • Progression-Free Survival\n\n",
    sep = ""
  )
  
  make_atlas("OS")
  make_atlas("DSS")
  make_atlas("PFS")
  
} else {
  
  stop("Invalid selection.")
}