# ==============================================================================
# Script Name: 08_Occupancy_Abundance_CoreTaxa_Analysis.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Fit logistic S-curve for mean relative abundance - occupancy relationship,
#          screen core taxa with occupancy threshold, generate scatter plot with 95% CI
# Input: Os_rb.xlsx, Os_eb.xlsx, Os_rf.xlsx, Os_ef.xlsx taxon abundance matrix stored in ./input
# Output: Occupancy-abundance scatter PNG/PDF, full & core taxa statistics Excel in ./output
# Dependencies: tidyverse, readxl, writexl, vegan, ggsci, minpack.lm, numDeriv
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified significance/filter threshold, easy batch modification across all scripts
global_alpha <- 0.05
core_occupancy_threshold <- 0.8 # Core taxa definition: detected in ≥80% samples

# -------------------------- Install & load required CRAN packages --------------------------
required_pkgs <- c("tidyverse", "readxl", "writexl", "vegan", "ggsci", "minpack.lm", "numDeriv")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# -------------------------- Global relative path configuration --------------------------
# Cross-platform compatible input & output folders, auto create missing directories
input_dir  <- "./input"
output_dir <- "./output"

if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created output directory: ./output")
}

# -------------------------- Custom Logistic S-curve fitting function with 95% CI --------------------------
# Model formula: occupancy = 1 / (1 + exp(-a*(log10_mean_RA - b)))
# Input: data frame containing log10_mean_RA and occupancy columns
# Output: fitted curve values + bounded 95% confidence interval (0 ≤ occupancy ≤ 1)
fit_logistic_s_curve <- function(occ_abun_df) {
  # Non-linear least squares fitting with robust LM algorithm
  nls_fit <- nlsLM(
    occupancy ~ 1 / (1 + exp(-a * (log10_mean_RA - b))),
    data = occ_abun_df,
    start = list(a = 2, b = -3),
    control = nls.lm.control(maxiter = 1024)
  )
  
  # Generate dense continuous x sequence for smooth prediction curve
  x_pred_seq <- seq(
    min(occ_abun_df$log10_mean_RA, na.rm = TRUE),
    max(occ_abun_df$log10_mean_RA, na.rm = TRUE),
    length.out = 1000
  )
  pred_input_df <- data.frame(log10_mean_RA = x_pred_seq)
  
  # Predict fitted occupancy values from logistic model
  pred_occupancy <- predict(nls_fit, newdata = pred_input_df)
  
  # Calculate 95% confidence interval via Jacobian & variance-covariance matrix
  coef_vec <- coef(nls_fit)
  cov_matrix <- vcov(nls_fit)
  jacobian_mat <- numDeriv::jacobian(
    function(par) {1 / (1 + exp(-par[1] * (x_pred_seq - par[2])))},
    c(coef_vec["a"], coef_vec["b"])
  )
  se_pred <- sqrt(diag(jacobian_mat %*% cov_matrix %*% t(jacobian_mat)))
  
  # Restrict CI bounds within valid occupancy logical range [0, 1]
  ci_upper <- pmin(pred_occupancy + 1.96 * se_pred, 1)
  ci_lower <- pmax(pred_occupancy - 1.96 * se_pred, 0)
  
  # Compile curve prediction result table with raw non-log transformed abundance
  curve_result <- data.frame(
    log10_mean_RA = x_pred_seq,
    mean_occupancy_pred = pred_occupancy,
    ci_upper = ci_upper,
    ci_lower = ci_lower,
    raw_RA = 10 ^ x_pred_seq
  )
  return(curve_result)
}

# -------------------------- Single community core taxa analysis workflow --------------------------
# Input args:
#   file_name: excel file name stored in ./input folder
#   plot_color: hex color code for core taxa scatter points
# Return list: ggplot object, full taxa abundance-occupancy table, filtered core taxa table
process_single_community <- function(file_name, plot_color) {
  full_file_path <- file.path(input_dir, file_name)
  if (!file.exists(full_file_path)) stop(paste("Input abundance file missing at path:", full_file_path))
  
  # Read OTU/genus abundance matrix, first column reserved for full taxonomic annotation
  raw_abundance <- read_excel(full_file_path, sheet = 1) %>%
    column_to_rownames(var = "#Taxonomy") %>%
    as.data.frame()
  
  if (nrow(raw_abundance) == 0) stop("Error: Input abundance table contains zero taxa entries")
  n_taxa <- nrow(raw_abundance)
  n_sample <- ncol(raw_abundance)
  message(paste("===== Start processing dataset:", file_name, "====="))
  message(paste("Taxa total count:", n_taxa, "| Sample total count:", n_sample))
  
  # Step1: Binary presence-absence matrix, calculate occupancy ratio per taxon
  pa_matrix <- 1 * (raw_abundance > 0)
  occupancy_vec <- rowSums(pa_matrix) / n_sample
  
  # Step2: Sample-wise total-sum relative abundance normalization
  rel_abundance_matrix <- decostand(raw_abundance, method = "total", MARGIN = 2)
  mean_rel_abundance <- rowMeans(rel_abundance_matrix, na.rm = TRUE)
  # Add tiny offset to avoid log10(0) infinite error
  log10_mean_ra <- log10(mean_rel_abundance + 1e-10)
  
  # Merge all taxon-level metrics + original raw abundance matrix
  full_taxon_meta <- data.frame(
    Taxonomy = rownames(raw_abundance),
    mean_relative_abundance = mean_rel_abundance,
    log10_mean_RA = log10_mean_ra,
    occupancy = occupancy_vec,
    is_core_taxon = ifelse(occupancy_vec >= core_occupancy_threshold, "Core", "Non-core"),
    stringsAsFactors = FALSE
  ) %>% left_join(raw_abundance %>% rownames_to_column("Taxonomy"), by = "Taxonomy")
  
  # Filter core taxa subset meeting occupancy ≥ predefined threshold
  core_taxon_subset <- full_taxon_meta %>% filter(is_core_taxon == "Core")
  message(paste("Identified core taxa count (≥", core_occupancy_threshold*100, "% occupancy):", nrow(core_taxon_subset)))
  
  # Fit logistic S-curve abundance-occupancy trend
  fitted_curve_data <- fit_logistic_s_curve(full_taxon_meta)
  
  # Generate standardized abundance-occupancy scatter plot
  occupancy_plot <- ggplot() +
    # Non-core taxa: low-transparency gray background points
    geom_point(
      data = full_taxon_meta %>% filter(is_core_taxon == "Non-core"),
      aes(x = log10_mean_RA, y = occupancy),
      color = "gray80", alpha = 0.6, size = 2
    ) +
    # Core taxa: highlighted custom color solid points
    geom_point(
      data = full_taxon_meta %>% filter(is_core_taxon == "Core"),
      aes(x = log10_mean_RA, y = occupancy),
      color = plot_color, size = 3, alpha = 0.9
    ) +
    # Fitted logistic main trend solid line
    geom_line(
      data = fitted_curve_data,
      aes(x = log10_mean_RA, y = mean_occupancy_pred),
      color = "black", linewidth = 1.2, alpha = 0.8
    ) +
    # 95% confidence interval dashed boundary lines
    geom_line(
      data = fitted_curve_data,
      aes(x = log10_mean_RA, y = ci_upper),
      color = "black", linetype = "dashed", linewidth = 1, alpha = 0.7
    ) +
    geom_line(
      data = fitted_curve_data,
      aes(x = log10_mean_RA, y = ci_lower),
      color = "black", linetype = "dashed", linewidth = 1, alpha = 0.7
    ) +
    # Horizontal threshold reference line for core taxa cutoff
    geom_hline(
      yintercept = core_occupancy_threshold,
      color = "black", linetype = "solid", linewidth = 0.8, alpha = 0.7
    ) +
    labs(
      x = "log10(Mean Relative Abundance)",
      y = "Occupancy (Detection frequency across samples)"
    ) +
    theme_light() +
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.margin = margin(10, 10, 10, 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  # Return complete analysis outputs for downstream export
  return(list(
    occupancy_scatter_plot = occupancy_plot,
    full_taxon_metadata = full_taxon_meta,
    core_taxon_dataset = core_taxon_subset,
    logistic_fit_curve = fitted_curve_data
  ))
}

# -------------------------- Batch analysis for four standard datasets (match script 03~07 input files) --------------------------
# Unified color palette consistent with prior community analysis scripts
community_config <- list(
  Rhizo_Bacteria = list(file = "Os_rb.xlsx", color = "#FF6347"),
  Endo_Bacteria  = list(file = "Os_eb.xlsx", color = "#32CD32"),
  Rhizo_Fungi    = list(file = "Os_rf.xlsx", color = "#9370DB"),
  Endo_Fungi     = list(file = "Os_ef.xlsx", color = "#FFD700")
)

# Batch loop run all four community datasets
all_community_result <- map(community_config, function(config) {
  process_single_community(file_name = config$file, plot_color = config$color)
})

# Export individual figures & statistics for each community
walk2(names(all_community_result), all_community_result, function(name, res) {
  fig_base <- paste0(name, "_Occupancy_Abundance_Scatter")
  excel_name <- paste0(name, "_Core_Taxa_Occupancy_Statistics.xlsx")
  
  # Export high resolution raster PNG & vector PDF
  ggsave(file.path(output_dir, paste0(fig_base, ".png")), res$occupancy_scatter_plot, width = 12, height = 8, dpi = 300)
  ggsave(file.path(output_dir, paste0(fig_base, ".pdf")), res$occupancy_scatter_plot, width = 12, height = 8)
  
  # Export multi-sheet statistical table
  write_xlsx(
    list(
      Full_Taxa_Metadata = res$full_taxon_metadata,
      Core_Taxa_Subset = res$core_taxon_dataset,
      Logistic_Fit_Curve_Data = res$logistic_fit_curve
    ),
    file.path(output_dir, excel_name)
  )
})

# Merge all core taxa from four habitats into one unified comparison table
all_core_combined <- map_dfr(names(all_community_result), function(name) {
  all_community_result[[name]]$core_taxon_dataset %>% mutate(Community = name)
})
write_xlsx(list(All_Combined_Core_Taxa = all_core_combined), file.path(output_dir, "All_Habitat_Core_Taxa_Merged_Summary.xlsx"))

# Print completion prompt
cat("\n===== All occupancy-core taxa analysis completed =====\n")
cat("All scatter figures, logistic curve data and taxon statistics saved to ./output folder\n")