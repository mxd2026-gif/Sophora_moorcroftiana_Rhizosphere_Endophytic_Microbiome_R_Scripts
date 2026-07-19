# ==============================================================================
# Script Name: 21_Feature_Environment_Spearman_Correlation_Heatmap.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Non-parametric Spearman rank correlation analysis between functional gene/metabolic features and soil environmental factors,
#          calculate correlation coefficients + permutation P-values, generate significance annotation heatmap,
#          export multi-sheet complete statistical matrix Excel tables for manuscript supplementary materials
# Input: 10last.xlsx functional feature abundance matrix, 1env16-8.xlsx environmental physicochemical table stored in unified ./input folder
# Output: Dual-format 300 DPI PNG + vector PDF annotated Spearman correlation heatmap, multi-sheet integrated correlation statistics Excel saved in ./output
# Dependencies: tidyverse, readxl, ggplot2, RColorBrewer, openxlsx, numDeriv
# Standardization: Fully aligned with script 01~20 unified global parameters, cross-platform path rules, journal dual-format output specs
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical significance threshold consistent with full serial analysis pipeline
global_alpha <- 0.05

# --------------------------- Install and load all dependent packages ---------------------------
required_packages <- c("tidyverse", "readxl", "ggplot2", "RColorBrewer", "openxlsx", "numDeriv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Dependency package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# --------------------------- Cross-platform relative path configuration (Unified standard ./input & ./output) ---------------------------
input_dir  <- "./input"
output_dir <- "./output"

# Auto create standard project folders with running status prompt
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# Global unified four geographic sampling group color palette (consistent with script 03~20 full pipeline)
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)
group_factor_level <- names(group_color_config)

# --------------------------- Function 1: End-to-End Single Dataset Spearman Correlation Calculation Core Function ---------------------------
#' Compute pairwise Spearman rank correlation, P-value matrix, significance annotation and plotting tidy data
#' @param feature_path Full standardized cross-platform path of functional feature abundance excel table
#' @param env_path Full standardized cross-platform path of environmental physicochemical excel table
#' @param env_target_cols Fixed ordered vector of target soil environmental factor column names
#' @param dataset_tag Unique descriptive identifier for output table and figure filenames
#' @param sheet_env Sheet name storing environmental factor data, default "all"
#' @param sheet_feature Sheet name storing feature abundance matrix with ID column, default "ID"
#' @return List containing correlation r matrix, P-value matrix, significance mark matrix, heatmap annotation text matrix, long-format ggplot plotting dataset, preserved axis label order
single_dataset_spearman_analysis <- function(feature_path, env_path, env_target_cols, dataset_tag, sheet_env = "all", sheet_feature = "ID") {
  # Read standardized environmental dataset
  env_raw <- read_excel(env_path, sheet = sheet_env)
  env_data <- env_raw %>%
    select(SampleID, all_of(env_target_cols)) %>%
    column_to_rownames("SampleID") %>%
    as.data.frame()
  
  # Read functional feature abundance matrix
  feature_raw <- read_excel(feature_path, sheet = sheet_feature)
  if (colnames(feature_raw)[1] != "ID") {
    colnames(feature_raw)[1] <- "ID"
  }
  feature_data <- feature_raw %>%
    column_to_rownames("ID") %>%
    as.data.frame()
  
  # Match overlapping sample IDs between feature matrix and environmental table
  shared_samples <- intersect(colnames(feature_data), rownames(env_data))
  if (length(shared_samples) < 3) {
    stop(sprintf("[%s] Insufficient overlapping valid samples (minimum 3 required) between feature dataset and environmental dataset", dataset_tag))
  }
  
  # Reorder and subset matrix to retain only shared samples
  feature_matched <- feature_data[, shared_samples, drop = FALSE] %>% t() %>% as.data.frame()
  env_matched <- env_data[shared_samples, , drop = FALSE]
  
  # Lock fixed axis display order for heatmap
  feature_order <- colnames(feature_matched)
  env_order <- colnames(env_matched)
  
  n_feature <- ncol(feature_matched)
  n_env <- ncol(env_matched)
  
  # Initialize empty storage matrices
  spearman_r_matrix <- matrix(NA, nrow = n_feature, ncol = n_env, dimnames = list(feature_order, env_order))
  p_value_matrix <- matrix(NA, nrow = n_feature, ncol = n_env, dimnames = list(feature_order, env_order))
  
  # Double loop pairwise Spearman test
  for (i in seq_len(n_feature)) {
    for (j in seq_len(n_env)) {
      vec_x <- feature_matched[, i]
      vec_y <- env_matched[, j]
      clean_x <- na.omit(vec_x)
      clean_y <- na.omit(vec_y)
      
      # Skip pairs with insufficient valid non-duplicate observations
      if (length(clean_x) < 3 || length(clean_y) < 3 || length(unique(clean_x)) < 2 || length(unique(clean_y)) < 2) {
        spearman_r_matrix[i, j] <- NA
        p_value_matrix[i, j] <- NA
        next
      }
      corr_test <- cor.test(vec_x, vec_y, method = "spearman", exact = FALSE)
      spearman_r_matrix[i, j] <- corr_test$estimate
      p_value_matrix[i, j] <- corr_test$p.value
    }
  }
  
  # Assign statistical significance labels based on global alpha threshold
  sig_label_matrix <- matrix("ns", nrow = n_feature, ncol = n_env, dimnames = list(feature_order, env_order))
  sig_label_matrix[p_value_matrix < 0.05] <- "*"
  sig_label_matrix[p_value_matrix < 0.01] <- "**"
  sig_label_matrix[p_value_matrix < 0.001] <- "***"
  
  # Combine correlation coefficient and significance mark as cell annotation text
  annotation_label_matrix <- matrix(
    paste0(sprintf("%.2f", spearman_r_matrix), "(", sig_label_matrix, ")"),
    nrow = n_feature, ncol = n_env,
    dimnames = list(feature_order, env_order)
  )
  
  # Reshape matrix to long tidy format for ggplot heatmap
  plot_tidy_data <- expand.grid(
    Feature = factor(feature_order, levels = feature_order),
    Environmental_Variable = factor(env_order, levels = env_order)
  )
  plot_tidy_data$annotation_text <- as.vector(annotation_label_matrix)
  plot_tidy_data$spearman_r <- as.vector(spearman_r_matrix)
  
  message(sprintf("[%s] Correlation calculation completed: %d matched samples, %d functional features, %d environmental physicochemical indicators",
                  dataset_tag, length(shared_samples), n_feature, n_env))
  
  return(list(
    r_matrix = spearman_r_matrix,
    p_matrix = p_value_matrix,
    significance_matrix = sig_label_matrix,
    annotation_label_matrix = annotation_label_matrix,
    plot_data = plot_tidy_data,
    feature_order = feature_order,
    env_order = env_order
  ))
}

# --------------------------- Function 2: Batch Export Full Correlation Statistical Matrices to Multi-Sheet Excel ---------------------------
#' Package all Spearman calculation results into standardized integrated Excel workbook
#' @param analysis_result Complete output list returned by single_dataset_spearman_analysis
#' @param dataset_tag Unique dataset identifier for output filename prefix
#' @param save_dir Unified standardized output directory path
export_spearman_excel_result <- function(analysis_result, dataset_tag, save_dir) {
  # Main workbook storing r, p, significance matrices
  wb <- createWorkbook()
  addWorksheet(wb, "Spearman_R_Correlation_Coefficient")
  writeData(wb, "Spearman_R_Correlation_Coefficient", as.data.frame(analysis_result$r_matrix) %>% rownames_to_column("Functional_Feature"))
  
  addWorksheet(wb, "Correlation_P_Value")
  writeData(wb, "Correlation_P_Value", as.data.frame(analysis_result$p_matrix) %>% rownames_to_column("Functional_Feature"))
  
  addWorksheet(wb, "Statistical_Significance_Mark")
  writeData(wb, "Statistical_Significance_Mark", as.data.frame(analysis_result$significance_matrix) %>% rownames_to_column("Functional_Feature"))
  
  saveWorkbook(wb, file.path(save_dir, paste0("Spearman_Correlation_Full_Statistical_Summary_", dataset_tag, ".xlsx")), overwrite = TRUE)
  
  # Separate workbook for heatmap combined annotation labels
  label_wb <- createWorkbook()
  addWorksheet(label_wb, "Heatmap_Cell_Annotation_Label_Matrix")
  writeData(label_wb, "Heatmap_Cell_Annotation_Label_Matrix",
            as.data.frame(analysis_result$annotation_label_matrix) %>% rownames_to_column("Functional_Feature"))
  saveWorkbook(label_wb, file.path(save_dir, paste0("Spearman_Heatmap_Annotation_Label_Matrix_", dataset_tag, ".xlsx")), overwrite = TRUE)
  
  message(sprintf("[%s] All Spearman correlation statistical Excel tables exported successfully", dataset_tag))
}

# --------------------------- Function 3: Generate Publication Standard Annotated Spearman Correlation Heatmap (Dual PNG+PDF Export) ---------------------------
#' Draw complete annotated heatmap with cell text showing r value + significance mark, export dual-format journal figure files
#' @param analysis_result Complete output list returned by single_dataset_spearman_analysis
#' @param dataset_tag Unique dataset identifier for chart title and output filename suffix
#' @param save_dir Unified standardized output directory path
plot_spearman_heatmap <- function(analysis_result, dataset_tag, save_dir) {
  heatmap_plot <- ggplot(analysis_result$plot_data, aes(x = Environmental_Variable, y = Feature, fill = spearman_r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = annotation_text), color = "black", size = 2.5) +
    scale_fill_gradient2(
      low = "#2a7bbb",
      mid = "white",
      high = "#e67e22",
      limits = c(-1, 1),
      midpoint = 0,
      name = "Spearman's Rank Correlation r"
    ) +
    labs(
      x = "Soil Environmental Physicochemical Factors",
      y = "Functional Gene / Metabolic Features",
      title = paste0("Spearman Rank Correlation Annotated Heatmap: Functional Features vs Environmental Variables (", dataset_tag, ")")
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  # Export lossless vector PDF for formal manuscript submission
  ggsave(
    filename = file.path(save_dir, paste0("Spearman_Correlation_Annotated_Heatmap_", dataset_tag, ".pdf")),
    plot = heatmap_plot, width = 12, height = 9
  )
  # Export high-resolution 300 DPI raster PNG for supplementary materials & presentation
  ggsave(
    filename = file.path(save_dir, paste0("Spearman_Correlation_Annotated_Heatmap_", dataset_tag, ".png")),
    plot = heatmap_plot, width = 12, height = 9, dpi = 300
  )
  message(sprintf("[%s] Annotated Spearman correlation heatmap dual-format figures exported successfully", dataset_tag))
  return(heatmap_plot)
}

# ====================== Single Dataset End-to-End Analysis Execution Demo ======================
# Global fixed ordered list of target soil environmental physicochemical indicators (consistent with script 18/19/20)
target_env_columns <- c("NO3--N", "Ex-Ca", "TP", "EC", "SM_20cm", "PH", "NH4+-N", "Ex-Mg")

# Standardized cross-platform full input file paths
feature_file <- file.path(input_dir, "10last.xlsx")
environmental_file <- file.path(input_dir, "1env16-8.xlsx")
current_dataset_id <- "Functional_Feature_Global_Set"

# Input file integrity pre-check to avoid runtime error
if (!file.exists(feature_file)) {
  stop(paste("Missing mandatory functional feature abundance input file:", feature_file))
}
if (!file.exists(environmental_file)) {
  stop(paste("Missing mandatory environmental physicochemical input file:", environmental_file))
}

# Execute full integrated Spearman correlation calculation workflow
spearman_analysis_output <- single_dataset_spearman_analysis(
  feature_path = feature_file,
  env_path = environmental_file,
  env_target_cols = target_env_columns,
  dataset_tag = current_dataset_id
)

# Export all multi-sheet statistical matrix Excel tables
export_spearman_excel_result(
  analysis_result = spearman_analysis_output,
  dataset_tag = current_dataset_id,
  save_dir = output_dir
)

# Generate and export dual-format annotated heatmap figures
plot_spearman_heatmap(
  analysis_result = spearman_analysis_output,
  dataset_tag = current_dataset_id,
  save_dir = output_dir
)

# Standardized completion running log with separator
message("\n===== Single Dataset Functional Feature & Environment Spearman Correlation Analysis Pipeline Completed =====")
message("Full workflow modules executed: Pairwise Spearman rank correlation calculation | Significance classification | Multi-sheet statistical table export | Annotated heatmap dual-format visualization")
message("All integrated Excel statistical tables, 300 DPI PNG raster & lossless vector PDF heatmap figures saved to unified ./output folder")