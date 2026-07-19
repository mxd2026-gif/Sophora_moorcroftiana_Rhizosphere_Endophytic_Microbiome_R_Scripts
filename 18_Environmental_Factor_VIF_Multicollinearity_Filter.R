# ==============================================================================
# Script Name: 18_Environmental_Factor_VIF_Multicollinearity_Filter.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Stepwise iterative VIF multicollinearity elimination for environmental factor matrix,
#          dual threshold screening (VIF=10 / VIF=5), generate VIF statistics tables & comparative bar plots,
#          output filtered standardized environmental matrix for subsequent RDA/Mantel correlation analysis
# Input: 1env16.xlsx environmental factor + sample grouping annotation table stored in unified ./input folder
# Output: Multi-sheet VIF iteration statistical Excel tables, dual-format PNG+PDF VIF bar charts,
#         filtered environmental dataset RDS intermediate file saved in ./output
# Dependencies: vegan, car, dplyr, readxl, ggplot2, stringr, tibble, writexl, tidyr, numDeriv
# Standardization: Fully aligned with script 01~17 unified global parameters, cross-platform path rules, journal dual-format output specs
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold consistent with full serial analysis pipeline
global_alpha <- 0.05

# --------------------------- Install and load all dependent packages ---------------------------
required_packages <- c("vegan", "car", "dplyr", "readxl", "ggplot2",
                       "stringr", "tibble", "writexl", "tidyr", "numDeriv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# --------------------------- Cross-platform relative path configuration ---------------------------
input_dir  <- "./input"
output_dir <- "./output"

# Auto create missing standard folders with running status prompt
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# --------------------------- Global fixed unified experimental parameters ---------------------------
# Fixed complete ordered environmental factor list (soil physicochemical indicators)
ENV_ORDER <- c("EC", "PH", "NO3--N", "NH4+-N", "TN", "TP", "AP", "TK", "AK",
               "OM", "OC", "AFe", "Ex-Ca", "Ex-Mg", "ST_20cm", "SM_20cm")
# Two-stage VIF screening thresholds (industry standard multicollinearity cutoff)
VIF_THRESHOLD_10 <- 10
VIF_THRESHOLD_5  <- 5
# Global unified four geographic sampling group color palette (consistent with script 03~17)
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)
group_factor_level <- names(group_color_config)

# --------------------------- Function 1: Standard environmental data preprocessing Z-score normalization ---------------------------
#' Normalize raw environmental data with Z-score scaling, separate sample metadata and numeric matrix
#' @param env_df Raw input data frame containing SampleID, Group and environmental factor columns
#' @return List contains scaled standardized matrix, sample grouping metadata, raw unstandardized numeric data, sample ID vector
preprocess_environment_data <- function(env_df) {
  sample_id_vector <- env_df$SampleID
  group_factor <- factor(env_df$Group, levels = group_factor_level)
  
  env_raw_numeric <- env_df[, ENV_ORDER, drop = FALSE]
  env_scaled_matrix <- as.data.frame(scale(env_raw_numeric))
  rownames(env_scaled_matrix) <- sample_id_vector
  
  sample_metadata <- data.frame(SampleID = sample_id_vector,
                                Group = group_factor,
                                row.names = sample_id_vector)
  
  return(list(
    env_scaled = env_scaled_matrix,
    env_meta   = sample_metadata,
    env_raw    = env_raw_numeric,
    sample_id  = sample_id_vector
  ))
}

# --------------------------- Function 2: Check completeness of required environmental columns ---------------------------
#' Verify all mandatory environmental variables exist in input table, terminate if missing
#' @param df Full input raw data frame
#' @param dataset_label Unique identification tag for current dataset console log
check_required_columns <- function(df, dataset_label) {
  missing_columns <- setdiff(ENV_ORDER, colnames(df))
  if (length(missing_columns) > 0) {
    stop(sprintf("Column validation failed for dataset [%s], missing mandatory environmental columns: %s",
                 dataset_label, paste(missing_columns, collapse = ",")))
  }
  message(sprintf("Mandatory environmental column validation passed for dataset: %s", dataset_label))
}

# --------------------------- Function 3: Calculate full pairwise VIF long-format table ---------------------------
#' Compute Variance Inflation Factor for each predictor variable in linear regression model
#' @param numeric_df Standardized complete numeric environmental matrix
#' @return Long-format table recording VIF value of each predictor under different response variables
calculate_vif_table <- function(numeric_df) {
  variable_list <- names(numeric_df)
  vif_result_list <- list()
  
  for (response_var in variable_list) {
    predictor_vars <- setdiff(variable_list, response_var)
    resp_quoted <- paste0("`", response_var, "`")
    pred_quoted <- paste0("`", predictor_vars, "`")
    formula_obj <- as.formula(paste(resp_quoted, "~", paste(pred_quoted, collapse = " + ")))
    
    lm_model <- lm(formula_obj, data = numeric_df)
    vif_values <- car::vif(lm_model)
    
    temp_df <- data.frame(
      response = response_var,
      predictor = names(vif_values),
      VIF = as.numeric(vif_values),
      stringsAsFactors = FALSE
    )
    vif_result_list[[response_var]] <- temp_df
  }
  return(bind_rows(vif_result_list))
}

# --------------------------- Function 4: Core stepwise iterative VIF elimination algorithm ---------------------------
#' Iteratively remove variables with highest maximum VIF until all factors meet dual threshold standards
#' @param scaled_env_df Z-score standardized environmental numeric matrix
#' @param th10 First-stage loose screening VIF threshold, default 10
#' @param th5 Second-stage strict screening VIF threshold, default 5
#' @param dataset_tag Unique tag for current dataset log and file naming
#' @return List storing iteration detail records, removed variable log, filtered environmental matrices under two thresholds
stepwise_vif_filter <- function(scaled_env_df, th10 = 10, th5 = 5, dataset_tag = "Dataset") {
  current_env_data <- scaled_env_df
  vif_iter_detail <- list()
  vif_summary_record <- list()
  variable_removed_log <- data.frame()
  
  iteration_step <- 1
  flag_vif10_pass <- FALSE
  flag_vif5_pass  <- FALSE
  env_vif_filter10 <- NULL
  env_vif_filter5  <- NULL
  
  while (TRUE) {
    vif_long_data <- calculate_vif_table(current_env_data)
    vif_summary_data <- vif_long_data %>%
      group_by(predictor) %>%
      summarise(max_VIF = max(VIF, na.rm = TRUE),
                mean_VIF = mean(VIF, na.rm = TRUE),
                .groups = "drop") %>%
      arrange(desc(max_VIF))
    
    vif_iter_detail[[paste0("Step_", iteration_step)]] <- vif_long_data
    vif_summary_record[[paste0("Step_", iteration_step)]] <- vif_summary_data
    
    max_vif_value <- vif_summary_data$max_VIF[1]
    max_vif_variable <- vif_summary_data$predictor[1]
    
    message(sprintf("[%s] Iteration %d | Current maximum VIF: %.3f | Variable to eliminate: %s",
                    dataset_tag, iteration_step, round(max_vif_value, 3), max_vif_variable))
    
    if (all(vif_summary_data$max_VIF < th10) && !flag_vif10_pass) {
      flag_vif10_pass <- TRUE
      env_vif_filter10 <- current_env_data
      message(sprintf("[%s] Iteration %d: All environmental factors satisfy VIF < %d loose threshold", dataset_tag, iteration_step, th10))
    }
    
    if (all(vif_summary_data$max_VIF < th5) && !flag_vif5_pass) {
      flag_vif5_pass <- TRUE
      env_vif_filter5 <- current_env_data
      message(sprintf("[%s] Iteration %d: All environmental factors satisfy strict VIF < %d threshold, iteration terminated",
                      dataset_tag, iteration_step, th5))
      break
    }
    
    current_env_data <- current_env_data %>% select(-all_of(max_vif_variable))
    variable_removed_log <- bind_rows(variable_removed_log, data.frame(
      step = iteration_step,
      removed_variable = max_vif_variable,
      max_vif = max_vif_value,
      stringsAsFactors = FALSE
    ))
    iteration_step <- iteration_step + 1
    
    if (ncol(current_env_data) <= 2) {
      warning(sprintf("[%s] Remaining environmental factors ≤ 2, forced termination of iterative filtering", dataset_tag))
      break
    }
  }
  
  final_vif_summary <- vif_summary_record[[length(vif_summary_record)]]
  
  return(list(
    vif_iteration_detail = vif_iter_detail,
    vif_summary_history = vif_summary_record,
    variable_removal_log = variable_removed_log,
    final_vif_statistics = final_vif_summary,
    env_data_vif10 = env_vif_filter10,
    env_data_vif5  = env_vif_filter5,
    total_iteration_steps = iteration_step - 1
  ))
}

# --------------------------- Function 5: VIF statistical bar chart batch export function (dual PNG+PDF) ---------------------------
#' Generate three types of standardized VIF publication charts: final VIF distribution, before-after comparison, iteration start/end comparison
#' @param vif_filter_result Output object returned by stepwise_vif_filter core function
#' @param raw_scaled_env Complete standardized environmental matrix before filtering
#' @param dataset_label Unique dataset identifier for chart title and output file prefix
#' @param save_dir Unified standardized output folder path
#' @param th10 Loose VIF threshold value
#' @param th5 Strict VIF threshold value
plot_vif_visual_results <- function(vif_filter_result, raw_scaled_env, dataset_label, save_dir, th10 = 10, th5 = 5) {
  vif_raw_all <- calculate_vif_table(raw_scaled_env)
  vif_before_filter <- vif_raw_all %>%
    group_by(predictor) %>%
    summarise(max_VIF = max(VIF, na.rm = TRUE), .groups = "drop") %>%
    mutate(filter_stage = "Before Filter")
  
  vif_after_filter <- vif_filter_result$final_vif_statistics %>%
    select(predictor, max_VIF) %>%
    mutate(filter_stage = "After Filter")
  
  # Chart 1: Final VIF distribution bar chart after strict filtering
  p_final_vif <- ggplot(vif_filter_result$final_vif_statistics,
                        aes(x = reorder(predictor, max_VIF), y = max_VIF)) +
    geom_col(fill = "#2C7FB8", width = 0.7) +
    geom_hline(yintercept = th10, linetype = "dashed", color = "red", linewidth = 1) +
    geom_hline(yintercept = th5,  linetype = "dashed", color = "orange", linewidth = 1) +
    coord_flip() +
    labs(
      title = paste0("Final Environmental Factor VIF Distribution - ", dataset_label),
      x = "Soil Environmental Indicators",
      y = "Maximum VIF Multicollinearity Coefficient"
    ) +
    theme_bw(base_size = 11)
  
  ggsave(file.path(save_dir, paste0("VIF_Final_Distribution_", dataset_label, ".png")),
         p_final_vif, width = 7, height = 6, dpi = 300)
  ggsave(file.path(save_dir, paste0("VIF_Final_Distribution_", dataset_label, ".pdf")),
         p_final_vif, width = 7, height = 6)
  
  # Chart 2: VIF before and after filtering grouped comparison bar chart
  vif_compare_dataset <- bind_rows(vif_before_filter, vif_after_filter)
  p_vif_compare <- ggplot(vif_compare_dataset,
                          aes(x = predictor, y = max_VIF, fill = filter_stage)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_hline(yintercept = th10, linetype = "dashed", color = "red", linewidth = 1) +
    geom_hline(yintercept = th5,  linetype = "dashed", color = "orange", linewidth = 1) +
    coord_flip() +
    labs(
      title = paste0("VIF Multicollinearity Comparison Before & After Filtering - ", dataset_label),
      x = "Soil Environmental Indicators",
      y = "Maximum VIF Multicollinearity Coefficient",
      fill = "Analysis Stage"
    ) +
    scale_fill_manual(values = c("Before Filter" = "#4472C4", "After Filter" = "#FF6347")) +
    theme_bw(base_size = 11)
  
  ggsave(file.path(save_dir, paste0("VIF_Before_After_Comparison_", dataset_label, ".png")),
         p_vif_compare, width = 8.5, height = 7, dpi = 300)
  ggsave(file.path(save_dir, paste0("VIF_Before_After_Comparison_", dataset_label, ".pdf")),
         p_vif_compare, width = 8.5, height = 7)
  
  # Chart3: VIF distribution of first iteration and last termination iteration
  iteration_names <- names(vif_filter_result$vif_summary_history)
  selected_iterations <- c(iteration_names[1], iteration_names[length(iteration_names)])
  
  for (iter in selected_iterations) {
    iter_num <- gsub("Step_", "", iter)
    iter_data <- vif_filter_result$vif_summary_history[[iter]]
    p_iter_vif <- ggplot(iter_data, aes(x = reorder(predictor, max_VIF), y = max_VIF)) +
      geom_col(fill = "#4472C4", width = 0.7) +
      geom_hline(yintercept = th10, linetype = "dashed", color = "red") +
      geom_hline(yintercept = th5,  linetype = "dashed", color = "orange") +
      coord_flip() +
      labs(
        title = paste0("VIF Iteration ", iter_num, " Distribution - ", dataset_label),
        x = "Soil Environmental Indicators",
        y = "Maximum VIF Multicollinearity Coefficient"
      ) +
      theme_bw(base_size = 10)
    
    ggsave(file.path(save_dir, paste0("VIF_Iteration_", iter_num, "_", dataset_label, ".png")),
           p_iter_vif, width = 7, height = 6, dpi = 300)
    ggsave(file.path(save_dir, paste0("VIF_Iteration_", iter_num, "_", dataset_label, ".pdf")),
           p_iter_vif, width = 7, height = 6)
  }
}

# --------------------------- Core single dataset integrated analysis wrapper function ---------------------------
#' End-to-end full VIF multicollinearity screening pipeline for single sheet environmental dataset
#' @param excel_file_path Full cross-platform standardized input excel file path
#' @param target_sheet Character name of target analysis sheet in excel
#' @param dataset_unique_tag Unique semantic identifier for all output table/figure filenames
#' @param output_folder Unified standardized output directory path
#' @return List containing complete VIF statistical results and standardized preprocessed environmental data
single_dataset_vif_analysis <- function(excel_file_path, target_sheet, dataset_unique_tag, output_folder) {
  raw_env_data <- read_excel(excel_file_path, sheet = target_sheet)
  check_required_columns(raw_env_data, dataset_unique_tag)
  
  processed_env <- preprocess_environment_data(raw_env_data)
  
  vif_result <- stepwise_vif_filter(
    scaled_env_df = processed_env$env_scaled,
    th10 = VIF_THRESHOLD_10,
    th5 = VIF_THRESHOLD_5,
    dataset_tag = dataset_unique_tag
  )
  
  # Export all multi-sheet VIF statistical tables to Excel
  write_xlsx(vif_result$vif_iteration_detail,
             file.path(output_folder, paste0("VIF_Iteration_Detail_", dataset_unique_tag, ".xlsx")))
  write_xlsx(vif_result$variable_removal_log,
             file.path(output_folder, paste0("VIF_Variable_Removal_Log_", dataset_unique_tag, ".xlsx")))
  write_xlsx(vif_result$final_vif_statistics,
             file.path(output_folder, paste0("VIF_Final_Statistics_", dataset_unique_tag, ".xlsx")))
  
  # Batch generate and export all standardized VIF bar chart figures
  plot_vif_visual_results(
    vif_filter_result = vif_result,
    raw_scaled_env = processed_env$env_scaled,
    dataset_label = dataset_unique_tag,
    save_dir = output_folder,
    th10 = VIF_THRESHOLD_10,
    th5 = VIF_THRESHOLD_5
  )
  
  return(list(
    vif_output = vif_result,
    env_processed_data = processed_env
  ))
}

# ====================== Single Sheet Task Execution Demo ======================
# Standardized full cross-platform input file path
input_excel <- file.path(input_dir, "1env16.xlsx")
# Terminate script if required environmental input table missing
if (!file.exists(input_excel)) {
  stop(paste("Missing mandatory environmental factor input file:", input_excel))
}

target_analysis_sheet <- "Rhizosphere"
dataset_id <- "Single_Rhizosphere_Environmental_Dataset"

# Execute complete end-to-end VIF multicollinearity filtering workflow
analysis_result <- single_dataset_vif_analysis(
  excel_file_path = input_excel,
  target_sheet = target_analysis_sheet,
  dataset_unique_tag = dataset_id,
  output_folder = output_dir
)

# Save full filtered standardized environmental dataset as intermediate RDS for downstream RDA/Mantel analysis
save(analysis_result, file = file.path(output_dir, "Filtered_Environmental_Dataset_VIF_Removed.rds"))

# Standardized completion running log with separator
message("\n===== Single Dataset Environmental Factor VIF Multicollinearity Screening Completed =====")
message(sprintf("Target analysis sheet: %s", target_analysis_sheet))
message("All VIF iteration statistical Excel tables, dual-format PNG/PDF bar charts and filtered environmental RDS file saved to unified ./output folder")