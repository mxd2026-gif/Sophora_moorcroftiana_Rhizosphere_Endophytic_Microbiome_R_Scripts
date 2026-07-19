# ==============================================================================
# Script Name: 14_Taxon_Contribution_Screening_Analysis.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Screen taxa by average contribution threshold, statistical classification & grouped bar plot visualization
# Input: OTU abundance matrix Excel table stored in unified ./input directory
# Output: High-res PNG raster & vector PDF bar chart, multi-sheet statistical Excel tables saved in ./output
# Dependencies: tidyverse, readxl, writexl, numDeriv
# Standardization: Fully aligned with script 01~13 unified path, global parameters, palette & output naming rules
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold variable (consistent with full pipeline scripts)
global_alpha <- 0.05

# -------------------------- Install & Load All Required Packages --------------------------
required_packages <- c("tidyverse", "readxl", "writexl", "numDeriv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# -------------------------- Cross-platform Relative Path Configuration --------------------------
input_dir  <- "./input"
output_dir <- "./output"

# Auto create missing folders with status prompt log
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# -------------------------- Global Fixed Parameter Setting --------------------------
# Taxon contribution screening threshold (single modification global effective)
contribution_threshold <- 95
# Unified global color palette consistent with script 03~13
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)

# -------------------------- Reusable Independent Core Function --------------------------
#' Process single community dataset for contribution screening, statistics and visualization
#' @param input_file_path Full path of input excel file
#' @param sheet_name Target sheet name in input excel
#' @param target_group_cols Character vector of replicate column names for calculating average contribution
#' @param community_id Unique semantic identifier of target community for standardized output naming
#' @param fill_color Custom hex fill color for high contribution group bar
#' @return List contains ggplot barplot object, multi-sheet excel list and classification count summary table
process_single_community <- function(input_file_path,
                                     sheet_name,
                                     target_group_cols,
                                     community_id,
                                     fill_color) {
  # Read full raw OTU dataset and retain all original taxon & sample columns
  raw_otu_data <- read_excel(input_file_path, sheet = sheet_name)
  
  # Check missing replicate columns and terminate with clear error prompt
  missing_columns <- setdiff(target_group_cols, colnames(raw_otu_data))
  if (length(missing_columns) > 0) {
    stop(sprintf("Missing required replicate columns in sheet [%s]: %s", sheet_name, paste(missing_columns, collapse = ", ")))
  }
  
  # Extract replicate abundance matrix, replace NA missing values with 0
  replicate_data <- raw_otu_data %>% select(all_of(target_group_cols))
  replicate_data <- replicate_data %>% mutate(across(everything(), ~replace_na(.x, 0)))
  
  # Calculate mean contribution across all biological replicates per single taxon
  average_contribution_value <- rowMeans(replicate_data, na.rm = TRUE)
  
  # Two-category classification based on preset global contribution threshold
  classification_label <- ifelse(
    average_contribution_value >= contribution_threshold,
    paste0("Average_Contribution_≥", contribution_threshold),
    paste0("Average_Contribution_<", contribution_threshold)
  )
  
  # Merge mean contribution index & classification tag back to full raw dataset
  total_analysis_dataset <- raw_otu_data %>%
    mutate(
      Average_Contribution = average_contribution_value,
      Classification_Group = classification_label
    )
  
  # Split full dataset into low / high contribution subgroups
  low_contribution_dataset <- total_analysis_dataset %>%
    filter(Classification_Group == paste0("Average_Contribution_<", contribution_threshold))
  high_contribution_dataset <- total_analysis_dataset %>%
    filter(Classification_Group == paste0("Average_Contribution_≥", contribution_threshold))
  
  # Count taxon quantity for each classification group, lock fixed axis display order
  classification_count_table <- total_analysis_dataset %>%
    count(Classification_Group, name = "Taxon_Count") %>%
    mutate(Classification_Group = factor(
      Classification_Group,
      levels = c(
        paste0("Average_Contribution_<", contribution_threshold),
        paste0("Average_Contribution_≥", contribution_threshold)
      )
    ))
  
  # Generate journal-standard grouped bar plot for taxon count statistics
  contribution_bar_plot <- ggplot(
    classification_count_table,
    aes(x = Classification_Group, y = Taxon_Count, fill = Classification_Group)
  ) +
    geom_col(width = 0.6, alpha = 0.9) +
    geom_text(aes(label = Taxon_Count), vjust = -0.6, size = 5, fontface = "bold") +
    scale_fill_manual(values = c("gray80", fill_color)) +
    labs(
      title = paste0(community_id, " | Taxon Average Contribution Classification Statistics"),
      x = "Contribution Classification Group",
      y = "Total Number of Taxa / OTUs"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    ) +
    ylim(0, max(classification_count_table$Taxon_Count) * 1.2)
  
  # Organize three standardized sheets for integrated output Excel
  excel_result_sheets <- list(
    Total_Full_Analysis_Dataset = total_analysis_dataset,
    Low_Contribution_Taxa_Subset = low_contribution_dataset,
    High_Contribution_Taxa_Subset = high_contribution_dataset
  )
  
  return(list(
    plot_object = contribution_bar_plot,
    excel_output = excel_result_sheets,
    count_result = classification_count_table
  ))
}

# -------------------------- Input File Pre-Check & Main Execution --------------------------
# Standardized cross-platform input file full path
input_excel_path <- file.path(input_dir, "3vs1REContribute611+9all.xlsx")
# Stop script if required input table missing
if (!file.exists(input_excel_path)) {
  stop(paste("Missing required input OTU abundance table:", input_excel_path))
}

# Run full contribution screening pipeline for rhizosphere microbial community
analysis_result <- process_single_community(
  input_file_path = input_excel_path,
  sheet_name = "3vs1R原9",
  target_group_cols = c("Group_SCNR", "Group_SRKR", "Group_SSNR"),
  community_id = "Rhizosphere_Microbial_Community",
  fill_color = group_color_config[["Lhasa"]]
)

# Semantic standardized output file prefix (no hard-coded serial numbers)
output_file_prefix = "Rhizosphere_Taxon_Average_Contribution_Screening_Result"

# Export high-resolution raster PNG & lossless vector PDF figures
ggsave(
  filename = file.path(output_dir, paste0(output_file_prefix, "_Classification_Barplot.png")),
  plot = analysis_result$plot_object,
  width = 10,
  height = 7,
  dpi = 300
)
ggsave(
  filename = file.path(output_dir, paste0(output_file_prefix, "_Classification_Barplot.pdf")),
  plot = analysis_result$plot_object,
  width = 10,
  height = 7
)

# Export multi-sheet integrated statistical Excel table
write_xlsx(
  x = analysis_result$excel_output,
  path = file.path(output_dir, paste0(output_file_prefix, "_Full_Statistical_Dataset.xlsx"))
)

# Console standardized running log output
message("\n===== Taxon contribution screening analysis completed successfully =====")
message(sprintf("Target community: %s", "Rhizosphere_Microbial_Community"))
message(sprintf("Taxa count below 95%% contribution threshold: %d", analysis_result$count_result$Taxon_Count[1]))
message(sprintf("Taxa count ≥ 95%% contribution threshold: %d", analysis_result$count_result$Taxon_Count[2]))
message("All bar chart figures & multi-sheet Excel statistical tables saved to unified ./output folder")