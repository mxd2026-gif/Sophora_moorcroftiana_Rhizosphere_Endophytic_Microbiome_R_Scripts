# ==============================================================================
# Script Name: 11_Differential_Genus_WelchT_StackedBar.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Pairwise two-sided Welch t-test differential genus abundance screening, BH-FDR multiple test correction,
#          generate stacked bar visualization of core / non-core differential taxa count across geographic groups
# Input: Four genus relative abundance Excel & core taxa occupancy annotation table stored in ./input
# Output: Per-community full pairwise statistics, FDR significant taxa table, differential count stacked bar PNG/PDF,
#         cross-community merged master summary Excel tables saved in unified ./output folder
# Dependencies: tidyverse, readxl, ggplot2, writexl, numDeriv
# Statistical standard: global unified significance threshold consistent with script 03~10
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold parameters (single modification global effective)
global_alpha <- 0.05
cutoff_log2FC <- 1
max_log2FC_truncate <- 15

# --------------------------- Auto install & load all required packages ---------------------------
required_packages <- c("tidyverse", "readxl", "ggplot2", "writexl", "numDeriv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# --------------------------- Cross-platform relative directory initialization ---------------------------
input_dir  <- "./input"
output_dir <- "./output"

# Auto create missing folders with status message
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# --------------------------- Fixed global sample & grouping metadata ---------------------------
# All sample IDs consistent with script 03~10
group_sample_mapping <- list(
  Shigatse = c("RRK1", "RRK2", "RRK3", "RRK4", "RRK5", "ERK1", "ERK2", "ERK3", "ERK4", "ERK5"),
  Lhasa    = c("RCN1", "RCN2", "RCN3", "RCN4", "RCN5", "ECN1", "ECN2", "ECN3", "ECN4", "ECN5"),
  Shannan  = c("RSN1", "RSN2", "RSN3", "RSN4", "RSN5", "ESN1", "ESN2", "ESN3", "ESN4", "ESN5"),
  Nyingchi = c("RML1", "RML2", "RML3", "RML4", "RML5", "EML1", "EML2", "EML3", "EML4", "EML5")
)
# Modular optimization suggestion: extract group_sample_mapping to shared global script to eliminate duplicate code

# Unified color palette consistent with all prior analysis scripts
group_fill_dark <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)
group_fill_light <- c(
  Shigatse = "#FFA07A",
  Lhasa    = "#90EE90",
  Shannan  = "#D8BFD8",
  Nyingchi = "#FFFFE0"
)

# Define all 6 pairwise geographic comparisons
comparison_labels <- c(
  "Shigatse vs Lhasa", "Shigatse vs Shannan", "Shigatse vs Nyingchi",
  "Lhasa vs Shannan", "Lhasa vs Nyingchi", "Shannan vs Nyingchi"
)
comparison_pairs <- lapply(comparison_labels, function(x) strsplit(x, " vs ")[[1]])
names(comparison_pairs) <- comparison_labels
elevation_order <- c("Shigatse","Lhasa","Shannan","Nyingchi")

# --------------------------- Standardized input file paths ---------------------------
abundance_input_files <- list(
  g_rb = file.path(input_dir, "Rhizosphere_Bacteria_Genus_Relative_Abundance.xlsx"),
  g_eb = file.path(input_dir, "Root_Endophytic_Bacteria_Genus_Relative_Abundance.xlsx"),
  g_rf = file.path(input_dir, "Rhizosphere_Fungi_Genus_Relative_Abundance.xlsx"),
  g_ef = file.path(input_dir, "Root_Endophytic_Fungi_Genus_Relative_Abundance.xlsx")
)
core_taxon_sheet_names <- c("g_rbhexin", "g_ebhexin", "g_rfhexin", "g_efhexin")
names(core_taxon_sheet_names) <- names(abundance_input_files)
core_taxon_dataset_path <- file.path(input_dir, "Core_Taxa_80Percent_Occupancy_Annotation.xlsx")

# Pre-check input file integrity
all_input_files <- c(unlist(abundance_input_files), core_taxon_dataset_path)
missing_files <- all_input_files[!file.exists(all_input_files)]
if (length(missing_files) > 0) {
  stop(paste0("Missing required input files:\n", paste(missing_files, collapse = "\n")))
}

# --------------------------- Core Function 1: Two-group Welch t-test statistical calculation ---------------------------
#' Calculate abundance mean, SD, log2FC, raw P-value via Welch t-test for pairwise group comparison
#' @param taxa_long_df Long-format abundance subset data of target two groups
#' @param groupA Name of primary comparison group
#' @param groupB Name of secondary comparison group
#' @return Tibble containing complete pairwise statistical indicators
single_pair_taxon_stat <- function(taxa_long_df, groupA, groupB) {
  grpA_abund <- taxa_long_df %>% filter(Group == groupA) %>% pull(Relative_Abundance)
  grpB_abund <- taxa_long_df %>% filter(Group == groupB) %>% pull(Relative_Abundance)
  
  meanA <- mean(grpA_abund, na.rm = TRUE)
  sdA <- sd(grpA_abund, na.rm = TRUE)
  meanB <- mean(grpB_abund, na.rm = TRUE)
  sdB <- sd(grpB_abund, na.rm = TRUE)
  
  mean_sd_A <- sprintf("%.2f±%.2f", meanA, sdA)
  mean_sd_B <- sprintf("%.2f±%.2f", meanB, sdB)
  
  # Avoid zero division error during fold change calculation
  meanA_safe <- ifelse(meanA == 0, 1e-8, meanA)
  meanB_safe <- ifelse(meanB == 0, 1e-8, meanB)
  fold_change <- meanA_safe / meanB_safe
  log2_fc <- log2(fold_change)
  log2_fc <- pmax(pmin(log2_fc, max_log2FC_truncate), -max_log2FC_truncate)
  
  # Two-sided Welch t-test with full error & warning capture
  t_test_result <- tryCatch(
    t.test(grpA_abund, grpB_abund, var.equal = FALSE),
    error = function(e) list(p.value = 1),
    warning = function(w) list(p.value = 1)
  )
  p_raw <- ifelse(is.na(t_test_result$p.value), 1, t_test_result$p.value)
  neg_log10_p <- -log10(p_raw + 1e-10)
  
  significance_symbol <- case_when(
    p_raw < 0.001 ~ "***",
    p_raw < 0.01 ~ "**",
    p_raw < global_alpha ~ "*",
    TRUE ~ "ns"
  )
  
  flag_p_sig <- p_raw < global_alpha
  flag_threshold_sig <- flag_p_sig & abs(log2_fc) > cutoff_log2FC
  flag_subthreshold_sig <- flag_p_sig & abs(log2_fc) > 0 & abs(log2_fc) <= cutoff_log2FC
  dominant_group_tag <- case_when(
    log2_fc > 0 ~ groupA,
    log2_fc < 0 ~ groupB,
    TRUE ~ "None"
  )
  abs_fc_over1_flag <- abs(log2_fc) > cutoff_log2FC
  fc_direction_tag <- case_when(
    log2_fc > cutoff_log2FC ~ "log2FC > 1",
    log2_fc < -cutoff_log2FC ~ "log2FC < -1",
    TRUE ~ "|log2FC| ≤ 1"
  )
  
  output_table <- tibble(
    Mean_GroupA = meanA,
    SD_GroupA = sdA,
    Mean_SD_GroupA = mean_sd_A,
    Mean_GroupB = meanB,
    SD_GroupB = sdB,
    Mean_SD_GroupB = mean_sd_B,
    Fold_Change = fold_change,
    log2FC = log2_fc,
    Abs_log2FC_Over_1 = abs_fc_over1_flag,
    log2FC_Direction = fc_direction_tag,
    P_Value = p_raw,
    neg_log10_P = neg_log10_p,
    Significance_Label = significance_symbol,
    Is_Significant_P = flag_p_sig,
    Is_Significant_FC_Threshold = flag_threshold_sig,
    Is_Significant_Subthreshold = flag_subthreshold_sig,
    Dominant_Habitat_Group = dominant_group_tag
  )
  return(output_table)
}

# --------------------------- Core Function 2: Full single community differential analysis pipeline ---------------------------
#' Complete differential abundance analysis for one single microbial community
#' @param community_id Short identifier of target community
#' @return List object storing all statistical results and plotting data
process_single_community <- function(community_id) {
  message(sprintf("\n===== Start processing community dataset: %s =====", community_id))
  # Step 1: Import genus abundance matrix, replace zero with tiny offset to avoid log crash
  abundance_raw <- read_excel(abundance_input_files[[community_id]]) %>%
    rename(Taxon_Full_ID = 1) %>%
    mutate(Taxon_Full_ID = as.character(Taxon_Full_ID)) %>%
    distinct(Taxon_Full_ID, .keep_all = TRUE)
  sample_column_names <- setdiff(colnames(abundance_raw), "Taxon_Full_ID")
  abundance_wide_clean <- abundance_raw %>% mutate(across(all_of(sample_column_names), ~ifelse(. == 0, 1e-8, .)))
  abundance_long <- abundance_wide_clean %>%
    pivot_longer(cols = all_of(sample_column_names), names_to = "SampleID", values_to = "Relative_Abundance") %>%
    mutate(Group = case_when(
      SampleID %in% group_sample_mapping$Shigatse ~ "Shigatse",
      SampleID %in% group_sample_mapping$Lhasa ~ "Lhasa",
      SampleID %in% group_sample_mapping$Shannan ~ "Shannan",
      SampleID %in% group_sample_mapping$Nyingchi ~ "Nyingchi",
      TRUE ~ "Exclude"
    )) %>% filter(Group != "Exclude")
  
  # Step 2: Import core taxon annotation sheet (80% occupancy threshold core taxa)
  core_taxon_table <- read_excel(core_taxon_dataset_path, sheet = core_taxon_sheet_names[[community_id]]) %>%
    rename(Taxon_Full_ID = 1, Taxon_Simple_Label = 2) %>%
    mutate(
      Taxon_Full_ID = as.character(Taxon_Full_ID),
      Taxon_Simple_Label = as.character(Taxon_Simple_Label),
      Is_Core_Taxon = TRUE
    )
  
  # Step3: Iterate all 6 pairwise group comparisons
  comparison_total_list <- list()
  full_detail_sheet_collection <- list()
  p_significant_subset_collection <- list()
  for (comp_tag in comparison_labels) {
    group_pair <- comparison_pairs[[comp_tag]]
    group_A <- group_pair[1]
    group_B <- group_pair[2]
    message(sprintf("  Running pairwise comparison: %s", comp_tag))
    
    # Calculate pairwise statistics per single taxon
    single_comp_result <- abundance_long %>%
      group_by(Taxon_Full_ID) %>%
      do(single_pair_taxon_stat(., groupA = group_A, groupB = group_B)) %>% ungroup() %>%
      left_join(core_taxon_table %>% select(Taxon_Full_ID, Taxon_Simple_Label, Is_Core_Taxon), by = "Taxon_Full_ID") %>%
      replace_na(list(Is_Core_Taxon = FALSE, Taxon_Simple_Label = "Unclassified")) %>%
      mutate(
        Genus_Short_Label = str_extract(Taxon_Full_ID, "(?<=g__)[^;]+"),
        Genus_Short_Label = ifelse(is.na(Genus_Short_Label) | Genus_Short_Label == "", Taxon_Simple_Label, Genus_Short_Label),
        Taxon_Label_Annotated = paste0(Genus_Short_Label, " (", Significance_Label, ")"),
        Is_Core_Significant = Is_Core_Taxon & Is_Significant_FC_Threshold,
        Comparison = comp_tag
      )
    # Independent BH-FDR correction for each pairwise comparison
    single_comp_result$FDR_BH_Corrected_P <- p.adjust(single_comp_result$P_Value, method = "BH")
    single_comp_result$Is_Significant_FDR <- single_comp_result$FDR_BH_Corrected_P < global_alpha
    
    p_significant_subset_collection[[comp_tag]] <- single_comp_result %>% filter(Is_Significant_P)
    # Merge raw sample abundance matrix into complete detail table
    full_detail_table <- left_join(single_comp_result, abundance_wide_clean, by = "Taxon_Full_ID")
    # Standard fixed column display order for clean output Excel
    standard_column_sequence <- c(
      "Taxon_Full_ID", "Taxon_Simple_Label", "Genus_Short_Label", "Is_Core_Taxon",
      sample_column_names,
      "Mean_GroupA", "SD_GroupA", "Mean_SD_GroupA",
      "Mean_GroupB", "SD_GroupB", "Mean_SD_GroupB",
      "Fold_Change", "log2FC", "Abs_log2FC_Over_1", "log2FC_Direction", "Dominant_Habitat_Group",
      "P_Value", "Is_Significant_P", "FDR_BH_Corrected_P", "Is_Significant_FDR",
      "neg_log10_P", "Significance_Label", "Is_Significant_FC_Threshold", "Is_Significant_Subthreshold",
      "Is_Core_Significant", "Taxon_Label_Annotated", "Comparison"
    )
    full_detail_table <- full_detail_table[, intersect(standard_column_sequence, colnames(full_detail_table))]
    full_detail_sheet_collection[[comp_tag]] <- full_detail_table
    comparison_total_list[[comp_tag]] <- single_comp_result
  }
  
  all_combined_comparison_data <- bind_rows(comparison_total_list)
  core_significant_taxon_dataset <- all_combined_comparison_data %>% filter(Is_Core_Significant)
  threshold_significant_taxon_dataset <- all_combined_comparison_data %>% filter(Is_Significant_FC_Threshold)
  fdr_significant_taxon_dataset <- all_combined_comparison_data %>% filter(Is_Significant_FDR)
  
  # Step4: Organize stacked bar plot input data (count core / non-core differential genera)
  bar_plot_raw_data <- threshold_significant_taxon_dataset %>%
    filter(!is.na(Dominant_Habitat_Group), Dominant_Habitat_Group != "None") %>%
    mutate(
      Dominant_Habitat_Group = factor(Dominant_Habitat_Group, levels = elevation_order),
      Bar_Sort_Index = case_when(
        Dominant_Habitat_Group == "Shigatse" & Comparison == "Shigatse vs Lhasa" ~ 1,
        Dominant_Habitat_Group == "Shigatse" & Comparison == "Shigatse vs Shannan" ~ 2,
        Dominant_Habitat_Group == "Shigatse" & Comparison == "Shigatse vs Nyingchi" ~ 3,
        Dominant_Habitat_Group == "Lhasa" & Comparison == "Shigatse vs Lhasa" ~ 4,
        Dominant_Habitat_Group == "Lhasa" & Comparison == "Lhasa vs Shannan" ~ 5,
        Dominant_Habitat_Group == "Lhasa" & Comparison == "Lhasa vs Nyingchi" ~ 6,
        Dominant_Habitat_Group == "Shannan" & Comparison == "Shigatse vs Shannan" ~ 7,
        Dominant_Habitat_Group == "Shannan" & Comparison == "Lhasa vs Shannan" ~ 8,
        Dominant_Habitat_Group == "Shannan" & Comparison == "Shannan vs Nyingchi" ~ 9,
        Dominant_Habitat_Group == "Nyingchi" & Comparison == "Shigatse vs Nyingchi" ~ 10,
        Dominant_Habitat_Group == "Nyingchi" & Comparison == "Lhasa vs Nyingchi" ~ 11,
        Dominant_Habitat_Group == "Nyingchi" & Comparison == "Shannan vs Nyingchi" ~ 12
      ),
      Bar_X_Axis_Label = paste0(Dominant_Habitat_Group, " - ", Comparison)
    ) %>%
    group_by(Bar_Sort_Index, Bar_X_Axis_Label, Dominant_Habitat_Group, Comparison) %>%
    summarise(
      Non_Core_Significant_Count = sum(!Is_Core_Taxon),
      Core_Significant_Count = sum(Is_Core_Taxon),
      Total_Significant_Count = n(),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Non_Core_Significant_Count, Core_Significant_Count), names_to = "Taxon_Category", values_to = "Taxon_Count") %>%
    mutate(
      Taxon_Category = factor(Taxon_Category, levels = c("Non_Core_Significant_Count", "Core_Significant_Count")),
      Fill_Color_Code = case_when(
        Taxon_Category == "Non_Core_Significant_Count" ~ group_fill_light[Dominant_Habitat_Group],
        Taxon_Category == "Core_Significant_Count" ~ group_fill_dark[Dominant_Habitat_Group]
      )
    ) %>% arrange(Bar_Sort_Index)
  bar_count_summary_table <- threshold_significant_taxon_dataset %>%
    filter(!is.na(Dominant_Habitat_Group), Dominant_Habitat_Group != "None") %>%
    group_by(Comparison, Dominant_Habitat_Group) %>%
    summarise(
      Non_Core_Significant_Count = sum(!Is_Core_Taxon),
      Core_Significant_Count = sum(Is_Core_Taxon),
      Total_Significant_Count = n(),
      .groups = "drop"
    ) %>% mutate(Community = community_id)
  
  return(list(
    community_id = community_id,
    full_pairwise_statistics = all_combined_comparison_data,
    core_significant_taxa = core_significant_taxon_dataset,
    threshold_significant_taxa = threshold_significant_taxon_dataset,
    fdr_significant_taxa = fdr_significant_taxon_dataset,
    p_significant_subset_collection = p_significant_subset_collection,
    full_detail_sheets = full_detail_sheet_collection,
    stacked_bar_plot_input = bar_plot_raw_data,
    bar_count_statistics = bar_count_summary_table
  ))
}

# --------------------------- Core Function3: Generate stacked bar chart ---------------------------
#' Draw stacked bar plot for differential taxon count distribution, semantic naming only
#' @param analysis_result Output list object from single community analysis function
#' @param output_folder Target directory to save figure files
plot_differential_taxa_stacked_bar <- function(analysis_result, output_folder) {
  comm_id <- analysis_result$community_id
  plot_input_df <- analysis_result$stacked_bar_plot_input
  bar_total_label_summary <- plot_input_df %>% group_by(Bar_X_Axis_Label, Bar_Sort_Index) %>% summarise(Total_Taxa = sum(Taxon_Count), .groups = "drop")
  
  target_plot <- ggplot(plot_input_df, aes(x = reorder(Bar_X_Axis_Label, Bar_Sort_Index), y = Taxon_Count)) +
    geom_col(aes(fill = Fill_Color_Code), position = "stack", color = "black", linewidth = 0.3) +
    geom_text(aes(label = Taxon_Count), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold", color = "white") +
    geom_text(data = bar_total_label_summary, aes(y = Total_Taxa + 0.3, label = Total_Taxa), size = 4, fontface = "bold", color = "black", inherit.aes = FALSE) +
    scale_fill_identity() +
    labs(
      x = "Geographic Elevation Group - Pairwise Comparison",
      y = "Number of Differentially Abundant Genera",
      title = paste0("Differential Genus Count Stacked Bar Visualization: ", comm_id)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  # Export figure with semantic community name only
  ggsave(
    filename = file.path(output_folder, paste0("Differential_Taxa_StackedBar_", comm_id, ".pdf")),
    plot = target_plot, width = 14, height = 7, dpi = 600
  )
  ggsave(
    filename = file.path(output_folder, paste0("Differential_Taxa_StackedBar_", comm_id, ".png")),
    plot = target_plot, width = 14, height = 7, dpi = 300
  )
  message(sprintf("Stacked bar plot exported successfully for community: %s", comm_id))
  return(target_plot)
}

# --------------------------- Core Function4: Export statistical Excel tables ---------------------------
#' Export all statistical output tables with descriptive semantic filenames
#' @param analysis_result Single community full analysis result list
#' @param output_folder Target directory for table exports
export_community_statistic_tables <- function(analysis_result, output_folder) {
  comm_id <- analysis_result$community_id
  # 1. Complete detail dataset with raw sample abundance and full FDR corrected indicators
  write_xlsx(
    analysis_result$full_detail_sheets,
    file.path(output_folder, paste0("Differential_Taxa_Full_Detail_Dataset_", comm_id, ".xlsx"))
  )
  # 2. Raw P-value significant taxa subset summary
  write_xlsx(
    analysis_result$p_significant_subset_collection,
    file.path(output_folder, paste0("Raw_P_Value_Significant_Taxa_Summary_", comm_id, ".xlsx"))
  )
  # 3. FDR BH corrected significant taxa table
  write_xlsx(
    list(FDR_Significant_Taxa_Dataset = analysis_result$fdr_significant_taxa),
    file.path(output_folder, paste0("FDR_BH_Corrected_Significant_Taxa_", comm_id, ".xlsx"))
  )
  # 4. Separate sheets for all pairwise / core sig / threshold sig taxa
  sheet_all_pair <- list()
  sheet_core_sig <- list()
  sheet_threshold_sig <- list()
  for (cmp_tag in comparison_labels) {
    sheet_all_pair[[cmp_tag]] <- analysis_result$full_pairwise_statistics %>% filter(Comparison == cmp_tag)
    sheet_core_sig[[cmp_tag]] <- analysis_result$core_significant_taxa %>% filter(Comparison == cmp_tag)
    sheet_threshold_sig[[cmp_tag]] <- analysis_result$threshold_significant_taxa %>% filter(Comparison == cmp_tag)
  }
  write_xlsx(sheet_all_pair, file.path(output_folder, paste0("All_Pairwise_Taxon_Statistical_Results_", comm_id, ".xlsx")))
  write_xlsx(sheet_core_sig, file.path(output_folder, paste0("Core_Differential_Significant_Taxa_", comm_id, ".xlsx")))
  write_xlsx(sheet_threshold_sig, file.path(output_folder, paste0("Threshold_Filtered_Significant_Taxa_", comm_id, ".xlsx")))
  message(sprintf("All statistical tables exported for community: %s", comm_id))
}

# ======================== Main Execution Pipeline ========================
# Global storage list for cross-community merged master summary tables
global_p_sig_merged <- list()
global_fdr_sig_merged <- list()
global_bar_stat_master <- tibble()

target_community_list <- names(abundance_input_files)
for (target_comm in target_community_list) {
  single_comm_output <- process_single_community(community_id = target_comm)
  # Generate stacked bar visualization
  plot_differential_taxa_stacked_bar(analysis_result = single_comm_output, output_folder = output_dir)
  # Export all semantic statistical tables
  export_community_statistic_tables(analysis_result = single_comm_output, output_folder = output_dir)
  # Merge cross-community summary datasets
  global_p_sig_merged[[target_comm]] <- bind_rows(single_comm_output$p_significant_subset_collection)
  global_fdr_sig_merged[[target_comm]] <- single_comm_output$fdr_significant_taxa
  global_bar_stat_master <- bind_rows(global_bar_stat_master, single_comm_output$bar_count_statistics)
}

# Export cross-community merged master summary tables
write_xlsx(global_p_sig_merged, file.path(output_dir, "Master_Dataset_All_Communities_RawP_Significant_Taxa.xlsx"))
write_xlsx(global_fdr_sig_merged, file.path(output_dir, "Master_Dataset_All_Communities_FDR_BH_Significant_Taxa.xlsx"))
write_xlsx(global_bar_stat_master, file.path(output_dir, "Differential_Taxa_Count_Comparison_Master_Statistics.xlsx"))

message("\n==================== All Differential Abundance Analysis Tasks Completed ====================")
message("All stacked bar visualizations, full pairwise statistics, FDR filtered tables saved in ./output/")
message("All filenames use pure semantic English labels, no hard-coded figure / table serial numbers")