# ==============================================================================
# Script Name: 04_WithinGroup_Distance_Boxplot.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Bray-Curtis within-group distance boxplot + PERMANOVA global test + pairwise R² calculation
# Input: Os_rb.xlsx, Os_eb.xlsx, Os_rf.xlsx, Os_ef.xlsx (taxon abundance matrix in ./input)
# Output: Within-group distance boxplot PNG/PDF, PERMANOVA & Kruskal-Wallis Excel result in ./output
# Dependencies: tidyverse, vegan, readxl, agricolae, writexl, gridExtra, picante, ape, ggsignif, ggpubr
# ==============================================================================
options(scipen = 999, digits = 4)
# Global significance threshold, easy to adjust uniformly
global_alpha <- 0.05

# Install and load required packages automatically
if (!requireNamespace("picante", quietly = TRUE)) install.packages("picante")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("ggsignif", quietly = TRUE)) install.packages("ggsignif")
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
if (!requireNamespace("agricolae", quietly = TRUE)) install.packages("agricolae")

library(tidyverse)
library(vegan)
library(readxl)
library(broom)
library(agricolae)
library(writexl)
library(gridExtra)
library(grid)
library(picante)
library(ape)
library(ggsignif)
library(ggpubr)

# Create standard input & output folders
input_dir <- "./input"
output_dir <- "./output"
if (!dir.exists(input_dir)) dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

process_community_analysis <- function(data_path, groups_list, plot_title, output_prefix, habitat_name) {
  # Import raw abundance matrix, clean negative / NA values to 0
  data <- read_excel(data_path) %>% rename(Taxon = 1)
  data_clean <- data %>% mutate(across(-Taxon, ~ ifelse(. < 0 | is.na(.), 0, .)))
  
  # Filter out taxa with all-zero abundance across all samples
  comm_matrix_raw <- data_clean %>%
    filter(if_any(-Taxon, ~ . > 0)) %>%
    column_to_rownames("Taxon") %>%
    t()
  comm_matrix <- comm_matrix_raw[rowSums(comm_matrix_raw) > 0, ]
  
  # Match sample grouping information
  sample_info <- tibble(Sample = rownames(comm_matrix)) %>%
    mutate(
      Group = case_when(
        Sample %in% groups_list$Shigatse ~ "Shigatse",
        Sample %in% groups_list$Lhasa ~ "Lhasa",
        Sample %in% groups_list$Shannan ~ "Shannan",
        Sample %in% groups_list$Nyingchi ~ "Nyingchi",
        TRUE ~ "Other"
      )
    ) %>%
    filter(Group != "Other") %>%
    column_to_rownames("Sample")
  sample_info$Group <- factor(sample_info$Group, levels = c("Shigatse", "Lhasa", "Shannan", "Nyingchi"))
  
  # Calculate Bray-Curtis distance matrix
  dist_matrix <- vegdist(comm_matrix, method = "bray")
  # PERMANOVA global test; permutations = 999 is standard, change to 9999 for higher precision
  permanova_result <- adonis2(dist_matrix ~ Group, data = sample_info, permutations = 999)
  global_R2 <- round(permanova_result$R2[1], 3)
  global_P <- round(permanova_result$Pr[F][1], 4)
  global_F <- round(permanova_result$F[1], 2)
  
  # Reshape distance matrix to long format for within-group boxplot
  dist_df <- as.matrix(dist_matrix) %>%
    as.data.frame() %>%
    rownames_to_column("Sample_Compare") %>%
    pivot_longer(cols = -Sample_Compare, names_to = "Sample", values_to = "BrayCurtis_Distance") %>%
    left_join(sample_info %>% rownames_to_column("Sample"), by = "Sample") %>%
    left_join(sample_info %>% rownames_to_column("Sample_Compare") %>% rename(Group2 = Group), by = "Sample_Compare") %>%
    filter(Sample == Group, BrayCurtis_Distance) %>%
    select(Group, BrayCurtis_Distance)
  
  # Kruskal test for within-group distance differences
  tukey_data <- dist_df %>% select(Group, BrayCurtis_Distance)
  sig_letters <- kruskal(tukey_data$BrayCurtis_Distance, trt = tukey_data$Group, alpha = global_alpha)
  sig_letters <- sig_letters$groups %>% rownames_to_column("Group") %>% rename(letters = groups)
  
  # Draw within-group distance boxplot
  boxplot_data <- dist_df %>%
    group_by(Group) %>%
    summarise(max_box = max(BrayCurtis_Distance[BrayCurtis_Distance <= boxplot.stats(BrayCurtis_Distance)$stats[5] * 1.5]), .groups = "drop")
  boxplot <- ggplot(dist_df, aes(x = Group, y = BrayCurtis_Distance, fill = Group)) +
    geom_boxplot(width = 0.2, outlier.shape = 16, outlier.size = 2, alpha = 0.7, size = 1.8) +
    scale_fill_manual(values = c("Shigatse" = "#FF6347", "Lhasa" = "#32CD32", "Shannan" = "#9370DB", "Nyingchi" = "#FFD700")) +
    ggtitle(plot_title) +
    labs(x = "Group", y = "Bray-Curtis Distance (Within Group)") +
    geom_text(data = sig_letters %>% left_join(boxplot_data, by = "Group"),
              aes(x = Group, y = max_box + 0.02, label = letters),
              size = 5, fontface = "bold") +
    annotate("text", x = 2.5, y = max(dist_df$BrayCurtis_Distance) * 1.15,
             label = paste0("PERMANOVA R² = ", global_R2, "\nF = ", global_F, "\nPERMANOVA P = ", global_P),
             size = 4.5, fontface = "bold") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.ticks.x = element_line()
    )
  
  # Pairwise group comparison: calculate R² between every two geographic groups
  # Two output tables separated: 1. Within-group distance data; 2. Pairwise comparison R²
  group_pairs <- combn(unique(sample_info$Group), 2, simplify = FALSE)
  between_r2 <- map_dfr(group_pairs, function(pair) {
    sub_samples <- rownames(sample_info)[sample_info$Group %in% pair]
    sub_comm <- comm_matrix[sub_samples, ]
    sub_dist <- vegdist(sub_comm, method = "bray")
    sub_perm <- adonis2(sub_dist ~ sub_group, permutations = 999)
    tibble(
      Comparison = paste(sort(pair), collapse = " vs "),
      R2 = round(sub_perm$R2[1], 3),
      F = round(sub_perm$F[1], 2),
      P = round(sub_perm$Pr[F][1], 4)
    )
  }) %>% mutate(Habitat = habitat_name)
  between_r2 <- between_r2 %>% add_row(
    Comparison = "Mean_R2",
    R2 = round(mean(between_r2$R2), 3),
    F = NA,
    P = NA,
    Habitat = habitat_name
  )
  
  # Export integrated statistical table to Excel
  excel_name <- paste0(output_prefix, "_PERMANOVA_Kruskal_Result.xlsx")
  write_xlsx(list(
    "PERMANOVA_Result" = permanova_result,
    "Kruskal_Letters" = sig_letters,
    "Between_Group_R2" = between_r2,
    "WithinGroup_Distance_Data" = dist_df
  ), file.path(output_dir, excel_name))
  
  # Export boxplot figures
  ggsave(file.path(output_dir, paste0(output_prefix, "_WithinGroup_Distance_Boxplot.png")), boxplot, dpi = 300, width = 6, height = 6)
  ggsave(file.path(output_dir, paste0(output_prefix, "_WithinGroup_Distance_Boxplot.pdf")), boxplot, width = 6, height = 6)
  
  return(list(boxplot = boxplot, between_r2 = between_r2, boxplot_data = dist_df))
}

# Fixed sample grouping consistent with previous scripts
groups_rhizo <- list(
  Shigatse = c("RRK1", "RRK2", "RRK3", "RRK4", "RRK5"),
  Lhasa = c("RCN1", "RCN2", "RCN3", "RCN4", "RCNS"),
  Shannan = c("RSN1", "RSN2", "RSN3", "RSN4", "RSN5"),
  Nyingchi = c("RML1", "RML2", "RML3", "RML4", "RML5")
)
groups_endo <- list(
  Shigatse = c("ERK1", "ERK2", "ERK3", "ERK4", "ERK5"),
  Lhasa = c("ECN1", "ECN2", "ECN3", "ECN4", "ECN5"),
  Shannan = c("ESN1", "ESN2", "ESN3", "ESN4", "ESN5"),
  Nyingchi = c("EML1", "EML2", "EML3", "EML4", "EML5")
)

# Run analysis for four datasets: Rhizosphere Bacteria / Endophytic Bacteria / Rhizosphere Fungi / Endophytic Fungi
bac_rhizo <- process_community_analysis(
  data_path = file.path(input_dir, "Os_rb.xlsx"), groups_list = groups_rhizo,
  plot_title = "Rhizosphere Bacterial Community",
  output_prefix = "Rhizo_Bacteria", habitat_name = "Rhizosphere"
)
bac_endo <- process_community_analysis(
  data_path = file.path(input_dir, "Os_eb.xlsx"), groups_list = groups_endo,
  plot_title = "Root Endophytic Bacterial Community",
  output_prefix = "Endo_Bacteria", habitat_name = "Root Endophytic"
)
fun_rhizo <- process_community_analysis(
  data_path = file.path(input_dir, "Os_rf.xlsx"), groups_list = groups_rhizo,
  plot_title = "Rhizosphere Fungal Community",
  output_prefix = "Rhizo_Fungi", habitat_name = "Rhizosphere"
)
fun_endo <- process_community_analysis(
  data_path = file.path(input_dir, "Os_ef.xlsx"), groups_list = groups_endo,
  plot_title = "Root Endophytic Fungal Community",
  output_prefix = "Endo_Fungi", habitat_name = "Root Endophytic"
)

# Merge all within-group distance raw data into one unified Excel file
write_xlsx(list(
  Rhizosphere_Bacteria = bac_rhizo$boxplot_data,
  Root_Endophytic_Bacteria = bac_endo$boxplot_data,
  Rhizosphere_Fungi = fun_rhizo$boxplot_data,
  Root_Endophytic_Fungi = fun_endo$boxplot_data
), file.path(output_dir, "All_WithinGroup_Distance_Dataset.xlsx"))

# Merge all pairwise R² results for cross-habitat comparison
all_between <- bind_rows(
  bac_rhizo$between_r2 %>% mutate(Type = "Bacteria"),
  bac_endo$between_r2 %>% mutate(Type = "Bacteria"),
  fun_rhizo$between_r2 %>% mutate(Type = "Fungi"),
  fun_endo$between_r2 %>% mutate(Type = "Fungi")
)
write_xlsx(list("All_Pairwise_PERMANOVA_R2_Data" = all_between), file.path(output_dir, "All_Pairwise_PERMANOVA_R2_Summary.xlsx"))

cat("\nAll within-group distance analysis finished, all outputs saved to ./output folder.\n")