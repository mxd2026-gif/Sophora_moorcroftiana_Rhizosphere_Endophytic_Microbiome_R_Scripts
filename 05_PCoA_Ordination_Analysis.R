# ==============================================================================
# Script Name: 05_PCoA_Ordination_Analysis.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: PCoA ordination based on Bray-Curtis distance, annotate PERMANOVA result on plot
# Input: Os_rb.xlsx, Os_eb.xlsx, Os_rf.xlsx, Os_ef.xlsx (taxon abundance matrix in ./input)
# Output: PCoA ordination PNG/PDF figure, coordinate & PERMANOVA result Excel in ./output
# Dependencies: tidyverse, vegan, readxl, picante, ggsignif, ggpubr, agricolae
# ==============================================================================
options(scipen = 999, digits = 4)
global_alpha <- 0.05
# Permutation times for PERMANOVA, set to 9999 for higher precision
perm_num <- 999

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
  # Import raw abundance matrix, replace negative / NA values with zero
  data <- read_excel(data_path) %>% rename(Taxon = 1)
  data_clean <- data %>% mutate(across(-Taxon, ~ ifelse(. < 0 | is.na(.), 0, .)))
  
  comm_matrix_raw <- data_clean %>%
    filter(if_any(-Taxon, ~ . > 0)) %>%
    column_to_rownames("Taxon") %>%
    t()
  # Filter samples with total zero abundance; samples with all 0 reads will be removed
  comm_matrix <- comm_matrix_raw[rowSums(comm_matrix_raw) > 0, ]
  
  # Match sample geographic grouping
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
  
  # Bray-Curtis distance matrix
  dist_matrix <- vegdist(comm_matrix, method = "bray")
  permanova_result <- adonis2(dist_matrix ~ Group, data = sample_info, permutations = perm_num)
  global_R2 <- round(permanova_result$R2[1], 3)
  global_P <- round(permanova_result$Pr[F][1], 4)
  global_F <- round(permanova_result$F[1], 2)
  
  # PCoA ordination; k=2 only extracts first 2 axes, adjust k if multi-axis visualization needed
  pcoa <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  pcoa_df <- as.data.frame(pcoa$points)
  colnames(pcoa_df) <- c("Axis1", "Axis2")
  pcoa_df$Group = sample_info$Group
  pcoa_df$Sample = rownames(sample_info)
  
  # Calculate axis variance explained percentage
  eig <- pcoa$eig
  explained <- eig / sum(eig) * 100
  
  # Draw PCoA ordination plot
  pcoa_plot <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Group, shape = Group)) +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 2) +
    scale_color_manual(values = c("Shigatse" = "#FF6347", "Lhasa" = "#32CD32", "Shannan" = "#9370DB", "Nyingchi" = "#FFD700")) +
    scale_shape_manual(values = c(16, 17, 18, 15)) +
    labs(
      x = paste0("PCoA Axis 1 (", round(explained[1], 1), "%)"),
      y = paste0("PCoA Axis 2 (", round(explained[2], 1), "%)"),
      color = "Group", shape = "Group"
    ) +
    ggtitle(plot_title) +
    annotate(
      "text",
      x = max(pcoa_df$Axis1) * 0.9,
      y = max(pcoa_df$Axis2) * 0.9,
      label = paste0("PERMANOVA R² = ", global_R2, "\nF = ", global_F, "\nPERMANOVA P = ", global_P),
      size = 4, fontface = "bold"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  
  # Pairwise PERMANOVA R² between every two geographic groups
  group_pairs <- combn(unique(sample_info$Group), 2, simplify = FALSE)
  between_r2 <- map_dfr(group_pairs, function(pair) {
    sub_samples <- rownames(sample_info)[sample_info$Group %in% pair]
    sub_comm <- comm_matrix[sub_samples, ]
    sub_dist <- vegdist(sub_comm, method = "bray")
    sub_perm <- adonis2(sub_dist ~ sub_group, permutations = perm_num)
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
  
  # Export integrated statistical Excel
  excel_name <- paste0(output_prefix, "_PCoA_PERMANOVA_Result.xlsx")
  write_xlsx(list(
    "PERMANOVA_Result" = permanova_result,
    "Between_Group_R2" = between_r2,
    "PCoA_Coordinate" = pcoa_df,
    "Axis_Variance_Explained" = data.frame(Axis = c("Axis1", "Axis2"), Variance_Percent = explained[1:2])
  ), file.path(output_dir, excel_name))
  
  # Export high-resolution vector figure
  ggsave(file.path(output_dir, paste0(output_prefix, "_PCoA_Ordination.png")), pcoa_plot, dpi = 300, width = 8, height = 6)
  ggsave(file.path(output_dir, paste0(output_prefix, "_PCoA_Ordination.pdf")), pcoa_plot, width = 8, height = 6)
  
  return(list(pcoa_plot = pcoa_plot, between_r2 = between_r2))
}

# Fixed sample grouping list, identical with script 03 & 04
# Suggestion: extract groups_rhizo / groups_endo to separate script 00_Global_Sample_Group_List.R for modular reuse
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
  plot_title = "Rhizosphere Bacterial Community", output_prefix = "Rhizo_Bacteria", habitat_name = "Rhizosphere"
)
bac_endo <- process_community_analysis(
  data_path = file.path(input_dir, "Os_eb.xlsx"), groups_list = groups_endo,
  plot_title = "Root Endophytic Bacterial Community", output_prefix = "Endo_Bacteria", habitat_name = "Root Endophytic"
)
fun_rhizo <- process_community_analysis(
  data_path = file.path(input_dir, "Os_rf.xlsx"), groups_list = groups_rhizo,
  plot_title = "Rhizosphere Fungal Community", output_prefix = "Rhizo_Fungi", habitat_name = "Rhizosphere"
)
fun_endo <- process_community_analysis(
  data_path = file.path(input_dir, "Os_ef.xlsx"), groups_list = groups_endo,
  plot_title = "Root Endophytic Fungal Community", output_prefix = "Endo_Fungi", habitat_name = "Root Endophytic"
)

# Merge all pairwise PERMANOVA R² results for unified comparison
all_between <- bind_rows(
  bac_rhizo$between_r2 %>% mutate(Type = "Bacteria"),
  bac_endo$between_r2 %>% mutate(Type = "Bacteria"),
  fun_rhizo$between_r2 %>% mutate(Type = "Fungi"),
  fun_endo$between_r2 %>% mutate(Type = "Fungi")
)
write_xlsx(list("All_Pairwise_PERMANOVA_R2_Data" = all_between), file.path(output_dir, "All_Pairwise_PERMANOVA_R2_Summary.xlsx"))

cat("\nAll PCoA ordination analysis finished, figures & statistics saved to ./output\n")