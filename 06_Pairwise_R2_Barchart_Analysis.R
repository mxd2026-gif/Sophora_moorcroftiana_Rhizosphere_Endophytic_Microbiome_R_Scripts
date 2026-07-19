# ==============================================================================
# Script Name: 06_Pairwise_R2_Barchart_Analysis.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Visualize pairwise PERMANOVA R² with bar plot, Wilcoxon test cross habitat/kingdom comparison
# Input: Os_rb.xlsx, Os_eb.xlsx, Os_rf.xlsx, Os_ef.xlsx (taxon abundance matrix stored in ./input)
# Output: R² comparison bar plot PNG/PDF, integrated Wilcoxon test statistical Excel in ./output
# Dependencies: tidyverse, vegan, readxl, picante, ggsignif, ggpubr, agricolae, writexl, ape
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified parameters, easy batch modification
global_alpha <- 0.05
perm_num <- 999 # Standard permutation for PERMANOVA; set to 9999 for higher precision

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

process_community_analysis <- function(
    data_path, groups_list, plot_title, output_prefix, habitat_name
) {
  # Import raw abundance matrix, replace negative / NA values with zero
  data <- read_excel(data_path) %>% rename(Taxon = 1)
  data_clean <- data %>% mutate(across(-Taxon, ~ ifelse(. < 0 | is.na(.), 0, .)))
  
  comm_matrix_raw <- data_clean %>%
    filter(if_any(-Taxon, ~ . > 0)) %>%
    column_to_rownames("Taxon") %>%
    t()
  # Filter out samples with total zero abundance across all taxa (will be removed automatically)
  comm_matrix <- comm_matrix_raw[rowSums(comm_matrix_raw) > 0, ]
  
  # Match sample geographic grouping information
  sample_info <- tibble(
    Sample = rownames(comm_matrix),
    Group = case_when(
      Sample %in% groups_list$Shigatse ~ "Shigatse",
      Sample %in% groups_list$Lhasa ~ "Lhasa",
      Sample %in% groups_list$Shannan ~ "Shannan",
      Sample %in% groups_list$Nyingchi ~ "Nyingchi",
      TRUE ~ "Other"
    )
  ) %>% filter(Group != "Other") %>% column_to_rownames("Sample")
  
  sample_info$Group <- factor(
    sample_info$Group,
    levels = c("Shigatse","Lhasa","Shannan","Nyingchi")
  )
  
  # Bray-Curtis distance & global PERMANOVA
  dist_matrix <- vegdist(comm_matrix, method = "bray")
  adonis2(dist_matrix ~ Group, data = sample_info, permutations = perm_num)
  
  # Calculate pairwise PERMANOVA R² between every two geographic site groups
  group_pairs <- combn(unique(sample_info$Group), 2, simplify = FALSE)
  between_r2 <- map_dfr(group_pairs, function(pair) {
    sub_samples <- rownames(sample_info)[sample_info$Group %in% pair]
    sub_comm <- comm_matrix[sub_samples, ]
    sub_dist <- vegdist(sub_comm, method = "bray")
    sub_group <- sample_info[sub_samples, "Group", drop = TRUE]
    sub_perm <- adonis2(sub_dist ~ sub_group, permutations = perm_num)
    tibble(
      Comparison = paste(sort(pair), collapse = " vs "),
      R2 = round(sub_perm$R2[1], 3),
      F = round(sub_perm$F[1], 2),
      P = round(sub_perm$`Pr(>F)`[1], 4)
    )
  }) %>% mutate(Habitat = habitat_name)
  
  # Append average R² row for summary statistics
  between_r2 <- between_r2 %>%
    add_row(
      Comparison = "Mean_R2",
      R2 = round(mean(between_r2$R2), 3),
      F = NA,
      P = NA,
      Habitat = habitat_name
    )
  return(list(between_r2 = between_r2))
}

# Fixed sample grouping list, fully consistent with script 03/04/05
# Modular optimization suggestion: extract groups_rhizo / groups_endo to 00_Global_Sample_Group_List.R to reduce duplicate code
groups_rhizo <- list(
  Shigatse = c("RRK1", "RRK2", "RRK3", "RRK4", "RRK5"),
  Lhasa = c("RCN1", "RCN2", "RCN3", "RCN4", "RCN5"),
  Shannan = c("RSN1", "RSN2", "RSN3", "RSN4", "RSN5"),
  Nyingchi = c("RML1", "RML2", "RML3", "RML4", "RML5")
)
groups_endo <- list(
  Shigatse = c("ERK1", "ERK2", "ERK3", "ERK4", "ERK5"),
  Lhasa = c("ECN1", "ECN2", "ECN3", "ECN4", "ECN5"),
  Shannan = c("ESN1", "ESN2", "ESN3", "ESN4", "ESN5"),
  Nyingchi = c("EML1", "EML2", "EML3", "EML4", "EML5")
)

# Run pairwise PERMANOVA R² calculation for four datasets
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

# Merge all pairwise R² data, remove average summary row for plotting
all_data <- bind_rows(
  bac_rhizo$between_r2 %>% filter(Comparison != "Mean_R2") %>% mutate(Habitat="Rhizosphere", Type="Bacteria"),
  bac_endo$between_r2 %>% filter(Comparison != "Mean_R2") %>% mutate(Habitat="Root Endophytic", Type="Bacteria"),
  fun_rhizo$between_r2 %>% filter(Comparison != "Mean_R2") %>% mutate(Habitat="Rhizosphere", Type="Fungi"),
  fun_endo$between_r2 %>% filter(Comparison != "Mean_R2") %>% mutate(Habitat="Root Endophytic", Type="Fungi")
) %>%
  mutate(
    sig = case_when(
      P < 0.001 ~ "***", P < 0.01 ~ "**", P < global_alpha ~ "*", TRUE ~ "ns"
    ),
    Habitat = factor(Habitat, levels = c("Rhizosphere", "Root Endophytic")),
    Type = factor(Type, levels = c("Bacteria", "Fungi"))
  )

# Calculate mean, SD, SE for bar plot error bars
bar_data <- all_data %>%
  group_by(Comparison, Habitat, Type) %>%
  summarise(
    Mean_R2 = mean(R2), SD_R2 = sd(R2), SE_R2 = sd(R2)/sqrt(n()), .groups = "drop"
  ) %>%
  mutate(
    Comparison = factor(Comparison,
                        levels = c("Shigatse vs Lhasa", "Shigatse vs Shannan", "Shigatse vs Nyingchi",
                                   "Lhasa vs Shannan", "Lhasa vs Nyingchi", "Shannan vs Nyingchi"))
  )

# Within-habitat Wilcoxon test: Bacteria vs Fungi per geographic comparison
group_tests <- map_dfr(levels(bar_data$Comparison), function(comp) {
  rhizo_bac <- all_data %>% filter(Comparison == comp, Habitat == "Rhizosphere", Type == "Bacteria") %>% pull(R2)
  rhizo_fun <- all_data %>% filter(Comparison == comp, Habitat == "Rhizosphere", Type == "Fungi") %>% pull(R2)
  endo_bac <- all_data %>% filter(Comparison == comp, Habitat == "Root Endophytic", Type == "Bacteria") %>% pull(R2)
  endo_fun <- all_data %>% filter(Comparison == comp, Habitat == "Root Endophytic", Type == "Fungi") %>% pull(R2)
  
  test_rhizo <- wilcox.test(rhizo_bac, rhizo_fun, paired = FALSE, exact = FALSE)
  test_endo <- wilcox.test(endo_bac, endo_fun, paired = FALSE, exact = FALSE)
  
  tibble(
    Comparison = comp, Habitat = c("Rhizosphere", "Root Endophytic"),
    Test = "Bacteria vs Fungi", Method = "Wilcoxon independent",
    Statistic = round(c(test_rhizo$statistic, test_endo$statistic), 3),
    P_value = round(c(test_rhizo$p.value, test_endo$p.value), 4),
    Sig = case_when(
      c(test_rhizo$p.value, test_endo$p.value) < 0.001 ~ "***",
      c(test_rhizo$p.value, test_endo$p.value) < 0.01 ~ "**",
      c(test_rhizo$p.value, test_endo$p.value) < global_alpha ~ "*", TRUE ~ "ns"
    )
  )
})

# Helper function for significance marker
get_sig <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < global_alpha ~ "*", TRUE ~ "ns")
}

# Cross-group overall Wilcoxon comparison
overall_tests <- tibble(
  Comparison = c("Rhizosphere Bacteria vs Rhizosphere Fungi",
                 "Root Endophytic Bacteria vs Root Endophytic Fungi",
                 "Rhizosphere Bacteria vs Root Endophytic Bacteria",
                 "Rhizosphere Fungi vs Root Endophytic Fungi",
                 "Rhizosphere All vs Root Endophytic All"),
  Method = rep("Wilcoxon independent", 5),
  Statistic = round(c(
    wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2),paired=F)$statistic,
    wilcox.test(all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2),paired=F)$statistic,
    wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2),paired=F)$statistic,
    wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2),paired=F)$statistic,
    wilcox.test(rowMeans(cbind(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2))),
                rowMeans(cbind(all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2))),paired=F)$statistic
  ),3),
  P_value = round(c(
    wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2),paired=F)$p.value,
    wilcox.test(all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2),paired=F)$p.value,
    wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2),paired=F)$p.value,
    wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2),paired=F)$p.value,
    wilcox.test(rowMeans(cbind(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2))),
                rowMeans(cbind(all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2))),paired=F)$p.value
  ),4),
  Sig = c(
    get_sig(wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2),paired=F)$p.value),
    get_sig(wilcox.test(all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2),paired=F)$p.value),
    get_sig(wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2),paired=F)$p.value),
    get_sig(wilcox.test(all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2),paired=F)$p.value),
    get_sig(wilcox.test(rowMeans(cbind(all_data%>%filter(Habitat=="Rhizosphere",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Rhizosphere",Type=="Fungi")%>%pull(R2))),
                        rowMeans(cbind(all_data%>%filter(Habitat=="Root Endophytic",Type=="Bacteria")%>%pull(R2), all_data%>%filter(Habitat=="Root Endophytic",Type=="Fungi")%>%pull(R2))),paired=F)$p.value)
  )
)

# Rhizosphere pairwise R² bar plot
p_rhizo <- bar_data %>% filter(Habitat == "Rhizosphere") %>%
  ggplot(aes(x = Comparison, y = Mean_R2, fill = Type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = Mean_R2 - SD_R2, ymax = Mean_R2 + SD_R2), width = 0.15, linewidth = 1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("Bacteria" = "#4a6fe3", "Fungi" = "#f2a65a")) +
  labs(title = "Rhizosphere", x = "Geographic Site Comparison", y = "Mean PERMANOVA R²") +
  theme_bw() + theme(plot.title = element_text(hjust=0.5,size=14,face="bold"), axis.text.x = element_text(angle=45,hjust=1,size=10), legend.position="bottom")

# Root Endophytic pairwise R² bar plot
p_endo <- bar_data %>% filter(Habitat == "Root Endophytic") %>%
  ggplot(aes(x = Comparison, y = Mean_R2, fill = Type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = Mean_R2 - SD_R2, ymax = Mean_R2 + SD_R2), width = 0.15, linewidth = 1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("Bacteria" = "#4a6fe3", "Fungi" = "#f2a65a")) +
  labs(title = "Root Endophytic", x = "Geographic Site Comparison", y = "Mean PERMAN