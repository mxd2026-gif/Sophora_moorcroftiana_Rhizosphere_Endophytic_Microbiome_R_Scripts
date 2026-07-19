# ==============================================================================
# Script Name: 03_AlphaShannon_Diversity_Analysis.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Calculate Shannon alpha diversity, Kruskal-Wallis test & post-hoc Tukey comparison, visualize boxplot
# Input: Os_rb.xlsx, Os_eb.xlsx, Os_rf.xlsx, Os_ef.xlsx (taxon abundance matrix stored in ./input)
# Output: Shannon diversity boxplot PNG/PDF, statistical result Excel table stored in ./output
# Dependencies: tidyverse, vegan, readxl, agricolae, gridExtra
# ==============================================================================
options(scipen = 999, digits = 4)

# Install and load required packages automatically
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("agricolae", quietly = TRUE)) install.packages("agricolae")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("grid")

library(tidyverse)
library(vegan)
library(readxl)
library(agricolae)
library(writexl)
library(gridExtra)
library(grid)

# Create output directory under current working directory
input_dir <- "./input"
output_dir <- "./output"
if (!dir.exists(input_dir)) dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Core function: filter taxon with total zero abundance across all samples, calculate Shannon index and statistical analysis
process_shannon <- function(data_path, groups_list, plot_title, output_name, comm_name) {
  # Import raw table and replace invalid negative / NA values with zero
  data <- read_excel(data_path) %>% 
    rename(Taxon = 1) %>% 
    mutate(across(-Taxon, ~ ifelse(. < 0 | is.na(.), 0, .)))
  
  # Reshape to long format and assign geographic grouping information
  long_data <- data %>% 
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Abundance") %>% 
    mutate(
      Abundance = as.numeric(Abundance),
      Group = factor(
        case_when(
          Sample %in% groups_list$Shigatse ~ "Shigatse",
          Sample %in% groups_list$Lhasa ~ "Lhasa",
          Sample %in% groups_list$Shannan ~ "Shannan",
          Sample %in% groups_list$Nyingchi ~ "Nyingchi",
          TRUE ~ "Other"
        ),
        levels = c("Shigatse", "Lhasa", "Shannan", "Nyingchi", "Other")
      )
    ) %>% 
    filter(Group != "Other") %>% 
    arrange(Group, Sample)
  
  # Count valid samples per group
  sample_check <- long_data %>% distinct(Sample, Group) %>% count(Group, name = "Sample_Count")
  cat("\nSample quantity statistics for ", comm_name, "\n")
  print(sample_check)
  
  # Transform to wide matrix and remove taxa whose total abundance equals zero across all samples
  all_wide <- long_data %>% 
    pivot_wider(id_cols = Taxon, names_from = Sample, values_from = Abundance, values_fill = 0) %>% 
    column_to_rownames("Taxon")
  all_mat <- as.matrix(all_wide)
  
  # Global filter: delete columns where total abundance of the whole taxon is zero
  all_mat_filter <- all_mat[rowSums(all_mat) > 0, , drop = FALSE]
  
  # Split filtered abundance matrix by geographic group and store in separated Excel sheets
  group_names <- c("Shigatse", "Lhasa", "Shannan", "Nyingchi")
  keep_sp <- colnames(all_mat_filter)
  filter_sheet_list <- list()
  for (g in group_names){
    sub_sample <- long_data %>% filter(Group == g) %>% pull(Sample) %>% unique()
    sub_df <- as.data.frame(t(all_mat_filter[sub_sample, drop=F]))
    filter_sheet_list[[g]] <- sub_df
  }
  filtered_table_file <- file.path(output_dir, paste0(output_name, "_Global_Filtered_Species_Abundance_Table.xlsx"))
  write_xlsx(filter_sheet_list, filtered_table_file)
  cat("\nFiltered species table saved at: ", filtered_table_file, "\n")
  
  comm_matrix <- all_mat_filter
  mode(comm_matrix) <- "numeric"
  
  # Calculate Shannon-Wiener diversity index with original zero values reserved
  shannon_index <- diversity(comm_matrix, index = "shannon")
  shannon_df <- data.frame(
    Sample = names(shannon_index),
    ShannonIndex = shannon_index
  )
  long_data %>% distinct(Sample, Group) %>% filter(Sample %in% names(shannon_index)) %>% pull(Group) -> group_vec
  shannon_df <- mutate(shannon_df, Group = factor(group_vec, levels = c("Shigatse","Lhasa","Shannan","Nyingchi")))
  
  # Kruskal-Wallis overall test and pseudo R² computation
  kruskal_result <- kruskal.test(shannonIndex ~ Group, data = shannon_df)
  k <- length(unique(shannon_df$Group))
  kruskal_R2 <- max(0, (H - k + 1) / (n - k))
  
  # Post-hoc multiple comparison and generate significance letter markers
  kruskal_mult_result <- kruskal(shannon_df$shannonIndex, shannon_df$Group, console = TRUE)
  letter_map <- kruskal_mult_result$groups %>% rownames_to_column("Group") %>% rename(letter = groups)
  shannon_df <- left_join(shannon_df, letter_map, by = "Group")
  
  # Generate single boxplot for alpha diversity distribution
  box_plot <- ggplot(shannon_df, aes(x = Group, y = ShannonIndex, fill = Group)) +
    geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 2, alpha = 0.7) +
    scale_fill_manual(values = c("Shigatse" = "#FF6347", "Lhasa" = "#32CD32", "Shannan" = "#9370DB", "Nyingchi" = "#FFD700")) +
    stat_summary(fun = max, geom = "text", aes(label = letter), vjust = -0.8, size = 5, fontface = "bold") +
    labs(x = NULL, y = "Shannon-Wiener index") +
    annotate("text", x = 2.5, y = max(shannon_df$shannonIndex) * 1.15,
             label = paste0("Kruskal-Wallis R² = ", round(kruskal_R2, 3), "\nP = ", formatC(kruskal_result$p.value, format = "f", digits = 4)),
             size = 5, fontface = "bold") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"), legend.position = "none",
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
  
  # Export single boxplot figure
  fig_base <- str_remove(output_name, ".png")
  ggsave(file.path(output_dir, paste0(fig_base, "_Shannon_AlphaBoxplot.png")), box_plot, dpi = 300, width = 6, height=6)
  ggsave(file.path(output_dir, paste0(fig_base, "_Shannon_AlphaBoxplot.pdf")), box_plot, width = 6, height=6)
  
  # Summarize descriptive statistics and test results for excel output
  stats_detail <- shannon_df %>% group_by(Group, .drop = F) %>%
    summarise(
      Shannon_Mean=mean(ShannonIndex), Shannon_SD=sd(ShannonIndex), Shannon_Median=median(ShannonIndex),
      Shannon_Min=min(ShannonIndex), Shannon_Max=max(ShannonIndex), .groups = "drop"
    )
  stats_detail <- data.frame(stats_detail, H_statistic=as.numeric(kruskal_result$statistic), df="n_total", k=length(unique(shannon_df$Group)),
                             R2=round(kruskal_R2,3), P.value=round(kruskal_result$p.value,4))
  multiple_comp <- kruskal_mult_result$groups %>% rownames_to_column("Comparison")
  res_list <- list(
    "Shannon_Index_Data" = shannon_df,
    "Group_Summary_Stats" = stats_detail,
    "Multiple_Comparison" = multiple_comp
  )
  excel_file <- paste0(fig_base, "_Shannon_AlphaDiversity_Full_Result.xlsx")
  write_xlsx(res_list, file.path(output_dir, excel_file))
  
  return(list(plot=box_plot,R2=kruskal_R2,P=kruskal_result$p.value,shannon_df=shannon_df))
}

# Fixed sample grouping information consistent with previous analysis scripts
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

# Run analysis for 4 datasets separately: Rhizosphere Bacteria / Endophytic Bacteria / Rhizosphere Fungi / Endophytic Fungi
res_rhizo <- process_shannon(
  data_path = file.path(input_dir, "Os_rb.xlsx"), groups_list = groups_rhizo,
  plot_title = "Rhizosphere Bacterial Community", output_prefix = "Rhizo_Bacteria", comm_name = "Rhizosphere_Bacteria"
)
res_endo <- process_shannon(
  data_path = file.path(input_dir, "Os_eb.xlsx"), groups_list = groups_endo,
  plot_title = "Root Endophytic Bacterial Community", output_prefix = "Endo_Bacteria", comm_name = "Endophytic_Bacteria"
)
fun_rhizo <- process_shannon(
  data_path = file.path(input_dir, "Os_rf.xlsx"), groups_list = groups_rhizo,
  plot_title = "Rhizosphere Fungal Community", output_prefix = "Rhizo_Fungi", comm_name = "Rhizosphere_Fungi"
)
fun_endo <- process_shannon(
  data_path = file.path(input_dir, "Os_ef.xlsx"), groups_list = groups_endo,
  plot_title = "Root Endophytic Fungal Community", output_prefix = "Endo_Fungi", comm_name = "Root_Endophytic_Fungi"
)

# Merge two single plots into one combined horizontal figure
combined_plot <- grid.arrange(res_rhizo$plot, res_endo$plot, ncol=2, widths=c(5.8,5.3))

# Save combined composite plot
ggsave(file.path(output_dir, "Combined_Rhizo_Endo_Shannon_Boxplot.png"), combined_plot, dpi=300,width=12,height=6)
ggsave(file.path(output_dir, "Combined_Rhizo_Endo_Shannon_Boxplot.pdf"), combined_plot, width=12, height=6)

# Output final statistical information in console
cat("\n===== All analytical outputs stored in ./output =====\n")
cat("Rhizosphere Bacteria | R² = ",res_rhizo$R2,"\nP = ",res_rhizo$P,"\n")
cat("Endophytic Bacteria | R² = ",res_endo$R2,"\nP = ",res_endo$P,"\n")
cat("Rhizosphere Fungi | R² = ",fun_rhizo$R2,"\nP = ",fun_rhizo$P,"\n")
cat("Endophytic Fungi | R² = ",fun_endo$R2,"\nP = ",fun_endo$P,"\n")