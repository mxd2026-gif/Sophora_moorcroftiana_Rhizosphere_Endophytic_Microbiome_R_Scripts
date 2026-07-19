# ==============================================================================
# Script Name: 10_KEGG_Level3_Top10_Barplot_ANOVA_Tukey.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: KEGG / Level3 functional abundance preprocessing, one-way ANOVA + Tukey HSD significance test,
#          horizontal mean abundance bar plot with error bar & compact significance letter annotation
# Input: 00kegg_80hexin.xlsx, 00level3_80hexin.xlsx functional abundance matrix stored in ./input
# Output: Merged full/top10/long-format Excel tables, per-dataset Tukey statistics Excel,
#         high-res PNG/PDF horizontal bar charts unified under ./output folder
# Dependencies: tidyverse, readxl, writexl, car
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified significance threshold consistent with script 03~09
global_alpha <- 0.05

# -------------------------- Install and load required packages --------------------------
required_pkgs <- c("tidyverse", "readxl", "writexl", "car")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# -------------------------- Unified standard relative path configuration --------------------------
# Align input/output directory with all prior scripts (01~09)
input_dir  <- "./input"
output_dir <- "./output"
# Auto create missing folders with prompt message
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# Standard input file paths
file_1 <- file.path(input_dir, "00kegg_80hexin.xlsx")
file_2 <- file.path(input_dir, "00level3_80hexin.xlsx")

# Pre-check input file existence
if (!file.exists(file_1) || !file.exists(file_2)) {
  stop(paste0("Missing input file! Required files:\n", file_1, "\n", file_2))
}

# -------------------------- Fixed geographic grouping & sample replicate definition --------------------------
# Fixed display order for plotting consistency
group_levels <- c("Shigatse", "Lhasa", "Shannan", "Nyingchi")
# Unified color palette fully consistent with scripts 03~09
group_palette <- c(
  "Shigatse" = "#FF6347",
  "Lhasa"    = "#32CD32",
  "Shannan"  = "#9370DB",
  "Nyingchi" = "#FFD700"
)

# Rhizosphere sample replicates (consistent with all community analysis scripts)
sample_set_1 <- list(
  Shigatse = c("RRK1", "RRK2", "RRK3", "RRK4", "RRK5"),
  Lhasa    = c("RCN1", "RCN2", "RCN3", "RCN4", "RCN5"),
  Shannan  = c("RSN1", "RSN2", "RSN3", "RSN4", "RSN5"),
  Nyingchi = c("RML1", "RML2", "RML3", "RML4", "RML5")
)

# Root endophytic sample replicates
sample_set_2 <- list(
  Shigatse = c("ERK1", "ERK2", "ERK3", "ERK4", "ERK5"),
  Lhasa    = c("ECN1", "ECN2", "ECN3", "ECN4", "ECN5"),
  Shannan  = c("ESN1", "ESN2", "ESN3", "ESN4", "ESN5"),
  Nyingchi = c("EML1", "EML2", "EML3", "EML4", "EML5")
)
# Modular optimization note: extract sample replicate lists to shared global script to eliminate duplicate code

# Dataset mapping: source file + sheet + habitat label + sample grouping + functional type flag
dataset_list <- list(
  ds_1 = list(src = file_1, sheet = "Rhexin", label = "Rhizosphere_KEGG", samples = sample_set_1, flag = "type1"),
  ds_2 = list(src = file_1, sheet = "Ehexin", label = "Endophytic_KEGG", samples = sample_set_2, flag = "type1"),
  ds_3 = list(src = file_2, sheet = "Rhexin", label = "Rhizosphere_Level3", samples = sample_set_1, flag = "type2"),
  ds_4 = list(src = file_2, sheet = "Ehexin", label = "Endophytic_Level3", samples = sample_set_2, flag = "type2")
)

# Alternating background row color for horizontal bar plot
row_color <- c("#F0F0F0", "#FFFFFF")

# -------------------------- Initialize storage list for merged final tables --------------------------
data_full_list  <- list()
data_top10_list <- list()
data_long_list  <- list()

# -------------------------- Batch processing loop for 4 functional datasets --------------------------
for (name in names(dataset_list)) {
  config <- dataset_list[[name]]
  cat("\n===== Start processing functional dataset:", config$label, "=====\n")
  
  raw_data <- read_excel(config$src, sheet = config$sheet)
  col_names <- colnames(raw_data)
  all_samples <- unlist(config$samples)
  match_samples <- intersect(all_samples, col_names)
  
  if (length(match_samples) == 0) {
    stop(paste("Error: No matching sample replicate columns found in dataset", config$label))
  }
  
  # Step1: Full raw table clean, calculate total summed abundance & global rank
  data_processed <- raw_data %>%
    select(1:3, all_of(match_samples)) %>%
    mutate(across(all_of(match_samples), as.numeric)) %>%
    rowwise() %>%
    mutate(total = sum(c_across(all_of(match_samples)), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(total)) %>%
    mutate(rank = row_number(), .after = total)
  
  data_full_list[[config$label]] <- data_processed
  
  # Step2: Extract top10 functional entries sorted by total cumulative abundance
  if (config$flag == "type1") {
    top_names <- data_processed %>% pull(KEGG_Name) %>% unique() %>% head(10)
    data_top10 <- data_processed %>% filter(KEGG_Name %in% top_names)
  } else {
    top_names <- data_processed %>% pull(Level3) %>% unique() %>% head(10)
    data_top10 <- data_processed %>% filter(Level3 %in% top_names)
  }
  data_top10_list[[config$label]] <- data_top10
  
  # Step3: Reshape wide abundance matrix to long format, log10 transform for ANOVA
  data_long <- data_top10 %>%
    pivot_longer(cols = all_of(match_samples), names_to = "sample", values_to = "abundance") %>%
    mutate(
      group = factor(case_when(
        sample %in% config$samples[["Shigatse"]] ~ "Shigatse",
        sample %in% config$samples[["Lhasa"]] ~ "Lhasa",
        sample %in% config$samples[["Shannan"]] ~ "Shannan",
        sample %in% config$samples[["Nyingchi"]] ~ "Nyingchi"
      ), levels = group_levels),
      log_abundance = log10(abundance + 1e-10) # Tiny offset to avoid log10(0) infinite value
    )
  
  # Lock vertical order of functional terms sorted by total abundance descending
  if (config$flag == "type1") {
    order_vec <- data_top10 %>% arrange(desc(total)) %>% pull(KEGG_Name) %>% unique()
    data_long$id <- factor(data_long$KEGG_Name, levels = order_vec)
  } else {
    order_vec <- data_top10 %>% arrange(desc(total)) %>% pull(Level3) %>% unique()
    data_long$id <- factor(data_long$Level3, levels = order_vec)
  }
  data_long_list[[config$label]] <- data_long
  
  # Step4: Calculate descriptive statistics (mean / SD raw & log abundance, max log value for label placement)
  stat_summary <- data_long %>%
    group_by(id, group) %>%
    summarise(
      mean_val = mean(abundance, na.rm = TRUE),
      sd_val = sd(abundance, na.rm = TRUE),
      mean_log = mean(log_abundance, na.rm = TRUE),
      sd_log = sd(log_abundance, na.rm = TRUE),
      max_log = max(log_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Step5: One-way ANOVA + Tukey HSD post-hoc test, generate compact significance letters
  letter_result <- data_long %>%
    group_by(id) %>%
    do({
      df <- .
      aov_model <- aov(log_abundance ~ group, data = df)
      tukey_res <- TukeyHSD(aov_model)
      p_values <- tukey_res$group[, "p adj"]
      
      group_mean <- df %>%
        group_by(group) %>%
        summarise(mean_log_val = mean(log_abundance, na.rm = TRUE)) %>%
        arrange(desc(mean_log_val))
      
      letter_tbl <- group_mean %>% mutate(letter = case_when(row_number() == 1 ~ "a", TRUE ~ "a"))
      
      # Assign distinct letter if pairwise comparison P < global significance threshold
      for (i in 2:nrow(letter_tbl)) {
        for (j in 1:(i - 1)) {
          pair_1 <- paste0(letter_tbl$group[j], "-", letter_tbl$group[i])
          pair_2 <- paste0(letter_tbl$group[i], "-", letter_tbl$group[j])
          p_cut <- ifelse(!is.na(p_values[pair_1]), p_values[pair_1], p_values[pair_2])
          if (!is.na(p_cut) && p_cut < global_alpha) {
            letter_tbl$letter[i] <- letters[i]
            break
          }
        }
      }
      letter_tbl %>% mutate(id = df$id[1])
    }) %>%
    left_join(stat_summary %>% select(id, group, max_log, mean_val, sd_val),
              by = c("id", "group"))
  
  # Export per-dataset statistical table with significance letters
  write_xlsx(
    list(Descriptive_Stats = stat_summary, Significance_Letter_Table = letter_result),
    path = file.path(output_dir, paste0("41_stats_", config$label, ".xlsx"))
  )
  
  # Alternating light gray / white background strip for horizontal bar plot rows
  bg_tbl <- data.frame(
    id = order_vec,
    pos = seq_along(order_vec),
    fill_col = rep(row_color, length.out = length(order_vec))
  ) %>% mutate(ymin = pos - 0.5, ymax = pos + 0.5)
  
  # Step6: Build standardized horizontal mean abundance bar plot with error bars & significance text
  plot_bar <- ggplot() +
    geom_rect(
      data = bg_tbl,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = fill_col),
      inherit.aes = FALSE, alpha = 0.3
    ) +
    scale_fill_identity() +
    geom_col(
      data = stat_summary,
      aes(x = mean_val, y = id, fill = group),
      position = position_dodge(width = 0.7),
      width = 0.6, color = NA
    ) +
    geom_errorbar(
      data = stat_summary,
      aes(x = mean_val, y = id, xmin = mean_val, xmax = mean_val + sd_val, group = group),
      position = position_dodge(width = 0.7),
      width = 0.2, linewidth = 0.7, orientation = "y"
    ) +
    geom_text(
      data = letter_result,
      aes(x = mean_val + sd_val + (max(stat_summary$mean_val, na.rm = TRUE) * 0.01),
          y = id, label = letter, color = group, group = group),
      position = position_dodge(width = 0.7),
      hjust = 0, size = 3.5, fontface = "bold", show.legend = FALSE
    ) +
    scale_fill_manual(values = group_palette) +
    scale_color_manual(values = group_palette) +
    labs(x = "Mean Functional Abundance ± SD", y = "") +
    theme_classic() +
    theme(
      panel.grid.major.x = element_line(linetype = "dashed", color = "gray90"),
      legend.position = "top",
      text = element_text(size = 14),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 11)
    )
  
  # Export high resolution raster PNG & vector PDF figures
  ggsave(
    filename = paste0("41_bar_", config$label, ".png"),
    plot = plot_bar,
    path = output_dir,
    width = 19, height = 25, units = "cm", dpi = 300
  )
  ggsave(
    filename = paste0("41_bar_", config$label, ".pdf"),
    plot = plot_bar,
    path = output_dir,
    width = 19, height = 25, units = "cm"
  )
  
  cat("Dataset analysis completed successfully:", config$label, "\n")
}

# Batch export integrated merged Excel tables combining all four functional datasets
write_xlsx(data_full_list,  file.path(output_dir, "41_1_Full_Functional_Abundance_Data.xlsx"))
write_xlsx(data_top10_list, file.path(output_dir, "41_2_Top10_Functional_Term_Summary.xlsx"))
write_xlsx(data_long_list,  file.path(output_dir, "41_3_Long_Format_Stat_Data.xlsx"))

cat("\n===== All KEGG & Level3 functional barplot ANOVA analysis finished =====\n")
cat("All statistical tables, horizontal bar charts saved in ./output folder\n")