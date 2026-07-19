# ==============================================================================
# Script Name: 09_Top10_Genus_Boxplot_ANOVA_Tukey.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Screen top10 dominant genera by total abundance, ANOVA + Tukey HSD significance test,
#          generate horizontal log10 abundance boxplot with compact letter display for geographic groups
# Input: 0g_80hexin36.xlsx genus-level taxonomic abundance matrix stored in ./input
# Output: Cleaned full genus raw table, Top10 genus statistical summary, long format data,
#         Tukey significance letter Excel, high-res PNG/PDF genus boxplots in ./output
# Dependencies: tidyverse, readxl, writexl, multcompView, car
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified significance threshold for all statistical tests, consistent with script 03~08
global_alpha <- 0.05

# -------------------------- Install and load required analytical packages --------------------------
required_pkgs <- c("tidyverse", "readxl", "writexl", "multcompView", "car")
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

# -------------------------- Fixed geographic grouping & sample replicate definition --------------------------
# Geographic display order for consistent plotting
group_order <- c("Shigatse", "Lhasa", "Shannan", "Nyingchi")

# Rhizosphere sample replicates, consistent with all previous analysis scripts
rhizo_replicates <- list(
  Shigatse = c("RRK1", "RRK2", "RRK3", "RRK4", "RRK5"),
  Lhasa    = c("RCN1", "RCN2", "RCN3", "RCN4", "RCN5"),
  Shannan  = c("RSN1", "RSN2", "RSN3", "RSN4", "RSN5"),
  Nyingchi = c("RML1", "RML2", "RML3", "RML4", "RML5")
)

# Root endophytic sample replicates, consistent with all previous analysis scripts
endophyte_replicates <- list(
  Shigatse = c("ERK1", "ERK2", "ERK3", "ERK4", "ERK5"),
  Lhasa    = c("ECN1", "ECN2", "ECN3", "ECN4", "ECN5"),
  Shannan  = c("ESN1", "ESN2", "ESN3", "ESN4", "ESN5"),
  Nyingchi = c("EML1", "EML2", "EML3", "EML4", "EML5")
)
# Modular optimization suggestion: extract rhizo/endo replicate lists to global shared script to reduce duplicate code

# -------------------------- Community metadata mapping --------------------------
# Match excel sheet name, community label and sample replicate grouping
communities <- list(
  g_rbhexin = list(sheet = "g_rbhexin", name = "Rhizosphere_Bacteria", reps = rhizo_replicates),
  g_ebhexin = list(sheet = "g_ebhexin", name = "Endosphere_Bacteria", reps = endophyte_replicates),
  g_rfhexin = list(sheet = "g_rfhexin", name = "Rhizosphere_Fungi", reps = rhizo_replicates),
  g_efhexin = list(sheet = "g_efhexin", name = "Endosphere_Fungi", reps = endophyte_replicates)
)

# Unified fixed color palette consistent with script 03~08
group_colors <- c(
  "Shigatse" = "#FF6347",
  "Lhasa"    = "#32CD32",
  "Shannan"  = "#9370DB",
  "Nyingchi" = "#FFD700"
)
# Alternating background fill color for genus rows in boxplot
bg_colors <- c("#F0F0F0", "#FFFFFF")

# Initialize empty list container to store intermediate tables for batch export
all_full_data      <- list()
all_top10_data     <- list()
all_long_data      <- list()

# -------------------------- Batch processing loop for four microbial communities --------------------------
for (comm in names(communities)) {
  current_comm <- communities[[comm]]
  cat("\n===== Start processing community dataset:", current_comm$name, "=====\n")
  
  # Read target sheet from unified input excel file
  input_file_path <- file.path(input_dir, "0g_80hexin36.xlsx")
  if (!file.exists(input_file_path)) stop(paste("Input abundance file missing:", input_file_path))
  data <- read_excel(input_file_path, sheet = current_comm$sheet)
  
  all_cols <- colnames(data)
  group_replicates <- current_comm$reps
  all_replicates <- unlist(group_replicates)
  existing_replicates <- intersect(all_replicates, all_cols)
  
  if (length(existing_replicates) == 0) {
    stop(paste("Error: No matching sample replicate columns found in community", current_comm$name))
  }
  
  # Step1: Clean raw genus abundance table, extract genus name from g__ taxonomy tag, sort by total reads
  data_full <- data %>%
    select(1, all_of(existing_replicates)) %>%
    mutate(across(1, as.character), across(all_of(existing_replicates), as.numeric)) %>%
    rename(Taxonomy = 1) %>%
    mutate(genus = str_extract(Taxonomy, "(?<=g__)[^;]+$"), .before = 1) %>%
    filter(!is.na(genus)) %>%
    rowwise() %>%
    mutate(total_abundance = sum(c_across(all_of(existing_replicates)), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(total_abundance)) %>%
    mutate(abundance_rank = row_number(), .after = total_abundance) %>%
    rename(Top30_genus = genus)
  
  all_full_data[[current_comm$name]] <- data_full
  
  # Step2: Extract top10 genera ranked by cumulative total abundance across all samples
  top_genus <- data_full %>% slice(1:10) %>% pull(Top30_genus)
  data_top10 <- data_full %>% filter(Top30_genus %in% top_genus)
  
  # Step3: Reshape wide abundance matrix to long format, log10 transform for statistical test
  long_data <- data_top10 %>%
    pivot_longer(cols = all_of(existing_replicates), names_to = "replicate", values_to = "abundance") %>%
    mutate(
      group = case_when(
        replicate %in% group_replicates[["Shigatse"]] ~ "Shigatse",
        replicate %in% group_replicates[["Lhasa"]] ~ "Lhasa",
        replicate %in% group_replicates[["Shannan"]] ~ "Shannan",
        replicate %in% group_replicates[["Nyingchi"]] ~ "Nyingchi"
      ),
      group = factor(group, levels = group_order),
      abundance_log = log10(abundance + 1e-10) # Add tiny offset to avoid log10(0) infinite value
    )
  
  # Lock genus vertical order sorted by total abundance descending
  genus_order <- data_top10 %>% arrange(desc(total_abundance)) %>% pull(Top30_genus)
  long_data <- long_data %>% mutate(Top30_genus = factor(Top30_genus, levels = genus_order))
  all_long_data[[current_comm$name]] <- long_data
  
  # Step4: Calculate mean & standard deviation abundance per genus × geographic group
  stats_mean_sd <- long_data %>%
    group_by(Top30_genus, group) %>%
    summarise(
      Mean = mean(abundance, na.rm = TRUE),
      SD   = sd(abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = group, values_from = c(Mean,SD), names_glue = "{.value}_{group}")
  
  # Step5: One-way ANOVA + Tukey HSD post-hoc test, extract minimum adjusted p-value per genus
  stats_p <- long_data %>%
    group_by(Top30_genus) %>%
    do({
      df <- .
      fit <- aov(abundance_log ~ group, data = df)
      tuk <- TukeyHSD(fit)$group
      p_val <- min(tuk[,"p adj"], na.rm = TRUE)
      tibble(P_Shigatse=p_val, P_Lhasa=p_val, P_Shannan=p_val, P_Nyingchi=p_val)
    })
  
  # Merge descriptive statistics + significance p-values into integrated top10 genus table
  top10_final <- data_top10 %>%
    left_join(stats_mean_sd, by = "Top30_genus") %>%
    left_join(stats_p, by = "Top30_genus")
  
  # Standardize column display order for clean output Excel
  fix_col <- c(
    "Top30_genus","Taxonomy","abundance_rank","total_abundance",
    "Mean_Shigatse","SD_Shigatse","P_Shigatse",
    "Mean_Lhasa","SD_Lhasa","P_Lhasa",
    "Mean_Shannan","SD_Shannan","P_Shannan",
    "Mean_Nyingchi","SD_Nyingchi","P_Nyingchi",
    existing_replicates
  )
  top10_final <- top10_final[,fix_col]
  all_top10_data[[current_comm$name]] <- top10_final
  
  # Step6: Generate compact significance letters via Tukey HSD for boxplot annotation
  letter_data <- long_data %>%
    group_by(Top30_genus) %>%
    do({
      df <- .
      fit <- aov(abundance_log ~ group, data = df)
      tuk <- TukeyHSD(fit)
      let_res <- multcompLetters4(fit, tuk, threshold = global_alpha)
      tibble(group = names(let_res$group$Letters), letter = let_res$group$Letters)
    }) %>% ungroup()
  
  # Extract maximum log abundance value for placing significance text above boxplots
  sum_log_max <- long_data %>% group_by(Top30_genus,group) %>% summarise(max_abun_log = max(abundance_log,na.rm=T), .groups="drop")
  letter_data <- left_join(letter_data, sum_log_max, by=c("Top30_genus","group"))
  letter_data$Top30_genus <- factor(letter_data$Top30_genus, levels=genus_order)
  letter_data$group <- factor(letter_data$group, levels=group_order)
  
  # Export Tukey significance letter & max abundance coordinate table
  stats_file <- paste0(current_comm$name, "_Genus_Tukey_Significance_Letter_Stats.xlsx")
  write_xlsx(
    list(Summary_Max_Log_Abundance = sum_log_max, Significance_Letter_Map = letter_data),
    path = file.path(output_dir, stats_file)
  )
  
  # Create alternating light gray/white background strip for genus rows
  genus_bg <- data.frame(
    Top30_genus = genus_order, y_pos = seq_along(genus_order), fill = rep(bg_colors, length.out = 10)
  ) %>% mutate(ymin = y_pos - 0.5, ymax = y_pos + 0.5)
  
  # Step7: Build standardized horizontal genus abundance boxplot with significance letters
  p_box <- ggplot() +
    geom_rect(data = genus_bg, aes(xmin=-Inf, xmax=Inf, ymin=ymin, ymax=ymax, fill=fill), inherit.aes=F, alpha=0.3) +
    scale_fill_identity() +
    geom_boxplot(data=long_data, aes(x=abundance_log, y=Top30_genus, fill=group),
                 position=position_dodge(0.7), width=0.6, outlier.size=2, linewidth=0.8) +
    geom_text(data=letter_data, aes(x=max_abun_log+0.18, y=Top30_genus, label=letter, color=group),
              position=position_dodge(0.7), hjust=0, size=3.5, fontface="bold", show.legend=F) +
    scale_fill_manual(values=group_colors) + scale_color_manual(values=group_colors) +
    xlab("Relative Abundance (log10 scale)") + ylab(current_comm$name) + theme_classic() +
    theme(panel.grid.major.x = element_line(linetype="dashed", color="gray90"), legend.position="top", text=element_text(size=14))
  
  # Export high-resolution raster PNG & vector PDF figures
  fig_png <- paste0(current_comm$name, "_Top10_Genus_Boxplot.png")
  fig_pdf <- paste0(current_comm$name, "_Top10_Genus_Boxplot.pdf")
  ggsave(file.path(output_dir, fig_png), p_box, width=19, height=25, units="cm", dpi=300)
  ggsave(file.path(output_dir, fig_pdf), p_box, width=19, height=25, units="cm")
  
  cat("Community analysis finished successfully:", current_comm$name, "\n")
}

# Batch export integrated Excel tables combining all four community datasets
write_xlsx(all_full_data,  path = file.path(output_dir, "All_Communities_Raw_Cleaned_Genus_Data.xlsx"))
write_xlsx(all_top10_data, path = file.path(output_dir, "All_Communities_Top10_Genus_Mean_SD_P_Result.xlsx"))
write_xlsx(all_long_data,  path = file.path(output_dir, "All_Communities_Long_Format_Abundance_Data.xlsx"))

cat("\n===== All genus-level ANOVA & boxplot analysis completed =====\n")
cat("All figures, raw tables and Tukey statistical results saved in ./output folder\n")