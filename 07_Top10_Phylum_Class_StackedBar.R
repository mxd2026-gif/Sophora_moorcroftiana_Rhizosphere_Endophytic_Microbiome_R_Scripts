# ==============================================================================
# Script Name: 07_Top10_Phylum_Class_StackedBar.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Extract top10 dominant phyla, generate relative abundance stacked bar chart at phylum level
# Input: Os_rb.xlsx, Os_eb.xlsx, Os_rf.xlsx, Os_ef.xlsx full taxonomic abundance matrix stored in ./input
# Output: Phylum-level stacked bar PNG/PDF, top10 phylum abundance statistical CSV & integrated Excel table in ./output
# Dependencies: tidyverse, readxl, stringr, tidyr, writexl, scales, purrr
# ==============================================================================
options(scipen = 999, digits = 4)
global_alpha <- 0.05

# Auto install & load all required packages
pkg_list <- c("readxl","dplyr","ggplot2","stringr","tidyr","scales","writexl","purrr")
for(pkg in pkg_list){
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
}

# Standard unified input & output relative path, auto create folders
input_dir <- "./input"
output_dir <- "./output"
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created output directory: ./output")
}

# Standardized input file paths, consistent with script 03~06
rhizo_bac_path <- file.path(input_dir, "Os_rb.xlsx")
endo_bac_path  <- file.path(input_dir, "Os_eb.xlsx")
rhizo_fun_path <- file.path(input_dir, "Os_rf.xlsx")
endo_fun_path  <- file.path(input_dir, "Os_ef.xlsx")

# Function1: Read raw taxonomy abundance matrix, extract phylum level annotation
# Taxonomy format require: contain p__ prefix for phylum rank; unassigned taxa marked as "Unclassified"
read_raw_full <- function(path, Group) {
  data <- read_excel(path)
  # Extract phylum string after "p__" taxonomy tag
  data$Phylum <- str_extract(data[[1]], "(?<=p__)[^;]+")
  data$Phylum[is.na(data$Phylum)] <- "Unclassified"
  data$Group <- Group
  colnames(data)[1] <- "Taxonomy"
  return(data)
}

# Combine four raw community abundance datasets
raw_full_list <- list(
  rhizo_bac = read_raw_full(rhizo_bac_path, "rhizo_bac"),
  endo_bac  = read_raw_full(endo_bac_path, "endo_bac"),
  rhizo_fun = read_raw_full(rhizo_fun_path, "rhizo_fun"),
  endo_fun  = read_raw_full(endo_fun_path, "endo_fun")
)
all_raw_full <- bind_rows(raw_full_list)
sample_cols <- setdiff(colnames(all_raw_full), c("Taxonomy","Phylum","Group"))
sample_cols_11 <- head(sample_cols,11)

# Function2: Reshape wide abundance table to long format, tag habitat & kingdom grouping
read_and_process <- function(path, Group, Community, Habitat) {
  data <- read_excel(path)
  data$Phylum <- str_extract(data[[1]], "(?<=p__)[^;]+")
  data$Phylum[is.na(data$Phylum)] <- "Unclassified"
  sample_cols <- colnames(data)[-c(1,ncol(data))]
  data_long <- pivot_longer(data, cols = all_of(sample_cols), names_to = "Sample", values_to = "Raw_Abundance")
  data_long <- data_long %>%
    mutate(
      Group = Group,
      Primary_Group = Community,  # Bacteria / Fungi
      Secondary_Group = Habitat,   # Rhizosphere / Endosphere
      Sample = as.character(Sample)
    )
  return(data_long)
}

# Generate unified long-format dataset for statistical calculation
rhizo_bac <- read_and_process(rhizo_bac_path, "rhizo_bac", "Bacteria", "Rhizosphere")
endo_bac  <- read_and_process(endo_bac_path,  "endo_bac",  "Bacteria", "Endosphere")
rhizo_fun <- read_and_process(rhizo_fun_path, "rhizo_fun", "Fungi",   "Rhizosphere")
endo_fun  <- read_and_process(endo_fun_path,  "endo_fun",  "Fungi",   "Endosphere")
all_raw_data <- bind_rows(rhizo_bac, endo_bac, rhizo_fun, endo_fun)

# Match sample ID to geographic location grouping, consistent with 03/04/05 scripts
location_mapping <- tibble(
  Sample = c("RRK1","RRK2","RRK3","RRK4","RRK5","RCN1","RCN2","RCN3","RCN4","RCN5","RSN1","RSN2","RSN3","RSN4","RSN5","RML1","RML2","RML3","RML4","RML5",
             "ERK1","ERK2","ERK3","ERK4","ERK5","ECN1","ECN2","ECN3","ECN4","ECN5","ESN1","ESN2","ESN3","ESN4","ESN5","EML1","EML2","EML3","EML4","EML5"),
  Location = factor(c(rep("Shigatse",5),rep("Lhasa",5),rep("Shannan",5),rep("Nyingchi",5),
                      rep("Shigatse",5),rep("Lhasa",5),rep("Shannan",5),rep("Nyingchi",5)),
                    levels = c("Shigatse","Lhasa","Shannan","Nyingchi"))
)
all_raw_data <- all_raw_data %>% left_join(location_mapping,by="Sample") %>% drop_na(Location)

# Normalize to relative abundance (%) per single sample total reads
all_raw_data <- all_raw_data %>%
  group_by(Group, Location, Sample) %>%
  mutate(Norm_Abundance = Raw_Abundance / sum(Raw_Abundance, na.rm = TRUE) * 100) %>%
  ungroup()

# Calculate mean & SD raw/relative abundance for each phylum per community group
group_class_stats <- all_raw_data %>%
  group_by(Group, Primary_Group, Secondary_Group, Phylum) %>%
  summarise(
    Mean_Raw = mean(Raw_Abundance, na.rm = TRUE),
    SD_Raw = sd(Raw_Abundance, na.rm = TRUE),
    Mean_Norm = mean(Norm_Abundance, na.rm = TRUE),
    SD_Norm = sd(Norm_Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(Rank = rank(desc(Mean_Raw), ties.method = "min")) %>%
  ungroup()

# Screen Top10 abundant phyla sorted by average raw read count
top10_by_group <- group_class_stats %>%
  group_by(Group) %>%
  arrange(desc(Mean_Raw), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# Merge all non-Top10 low-abundance taxa into unified "Other" category for plotting & statistics
## Step 1: Sum raw abundance for non-top10 phyla
other_raw_sum <- all_raw_full %>%
  anti_join(top10_by_group[,c("Group","Phylum")], by = c("Group","Phylum")) %>%
  group_by(Group) %>%
  summarise(across(all_of(sample_cols_11), sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(Phylum = "Other")

## Step 2: Calculate summary stats for merged "Other" group
other_by_group <- all_raw_data %>%
  anti_join(top10_by_group[,c("Group","Phylum")], by = c("Group","Phylum")) %>%
  group_by(Group, Primary_Group, Secondary_Group) %>%
  summarise(
    Phylum = "Other",
    Mean_Raw = mean(Raw_Abundance, na.rm = TRUE),
    SD_Raw = sd(Raw_Abundance, na.rm = TRUE),
    Mean_Norm = mean(Norm_Abundance, na.rm = TRUE),
    SD_Norm = sd(Norm_Abundance, na.rm = TRUE),
    Rank = NA,
    .groups = "drop"
  ) %>%
  left_join(other_raw_sum, by = c("Group","Phylum"))

# Average abundance aggregated by geographic sampling site
loc_mean_df <- all_raw_data %>%
  group_by(Group, Primary_Group, Secondary_Group, Location, Phylum) %>%
  summarise(Mean_Loc = mean(Raw_Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Location, values_from = Mean_Loc)

# Full group-phylum combination table, fill missing abundance with zero
group_class_full <- expand.grid(
  Group = unique(all_raw_full$Group),
  Phylum = unique(all_raw_full$Phylum),
  stringsAsFactors = FALSE
)
raw_wide_all <- all_raw_full %>%
  select(Group, Phylum, all_of(sample_cols_11)) %>%
  full_join(group_class_full, by = c("Group","Phylum")) %>%
  mutate(across(all_of(sample_cols_11), ~replace_na(.x, 0)))

# Auxiliary function: attach original sample raw abundance columns to statistical table
add_sample_col <- function(df) {
  df %>% left_join(raw_wide_all, by = c("Group","Phylum"))
}

# Prepare 3 integrated sheets for master output Excel
## Sheet1: Top10 phyla + Other group, location averaged abundance + raw sample data
sheet1_top10_loc <- bind_rows(top10_by_group, other_by_group) %>%
  left_join(loc_mean_df, by = c("Group","Primary_Group","Secondary_Group","Phylum")) %>%
  select(Group, Primary_Group, Secondary_Group, Phylum, Rank, Mean_Raw, SD_Raw, Mean_Norm, SD_Norm, everything()) %>%
  group_by(Group) %>% arrange(desc(Mean_Raw), .by_group = TRUE) %>% ungroup() %>%
  add_sample_col()

## Sheet2: Full all-phylum statistical ranking table (including rare low-abundance taxa)
all_class_all_stats <- all_raw_data %>%
  group_by(Group, Primary_Group, Secondary_Group, Location, Phylum) %>%
  summarise(Mean_Loc = mean(Raw_Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Location, values_from = Mean_Loc) %>%
  left_join(group_class_stats, by = c("Group","Primary_Group","Secondary_Group","Phylum")) %>%
  select(Group, Primary_Group, Secondary_Group, Phylum, Rank, Mean_Raw, SD_Raw, Mean_Norm, SD_Norm, everything()) %>%
  add_sample_col()

## Sheet3: Raw sample abundance matrix only for Top10 + Other taxa
top10other_key <- bind_rows(top10_by_group, other_by_group)[,c("Group","Phylum")]
raw_sub <- all_raw_data %>% inner_join(top10other_key, by = c("Group","Phylum"))
raw_wide <- raw_sub %>% pivot_wider(names_from = Sample, values_from = Raw_Abundance)
sheet3_raw <- raw_wide %>%
  left_join(group_class_stats[,c("Group","Phylum","Rank","Mean_Raw","SD_Raw","Mean_Norm","SD_Norm")], by = c("Group","Phylum")) %>%
  select(Group, Primary_Group, Secondary_Group, Phylum, Rank, Mean_Raw, SD_Raw, Mean_Norm, SD_Norm, everything()) %>%
  group_by(Group) %>% arrange(desc(Mean_Raw), .by_group = TRUE) %>% ungroup() %>%
  add_sample_col()

# Export integrated master statistical Excel file
write_xlsx(
  list(
    "Top10_Other_Location_Average_Abundance" = sheet1_top10_loc,
    "All_Phylum_Abundance_Statistics_Ranking" = all_class_all_stats,
    "Top10_Other_Raw_Sample_Abundance" = sheet3_raw
  ),
  file.path(output_dir, "Phylum_Abundance_Statistics_Result.xlsx")
)

# Build plotting dataset for stacked relative abundance bar chart
plot_raw <- all_raw_data %>%
  inner_join(top10other_key, by = c("Group","Phylum")) %>%
  group_by(Primary_Group, Secondary_Group, Location, Phylum) %>%
  summarise(Mean_Loc = mean(Raw_Abundance, na.rm = TRUE), .groups = "drop")

# Aggregate non-Top10 taxa into single "Other" series for plotting
other_plot <- all_raw_data %>%
  anti_join(top10_by_group[,c("Group","Phylum")], by = c("Group","Phylum")) %>%
  group_by(Primary_Group, Secondary_Group, Location) %>%
  summarise(Mean_Loc = sum(Raw_Abundance, na.rm = TRUE), Phylum = "Other", .groups = "drop")

# Merge plotting data & recalculate percentage relative abundance per site
plot_all <- bind_rows(plot_raw, other_plot) %>%
  group_by(Primary_Group, Secondary_Group, Location) %>%
  mutate(Norm_Perc = Mean_Loc / sum(Mean_Loc) * 100) %>%
  ungroup()

# Split dataset into Bacteria / Fungi subplots
bac_plot <- plot_all %>% filter(Primary_Group == "Bacteria")
fun_plot <- plot_all %>% filter(Primary_Group == "Fungi")

# Export separate CSV statistical tables for bacterial & fungal Top10 phyla
bac_top10_csv <- top10_by_group %>%
  filter(Primary_Group == "Bacteria") %>%
  select(Group, Primary_Group, Secondary_Group, Phylum, Rank, Mean_Raw, SD_Raw, Mean_Norm, SD_Norm) %>%
  arrange(Group, Rank)
write.csv(bac_top10_csv, file.path(output_dir, "Bacteria_Phylum_Top10_Statistics.csv"), row.names = FALSE)

fun_top10_csv <- top10_by_group %>%
  filter(Primary_Group == "Fungi") %>%
  select(Group, Primary_Group, Secondary_Group, Phylum, Rank, Mean_Raw, SD_Raw, Mean_Norm, SD_Norm) %>%
  arrange(Group, Rank)
write.csv(fun_top10_csv, file.path(output_dir, "Fungi_Phylum_Top10_Statistics.csv"), row.names = FALSE)

# Fixed color palette for 10 dominant phyla, grey for merged "Other" group
col10 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#1B9E77")

# Set phylum factor levels sorted by average abundance (match legend order)
bac_order <- top10_by_group %>% filter(Primary_Group == "Bacteria") %>% arrange(desc(Mean_Raw)) %>% pull(Phylum) %>% unique() %>% c("Other")
bac_color <- setNames(c(col10[1:(length(bac_order)-1)], "#D3D3D3"), bac_order)
bac_plot$Phylum <- factor(bac_plot$Phylum, levels = bac_order)

fun_order <- top10_by_group %>% filter(Primary_Group == "Fungi") %>% arrange(desc(Mean_Raw)) %>% pull(Phylum) %>% unique() %>% c("Other")
fun_color <- setNames(c(col10[1:(length(fun_order)-1)], "#D3D3D3"), fun_order)
fun_plot$Phylum <- factor(fun_plot$Phylum, levels = fun_order)

# Lock geographic location display order
loc_lev <- c("Shigatse","Lhasa","Shannan","Nyingchi")
bac_plot$Location <- factor(bac_plot$Location, levels = loc_lev)
fun_plot$Location <- factor(fun_plot$Location, levels = loc_lev)

# Universal stacked bar plotting function, unified journal-style theme
plot_community <- function(data, color_map, title) {
  ggplot(data, aes(x = Location, y = Norm_Perc, fill = Phylum)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(values = color_map, name = "Phylum") +
    scale_y_continuous(labels = percent_format(scale = 1), expand = c(0,0), limits = c(0,100)) +
    facet_wrap(~Secondary_Group, nrow = 1) +
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      strip.background = element_rect(color = "black", fill = "white"),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(y = "Relative Abundance (%)")
}

# Generate final stacked bar plots
p_bac <- plot_community(bac_plot, bac_color, "Bacterial Community (Phylum Level)")
p_fun <- plot_community(fun_plot, fun_color, "Fungal Community (Phylum Level)")

# Export high-resolution raster PNG & vector PDF figures to ./output
ggsave(file.path(output_dir, "Bacteria_Phylum_Stacked_Barplot.png"), p_bac, width = 14, height = 8, dpi = 300)
ggsave(file.path(output_dir, "Bacteria_Phylum_Stacked_Barplot.pdf"), p_bac, width = 14, height = 8)
ggsave(file.path(output_dir, "Fungi_Phylum_Stacked_Barplot.png"), p_fun, width = 14, height = 8, dpi = 300)
ggsave(file.path(output_dir, "Fungi_Phylum_Stacked_Barplot.pdf"), p_fun, width = 14, height = 8)

# Render plots in RStudio viewer
print(p_bac)
print(p_fun)

cat("\n===== Phylum composition stacked bar analysis finished =====\n")
cat("All Excel statistics, CSV tables and figures saved to ./output folder\n")