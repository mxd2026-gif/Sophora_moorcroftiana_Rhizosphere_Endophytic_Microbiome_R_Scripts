# ==============================================================================
# Script Name: 16_Taxon_Gene_Abundance_Heatmap_Visualization.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Microbial taxon & functional KEGG gene abundance heatmap visualization with independent color scale legends,
#          support custom five-segment gradient palette for taxa & continuous blue-yellow-red palette for genes
# Input: 0sig+mean.xlsx (mean abundance + significance annotation) stored in unified ./input folder
# Output: Separate color scale legend PDFs, per-sheet high-res PNG + vector PDF heatmaps saved in ./output
# Dependencies: readxl, pheatmap, dplyr, ggplot2, scales, numDeriv
# Standardization: Fully aligned with script 01~15 unified path rules, global parameters, dual-format output & journal figure specs
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold consistent with full serial analysis pipeline
global_alpha <- 0.05

# -------------------------- Install and load all required packages --------------------------
required_packages <- c("readxl", "pheatmap", "dplyr", "ggplot2", "scales", "numDeriv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# -------------------------- Cross-platform relative path configuration --------------------------
input_dir  <- "./input"
output_dir <- "./output"

# Auto create missing directories with running status prompt
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# -------------------------- Global unified experimental configuration --------------------------
# Fixed geographic group column names for mean abundance (fully consistent with script 09~11 statistical outputs)
group_order <- c("Mean_Group_Shigatse", "Mean_Group_Lhasa", "Mean_Group_Shannan", "Mean_Group_Nyingchi")
# Significance letter annotation columns (Tukey HSD test results from prior ANOVA analysis)
sig_cols    <- c("Sig_Shigatse", "Sig_Lhasa", "Sig_Shannan", "Sig_Nyingchi")
# Standardized display labels for four elevation groups (consistent with full pipeline 01~15)
group_labels <- c("Shigatse", "Lhasa", "Shannan", "Nyingchi")
# Global unified community color palette (consistent structure with all prior scripts)
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)

# Target sheet name lists for gene functional datasets and taxon genus datasets
gene_sheet_list <- c("KEGG_Name_E80", "KEGG_Name_R80", "1vs3_R")
taxon_sheet_list <- c("Genus_eb", "Genus_rb", "1vs3_rb")

# -------------------------- Global abundance range calculation for unified color scale --------------------------
#' Calculate global log10 transformed abundance range across multiple sheets for consistent color scaling
#' @param input_path Full standardized path of input abundance excel file
#' @param sheet_vector Character vector of target sheet names to include in range calculation
#' @return Numeric vector containing all log10(abundance+1) values across all target sheets
get_global_abundance_range <- function(input_path, sheet_vector) {
  all_log_values <- c()
  for (sheet_name in sheet_vector) {
    raw_df <- read_excel(input_path, sheet = sheet_name)
    abundance_matrix <- as.matrix(raw_df[, group_order])
    abundance_matrix[is.na(abundance_matrix)] <- 0
    log_matrix <- log10(abundance_matrix + 1)
    all_log_values <- c(all_log_values, as.vector(log_matrix))
  }
  return(all_log_values)
}

# -------------------------- Custom color palette definition --------------------------
# Custom five-segment gradient color palette for taxon datasets (top journal publication standard)
taxon_color_segment1 <- colorRampPalette(c("#08306b", "#2171b5", "#6baed6", "#fd8d3c", "#e6550d"))(100)
taxon_color_segment2 <- colorRampPalette(c("#ffffcc", "#ffeda0", "#fed976", "#feb24c"))(100)
taxon_color_segment3 <- colorRampPalette(c("#e5f5e0", "#a1d99b", "#74c476", "#41ab5d", "#238b45"))(100)
taxon_color_segment4 <- colorRampPalette(c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#756bb1", "#6a51a3"))(50)
taxon_color_segment5 <- colorRampPalette(c("#fee0d2", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15"))(20)
taxon_full_color_palette <- c(taxon_color_segment1, taxon_color_segment2, taxon_color_segment3, taxon_color_segment4, taxon_color_segment5)

# Break points for five-segment taxon color scale (non-linear segmentation for low abundance differentiation)
break_1 <- seq(0, 0.0001, length.out = 101)
break_2 <- seq(0.0001, 0.001, length.out = 101)[-1]
break_3 <- seq(0.001, 0.01, length.out = 101)[-1]
break_4 <- seq(0.01, 0.05, length.out = 51)[-1]
break_5 <- seq(0.2, 0.26, length.out = 21)[-1]
taxon_color_breaks <- c(break_1, break_2, break_3, break_4, break_5)

# Continuous blue-yellow-red gradient palette for KEGG functional gene datasets
gene_color_palette <- colorRampPalette(
  c("#081d58","#253494","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4",
    "#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#b10026")
)(200)

# -------------------------- Independent color scale legend generation functions --------------------------
#' Generate & export five separate segment color scale legends for taxon heatmap (independent layout for manuscript)
#' @return No return value, PDF legend files saved directly to output directory
generate_taxon_color_legend <- function() {
  # Segment 1: Deep blue - orange low abundance range
  seg1_data <- data.frame(x = 1, y = seq(0, 0.0001, length.out = 100))
  plot_seg1 <- ggplot(seg1_data, aes(x, y)) +
    geom_tile(aes(fill = y), color = NA) +
    scale_fill_gradientn(colours = taxon_color_segment1, limits = c(0, 0.0001)) +
    labs(fill = "log10(Abundance+1)") +
    theme_void() +
    theme(legend.position = "right", legend.key.height = unit(10, "cm"))
  
  # Segment 2: Light yellow - orange low-mid abundance range
  seg2_data <- data.frame(x = 1, y = seq(0.0001, 0.001, length.out = 100))
  plot_seg2 <- ggplot(seg2_data, aes(x, y)) +
    geom_tile(aes(fill = y), color = NA) +
    scale_fill_gradientn(colours = taxon_color_segment2, limits = c(0.0001, 0.001)) +
    labs(fill = "") +
    theme_void() +
    theme(legend.position = "right", legend.key.height = unit(10, "cm"))
  
  # Segment 3: Green mid abundance range
  seg3_data <- data.frame(x = 1, y = seq(0.001, 0.01, length.out = 100))
  plot_seg3 <- ggplot(seg3_data, aes(x, y)) +
    geom_tile(aes(fill = y), color = NA) +
    scale_fill_gradientn(colours = taxon_color_segment3, limits = c(0.001, 0.01)) +
    labs(fill = "") +
    theme_void() +
    theme(legend.position = "right", legend.key.height = unit(10, "cm"))
  
  # Segment 4: Purple mid-high abundance range
  seg4_data <- data.frame(x = 1, y = seq(0.01, 0.05, length.out = 50))
  plot_seg4 <- ggplot(seg4_data, aes(x, y)) +
    geom_tile(aes(fill = y), color = NA) +
    scale_fill_gradientn(colours = taxon_color_segment4, limits = c(0.01, 0.05)) +
    labs(fill = "") +
    theme_void() +
    theme(legend.position = "right", legend.key.height = unit(5, "cm"))
  
  # Segment 5: Red high abundance range
  seg5_data <- data.frame(x = 1, y = seq(0.2, 0.26, length.out = 20))
  plot_seg5 <- ggplot(seg5_data, aes(x, y)) +
    geom_tile(aes(fill = y), color = NA) +
    scale_fill_gradientn(colours = taxon_color_segment5, limits = c(0.2, 0.26)) +
    labs(fill = "") +
    theme_void() +
    theme(legend.position = "right", legend.key.height = unit(3, "cm"))
  
  # Export all five segment legend vector PDF files
  ggsave(file.path(output_dir, "Taxon_Color_Scale_Segment_1.pdf"), plot_seg1, width = 2.5, height = 12, dpi = 300)
  ggsave(file.path(output_dir, "Taxon_Color_Scale_Segment_2.pdf"), plot_seg2, width = 2.5, height = 12, dpi = 300)
  ggsave(file.path(output_dir, "Taxon_Color_Scale_Segment_3.pdf"), plot_seg3, width = 2.5, height = 12, dpi = 300)
  ggsave(file.path(output_dir, "Taxon_Color_Scale_Segment_4.pdf"), plot_seg4, width = 2.5, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "Taxon_Color_Scale_Segment_5.pdf"), plot_seg5, width = 2.5, height = 4, dpi = 300)
  message("Taxon five-segment color scale legends exported successfully")
}

#' Generate & export continuous unified color scale legend for KEGG gene heatmap
#' @param g_min Minimum value of global log abundance range
#' @param g_max Maximum value of global log abundance range
#' @return No return value, PDF legend file saved directly to output directory
generate_gene_color_legend <- function(g_min, g_max) {
  legend_data <- data.frame(x = 1, y = seq(g_min, g_max, length.out = 200))
  legend_plot <- ggplot(legend_data, aes(x, y)) +
    geom_tile(aes(fill = y), color = NA) +
    scale_fill_gradientn(colours = gene_color_palette, limits = c(g_min, g_max)) +
    labs(fill = "log10(Abundance+1)") +
    theme_void() +
    theme(legend.position = "right", legend.key.height = unit(16, "cm"))
  
  ggsave(file.path(output_dir, "Gene_Continuous_Color_Scale.pdf"), legend_plot, width = 2.5, height = 18, dpi = 300)
  message("Gene continuous color scale legend exported successfully")
}

# -------------------------- Core reusable single-sheet heatmap rendering function --------------------------
#' Generate non-clustered abundance heatmap with significance annotation, export dual PNG+PDF format
#' @param input_path Full standardized path of input abundance excel file
#' @param target_sheet Character name of target analysis sheet
#' @param color_break Numeric vector of color scale break points
#' @param color_palette Custom color vector for heatmap rendering
#' @param output_file_prefix Semantic prefix for output heatmap filenames
#' @return No return value, heatmap files saved directly to output directory
render_single_heatmap <- function(input_path, target_sheet, color_break, color_palette, output_file_prefix) {
  raw_data <- read_excel(input_path, sheet = target_sheet)
  id_column <- colnames(raw_data)[1]
  # Deduplicate rows by first ID column to avoid duplicate row name error
  clean_data <- distinct(raw_data, !!sym(id_column), .keep_all = TRUE)
  
  # Construct log10 transformed abundance matrix
  abundance_matrix <- as.matrix(clean_data[, group_order])
  rownames(abundance_matrix) <- clean_data[[id_column]]
  abundance_matrix[is.na(abundance_matrix)] <- 0
  log_abundance_matrix <- log10(abundance_matrix + 1)
  
  # Construct significance letter annotation matrix, replace NA with empty string
  sig_annotation_matrix <- as.matrix(clean_data[, sig_cols])
  sig_annotation_matrix[is.na(sig_annotation_matrix)] <- ""
  
  # Export vector PDF heatmap for manuscript submission
  pheatmap(
    mat = log_abundance_matrix,
    color = color_palette,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_col = group_labels,
    fontsize = 8.5,
    display_numbers = sig_annotation_matrix,
    number_color = "black",
    breaks = color_break,
    legend = FALSE,
    filename = file.path(output_dir, paste0(output_file_prefix, target_sheet, ".pdf")),
    width = 8,
    height = max(3, nrow(log_abundance_matrix) / 3)
  )
  
  # Export high-resolution 300 DPI PNG heatmap for preview / presentation
  pheatmap(
    mat = log_abundance_matrix,
    color = color_palette,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_col = group_labels,
    fontsize = 8.5,
    display_numbers = sig_annotation_matrix,
    number_color = "black",
    breaks = color_break,
    legend = FALSE,
    filename = file.path(output_dir, paste0(output_file_prefix, target_sheet, ".png")),
    width = 8,
    height = max(3, nrow(log_abundance_matrix) / 3),
    res = 300
  )
  
  message(sprintf("Dual-format heatmap exported successfully for sheet: %s", target_sheet))
}

# -------------------------- Input integrity pre-check & main workflow execution --------------------------
# Standardized cross-platform full input file path
input_file_path <- file.path(input_dir, "0sig+mean.xlsx")
# Terminate script if required input abundance file missing
if (!file.exists(input_file_path)) {
  stop(paste("Missing required abundance + significance input file:", input_file_path))
}

# Calculate global log abundance range for unified gene color scale across all gene sheets
gene_log_values <- get_global_abundance_range(input_file_path, gene_sheet_list)
gene_min <- min(gene_log_values)
gene_max <- max(gene_log_values)
gene_color_breaks <- seq(gene_min, gene_max, length.out = 201)

# Generate independent color scale legend files
generate_taxon_color_legend()
generate_gene_color_legend(gene_min, gene_max)

# Example single-sheet heatmap rendering (wrap with loop for batch analysis of all sheets)
# Taxon genus dataset example
render_single_heatmap(
  input_path = input_file_path,
  target_sheet = "Genus_eb",
  color_break = taxon_color_breaks,
  color_palette = taxon_full_color_palette,
  output_file_prefix = "Taxon_Abundance_Heatmap_"
)

# KEGG functional gene dataset example
render_single_heatmap(
  input_path = input_file_path,
  target_sheet = "KEGG_Name_E80",
  color_break = gene_color_breaks,
  color_palette = gene_color_palette,
  output_file_prefix = "Gene_Abundance_Heatmap_"
)

# Completion console log
message("\n===== All taxon & gene abundance heatmap visualization tasks completed =====")
message("Custom five-segment gradient palette applied to taxon datasets, continuous palette applied to gene datasets")
message("All color scale legends & dual-format heatmap files saved to unified ./output folder")