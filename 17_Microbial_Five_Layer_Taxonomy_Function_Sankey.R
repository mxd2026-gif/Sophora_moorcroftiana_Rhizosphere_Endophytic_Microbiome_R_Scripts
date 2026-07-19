# ==============================================================================
# Script Name: 17_Microbial_Five_Layer_Taxonomy_Function_Sankey.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Five-layer alluvial/Sankey diagram visualization linking microbial taxonomic classification (Phylum-Class-Genus)
#          and KEGG functional hierarchy (KEGG_Name-Level3), unified color mapping for journal publication
# Input: 00sangji.xlsx taxonomic + KEGG full annotation table stored in unified ./input folder
# Output: High-res 300 DPI PNG raster & lossless vector PDF Sankey figures saved in ./output
# Dependencies: readxl, ggplot2, ggalluvial, dplyr, numDeriv
# Standardization: Fully aligned with script 01~16 unified global parameters, cross-platform paths, dual-format output rules
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold consistent with full serial analysis pipeline
global_alpha <- 0.05

# -------------------------- Install and load all required packages --------------------------
required_packages <- c("readxl", "ggplot2", "ggalluvial", "dplyr", "numDeriv")
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

# Auto create missing standard folders with running status prompt
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# -------------------------- Global unified experimental configuration --------------------------
# Fixed five sequential axis order for 5-layer Sankey: Phylum > Class > Genus > KEGG_Name > Level3
sankey_column_order <- c("Phylum", "Class", "Genus", "KEGG_Name", "Level3")
# Global unified community geographic color palette (consistent with script 03~16 full pipeline)
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)

# Custom distinct color palette for top phylum taxonomic strata
phylum_color_palette <- c(
  "#1f77b4",
  "#ff7f0e",
  "#2ca02c",
  "#d62728",
  "#9467bd",
  "#8c564b",
  "#e377c2",
  "#bcbd22",
  "#17becf"
)

# Fixed gray color mapping for unknown/unclassified KEGG Level3 functional entries
special_level3_color_map <- c(
  "99997 Function unknown" = "#d3d3d3",
  "99994 Others"           = "#a9a9a9"
)

# Low-saturation neutral color series for regular annotated Level3 functional pathways
base_low_saturation_colors <- c(
  "#4a6c9b", "#6b8baf", "#8da0c2", "#b0b9d4",
  "#c98a8a", "#d9a3a3", "#e4b8b8", "#e9c5c5",
  "#c9c28a", "#d9d3a3", "#e4e0b8", "#e9e6c5",
  "#c9a36b", "#d9b48c", "#e4c5a3", "#e9d1b8"
)

# -------------------------- Core reusable single-sheet Sankey dataset processing function --------------------------
#' Process single excel annotation sheet, generate full color mapping and complete five-layer Sankey ggplot object
#' @param input_file_path Full standardized cross-platform path of input annotation excel file
#' @param target_sheet Character name of target analysis sheet
#' @return List object containing cleaned frequency dataset, fully assembled Sankey plot ggplot object
process_single_sankey_dataset <- function(input_file_path, target_sheet) {
  # Read raw taxonomic + functional full annotation table
  raw_data <- read_excel(input_file_path, sheet = target_sheet)
  
  # Standard data preprocessing: lock fixed five-layer columns, remove NA, calculate path frequency count
  clean_dataset <- raw_data %>%
    select(all_of(sankey_column_order)) %>%
    drop_na() %>%
    group_by(across(all_of(sankey_column_order))) %>%
    summarise(frequency = n(), .groups = "drop")
  
  # Build color lookup table for all KEGG Level3 functional categories
  all_level3_items <- unique(clean_dataset$Level3)
  regular_level3_items <- setdiff(all_level3_items, names(special_level3_color_map))
  regular_level3_color_vector <- rep(base_low_saturation_colors, length.out = length(regular_level3_items))
  names(regular_level3_color_vector) <- regular_level3_items
  level3_full_color_map <- c(regular_level3_color_vector, special_level3_color_map)
  
  # Build phylum-matched color lookup table for taxonomic strata
  all_phylum_items <- unique(clean_dataset$Phylum)
  phylum_color_vector <- rep(phylum_color_palette, length.out = length(all_phylum_items))
  names(phylum_color_vector) <- all_phylum_items
  
  # Initialize base alluvial/Sankey canvas framework
  sankey_plot <- ggplot(
    clean_dataset,
    aes(
      y = frequency,
      axis1 = Phylum,
      axis2 = Class,
      axis3 = Genus,
      axis4 = KEGG_Name,
      axis5 = Level3
    )
  ) +
    geom_alluvium(aes(fill = Level3), width = 1/12, alpha = 0.8) +
    scale_fill_manual(values = level3_full_color_map, drop = FALSE) +
    geom_stratum(fill = NA, color = "black", linewidth = 0.5, width = 1/8) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 2.2,
      color = "black",
      fontface = "bold",
      nudge_x = 0.12
    ) +
    scale_x_discrete(limits = sankey_column_order, expand = c(0.15, 0.15)) +
    scale_y_continuous(breaks = NULL) +
    labs(
      title = paste0("Five-Layer Taxonomy-Function Sankey Diagram | ", target_sheet),
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(size = 9, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  stratum_fixed_width <- 1 / 8
  
  # Fill solid color blocks for Phylum (1st axis) strata
  phylum_stat_data <- clean_dataset %>%
    group_by(Phylum) %>%
    summarise(total_frequency = sum(frequency), .groups = "drop") %>%
    arrange(desc(total_frequency)) %>%
    mutate(
      y_max = cumsum(total_frequency),
      y_min = lag(y_max, default = 0),
      x_min = 1 - stratum_fixed_width / 2,
      x_max = 1 + stratum_fixed_width / 2
    )
  
  for (row_idx in seq_len(nrow(phylum_stat_data))) {
    sankey_plot <- sankey_plot + annotate(
      "rect",
      xmin = phylum_stat_data$x_min[row_idx],
      xmax = phylum_stat_data$x_max[row_idx],
      ymin = phylum_stat_data$y_min[row_idx],
      ymax = phylum_stat_data$y_max[row_idx],
      fill = phylum_color_vector[phylum_stat_data$Phylum[row_idx]],
      color = NA
    )
  }
  
  # Fill solid color blocks for Class (2nd axis) strata (inherit phylum color)
  class_stat_data <- clean_dataset %>%
    group_by(Class, Phylum) %>%
    summarise(total_frequency = sum(frequency), .groups = "drop") %>%
    arrange(desc(total_frequency)) %>%
    mutate(
      y_max = cumsum(total_frequency),
      y_min = lag(y_max, default = 0),
      x_min = 2 - stratum_fixed_width / 2,
      x_max = 2 + stratum_fixed_width / 2
    )
  
  for (row_idx in seq_len(nrow(class_stat_data))) {
    sankey_plot <- sankey_plot + annotate(
      "rect",
      xmin = class_stat_data$x_min[row_idx],
      xmax = class_stat_data$x_max[row_idx],
      ymin = class_stat_data$y_min[row_idx],
      ymax = class_stat_data$y_max[row_idx],
      fill = phylum_color_vector[class_stat_data$Phylum[row_idx]],
      color = NA
    )
  }
  
  # Fill solid color blocks for Genus (3rd axis) strata (inherit phylum color)
  genus_stat_data <- clean_dataset %>%
    group_by(Genus, Phylum) %>%
    summarise(total_frequency = sum(frequency), .groups = "drop") %>%
    arrange(desc(total_frequency)) %>%
    mutate(
      y_max = cumsum(total_frequency),
      y_min = lag(y_max, default = 0),
      x_min = 3 - stratum_fixed_width / 2,
      x_max = 3 + stratum_fixed_width / 2
    )
  
  for (row_idx in seq_len(nrow(genus_stat_data))) {
    sankey_plot <- sankey_plot + annotate(
      "rect",
      xmin = genus_stat_data$x_min[row_idx],
      xmax = genus_stat_data$x_max[row_idx],
      ymin = genus_stat_data$y_min[row_idx],
      ymax = genus_stat_data$y_max[row_idx],
      fill = phylum_color_vector[genus_stat_data$Phylum[row_idx]],
      color = NA
    )
  }
  
  # Fill blank white background for intermediate KEGG_Name (4th axis) strata
  kegg_stat_data <- clean_dataset %>%
    group_by(KEGG_Name) %>%
    summarise(total_frequency = sum(frequency), .groups = "drop") %>%
    arrange(desc(total_frequency)) %>%
    mutate(
      y_max = cumsum(total_frequency),
      y_min = lag(y_max, default = 0),
      x_min = 4 - stratum_fixed_width / 2,
      x_max = 4 + stratum_fixed_width / 2
    )
  
  for (row_idx in seq_len(nrow(kegg_stat_data))) {
    sankey_plot <- sankey_plot + annotate(
      "rect",
      xmin = kegg_stat_data$x_min[row_idx],
      xmax = kegg_stat_data$x_max[row_idx],
      ymin = kegg_stat_data$y_min[row_idx],
      ymax = kegg_stat_data$y_max[row_idx],
      fill = "white",
      color = NA
    )
  }
  
  # Fill custom color blocks for Level3 functional pathway (5th axis) strata
  level3_stat_data <- clean_dataset %>%
    group_by(Level3) %>%
    summarise(total_frequency = sum(frequency), .groups = "drop") %>%
    arrange(desc(total_frequency)) %>%
    mutate(
      y_max = cumsum(total_frequency),
      y_min = lag(y_max, default = 0),
      x_min = 5 - stratum_fixed_width / 2,
      x_max = 5 + stratum_fixed_width / 2
    )
  
  for (row_idx in seq_len(nrow(level3_stat_data))) {
    sankey_plot <- sankey_plot + annotate(
      "rect",
      xmin = level3_stat_data$x_min[row_idx],
      xmax = level3_stat_data$x_max[row_idx],
      ymin = level3_stat_data$y_min[row_idx],
      ymax = level3_stat_data$y_max[row_idx],
      fill = level3_full_color_map[level3_stat_data$Level3[row_idx]],
      color = NA
    )
  }
  
  return(list(
    plot_object = sankey_plot,
    raw_clean_data = clean_dataset
  ))
}

# -------------------------- Independent standardized dual-format figure export function --------------------------
#' Export Sankey diagram to high-res 300 DPI PNG raster & lossless vector PDF (journal standard dual output)
#' @param plot_obj Complete assembled ggplot Sankey figure object
#' @param save_dir Unified standardized output directory path
#' @param sheet_id Target sheet identifier for semantic descriptive output filenames
export_sankey_figure <- function(plot_obj, save_dir, sheet_id) {
  # Export high-resolution 300 DPI PNG for preview, presentation & supplementary materials
  png_output_path <- file.path(save_dir, paste0("Taxonomy_Function_FiveLayer_Sankey_", sheet_id, ".png"))
  ggsave(
    filename = png_output_path,
    plot = plot_obj,
    width = 16,
    height = 10,
    dpi = 300
  )
  
  # Export lossless vector PDF for formal manuscript submission
  pdf_output_path <- file.path(save_dir, paste0("Taxonomy_Function_FiveLayer_Sankey_", sheet_id, ".pdf"))
  ggsave(
    filename = pdf_output_path,
    plot = plot_obj,
    width = 16,
    height = 10
  )
  
  message(sprintf("Five-layer taxonomy-function Sankey dual-format figures exported successfully for sheet: %s", sheet_id))
}

# -------------------------- Input integrity pre-check & single sheet main execution pipeline --------------------------
# Standardized cross-platform full input file path
input_excel_file <- file.path(input_dir, "00sangji.xlsx")
# Terminate script with clear prompt if required annotation input file missing
if (!file.exists(input_excel_file)) {
  stop(paste("Missing required taxonomic & functional annotation input file:", input_excel_file))
}

# Automatically detect first sheet as target analysis sheet (wrap with for loop for batch multi-sheet analysis)
target_analysis_sheet <- excel_sheets(input_excel_file)[1]

# Run full five-layer Sankey data cleaning, color mapping & plot assembly workflow
analysis_result <- process_single_sankey_dataset(input_excel_file, target_analysis_sheet)

# Export standardized dual-format Sankey visualization figures
export_sankey_figure(
  plot_obj = analysis_result$plot_object,
  save_dir = output_dir,
  sheet_id = target_analysis_sheet
)

# Standardized completion running log with separator
message("\n===== Single-Sheet Five-Layer Taxonomy-Function Sankey Visualization Analysis Completed =====")
message(sprintf("Processed annotation sheet: %s", target_analysis_sheet))
message("All high-res PNG raster & vector PDF Sankey figures saved to unified ./output folder")