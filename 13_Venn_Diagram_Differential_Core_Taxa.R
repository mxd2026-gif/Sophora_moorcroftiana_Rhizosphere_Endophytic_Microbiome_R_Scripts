# ==============================================================================
# Script Name: 13_Venn_Diagram_Differential_Core_Taxa.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Two-set Venn diagram visualization for unique & overlapping taxa between differential abundant genera and core occupancy taxa
# Input: Species_Statistics_Result.xlsx comparative taxon count table stored in ./input
# Output: High-res PNG raster & vector PDF Venn figures saved in unified ./output folder
# Dependencies: VennDiagram, grid, readxl, numDeriv
# Standardization: Fully aligned with script 01~12 global color palette, cross-platform relative paths, English semantic naming
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold consistent with all prior analysis scripts
global_alpha <- 0.05

# -------------------------- Install & Load Required Packages --------------------------
required_packages <- c("VennDiagram", "grid", "readxl", "numDeriv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# -------------------------- Cross-platform Relative Path Configuration --------------------------
input_dir  <- "./input"
output_dir <- "./output"

# Auto create missing directories with status prompt
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# -------------------------- Global unified color palette consistent with script 03~12 --------------------------
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)

# -------------------------- Reusable Two-Set Venn Plot Core Function --------------------------
#' Generate standardized publication-ready two-set Venn diagram (Differential Taxa vs Core Taxa)
#' @param diff_only Integer count of taxa unique to differential abundant set
#' @param core_only Integer count of taxa unique to core occupancy set
#' @param intersect_num Integer count of shared overlapping taxa between two sets
#' @param group_name Descriptive label for target comparison group
#' @param diff_color Hex fill color for differential taxa circle (left)
#' @param core_color Hex fill color for core taxa circle (right)
#' @param output_filename Semantic prefix name for exported figure files
#' @return Grid graphical object of Venn diagram for optional further editing
generate_venn_plot <- function(diff_only,
                               core_only,
                               intersect_num,
                               group_name,
                               diff_color,
                               core_color,
                               output_filename) {
  
  # Construct dummy element vectors to define set sizes for Venn rendering
  diff_set_elements <- c(rep("Differential_Unique", diff_only), rep("Overlap", intersect_num))
  core_set_elements <- c(rep("Core_Unique", core_only), rep("Overlap", intersect_num))
  
  # Named two-set list: Left = Differential Taxa, Right = Core Taxa
  venn_dataset <- list(
    "Differential Taxa" = diff_set_elements,
    "Core Taxa" = core_set_elements
  )
  
  # Build Venn graphical object with journal-standard aesthetic parameters
  venn_plot_object <- venn.diagram(
    x = venn_dataset,
    filename = NULL,
    output = TRUE,
    col = "black",
    fill = c(diff_color, core_color),
    alpha = 0.5,
    label.col = "black",
    cat.col = c("black", "black"),
    cat.cex = 1.1,
    cex = 1.1,
    cat.fontface = "bold",
    fontface = "bold",
    margin = 0.2,
    rotation.degree = 0,
    category.names = names(venn_dataset),
    height = 2000,
    width = 2000,
    resolution = 300
  )
  
  # Render Venn plot inside active R graphics device
  grid.draw(venn_plot_object)
  
  # Export high-resolution raster PNG figure
  png(file.path(output_dir, paste0(output_filename, ".png")), width = 2000, height = 2000, res = 300)
  grid.draw(venn_plot_object)
  dev.off()
  
  # Export lossless vector PDF figure
  pdf(file.path(output_dir, paste0(output_filename, ".pdf")), width = 8, height = 8)
  grid.draw(venn_plot_object)
  dev.off()
  
  message(sprintf("Two-set Venn diagram exported successfully for comparison group: %s", group_name))
  return(venn_plot_object)
}

# -------------------------- Input Data Import & Integrity Pre-Check --------------------------
# Standardized full input file path
input_stats_file <- file.path(input_dir, "Species_Statistics_Result.xlsx")
# Stop execution if input table missing
if (!file.exists(input_stats_file)) {
  stop(paste("Missing required input statistics table:", input_stats_file))
}

# Read taxon comparison count table
stats_data <- read_excel(input_stats_file)

# Target pairwise geographic comparison label
target_comparison <- "Shigatse vs Nyingchi"
target_row <- stats_data[stats_data$Comparison == target_comparison, ]

# Extract integer count values for Venn set construction
unique_diff_count  <- target_row$Diff_Only
unique_core_count  <- target_row$Core_Only
overlap_taxa_count <- target_row$Intersection

# Console log summary of taxon counts for quick validation
message("\n===== Taxon set count summary for ", target_comparison, " =====")
message(sprintf("Unique differential abundant taxa: %d", unique_diff_count))
message(sprintf("Unique core occupancy taxa: %d", unique_core_count))
message(sprintf("Shared overlapping taxa: %d", overlap_taxa_count))

# -------------------------- Main Visualization Execution Entry --------------------------
# Unified palette consistent with full pipeline: Nyingchi series yellow tone system
venn_result <- generate_venn_plot(
  diff_only = unique_diff_count,
  core_only = unique_core_count,
  intersect_num = overlap_taxa_count,
  group_name = target_comparison,
  diff_color = "#FFFACD",
  core_color = group_color_config[["Nyingchi"]],
  output_filename = "Venn_Differential_Core_Taxa_Shigatse_Nyingchi"
)

message("\n===== Full two-set Venn diagram analysis workflow completed =====")
message("Vector PDF & high-res PNG Venn figures saved to ./output folder")