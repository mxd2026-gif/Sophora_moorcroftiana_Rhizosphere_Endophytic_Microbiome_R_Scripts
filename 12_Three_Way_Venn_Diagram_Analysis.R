# ==============================================================================
# Script Name: 12_Three_Way_Venn_Diagram_Analysis.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Generate standardized 3-set Venn diagram to quantify unique & shared taxa across three geographic elevation groups
# Input: Three_Group_Taxon_Set.xlsx taxon ID membership table stored in ./input
# Output: High-res PNG raster & vector PDF Venn figures saved in unified ./output folder
# Dependencies: VennDiagram, grDevices, grid, numDeriv, readxl
# Standardization: Fully aligned with script 01~11 global palette, relative cross-platform paths, English semantic naming only
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified significance/plot configuration consistent with all prior scripts
global_alpha <- 0.05

# -------------------------- Install and load all required packages --------------------------
required_packages <- c("VennDiagram", "grDevices", "grid", "numDeriv", "readxl")
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

# Auto create missing directories with status prompt
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# -------------------------- Global unified plotting metadata (match script 03~11 palette) --------------------------
# Fixed group color palette fully consistent with all community analysis scripts
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)

# Fixed display order of three groups rendered in Venn plot
fixed_group_order <- c("Shigatse", "Nyingchi", "Lhasa")

# -------------------------- Reusable independent 3-way Venn plotting core function --------------------------
#' Generate standardized publication-ready 3-set Venn diagram from three unique taxon ID vectors
#' @param set_1 Character vector of unique taxon IDs belonging to Shigatse
#' @param set_2 Character vector of unique taxon IDs belonging to Nyingchi
#' @param set_3 Character vector of unique taxon IDs belonging to Lhasa
#' @param output_name Semantic prefix for exported figure filenames
#' @return Grid graphical object of Venn diagram for optional further editing
generate_three_way_venn <- function(set_1, set_2, set_3, output_name) {
  # Assemble named list of taxon sets for Venn engine
  venn_dataset <- list(
    Shigatse = set_1,
    Nyingchi = set_2,
    Lhasa = set_3
  )
  
  # Build Venn graphical object with journal-standard aesthetic parameters
  venn_graph <- venn.diagram(
    x = venn_dataset,
    filename = NULL,
    col = "black",
    fill = group_color_config[fixed_group_order],
    alpha = 0.7,
    label.col = "black",
    cex = 1.2,
    fontface = "bold",
    category.names = fixed_group_order,
    cat.col = "black",
    cat.cex = 1.2,
    cat.fontface = "bold",
    margin = 0.1,
    height = 2000,
    width = 2000,
    resolution = 300,
    group.order = fixed_group_order
  )
  
  # Render Venn plot inside R active graphics device
  grid.draw(venn_graph)
  
  # Export lossless vector PDF figure
  pdf(file.path(output_dir, paste0(output_name, ".pdf")), width = 8, height = 8)
  grid.draw(venn_graph)
  dev.off()
  
  # Export high-resolution raster PNG figure
  png(file.path(output_dir, paste0(output_name, ".png")), width = 2000, height = 2000, res = 300)
  grid.draw(venn_graph)
  dev.off()
  
  message(sprintf("3-way Venn diagram exported successfully: %s", output_name))
  return(venn_graph)
}

# -------------------------- Input dataset import & pre-check --------------------------
# Define standardized input file path
input_taxon_file <- file.path(input_dir, "Three_Group_Taxon_Set.xlsx")
# Stop execution if input file missing
if (!file.exists(input_taxon_file)) {
  stop(paste("Missing required input file:", input_taxon_file))
}

# Read taxon ID membership table
taxon_set_data <- read_excel(input_taxon_file)

# Extract deduplicated unique taxon ID vectors per geographic group
shigatse_taxon_set <- unique(taxon_set_data$Shigatse_Taxon_ID)
nyingchi_taxon_set <- unique(taxon_set_data$Nyingchi_Taxon_ID)
lhasa_taxon_set <- unique(taxon_set_data$Lhasa_Taxon_ID)

# Print set size summary for console log reference
message("\n===== Taxon set size summary =====")
message(sprintf("Shigatse total unique taxa: %d", length(shigatse_taxon_set)))
message(sprintf("Nyingchi total unique taxa: %d", length(nyingchi_taxon_set)))
message(sprintf("Lhasa total unique taxa: %d", length(lhasa_taxon_set)))

# -------------------------- Main visualization execution entry --------------------------
final_venn_result <- generate_three_way_venn(
  set_1 = shigatse_taxon_set,
  set_2 = nyingchi_taxon_set,
  set_3 = lhasa_taxon_set,
  output_name = "Three_Group_Taxa_Venn_Distribution"
)

message("\n===== Full 3-way Venn diagram analysis workflow completed =====")
message("Vector PDF & high-res PNG figures saved to ./output folder")