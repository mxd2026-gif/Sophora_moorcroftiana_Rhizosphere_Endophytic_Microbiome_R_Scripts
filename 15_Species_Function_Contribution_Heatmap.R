# ==============================================================================
# Script Name: 15_Species_Function_Contribution_Heatmap.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Species genus & KEGG pathway functional contribution data cleaning, total contribution statistics,
#          generate publication standardized non-clustered heatmap for four geographic elevation groups
# Input: 0Contribute%.xlsx genus-species & KEGG functional contribution matrix stored in unified ./input folder
# Output: Descriptive statistical Excel table, high-res PNG raster & vector PDF contribution heatmap saved in ./output
# Dependencies: readxl, writexl, dplyr, tibble, pheatmap, numDeriv
# Standardization: Fully aligned with script 01~14 unified path, global parameters, output naming & journal figure rules
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold consistent with full serial analysis pipeline
global_alpha <- 0.05

# -------------------------- Install and load all required packages --------------------------
required_packages <- c("readxl", "writexl", "dplyr", "tibble", "pheatmap", "numDeriv")
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

# -------------------------- Global unified experimental plotting configuration --------------------------
# Fixed consistent display order of four geographic elevation groups (match all scripts 03~14)
group_order <- c("Shigatse", "Lhasa", "Shannan", "Nyingchi")
# Blue-white-orange gradient heatmap color palette fixed for species & functional contribution visualization
heatmap_color_palette <- colorRampPalette(c("#1f77b4", "white", "#ff7f0e"))(100)
# Global unified community color palette consistent with full pipeline
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)

# -------------------------- Core reusable single-sheet contribution analysis function --------------------------
#' Process single excel sheet for species-genus & KEGG functional contribution statistics and heatmap matrix construction
#' @param input_file_path Full standardized file path of input contribution excel
#' @param target_sheet Character name of target analysis sheet
#' @return List object containing cleaned statistical table, heatmap numeric matrix, dynamic figure dimension parameters
process_single_contribution_sheet <- function(input_file_path, target_sheet) {
  # Read raw species genus + KEGG functional contribution dataset
  raw_data <- read_excel(input_file_path, sheet = target_sheet)
  
  # Branch split: distinguish rhizosphere / endophytic sample column mapping, unify group column names
  if (target_sheet %in% c("KEGG_Name_R80", "1vs3_R")) {
    clean_data <- raw_data %>%
      select(KEGG_Name, Contrib_Genus,
             Contrib_Group_SCNR, Contrib_Group_SSNR, Contrib_Group_SRKR, Contrib_Group_SMLR) %>%
      rename(Lhasa    = Contrib_Group_SCNR,
             Shannan  = Contrib_Group_SSNR,
             Shigatse = Contrib_Group_SRKR,
             Nyingchi = Contrib_Group_SMLR)
  } else {
    clean_data <- raw_data %>%
      select(KEGG_Name, Contrib_Genus,
             Contrib_Group_SCNE, Contrib_Group_SSNE, Contrib_Group_SRKE, Contrib_Group_SMLE) %>%
      rename(Lhasa    = Contrib_Group_SCNE,
             Shannan  = Contrib_Group_SSNE,
             Shigatse = Contrib_Group_SRKE,
             Nyingchi = Contrib_Group_SMLE)
  }
  
  # Standard data preprocessing: numeric conversion, NA fill, total contribution calculation, filter zero-contribution entries
  clean_data <- clean_data %>%
    mutate(across(all_of(group_order), ~as.numeric(.x)),
           across(all_of(group_order), ~replace_na(.x, 0))) %>%
    mutate(Total_Contribution = rowSums(across(all_of(group_order)))) %>%
    filter(Total_Contribution > 0) %>%
    arrange(desc(Total_Contribution))
  
  # Integrate standardized annotation metadata into full statistical output table
  result_table <- clean_data %>%
    mutate(Source_Sheet = target_sheet,
           Calculation_Method = "Contribution = Average species/genus functional contribution value per geographic group",
           Sample_Type = ifelse(target_sheet %in% c("KEGG_Name_R80", "1vs3_R"), "Rhizosphere", "Endophytic")) %>%
    select(Source_Sheet, KEGG_Name, Contrib_Genus,
           Shigatse, Lhasa, Shannan, Nyingchi,
           Total_Contribution, Calculation_Method, Sample_Type)
  
  # Construct heatmap matrix, combine KEGG pathway + species genus as unique row labels, lock fixed group column order
  heatmap_matrix <- clean_data %>%
    mutate(Row_Label = paste0(KEGG_Name, " | ", Contrib_Genus)) %>%
    select(Row_Label, all_of(group_order)) %>%
    column_to_rownames("Row_Label")
  heatmap_matrix <- heatmap_matrix[, group_order]
  
  row_number <- nrow(heatmap_matrix)
  
  # Dynamic auto-adjust figure size based on total species & functional entry count
  if (target_sheet == "KEGG_Name_E80") {
    png_w <- 2800
    png_h <- 10000
    pdf_w <- 10
    pdf_h <- 35
  } else {
    png_w <- 2800
    png_h <- max(3500, row_number * 110)
    pdf_w <- 10
    pdf_h <- max(15, row_number * 0.45)
  }
  
  return(list(
    stat_table = result_table,
    heat_matrix = heatmap_matrix,
    png_width = png_w,
    png_height = png_h,
    pdf_width = pdf_w,
    pdf_height = pdf_h,
    total_rows = row_number
  ))
}

# -------------------------- Independent standardized heatmap export function --------------------------
#' Generate & export non-clustered species-function contribution heatmap to high-res PNG and lossless vector PDF
#' @param heat_matrix Numeric contribution matrix for heatmap rendering
#' @param fig_width Pixel width of PNG figure
#' @param fig_height Pixel height of PNG figure
#' @param row_total Integer total number of species/functional rows for auto font scaling
#' @param sheet_name Sheet identifier for figure title & semantic output filename
#' @param save_dir Unified output folder path
export_contribution_heatmap <- function(heat_matrix, fig_width, fig_height, row_total, sheet_name, save_dir) {
  # Export high resolution raster PNG figure (300 DPI journal standard)
  png(file.path(save_dir, paste0("Species_Function_Contribution_Heatmap_", sheet_name, ".png")),
      width = fig_width, height = fig_height, res = 300)
  pheatmap(heat_matrix,
           color = heatmap_color_palette,
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           border_color = "gray90",
           fontsize_row = ifelse(row_total > 50, 6, 8),
           fontsize_col = 14,
           treeheight_row = 0,
           treeheight_col = 0,
           angle_col = 0,
           main = paste0("Species & Functional Contribution Heatmap | ", sheet_name))
  dev.off()
  
  # Export lossless vector PDF figure for manuscript submission
  pdf(file.path(save_dir, paste0("Species_Function_Contribution_Heatmap_", sheet_name, ".pdf")),
      width = fig_width / 280, height = fig_height / 280 * 10)
  pheatmap(heat_matrix,
           color = heatmap_color_palette,
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           border_color = "gray90",
           fontsize_row = ifelse(row_total > 50, 6, 8),
           fontsize_col = 14,
           treeheight_row = 0,
           treeheight_col = 0,
           angle_col = 0,
           main = paste0("Species & Functional Contribution Heatmap | ", sheet_name))
  dev.off()
  
  message(sprintf("Species & functional contribution heatmap exported successfully for sheet: %s", sheet_name))
}

# -------------------------- Input integrity pre-check & single sheet main execution --------------------------
# Standardized cross-platform full input file path
input_excel_path <- file.path(input_dir, "0Contribute%.xlsx")
# Terminate script if required input matrix file missing
if (!file.exists(input_excel_path)) {
  stop(paste("Missing required species & functional contribution input file:", input_excel_path))
}

# Target single analysis sheet (can loop all sheets for batch analysis)
target_analysis_sheet <- "KEGG_Name_R80"

# Run full species-function contribution data cleaning & matrix construction pipeline
analysis_result <- process_single_contribution_sheet(input_excel_path, target_analysis_sheet)

# Export multi-column descriptive statistical Excel table
write_xlsx(analysis_result$stat_table,
           file.path(output_dir, paste0("Species_Function_Contribution_Full_Statistical_Table_", target_analysis_sheet, ".xlsx")))

# Render and save standardized species&function contribution heatmap figure
export_contribution_heatmap(
  heat_matrix = analysis_result$heat_matrix,
  fig_width = analysis_result$png_width,
  fig_height = analysis_result$png_height,
  row_total = analysis_result$total_rows,
  sheet_name = target_analysis_sheet,
  save_dir = output_dir
)

# Completion console log
message("\n===== Single-sheet Species & Functional Contribution Heatmap Analysis Completed =====")
message(sprintf("Processed sheet: %s | Total valid species/functional entries: %d", target_analysis_sheet, analysis_result$total_rows))
message("All statistical Excel tables & contribution heatmap PNG/PDF figures saved to unified ./output folder")