# ==============================================================================
# Script Name: 19_Microbial_Community_Environment_RDA_Analysis.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Redundancy Analysis (RDA) linking genus-level microbial community composition and VIF-filtered environmental factors,
#          Hellinger transformation + permutation ANOVA significance test, generate ordination biplot with grouping ellipses,
#          export full coordinate statistics and axis interpretation tables for manuscript
# Input: Rhizosphere_Bacteria_Genus.xlsx genus abundance matrix in ./input;
#        Filtered_Environmental_Dataset_VIF_Removed.rds intermediate file output from script 18 VIF screening
# Output: Dual-format 300 DPI PNG + vector PDF RDA ordination plots, multi-sheet Excel RDA statistical result tables saved in ./output
# Dependencies: vegan, dplyr, readxl, ggplot2, stringr, tibble, writexl, tidyr, numDeriv
# Standardization: Fully aligned with script 01~18 unified global parameters, cross-platform path rules, journal dual-format output specs
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold consistent with full serial analysis pipeline
global_alpha <- 0.05

# --------------------------- Install & Load All Dependent Packages ---------------------------
required_packages <- c("vegan", "dplyr", "readxl", "ggplot2",
                       "stringr", "tibble", "writexl", "tidyr", "numDeriv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Dependency package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# --------------------------- Cross-platform Relative Path Configuration ---------------------------
input_dir  <- "./input"
output_dir <- "./output"

# Auto create standard project folders with running prompt
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# Global unified fixed experimental parameters
# Intermediate environmental RDS output from No.18 VIF filtering script
env_rds_path <- file.path(output_dir, "Filtered_Environmental_Dataset_VIF_Removed.rds")
permutation_times <- 999
# Global four geographic sampling group color palette (consistent with script 03~18 full pipeline)
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)
group_factor_level <- names(group_color_config)

# Pre-check intermediate environmental file from VIF screening
if (!file.exists(env_rds_path)) {
  stop(paste("Missing mandatory filtered environmental intermediate file:", env_rds_path,
             "\nPlease run script 18 VIF multicollinearity filtering first to generate this file."))
}
# Load full VIF analysis result object
load(env_rds_path)
# Extract core standardized environmental data for RDA
rhizo_analysis_res <- analysis_result

# --------------------------- Function 1: Microbial OTU/Genus Abundance Matrix Preprocessing ---------------------------
#' Extract genus-level abundance matrix, remove all-zero taxa, unify sample ID order matching environmental metadata
#' @param otu_excel_path Full standardized path of genus abundance excel table
#' @param env_sample_id Sample ID vector extracted from preprocessed environmental metadata
#' @return List containing filtered genus abundance matrix and ordered matched sample ID vector
preprocess_microbial_otu <- function(otu_excel_path, env_sample_id) {
  raw_otu <- read_excel(otu_excel_path)
  colnames(raw_otu)[1] <- "Taxonomy"
  
  # Extract genus annotation string from full taxonomic lineage
  clean_otu <- raw_otu %>%
    mutate(Genus_label = stringr::str_extract(Taxonomy, "(?<=g__)[^;]+"),
           Genus_label = ifelse(is.na(Genus_label), "Unclassified", Genus_label)) %>%
    filter(!is.na(Genus_label))
  
  # Resolve duplicate genus label names
  clean_otu$Genus_label <- make.unique(clean_otu$Genus_label, sep = "_")
  
  # Isolate sample abundance numeric matrix
  abundance_matrix <- clean_otu[, !colnames(clean_otu) %in% c("Taxonomy", "Genus_label")]
  rownames(abundance_matrix) <- clean_otu$Genus_label
  # Delete taxa with total abundance equal to zero across all samples
  abundance_matrix <- abundance_matrix[rowSums(abundance_matrix) > 0, ]
  
  # Match shared samples between community matrix and environmental dataset
  shared_samples <- intersect(colnames(abundance_matrix), env_sample_id)
  if (length(shared_samples) == 0) {
    stop("No overlapping sample IDs detected between microbial genus matrix and environmental dataset, please check sample naming consistency")
  }
  
  # Reorder matrix columns strictly follow environmental sample ID sequence
  abundance_matrix <- abundance_matrix[, shared_samples, drop = FALSE]
  abundance_matrix <- abundance_matrix[, match(env_sample_id, colnames(abundance_matrix)), drop = FALSE]
  
  return(list(
    community_matrix = abundance_matrix,
    matched_sample_id = colnames(abundance_matrix)
  ))
}

# --------------------------- Function 2: Reorder environmental matrix according to community sample sequence ---------------------------
#' Synchronize row order of standardized environmental matrix with microbial sample columns
#' @param env_scaled_matrix Z-score normalized environmental matrix from VIF filtered result
#' @param target_sample_vector Ordered sample ID vector from processed community matrix
#' @return Reordered environmental matrix aligned with microbial data
align_environment_matrix <- function(env_scaled_matrix, target_sample_vector) {
  env_scaled_matrix[target_sample_vector, , drop = FALSE]
}

# --------------------------- Function 3: Core RDA Model Fitting & Permutation Significance Test ---------------------------
#' Perform redundancy analysis with Hellinger transformation, full permutation ANOVA for global significance
#' @param community_matrix Genus abundance matrix (row = genus, column = sample)
#' @param filtered_env_matrix Standardized environmental matrix after VIF<5 filtering
#' @param sample_metadata Sample ID + geographic grouping metadata table
#' @param perm_times Permutation cycle number for Monte-Carlo significance test
#' @return Complete RDA result list including model object, ANOVA table, axis explained variance, site/environment loading coordinates
run_single_community_rda <- function(community_matrix, filtered_env_matrix, sample_metadata, perm_times = 999) {
  if (ncol(filtered_env_matrix) < 1) {
    warning("No valid environmental predictors retained after VIF filtering, terminate RDA calculation")
    return(NULL)
  }
  if (!identical(rownames(filtered_env_matrix), colnames(community_matrix))) {
    stop("Sample row order mismatch: environmental matrix rows inconsistent with community matrix columns")
  }
  
  # Hellinger transformation for species community data (standard for RDA ordination)
  species_hellinger <- vegan::decostand(t(community_matrix), method = "hellinger")
  env_data <- as.data.frame(filtered_env_matrix)
  
  # Automatically construct multi-predictor linear formula
  predictor_vars <- colnames(env_data)
  quoted_vars <- paste0("`", predictor_vars, "`")
  rda_formula <- as.formula(paste("species_hellinger ~", paste(quoted_vars, collapse = " + ")))
  
  # Fit RDA constrained ordination model + permutation ANOVA test
  rda_model <- rda(rda_formula, data = env_data)
  anova_result <- anova(rda_model, permutations = perm_times)
  
  # Extract percentage explained variance of top two RDA axes
  rda_summary <- summary(rda_model)
  axis1_explain <- round(rda_summary$cont$importance[2, 1] * 100, 2)
  axis2_explain <- round(rda_summary$cont$importance[2, 2] * 100, 2)
  
  # Extract sample (site) ordination coordinates and merge grouping metadata
  site_coords <- as.data.frame(scores(rda_model, display = "sites"))
  site_coords$SampleID <- rownames(site_coords)
  site_coords <- merge(site_coords, sample_metadata, by = "SampleID", all.x = TRUE)
  
  # Extract environmental factor loading vector coordinates
  env_coords <- as.data.frame(scores(rda_model, display = "bp"))
  env_coords$Environmental_Factor <- rownames(env_coords)
  
  # Convert permutation ANOVA output to standard data frame table
  anova_table <- as.data.frame(anova_result)
  anova_table$Variable_Term <- rownames(anova_table)
  
  # Organize single-row axis explanation summary table
  overall_explanation <- data.frame(
    RDA1_Explained_Rate = axis1_explain,
    RDA2_Explained_Rate = axis2_explain,
    stringsAsFactors = FALSE
  )
  
  return(list(
    rda_model = rda_model,
    permutation_anova = anova_result,
    anova_stat_table = anova_table,
    axis_explanation = overall_explanation,
    rda1_explain = axis1_explain,
    rda2_explain = axis2_explain,
    site_coordinate = site_coords,
    environmental_coordinate = env_coords
  ))
}

# --------------------------- Function 4: Standardized RDA Biplot Export Function (Dual PNG+PDF Output) ---------------------------
#' Generate publication-ready RDA ordination biplot with group confidence ellipses and environmental loading arrows
#' @param rda_result Full result list output from run_single_community_rda core function
#' @param plot_title Descriptive figure title for manuscript
#' @param file_suffix Unique semantic suffix for standardized output filenames
#' @param save_dir Unified output folder path
export_rda_ordination_plot <- function(rda_result, plot_title, file_suffix, save_dir) {
  if (is.null(rda_result)) return()
  if (dev.cur() > 1) dev.off()
  
  rda_plot <- ggplot() +
    geom_point(
      data = rda_result$site_coordinate,
      aes(x = RDA1, y = RDA2, color = Group),
      size = 3, alpha = 0.8
    ) +
    stat_ellipse(
      data = rda_result$site_coordinate,
      aes(x = RDA1, y = RDA2, color = Group),
      linewidth = 1
    ) +
    geom_segment(
      data = rda_result$environmental_coordinate,
      aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "blue", alpha = 0.8, linewidth = 0.8
    ) +
    geom_text(
      data = rda_result$environmental_coordinate,
      aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = Environmental_Factor),
      color = "darkblue", size = 3.5, check_overlap = TRUE
    ) +
    labs(
      title = plot_title,
      x = paste0("RDA1 (", rda_result$rda1_explain, "%)"),
      y = paste0("RDA2 (", rda_result$rda2_explain, "%)"),
      color = "Sampling Geographic Group"
    ) +
    scale_color_manual(values = group_color_config, breaks = group_factor_level) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80")
  
  # Export high-resolution 300 DPI raster PNG for supplementary materials
  ggsave(
    filename = file.path(save_dir, paste0("RDA_Ordination_Biplot_", file_suffix, ".png")),
    plot = rda_plot, width = 10, height = 8, dpi = 300
  )
  # Export lossless vector PDF for formal journal manuscript submission
  ggsave(
    filename = file.path(save_dir, paste0("RDA_Ordination_Biplot_", file_suffix, ".pdf")),
    plot = rda_plot, width = 10, height = 8
  )
  return(rda_plot)
}

# --------------------------- Function 5: End-to-End Single Community Integrated RDA Wrapper Function ---------------------------
#' Complete independent RDA analysis pipeline for one single microbial community dataset
#' @param otu_excel_path Full standardized path of genus abundance excel table
#' @param env_sample_id Sample ID vector extracted from preprocessed environmental metadata
#' @param filtered_env_matrix VIF<5 Z-score normalized environmental matrix
#' @param sample_meta Sample ID + geographic grouping metadata table
#' @param community_tag Unique descriptive identifier for chart title and output file prefix
#' @param perm_times Permutation iteration count for global significance ANOVA test
#' @param out_dir Unified standardized output directory path
#' @return Full RDA statistical result list containing coordinates and ANOVA tables
single_community_rda_analysis <- function(otu_excel_path,
                                          env_sample_id,
                                          filtered_env_matrix,
                                          sample_meta,
                                          community_tag,
                                          perm_times = 999,
                                          out_dir) {
  # Step1: Preprocess genus abundance matrix, filter zero-abundance taxa, match sample IDs
  otu_proc_res <- preprocess_microbial_otu(otu_excel_path, env_sample_id)
  # Step2: Synchronize environmental matrix row order with microbial sample columns
  aligned_env <- align_environment_matrix(filtered_env_matrix, otu_proc_res$matched_sample_id)
  # Step3: Fit core RDA ordination model + permutation ANOVA test
  rda_res <- run_single_community_rda(
    community_matrix = otu_proc_res$community_matrix,
    filtered_env_matrix = aligned_env,
    sample_metadata = sample_meta,
    perm_times = perm_times
  )
  # Step4: Generate and export standardized RDA biplot dual-format figures
  plot_title <- paste0("RDA Ordination Analysis: ", community_tag, " | All Predictors VIF < 5")
  export_rda_ordination_plot(
    rda_result = rda_res,
    plot_title = plot_title,
    file_suffix = gsub(" ", "_", community_tag),
    save_dir = out_dir
  )
  return(rda_res)
}

# ====================== Single Community Full Pipeline Execution Demo ======================
# Standardized cross-platform full input file path for rhizosphere bacterial genus table
demo_otu_file <- file.path(input_dir, "Rhizosphere_Bacteria_Genus.xlsx")
# Terminate script if microbial abundance input file missing
if (!file.exists(demo_otu_file)) {
  stop(paste("Missing required genus abundance input file:", demo_otu_file))
}

# Run complete end-to-end RDA analysis workflow for rhizosphere bacterial community
single_rda_result <- single_community_rda_analysis(
  otu_excel_path = demo_otu_file,
  env_sample_id = rhizo_analysis_res$env_processed_data$sample_id,
  filtered_env_matrix = rhizo_analysis_res$vif_output$env_data_vif5,
  sample_meta = rhizo_analysis_res$env_processed_data$env_meta,
  community_tag = "Rhizosphere Bacteria",
  perm_times = permutation_times,
  out_dir = output_dir
)

# Organize all multi-sheet statistical tables and export integrated Excel file
single_result_list <- list(
  Site_Coordinates = single_rda_result$site_coordinate,
  Environmental_Loading_Vectors = single_rda_result$environmental_coordinate,
  Permutation_ANOVA_Table = single_rda_result$anova_stat_table,
  Axis_Variance_Explained_Rate = single_rda_result$axis_explanation
)
writexl::write_xlsx(
  single_result_list,
  file.path(output_dir, paste0("RDA_Full_Statistical_Result_", "Rhizosphere_Bacteria", ".xlsx"))
)

# Standardized completion running log with separator
message("\n===== Single Microbial Community RDA Ordination Analysis Completed =====")
message("Target dataset: Rhizosphere Bacteria genus abundance matrix")
message("All dual-format RDA biplot figures and multi-sheet statistical Excel tables saved to unified ./output folder")