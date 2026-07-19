# ==============================================================================
# Script Name: 20_Mantel_Community_Environment_Correlation_Multicollinearity.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Comprehensive multi-module environmental-microbe linkage analysis:
#          1) Environmental factor pairwise Pearson correlation matrix & significance test
#          2) VIF multicollinearity diagnosis bar chart & statistics
#          3) Bray-Curtis community distance vs environmental Euclidean distance Mantel permutation test
#          4) Integrated combined plot: environmental correlation heatmap + Mantel linkage curve overlay via linkET
# Input: 1env16-8.xlsx environmental table, genus/KEGG functional abundance matrices stored in unified ./input folder
# Output: Multi-sheet integrated statistical Excel table, dual-format 300 DPI PNG + vector PDF VIF plot & Mantel combined heatmap saved in ./output
# Dependencies: tidyverse, RColorBrewer, ggnewscale, readxl, writexl, vegan, openxlsx, car, devtools, linkET, numDeriv
# Standardization: Fully aligned with script 01~19 unified global parameters, cross-platform path rules, journal dual-format output specs
# ==============================================================================
options(scipen = 999, digits = 4)
# Global unified statistical threshold consistent with full serial analysis pipeline
global_alpha <- 0.05

# --------------------------- Install and load all dependent packages ---------------------------
required_cran_pkgs <- c("tidyverse", "RColorBrewer", "ggnewscale", "readxl", "writexl",
                        "vegan", "openxlsx", "car", "devtools", "numDeriv")
for (pkg in required_cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("CRAN package %s installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# Install linkET from GitHub (required for Mantel combined correlation heatmap)
if (!requireNamespace("linkET", quietly = TRUE)) {
  devtools::install_github("Hy4m/linkET", force = TRUE)
  message("linkET visualization package installed from GitHub repository")
}
library(linkET)

# --------------------------- Cross-platform relative path configuration (Unified 01~19 standard ./input & ./output) ---------------------------
# Fix inconsistent raw_input/output_result directory to unified ./input / ./output matching serial scripts 01-19
input_dir  <- "./input"
output_dir <- "./output"

# Auto create standard project folders with running status prompt
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created standard output directory: ./output")
}

# Global unified four geographic sampling group color palette (consistent with script 03~19 full pipeline)
group_color_config <- c(
  Shigatse = "#FF6347",
  Lhasa    = "#32CD32",
  Shannan  = "#9370DB",
  Nyingchi = "#FFD700"
)
group_factor_level <- names(group_color_config)

# --------------------------- Global fixed standardized sample ID mapping ---------------------------
# Short unified environmental sample ID (consistent across VIF/RDA scripts 18/19)
env_sample_ids <- c("RK1", "RK2", "RK3", "RK4", "RK5",
                    "CN1", "CN2", "CN3", "CN4", "CN5",
                    "SN1", "SN2", "SN3", "SN4", "SN5",
                    "ML1", "ML2", "ML3", "ML4", "ML5")
# Rhizosphere sample ID prefix R-
rhizo_sample_ids <- c("RRK1", "RRK2", "RRK3", "RRK4", "RRK5",
                      "RCN1", "RCN2", "RCN3", "RCN4", "RCN5",
                      "RSN1", "RSN2", "RSN3", "RSN4", "RSN5",
                      "RML1", "RML2", "RML3", "RML4", "RML5")
# Endophytic sample ID prefix E-
endo_sample_ids <- c("ERK1", "ERK2", "ERK3", "ERK4", "ERK5",
                     "ECN1", "ECN2", "ECN3", "ECN4", "ECN5",
                     "ESN1", "ESN2", "ESN3", "ESN4", "ESN5",
                     "EML1", "EML2", "EML3", "EML4", "EML5")

# --------------------------- Function 1: Sample ID conversion standardization ---------------------------
#' Convert long rhizosphere/endophytic extended sample ID to unified short environmental sample ID
#' @param df Community or environmental data frame matrix
#' @param original_id_vector Raw full extended sample ID vector (rhizo/endo)
#' @param target_id_vector Standard unified short environmental sample ID vector
#' @return Data frame with standardized row/column sample ID matching environmental dataset
convert_sample_ids <- function(df, original_id_vector, target_id_vector) {
  if (all(rownames(df) %in% original_id_vector)) {
    idx <- match(rownames(df), original_id_vector)
    rownames(df) <- target_id_vector[idx]
  }
  if (all(colnames(df) %in% original_id_vector)) {
    idx <- match(colnames(df), original_id_vector)
    colnames(df) <- target_id_vector[idx]
  }
  return(df)
}

# --------------------------- Function 2: Single community preprocessing & Bray-Curtis distance calculation ---------------------------
#' Read genus/KEGG abundance table, aggregate total abundance, standardize sample ID, compute Bray-Curtis dissimilarity matrix
#' @param file_path Excel file full standardized path of community abundance matrix
#' @param id_type Sample ID tag: rhizo / endo
#' @return List: standardized community abundance dataframe, Bray-Curtis distance dissimilarity matrix
preprocess_single_community <- function(file_path, id_type = "rhizo") {
  raw_comm <- read_excel(file_path)
  tax_col <- if ("#Taxonomy" %in% colnames(raw_comm)) "#Taxonomy" else "KEGG_Name"
  
  comm_agg <- raw_comm %>%
    group_by(.data[[tax_col]]) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
    column_to_rownames(var = tax_col) %>%
    t() %>%
    as.data.frame()
  
  # Convert extended sample ID to unified short environmental ID
  if (id_type == "rhizo") {
    comm_agg <- convert_sample_ids(comm_agg, rhizo_sample_ids, env_sample_ids)
  } else if (id_type == "endo") {
    comm_agg <- convert_sample_ids(comm_agg, endo_sample_ids, env_sample_ids)
  }
  
  # Reorder rows strictly follow fixed environmental sample ID sequence
  comm_agg <- comm_agg[env_sample_ids, ]
  comm_dist <- vegdist(comm_agg, method = "bray", na.rm = TRUE)
  
  return(list(community_matrix = comm_agg, distance_matrix = comm_dist))
}

# --------------------------- Function 3: Single community Mantel permutation test core function ---------------------------
#' Perform Mantel test between community Bray-Curtis dissimilarity and standardized environmental Euclidean distance matrix
#' @param comm_dist Bray-Curtis dissimilarity matrix of target microbial community
#' @param env_dist Standardized environmental factor Euclidean distance matrix
#' @param community_name Unique descriptive community name tag for table title
#' @param perm_times Permutation iteration count for Monte-Carlo significance test, default 999
#' @return Single-row Mantel test statistical result table
run_single_mantel_test <- function(comm_dist, env_dist, community_name, perm_times = 999) {
  mantel_out <- mantel(comm_dist, env_dist, permutations = perm_times, na.rm = TRUE)
  mantel_df <- data.frame(
    community = community_name,
    mantel_r = mantel_out$statistic,
    mantel_p = mantel_out$signif,
    stringsAsFactors = FALSE
  )
  return(mantel_df)
}

# --------------------------- Step 1: Load and preprocess standardized environmental dataset ---------------------------
env_file_path <- file.path(input_dir, "1env16-8.xlsx")
# Input file integrity pre-check
if (!file.exists(env_file_path)) {
  stop(paste("Missing mandatory environmental factor input file:", env_file_path))
}
env_raw <- read_excel(env_file_path, sheet = "all")
env_target_cols <- c("NO3--N", "Ex-Ca", "TP", "EC", "SM_20cm", "PH", "NH4+-N", "Ex-Mg")

numeric_index <- which(sapply(env_raw, is.numeric))
colnames(env_raw)[numeric_index] <- env_target_cols

env_data <- env_raw %>%
  select(SampleID, all_of(env_target_cols)) %>%
  column_to_rownames("SampleID") %>%
  as.data.frame()

# Sample ID completeness validation
if (!all(env_sample_ids %in% rownames(env_data))) {
  stop("Partial standard sample IDs missing in environmental dataset, please check sample naming consistency")
}
env_data <- env_data[env_sample_ids, ]

# Numeric variable validation
if (!all(sapply(env_data, is.numeric))) {
  stop("Non-numeric categorical variables detected in environmental factor matrix")
}
env_variable_order <- colnames(env_data)
message("Standardized environmental factor dataset loaded and validated successfully")

# Z-score standardization + Euclidean distance matrix for Mantel test
env_standardized <- decostand(env_data, method = "standardize", na.rm = TRUE)
env_euclidean_dist <- vegdist(env_standardized, method = "euclidean", na.rm = TRUE)

# --------------------------- Step 2: Environmental pairwise Pearson correlation & significance matrix ---------------------------
pearson_corr_matrix <- cor(env_data, method = "pearson", use = "complete.obs")
p_value_matrix <- matrix(NA, nrow = ncol(env_data), ncol = ncol(env_data),
                         dimnames = list(colnames(env_data), colnames(env_data)))

# Calculate two-tailed Pearson correlation p-value matrix
for (i in seq_len(ncol(env_data))) {
  for (j in seq_len(ncol(env_data))) {
    if (i != j) {
      p_value_matrix[i, j] <- cor.test(env_data[, i], env_data[, j], method = "pearson")$p.value
    } else {
      p_value_matrix[i, j] <- 1
    }
  }
}

# Reshape correlation matrix to long format statistical table
env_correlation_result <- as.data.frame(pearson_corr_matrix) %>%
  rownames_to_column("Env_Variable_1") %>%
  pivot_longer(cols = -Env_Variable_1, names_to = "Env_Variable_2", values_to = "Pearson_r") %>%
  mutate(P_value = as.vector(p_value_matrix)) %>%
  mutate(
    Significance = ifelse(P_value < global_alpha, "Significant", "Non-significant"),
    Correlation_Strength = case_when(
      abs(Pearson_r) < 0.2 ~ "Weak",
      abs(Pearson_r) >= 0.2 & abs(Pearson_r) < 0.3 ~ "Moderate",
      TRUE ~ "Strong"
    )
  )
message("Environmental factor pairwise Pearson correlation & significance calculation completed")

# --------------------------- Step 3: VIF multicollinearity diagnosis & horizontal bar plot ---------------------------
# Dummy linear regression for VIF calculation
vif_linear_model <- lm(rnorm(nrow(env_data)) ~ ., data = env_data)
vif_result_vector <- car::vif(vif_linear_model)

vif_result_df <- data.frame(
  Environmental_Variable = names(vif_result_vector),
  VIF_Value = as.numeric(vif_result_vector),
  stringsAsFactors = FALSE
) %>%
  mutate(
    VIF_Plot_Value = ifelse(is.infinite(VIF_Value), 100, VIF_Value),
    Environmental_Variable = factor(Environmental_Variable, levels = rev(env_variable_order))
  )

# Warn high multicollinearity variables exceeding VIF=10 threshold
high_multicollinearity <- filter(vif_result_df, VIF_Plot_Value > 10)
if (nrow(high_multicollinearity) > 0) {
  message("WARNING: High multicollinearity detected (VIF > 10), variables listed below:")
  print(high_multicollinearity[, c("Environmental_Variable", "VIF_Value")])
}

# Standardized VIF horizontal bar plot (journal publication style)
vif_plot <- ggplot(vif_result_df, aes(x = VIF_Plot_Value, y = Environmental_Variable)) +
  geom_segment(aes(x = 0, xend = VIF_Plot_Value, y = Environmental_Variable, yend = Environmental_Variable),
               color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Variance Inflation Factor (VIF)", y = "", title = "Environmental Factor Multicollinearity VIF Diagnosis") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    panel.grid.major.y = element_blank()
  )
message("VIF multicollinearity diagnosis & plotting workflow completed")

# --------------------------- Step 4: Community dataset metadata list for batch Mantel analysis ---------------------------
community_metadata_list <- list(
  list(
    file_path = file.path(input_dir, "g_rb.xlsx"),
    sample_id_type = "rhizo",
    community_label = "Rhizosphere Bacterial Community (Genus-level)",
    sheet_tag = "Rhizosphere_Bacteria"
  ),
  list(
    file_path = file.path(input_dir, "g_rf.xlsx"),
    sample_id_type = "rhizo",
    community_label = "Rhizosphere Fungal Community (Genus-level)",
    sheet_tag = "Rhizosphere_Fungi"
  ),
  list(
    file_path = file.path(input_dir, "g_eb.xlsx"),
    sample_id_type = "endo",
    community_label = "Root Endophytic Bacterial Community (Genus-level)",
    sheet_tag = "Endophytic_Bacteria"
  ),
  list(
    file_path = file.path(input_dir, "g_ef.xlsx"),
    sample_id_type = "endo",
    community_label = "Root Endophytic Fungal Community (Genus-level)",
    sheet_tag = "Endophytic_Fungi"
  ),
  list(
    file_path = file.path(input_dir, "2Rkegg_Name.xlsx"),
    sample_id_type = "rhizo",
    community_label = "Rhizosphere Functional Genes (KEGG KO)",
    sheet_tag = "Rhizosphere_KEGG_Function"
  ),
  list(
    file_path = file.path(input_dir, "2Ekegg_Name.xlsx"),
    sample_id_type = "endo",
    community_label = "Root Endophytic Functional Genes (KEGG KO)",
    sheet_tag = "Endophytic_KEGG_Function"
  )
)

# Demo execution: Single community Mantel test example (wrap for full batch loop)
target_comm <- community_metadata_list[[1]]
# Input community file pre-check
if (!file.exists(target_comm$file_path)) {
  stop(paste("Missing community abundance input file:", target_comm$file_path))
}
comm_proc_res <- preprocess_single_community(target_comm$file_path, target_comm$sample_id_type)
single_mantel_res <- run_single_mantel_test(
  comm_dist = comm_proc_res$distance_matrix,
  env_dist = env_euclidean_dist,
  community_name = target_comm$community_label,
  perm_times = 999
)

# Initialize global Mantel result storage table
mantel_total_result <- single_mantel_res

# --------------------------- Step 5: Mantel test result categorical grouping for visualization ---------------------------
mantel_plot_data <- mantel_total_result %>%
  mutate(
    Mantel_R_Category = cut(mantel_r, breaks = c(-Inf, 0.2, 0.3, Inf),
                            labels = c("< 0.2", "0.2 - 0.3", ">= 0.3")),
    Mantel_P_Category = cut(mantel_p, breaks = c(-Inf, 0.005, 0.05, Inf),
                            labels = c("< 0.005", "0.01 - 0.05", ">= 0.05")),
    Significance = ifelse(mantel_p < global_alpha, "Significant", "Non-significant"),
    Correlation_Strength = case_when(
      abs(mantel_r) < 0.2 ~ "Weak",
      abs(mantel_r) >= 0.2 & abs(mantel_r) < 0.3 ~ "Moderate",
      TRUE ~ "Strong"
    )
  )

# --------------------------- Step 6: Integrated combined plot: Environmental correlation heatmap + overlay Mantel linkage curves ---------------------------
env_corr_plot_data <- correlate(env_data)
mantel_combined_plot <- qcorrplot(env_corr_plot_data,
                                  grid_col = "grey50",
                                  grid_size = 0.5,
                                  type = "upper",
                                  diag = FALSE) +
  geom_square() +
  geom_mark(size = 4.2, only_mark = TRUE, sig_level = c(0.05, 0.01, 0.001),
            sig_thres = global_alpha, colour = "white", stroke = 0.1) +
  geom_couple(data = mantel_plot_data,
              aes(color = Mantel_P_Category, size = Mantel_R_Category),
              label.size = 3.88, nudge_x = 0.2, curvature = nice_curvature(by = "from"),
              point_fill = "white", point_color = "gray50") +
  scale_fill_gradient2(limits = c(-0.8, 0.8), mid = "white", low = "#2a7bbb", high = "#e67e22") +
  scale_size_manual(values = c("< 0.2" = 0.3, "0.2 - 0.3" = 0.7, ">= 0.3" = 1.1)) +
  scale_color_manual(values = c("< 0.005" = "#2d5d7a", "0.01 - 0.05" = "#5a90b3", ">= 0.05" = "#f39c12")) +
  guides(
    size = guide_legend(title = "Mantel's r correlation strength", order = 2),
    colour = guide_legend(title = "Mantel permutation P-value", order = 1),
    fill = guide_colorbar(title = "Pearson's r environmental correlation", order = 3)
  ) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

# --------------------------- Step 7: Batch export integrated multi-sheet Excel statistical tables & dual-format figures ---------------------------
# Create unified integrated result workbook
workbook <- createWorkbook()
addWorksheet(workbook, "Environmental_Pearson_Correlation")
writeData(workbook, "Environmental_Pearson_Correlation", env_correlation_result)

addWorksheet(workbook, "VIF_Multicollinearity_Statistics")
writeData(workbook, "VIF_Multicollinearity_Statistics", select(vif_result_df, -VIF_Plot_Value))

addWorksheet(workbook, "Mantel_Total_Analysis_Result")
writeData(workbook, "Mantel_Total_Analysis_Result", mantel_plot_data)

# Independent sheet for single community Mantel test record
addWorksheet(workbook, target_comm$sheet_tag)
writeData(workbook, target_comm$sheet_tag, single_mantel_res)

# Save integrated statistical table to standardized output directory
saveWorkbook(workbook, file.path(output_dir, "Mantel_Env_Corr_VIF_Integrated_Statistical_Result.xlsx"), overwrite = TRUE)

# Export Mantel combined linkage heatmap dual-format PNG+PDF
ggsave(mantel_combined_plot,
       filename = file.path(output_dir, "Mantel_Community_Environmental_Integrated_Linkage_Plot.pdf"),
       width = 14, height = 9)
ggsave(mantel_combined_plot,
       filename = file.path(output_dir, "Mantel_Community_Environmental_Integrated_Linkage_Plot.png"),
       width = 14, height = 9, dpi = 300)

# Export VIF multicollinearity diagnosis bar plot dual-format PNG+PDF
ggsave(vif_plot,
       filename = file.path(output_dir, "VIF_Multicollinearity_Diagnosis_Bar_Plot.pdf"),
       width = 8, height = 6)
ggsave(vif_plot,
       filename = file.path(output_dir, "VIF_Multicollinearity_Diagnosis_Bar_Plot.png"),
       width = 8, height = 6, dpi = 300)

# Standardized completion running log with separator
message("\n===== Full Mantel Environmental-Microbe Linkage Multi-Module Analysis Completed =====")
message("Analysis modules executed: Pearson environmental correlation | VIF multicollinearity | Mantel permutation test | Integrated linkage visualization")
message("All integrated statistical Excel tables, 300 DPI PNG raster & lossless vector PDF figures saved to unified ./output folder")