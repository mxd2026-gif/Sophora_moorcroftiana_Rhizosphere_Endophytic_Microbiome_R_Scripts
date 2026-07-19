# Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
Full standardized serial R analysis pipeline for rhizosphere & root endophytic microbiome multi-omics of Sophora moorcroftiana.

## Repository Overview
Continuous numbered scripts 01 ~ 21, unified global parameters, cross-platform compatible (Windows/macOS/Linux).
All plotting functions support dual export: 300 DPI PNG raster + lossless vector PDF for journal submission.
Uniform environmental factor list, sampling site grouping color palette, global significance threshold α=0.05.

## Pipeline Workflow Sequence
1. Field sampling spatial mapping & terrain 3D visualization
2. Alpha diversity (Shannon) statistical analysis
3. Community dissimilarity boxplot & PCoA ordination
4. Inter-group explanatory R² comparison bar chart
5. Taxonomic composition stacked bar (Phylum / Class / Genus)
6. Core species occupancy-abundance screening
7. Differential taxa ANOVA + Tukey post-hoc test, Welch t-test
8. Three-set Venn diagram for shared & unique taxa
9. Taxonomic contribution screening
10. Five-layer Taxonomy-KO functional Sankey alluvial diagram
11. Environmental factor VIF multicollinearity filtering
12. RDA constrained ordination linking community & filtered environment
13. Mantel test + environmental pairwise Pearson correlation integrated heatmap
14. Spearman rank correlation between functional features & soil physicochemical factors

## Script Function Directory Table
| Script No. | Core Function Description |
| ---- | ---- |
| 01 | Sampling site spatial distribution map |
| 02 | Elevation 3D cone terrain visualization |
| 03 | Alpha (Shannon) diversity statistical analysis |
| 04 | Within-group Bray-Curtis distance boxplot |
| 05 | PCoA unconstrained ordination analysis |
| 06 | Pairwise permutation R² comparison bar chart |
| 07 | Top10 Phylum & Class stacked abundance bar |
| 08 | Core taxa occupancy & total abundance screening |
| 09 | Top10 Genus boxplot + ANOVA Tukey test |
| 10 | Top10 KEGG Level3 pathway bar + ANOVA Tukey |
| 11 | Differential Genus Welch t-test stacked bar |
| 12 | Three-way Venn diagram for community overlap |
| 13 | Venn diagram of differential / core taxa |
| 14 | Taxonomic contribution screening analysis |
| 15 | Species-Functional feature correlation heatmap |
| 16 | Taxon & KO gene abundance annotated heatmap |
| 17 | Five-layer Phylum-Class-Genus-KEGG-Level3 Sankey diagram |
| 18 | Environmental factor VIF multicollinearity filter |
| 19 | Microbial community & environment RDA ordination |
| 20 | Mantel test + environmental Pearson correlation combined plot |
| 21 | Functional feature & environment Spearman correlation heatmap |

## Running Prerequisite
1. Create two folders `./input` and `./output` under your working directory
2. Place all Excel raw data tables into `./input` folder
3. All dependent packages will be installed automatically when executing each script
4. All statistical figures & Excel result tables will be auto-saved to `./output`

## Unified Standardization Specifications
1. Fixed unified input directory `./input`, output directory `./output` across all 01~21 scripts
2. Global variable `global_alpha = 0.05` for all significance hypothesis tests
3. Fixed four geographic sampling group color palette consistent in all visualization scripts
4. Every core custom function equipped with complete roxygen comment, supports external batch loop invocation
5. Cross-platform path compatible, works normally on Windows / macOS / Linux
