# ==============================================================================
# Script Name: 01_Sampling_Site_Spatial_Map.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Generate spatial distribution map of field sampling sites in Xizang Autonomous Region
# Input: Sampling_Coordinate.csv stored in ./input (Longitude, Latitude, Site, Altitude, Soil_Type)
# Output: PNG raster image and PDF vector map stored in ./output folder
# Dependencies: sf, ggplot2, dplyr, rnaturalearth, rnaturalearthdata, ggspatial, numDeriv
# ==============================================================================

# ---------------- Install & load all required packages ----------------
required_packages <- c("sf", "ggplot2", "dplyr", "rnaturalearth",
                       "rnaturalearthdata", "ggspatial", "numDeriv")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(sprintf("Package %s has been installed successfully", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ---------------- Global relative path configuration ----------------
input_dir <- "./input"
output_dir <- "./output"

# Automatically create folders if directories do not exist
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created input directory: ./input")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Created output directory: ./output")
}

# ---------------- Custom projected coordinate system ----------------
custom_proj_crs <- "+proj=aegd +lat_0=35 +lon_0=105 +ellps=WGS84 +units=m +no_defs"

# ---------------- Core reusable mapping function ----------------
generate_sampling_site_map <- function() {
  # Load national administrative boundary data
  china_admin_sf <- ne_states(country = "China", returnclass = "sf")
  xizang_boundary_sf <- china_admin_sf %>%
    filter(name == "Xizang")
  
  if (nrow(xizang_boundary_sf) == 0) {
    stop("Error: Xizang administrative boundary spatial data cannot be retrieved")
  }
  
  # Coordinate transformation
  xizang_boundary_sf <- st_transform(xizang_boundary_sf, crs = custom_proj_crs)
  
  # Note: rnaturalearth/ggspatial may fail to download boundary data in mainland China; pre-download shapefile to ./input if needed
  # Download global river spatial data and filter target river
  global_river_sf <- ne_download(
    type = "rivers_lake_centerlines",
    category = "physical",
    returnclass = "sf"
  ) %>% st_transform(crs = custom_proj_crs)
  
  yarlung_zangbo_river_sf <- global_river_sf %>%
    filter(grepl("Brahmaputra", name, ignore.case = TRUE))
  
  # Sampling site attribute and coordinate dataset
  # Priority 1: Read external coordinate CSV from ./input
  csv_coord_path <- file.path(input_dir, "Sampling_Coordinate.csv")
  if (file.exists(csv_coord_path)) {
    sampling_site_df <- read.csv(csv_coord_path, header = TRUE)
  } else {
  # Priority 2: Built-in hard-coded coordinates (run directly without input CSV)
    sampling_site_df <- data.frame(
      Longitude = c(88.94396557, 90.95770474, 91.81402589, 94.39786388),
      Latitude  = c(29.24944748, 29.42912261, 29.28141954, 29.35133819),
      Site      = c("Shigatse", "Lhasa", "Shannan", "Nyingchi"),
      Altitude  = c(3801, 3572, 3572, 2944),
      Soil_Type = c("Gravel_soil", "Alpine_meadow", "Brown_soil", "Sandy_loam")
    )
  }
  
  # Convert data frame to spatial sf object
  sampling_point_sf <- st_as_sf(
    sampling_site_df,
    coords = c("Longitude", "Latitude"),
    crs = 4326
  ) %>% st_transform(crs = custom_proj_crs)
  
  # Prepare river annotation with leader line
  river_annotation_df <- data.frame()
  if (nrow(yarlung_zangbo_river_sf) > 0) {
    river_coords_mat <- st_coordinates(yarlung_zangbo_river_sf)
    river_annotation_df <- data.frame(
      label = "Yarlung Zangbo River",
      text_x = river_coords_mat[1, 1] + 50000,
      text_y = river_coords_mat[1, 2] + 50000,
      river_x = river_coords_mat[1, 1],
      river_y = river_coords_mat[1, 2]
    )
  } else {
    river_annotation_df <- data.frame(
      label = "Yarlung Zangbo River",
      text_x = 9200000,
      text_y = 3000000,
      river_x = 9150000,
      river_y = 2950000
    )
  }
  
  # Layered map drawing
  spatial_map <- ggplot() +
    geom_sf(
      data = xizang_boundary_sf,
      fill = "antiquewhite",
      color = "gray50",
      linewidth = 0.4
    ) +
    geom_sf(
      data = yarlung_zangbo_river_sf,
      color = "gray50",
      linewidth = 1,
      alpha = 0.7
    ) +
    geom_segment(
      data = river_annotation_df,
      aes(x = text_x, y = text_y, xend = river_x, yend = river_y),
      color = "gray50",
      linewidth = 0.4
    ) +
    geom_text(
      data = river_annotation_df,
      aes(x = text_x, y = text_y, label = label),
      color = "gray50",
      size = 3,
      vjust = -1,
      hjust = 0
    ) +
    geom_sf(
      data = sampling_point_sf,
      aes(color = Site),
      size = 5,
      alpha = 0.8
    ) +
    scale_color_manual(values = c(
      "Shigatse" = "#FF6347",
      "Lhasa" = "#32CD32",
      "Shannan" = "#9370DB",
      "Nyingchi" = "#FFD700"
    )) +
    annotation_north_arrow(
      location = "tl",
      style = north_arrow_fancy_orienteering
    ) +
    annotation_scale(location = "bl") +
    labs(
      title = "Sampling Sites in Xizang Autonomous Region",
      color = "Sampling Site"
    ) +
    theme(legend.position = "right")
  
  return(spatial_map)
}

# ---------------- Main execution entry ----------------
# Generate map object
sampling_distribution_map <- generate_sampling_site_map()
print(sampling_distribution_map)

# Export high-resolution map files
png_output_path <- file.path(output_dir, "Sampling_Sites_Distribution_Map.png")
pdf_output_path <- file.path(output_dir, "Sampling_Sites_Distribution_Map.pdf")

ggsave(
  filename = png_output_path,
  plot = sampling_distribution_map,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = pdf_output_path,
  plot = sampling_distribution_map,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

message("Map generation completed. Output files saved to ./output directory.")