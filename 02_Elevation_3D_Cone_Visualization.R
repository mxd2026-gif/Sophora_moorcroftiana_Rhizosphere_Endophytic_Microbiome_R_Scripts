# ==============================================================================
# Script Name: 02_Elevation_3D_Cone_Visualization.R
# Repository: https://github.com/mxd2026-gif/Sophora_moorcroftiana_Rhizosphere_Endophytic_Microbiome_R_Scripts
# Purpose: Generate 3D cone-shaped column plot for geographic elevation visualization
# Input: Elevation_Location_Data.csv stored in ./input, required columns: Location, Elevation, Color
# Output: PNG preview, EPS & PDF vector figures, interactive HTML widget saved in ./output folder
# Dependencies: rgl, htmlwidgets
# ==============================================================================

# ---------------- Install & load required packages ----------------
required_packages <- c("rgl", "htmlwidgets")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message(paste("Package", pkg, "installed successfully"))
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

# Input file name requirement: CSV with 3 columns: Location, Elevation, Color
input_file_name <- "Elevation_Location_Data.csv"
full_input_path <- file.path(input_dir, input_file_name)

# ---------------- Custom reusable functions ----------------
#' Read elevation metadata from CSV file
#' @param file_path Full path of input CSV file
#' @return Clean Location-elevation-color data frame
read_elevation_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("Error: Input file not found at path:", file_path,
               "\nPlease place 'Elevation_Location_Data.csv' inside ./input folder"))
  }
  raw_data <- read.csv(file_path, stringsAsFactors = FALSE)
  # Check required columns
  required_cols <- c("Location", "Elevation", "Color")
  if (!all(required_cols %in% colnames(raw_data))) {
    stop("Input CSV must contain these columns: Location, Elevation, Color")
  }
  raw_data$Location <- factor(raw_data$Location, levels = unique(raw_data$Location))
  return(raw_data)
}

#' Draw single 3D cone geometry for each sampling site
#' @param x X coordinate index
#' @param y Fixed Y axis position
#' @param z Elevation value (cone height)
#' @param fill_color Hex color code for cone surface
#' @param cone_base_radius Radius of cone bottom base
draw_single_3d_cone <- function(x, y, z, fill_color, cone_base_radius = 0.3) {
  base_center_point <- c(x, y, 0)
  # Generate circular base coordinates
  theta_sequence <- seq(from = 0, to = 2 * pi, length.out = 50)
  base_x_coords <- base_center_point[1] + cone_base_radius * cos(theta_sequence)
  base_y_coords <- base_center_point[2] + cone_base_radius * sin(theta_sequence)
  base_z_coords <- rep(0, length(theta_sequence))
  
  # Construct triangular facets to form cone surface
  for (i in seq_len(length(theta_sequence) - 1)) {
    triangle_vertices <- rbind(
      apex.point <- c(x, y, z),
      c(base_x_coords[i], base_y_coords[i], base_z_coords[i]),
      c(base_x_coords[i + 1], base_y_coords[i + 1], base_z_coords[i + 1])
    )
    polygon3d(triangle_vertices, col = fill_color, alpha = 1, add = TRUE)
  }
  # Add elevation value text annotation above cone apex
  text3d(x = x, y = y, z = z + 200, text = as.character(z), col = "black", cex = 1.2)
}

#' Fit spline smooth curve above all 3D cones to show elevation trend
#' @param x_vec Numeric X index of each sampling site
#' @param z_vec Elevation values for each location
#' @param y_fixed Fixed Y position for trend line
draw_elevation_trend_curve <- function(x_vec, z_vec, y_fixed = 1.2) {
  spline_fit_result <- spline(x = x_vec, y = z_vec + 300, n = 100)
  smooth_x <- spline_fit_result$x
  smooth_z <- spline_fit_result$y
  smooth_y <- rep(y_fixed, length(smooth_x))
  
  # Draw continuous black smooth line
  for (i in seq_len(length(smooth_x) - 1)) {
    segments3d(
      x = c(smooth_x[i], smooth_x[i + 1]),
      y = c(smooth_y[i], smooth_y[i + 1]),
      z = c(smooth_z[i], smooth_z[i + 1]),
      color = "black", lwd = 2
    )
  }
}

#' Main function to render full 3D elevation visualization
#' @param elevation_data Input data frame with Location, Elevation, Color columns
render_3d_elevation_plot <- function(elevation_data) {
  # Extract plotting variables
  site_index <- seq_len(nrow(elevation_data))
  site_labels <- as.character(elevation_data$Location)
  fixed_y_pos <- 1
  elevation_values <- elevation_data$Elevation
  site_colors <- elevation_data$Color
  
  # Set global axis range
  x_axis_limit <- c(0.5, nrow(elevation_data) + 0.5)
  y_axis_limit <- c(0.5, 1.5)
  z_axis_limit <- c(0, 4000)
  
  # Initialize empty 3D scene
  plot3d(
    x = 0, y = 0, z = 0,
    xlim = x_axis_limit,
    ylim = y_axis_limit,
    zlim = z_axis_limit,
    xlab = "", ylab = "", zlab = "",
    type = "n"
  )
  
  # Iterate to draw each 3D cone
  for (i in site_index) {
    draw_single_3d_cone(
      x = i,
      y = fixed_y_pos,
      z = elevation_values[i],
      fill_color = site_colors[i],
      cone_base_radius = 0.3
    )
  }
  
  # Add smooth elevation trend curve
  draw_elevation_trend_curve(x_vec = site_index, z_vec = elevation_values)
  
  # Customize X axis (Location labels)
  axis3d(
    edge = "x-",
    at = site_index,
    labels = site_labels,
    tick = FALSE,
    line = -0.5,
    cex.axis = 1,
    col.axis = "black",
    hasTicks = FALSE
  )
  
  # Customize Z axis (elevation ticks)
  axis3d(
    edge = "z-",
    at = seq(from = 0, to = 4000, by = 1000),
    labels = seq(from = 0, to = 4000, by = 1000),
    tick = FALSE,
    line = -0.5,
    cex.axis = 1,
    col.axis = "black",
    hasTicks = FALSE
  )
  
  # Hide Y axis completely
  axis3d(edge = "y-", labels = FALSE, tick = FALSE)
  
  # Remove default 3D bounding box
  rgl.clear(type = "bbox")
  
  # Optimize 3D viewing angle
  view3d(theta = -70, phi = 35, fov = 60, zoom = 0.9)
  
  # Return interactive 3D widget
  return(rglwidget())
}

# ---------------- Main execution workflow ----------------
# Step 1: Read input elevation data
elevation_dataset <- read_elevation_data(full_input_path)
message(paste("Successfully loaded", nrow(elevation_dataset), "geographic sampling sites"))

# Step 2: Generate 3D visualization
interactive_3d_widget <- render_3d_elevation_plot(elevation_dataset)

# Step 3: Export multiple format figures
png_output_path <- file.path(output_dir, "3D_Elevation_Cone_Plot.png")
eps_output_path <- file.path(output_dir, "3D_Elevation_Cone_Plot.eps")
pdf_output_path <- file.path(output_dir, "3D_Elevation_Cone_Plot.pdf")

rgl.snapshot(filename = png_output_path)
rgl.postscript(filename = eps_output_path, fmt = "eps")
rgl.postscript(filename = pdf_output_path, fmt = "pdf")

# Export interactive 3D HTML widget
htmlwidgets::saveWidget(interactive_3d_widget, file.path(output_dir, "Elevation_3D_Interactive.html"), selfcontained = TRUE)
message("All visualization files successfully saved to ./output folder")

# Display interactive 3D plot in RStudio viewer
print(interactive_3d_widget)

# Clear memory & release cache
rm(list = ls())
gc()