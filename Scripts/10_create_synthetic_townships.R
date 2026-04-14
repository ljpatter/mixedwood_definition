# ---
# title: "Create synthetic townships"
# author: "Leonard Patterson"
# created: "2025-08-11"
# description: This script creates synthetic intimate and segregated mixedwoods that will 
#             then be used for generating predictions using output from the multi scale BRT model
# ---

# Remove objects from environment
rm(list = ls())

# Load libraries
library(dplyr)
library(sf)
library(terra)
library(landscapemetrics)
library(readr)


# Define township center coordinates (change as needed)
township_center <- c(x = 	623760, y = 6145354)   # Example UTM coordinates

# Define township size (10 km x 10 km)
township_size <- 10000   # in meters
half_size <- township_size / 2

# Set CRS
target_crs <- "EPSG:3400"

# Create township polygon
township_polygon <- st_polygon(list(rbind(
  c(township_center["x"] - half_size, township_center["y"] - half_size),
  c(township_center["x"] + half_size, township_center["y"] - half_size),
  c(township_center["x"] + half_size, township_center["y"] + half_size),
  c(township_center["x"] - half_size, township_center["y"] + half_size),
  c(township_center["x"] - half_size, township_center["y"] - half_size)
))) %>% st_sfc(crs = target_crs)

township_sf <- st_sf(geometry = township_polygon)

# Export township polygon shapefile
dir.create("Output/Spatial Data/Township polygon", recursive = TRUE, showWarnings = FALSE)
st_write(township_sf, "Output/Spatial Data/Township polygon/township_polygon.shp", delete_layer = TRUE)

# Parameters
cell_size <- 5
township_length_m <- 10000
ncol <- nrow <- township_length_m / cell_size
block_size_m <- 500
block_size_pixels <- round(block_size_m / cell_size)
blocks_per_side <- ceiling(ncol / block_size_pixels)

# Convert matrix to raster
matrix_to_raster <- function(mat, crs_string, center_x, center_y, total_size) {
  half_total_size <- total_size / 2
  xmin_val <- center_x - half_total_size
  xmax_val <- center_x + half_total_size
  ymin_val <- center_y - half_total_size
  ymax_val <- center_y + half_total_size
  rast_obj <- rast(nrows = nrow(mat), ncols = ncol(mat),
                   xmin = xmin_val, xmax = xmax_val,
                   ymin = ymin_val, ymax = ymax_val,
                   crs = crs_string)
  values(rast_obj) <- as.vector(t(mat))
  return(rast_obj)
}

# Intimate: random 0/1
generate_intimate <- function() {
  matrix(sample(c(0,1), size = nrow * ncol, replace = TRUE), nrow = nrow, ncol = ncol)
}

# Segregated blocks with half-half config
generate_half_half_block <- function(size) {
  block <- matrix(0, nrow = size, ncol = size)
  side <- sample(c("top", "bottom", "left", "right"), 1)
  if (side %in% c("top", "bottom")) {
    rows_to_fill <- floor(size / 2)
    if (side == "top") {
      block[1:rows_to_fill, ] <- 1
    } else {
      block[(size - rows_to_fill + 1):size, ] <- 1
    }
  } else {
    cols_to_fill <- floor(size / 2)
    if (side == "left") {
      block[, 1:cols_to_fill] <- 1
    } else {
      block[, (size - cols_to_fill + 1):size] <- 1
    }
  }
  return(block)
}

# Stitch segregated blocks
generate_segregated <- function() {
  full_mat <- matrix(0, nrow = nrow, ncol = ncol)
  for (row_block in 0:(blocks_per_side - 1)) {
    for (col_block in 0:(blocks_per_side - 1)) {
      block <- generate_half_half_block(block_size_pixels)
      r_start <- row_block * block_size_pixels + 1
      r_end <- min((row_block + 1) * block_size_pixels, nrow)
      c_start <- col_block * block_size_pixels + 1
      c_end <- min((col_block + 1) * block_size_pixels, ncol)
      block_rows <- r_end - r_start + 1
      block_cols <- c_end - c_start + 1
      if (block_rows != block_size_pixels || block_cols != block_size_pixels) {
        block <- block[1:block_rows, 1:block_cols]
      }
      full_mat[r_start:r_end, c_start:c_end] <- block
    }
  }
  return(full_mat)
}

# Generate landscapes
intimate_mat <- generate_intimate()
segregated_mat <- generate_segregated()

# Convert to SpatRaster and set CRS
intimate_raster <- matrix_to_raster(intimate_mat, target_crs, township_center["x"], township_center["y"], township_size)
segregated_raster <- matrix_to_raster(segregated_mat, target_crs, township_center["x"], township_center["y"], township_size)

values(intimate_raster) <- as.integer(values(intimate_raster))
values(segregated_raster) <- as.integer(values(segregated_raster))

# Define black-white palette
bw_palette <- c("white", "black")

# --- Plotting section with no axis or title ---
old_par <- par(no.readonly = TRUE)
par(mar = c(0, 0, 0, 0))

# On-screen plot
plot(intimate_raster, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
plot(segregated_raster, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)

par(old_par)


# --- Save clean TIFF plots ---

tiff("Output/Figures/Synthetic landscapes/intimate_mixedwood_10m.tiff",
     width = 800, height = 800, res = 600)

par(mar = c(0, 0, 0, 0))
plot(intimate_raster, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
dev.off()

tiff("Output/Figures/Synthetic landscapes/segregated_mixedwood_10m.tiff",
     width = 800, height = 800, res = 600)

par(mar = c(0, 0, 0, 0))
plot(segregated_raster, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
dev.off()




# Calculate CLUMPY
intimate_clumpy <- lsm_c_clumpy(intimate_raster)
segregated_clumpy <- lsm_c_clumpy(segregated_raster)

print(intimate_clumpy)
print(segregated_clumpy)

# Save rasters
dir.create("Output/Spatial Data/Simulated landscapes", recursive = TRUE, showWarnings = FALSE)
writeRaster(intimate_raster, filename = "Output/Spatial Data/Simulated landscapes/intimate_mixedwood_5m.tif", overwrite = TRUE)
writeRaster(segregated_raster, filename = "Output/Spatial Data/Simulated landscapes/segregated_mixedwood_5m.tif", overwrite = TRUE)








# -------------------------------
# Grid of point-count locations
#  - spacing: 2 km
#  - inset:   1 km from the edge
# -------------------------------

# Assumes: `township_sf` already exists (EPSG:3400)
# from your earlier block where you created the 10x10 km polygon.

# 1) Inset polygon by 1 km to keep points 1 km from edge
inset_1km <- sf::st_buffer(township_sf, dist = -1000)

# Safety: if inset produces MULTIPOLYGON/GEOMETRYCOLLECTION, cast to POLYGON
inset_1km <- suppressWarnings(sf::st_cast(inset_1km, "POLYGON"))

# 2) Make a 2 km grid of **centers** over the inset polygon
#    (use what = "centers" so we get point locations directly)
grid_centers <- sf::st_make_grid(
  inset_1km,
  cellsize = 2000,      # 2 km spacing
  what     = "centers",
  square   = TRUE
)

# 3) Keep only points that fall inside the inset polygon
grid_centers <- grid_centers[sf::st_within(grid_centers, sf::st_geometry(inset_1km), sparse = FALSE)[,1]]

# 4) Package as an sf layer your downstream code expects
points_sf <- sf::st_sf(
  id = seq_along(grid_centers),
  geometry = grid_centers
)

# (Optional) quick check
cat("Created", nrow(points_sf), "point-count locations.\n")

# (Optional) save the grid
dir.create("Output/Spatial Data/Simulated landscapes", recursive = TRUE, showWarnings = FALSE)
sf::st_write(points_sf, "Output/Spatial Data/Simulated landscapes/point_grid_2km_inset1km.shp", delete_layer = TRUE)




# =========================
# Overlay points on rasters
# =========================

# Convert to terra vector for fast plotting
points_vect <- terra::vect(points_sf)

# --- On-screen side-by-side (no axes/legend) ---
old_par <- par(no.readonly = TRUE)
par(mfrow = c(1,2), mar = c(0,0,0,0))

# Intimate
plot(intimate_raster, col = c("white","black"), legend = FALSE, axes = FALSE, box = FALSE)
plot(points_vect, add = TRUE, pch = 21, cex = 0.6, col = "red", bg = "red")

# Segregated
plot(segregated_raster, col = c("white","black"), legend = FALSE, axes = FALSE, box = FALSE)
plot(points_vect, add = TRUE, pch = 21, cex = 0.6, col = "red", bg = "red")

par(old_par)

# --- Save clean PNGs with overlays ---
dir.create("Figures/Synthetic landscapes", recursive = TRUE, showWarnings = FALSE)

png("Figures/Synthetic landscapes/intimate_mixedwood_5m_points.png",
    width = 900, height = 900, res = 120)
par(mar = c(0,0,0,0))
plot(intimate_raster, col = c("white","black"), legend = FALSE, axes = FALSE, box = FALSE)
plot(points_vect, add = TRUE, pch = 21, cex = 0.7, col = "red", bg = "red")
dev.off()

png("Figures/Synthetic landscapes/segregated_mixedwood_5m_points.png",
    width = 900, height = 900, res = 120)
par(mar = c(0,0,0,0))
plot(segregated_raster, col = c("white","black"), legend = FALSE, axes = FALSE, box = FALSE)
plot(points_vect, add = TRUE, pch = 21, cex = 0.7, col = "red", bg = "red")
dev.off()






############# Proportion conifer (mean) and CLUMPY at each point & extent

# points_sf must already exist (your 2-km grid inset 1 km)

# Load rasters (0 = deciduous, 1 = conifer)
rasters <- list(
  Intimate   = rast("Output/Spatial Data/Simulated landscapes/intimate_mixedwood_5m.tif"),
  Segregated = rast("Output/Spatial Data/Simulated landscapes/segregated_mixedwood_5m.tif")
)

# Spatial extents in meters
extents <- c(150, 500, 1000)

# Output container
pt_vals <- vector("list", length = length(rasters) * length(extents) * nrow(points_sf))
k <- 0L

for (layout_name in names(rasters)) {
  r <- rasters[[layout_name]]
  
  # make sure values are integer 0/1 for landscapemetrics
  values(r) <- as.integer(round(values(r)))
  
  for (radius_m in extents) {
    message("Processing layout: ", layout_name, "  extent: ", radius_m, " m")
    
    for (i in seq_len(nrow(points_sf))) {
      pt <- points_sf[i, , drop = FALSE]
      
      # circular buffer around point
      buf_sf   <- st_buffer(pt, dist = radius_m)
      buf_vect <- vect(buf_sf)
      
      # crop/mask raster to buffer
      r_crop <- try(crop(r, buf_vect), silent = TRUE)
      if (inherits(r_crop, "try-error")) next
      r_mask <- try(mask(r_crop, buf_vect), silent = TRUE)
      if (inherits(r_mask, "try-error")) next
      
      # ---- proportion conifer (mean of 0/1) ----
      prop_con_val <- try({
        as.numeric(global(r_mask, mean, na.rm = TRUE)[1, 1])
      }, silent = TRUE)
      if (inherits(prop_con_val, "try-error")) prop_con_val <- NA_real_
      
      # ---- CLUMPY for class==1 (conifer) ----
      # use categorical raster for landscapemetrics
      r_cat <- as.factor(r_mask)
      clumpy_val <- try({
        met <- lsm_c_clumpy(r_cat)
        # class-level row for conifer (=1)
        row1 <- met[met$level == "class" & met$class == 1, , drop = FALSE]
        if (nrow(row1) > 0) as.numeric(row1$value[1]) else NA_real_
      }, silent = TRUE)
      if (inherits(clumpy_val, "try-error")) clumpy_val <- NA_real_
      
      # store
      k <- k + 1L
      pt_vals[[k]] <- tibble::tibble(
        x       = st_coordinates(pt)[1],
        y       = st_coordinates(pt)[2],
        layout  = layout_name,
        extent  = radius_m,
        prop_con = prop_con_val,
        clumpy   = clumpy_val
      )
    }
  }
}

# Combine all results
metrics_points_df <- dplyr::bind_rows(pt_vals)

# Optional rule: if CLUMPY is NA (e.g., only one class present), set to 1
metrics_points_df <- metrics_points_df %>%
  mutate(clumpy = if_else(is.na(clumpy), 1, clumpy))

# Save
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)
saveRDS(metrics_points_df, "Output/Tabular Data/propcon_clumpy_by_point_extent.rds")

# Quick peek
metrics_points_df %>% count(layout, extent) %>% print(n = Inf)


# metrics_points_df must already exist (one row per point, per extent, per layout)
# Summarize across points:
extent_means <- metrics_points_df %>%
  group_by(layout, extent) %>%
  summarise(
    mean_prop_con = mean(prop_con, na.rm = TRUE),
    mean_clumpy   = mean(clumpy,   na.rm = TRUE),
    sd_prop_con   = sd(prop_con,   na.rm = TRUE),
    sd_clumpy     = sd(clumpy,     na.rm = TRUE),
    n_points      = n(),
    .groups = "drop"
  ) %>%
  arrange(layout, extent)

print(extent_means, n = Inf)

# (optional) save
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)
write_csv(extent_means, "Output/Tabular Data/mean_propcon_clumpy_by_layout_extent.csv")













###################### FOR FIGURE CREATION ONLY

########## Extract a 2 by 2 km cell from township landscape for figure


# Define township center coordinates and size
township_center <- c(x = 528551, y = 6000000)
township_size <- 10000
half_size <- township_size / 2
cell_size <- 10
target_crs <- "EPSG:3400"

# Derived grid size in cells
ncol <- nrow <- township_size / cell_size

# Matrix to SpatRaster
matrix_to_raster <- function(mat, crs_string, center_x, center_y, total_size) {
  half_total_size <- total_size / 2
  rast_obj <- rast(
    nrows = nrow(mat), ncols = ncol(mat),
    xmin = center_x - half_total_size,
    xmax = center_x + half_total_size,
    ymin = center_y - half_total_size,
    ymax = center_y + half_total_size,
    crs = crs_string
  )
  values(rast_obj) <- as.vector(t(mat))
  return(rast_obj)
}

# Intimate = random 0/1
generate_intimate <- function() {
  matrix(sample(c(0, 1), size = nrow * ncol, replace = TRUE), nrow = nrow, ncol = ncol)
}

# Segregated half-half blocks
generate_half_half_block <- function(size) {
  block <- matrix(0, nrow = size, ncol = size)
  side <- sample(c("top", "bottom", "left", "right"), 1)
  if (side %in% c("top", "bottom")) {
    rows <- floor(size / 2)
    if (side == "top") {
      block[1:rows, ] <- 1
    } else {
      block[(size - rows + 1):size, ] <- 1
    }
  } else {
    cols <- floor(size / 2)
    if (side == "left") {
      block[, 1:cols] <- 1
    } else {
      block[, (size - cols + 1):size] <- 1
    }
  }
  return(block)
}

generate_segregated <- function() {
  full_mat <- matrix(0, nrow = nrow, ncol = ncol)
  block_size_pixels <- 500 / cell_size
  blocks_per_side <- ceiling(nrow / block_size_pixels)
  
  for (i in 0:(blocks_per_side - 1)) {
    for (j in 0:(blocks_per_side - 1)) {
      block <- generate_half_half_block(block_size_pixels)
      r1 <- i * block_size_pixels + 1
      r2 <- min((i + 1) * block_size_pixels, nrow)
      c1 <- j * block_size_pixels + 1
      c2 <- min((j + 1) * block_size_pixels, ncol)
      full_mat[r1:r2, c1:c2] <- block[1:(r2 - r1 + 1), 1:(c2 - c1 + 1)]
    }
  }
  return(full_mat)
}

# Generate and convert to raster
intimate_mat <- generate_intimate()
segregated_mat <- generate_segregated()

intimate_raster <- matrix_to_raster(intimate_mat, target_crs, township_center["x"], township_center["y"], township_size)
segregated_raster <- matrix_to_raster(segregated_mat, target_crs, township_center["x"], township_center["y"], township_size)

# Define bottom-right corner extent (2.5 km × 2.5 km block)
xmin_crop <- township_center["x"] + half_size - 2500
xmax_crop <- township_center["x"] + half_size
ymin_crop <- township_center["y"] - half_size
ymax_crop <- township_center["y"] - half_size + 2500
crop_extent <- ext(xmin_crop, xmax_crop, ymin_crop, ymax_crop)

# Crop both rasters
intimate_crop <- crop(intimate_raster, crop_extent)
segregated_crop <- crop(segregated_raster, crop_extent)

# ---------------------------------------
# Create 500 m grid from cropped raster extent
# ---------------------------------------
make_grid <- function(r, cellsize = 500, crs_string = "EPSG:3400") {
  e <- ext(r)
  
  bbox_poly <- st_as_sfc(
    st_bbox(
      c(
        xmin = xmin(r),
        ymin = ymin(r),
        xmax = xmax(r),
        ymax = ymax(r)
      ),
      crs = st_crs(crs_string)
    )
  )
  
  g <- st_make_grid(
    bbox_poly,
    cellsize = cellsize,
    square = TRUE,
    what = "polygons"
  )
  
  st_sf(geometry = g)
}

grid_overlay <- make_grid(intimate_crop, cellsize = 500)

# ---------------------------------------
# Save TIFFs WITH grid
# ---------------------------------------
dir.create("Output/Figures/Synthetic landscapes", recursive = TRUE, showWarnings = FALSE)

bw_palette <- c("white", "black")

tiff("Output/Figures/Synthetic landscapes/intimate_block_2.5km_grid.tiff",
     width = 800, height = 800, res = 600, compression = "lzw")

par(mar = c(0, 0, 0, 0))
plot(intimate_crop, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
plot(st_geometry(grid_overlay), add = TRUE, border = "black", lwd = 1)
dev.off()

tiff("Output/Figures/Synthetic landscapes/segregated_block_2.5km_grid.tiff",
     width = 800, height = 800, res = 600, compression = "lzw")

par(mar = c(0, 0, 0, 0))
plot(segregated_crop, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
plot(st_geometry(grid_overlay), add = TRUE, border = "black", lwd = 1)
dev.off()

# ---------------------------------------
# Optional: also save cropped rasters without grid
# ---------------------------------------
tiff("Output/Figures/Synthetic landscapes/intimate_block_2.5km.tiff",
     width = 800, height = 800, res = 600, compression = "lzw")

par(mar = c(0, 0, 0, 0))
plot(intimate_crop, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
dev.off()

tiff("Output/Figures/Synthetic landscapes/segregated_block_2.5km.tiff",
     width = 800, height = 800, res = 600, compression = "lzw")

par(mar = c(0, 0, 0, 0))
plot(segregated_crop, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
dev.off()

# Save cropped rasters
dir.create("Output/Spatial Data/Simulated landscapes", recursive = TRUE, showWarnings = FALSE)
writeRaster(intimate_crop, filename = "Output/Spatial Data/Simulated landscapes/intimate_block_2.5km.tif", overwrite = TRUE)
writeRaster(segregated_crop, filename = "Output/Spatial Data/Simulated landscapes/segregated_block_2.5km.tif", overwrite = TRUE)

# CLUMPY
print(lsm_c_clumpy(intimate_crop))
print(lsm_c_clumpy(segregated_crop))









################# Add 500 m x 500 m grid overlay for figure

# ---------------------------------------
# Create 500 m x 500 m grid for township
# ---------------------------------------
grid_500m <- st_make_grid(
  township_sf,
  cellsize = 500,
  what = "polygons",
  square = TRUE
)

grid_500m_sf <- st_sf(geometry = grid_500m)

# Clip grid to township boundary
grid_500m_sf <- st_intersection(grid_500m_sf, township_sf)

# Save grid shapefile
dir.create("Output/Spatial Data/Township polygon", recursive = TRUE, showWarnings = FALSE)
st_write(
  grid_500m_sf,
  "Output/Spatial Data/Township polygon/township_grid_500m.shp",
  delete_layer = TRUE
)

# ---------------------------------------
# On-screen plot with grid overlay
# ---------------------------------------
old_par <- par(no.readonly = TRUE)
par(mar = c(0, 0, 0, 0))

plot(intimate_raster, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
plot(st_geometry(grid_500m_sf), add = TRUE, border = "black", lwd = 2)

plot(segregated_raster, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
plot(st_geometry(grid_500m_sf), add = TRUE, border = "black", lwd = 2)

par(old_par)

# ---------------------------------------
# Save clean TIFF plots WITH grid overlay
# ---------------------------------------
dir.create("Output/Figures/Synthetic landscapes", recursive = TRUE, showWarnings = FALSE)

tiff("Output/Figures/Synthetic landscapes/intimate_mixedwood_5m_grid.tiff",
     width = 800, height = 800, res = 600, compression = "lzw")

par(mar = c(0, 0, 0, 0))
plot(intimate_raster, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
plot(st_geometry(grid_500m_sf), add = TRUE, border = "black", lwd = 2)
dev.off()

tiff("Output/Figures/Synthetic landscapes/segregated_mixedwood_5m_grid.tiff",
     width = 800, height = 800, res = 600, compression = "lzw")

par(mar = c(0, 0, 0, 0))
plot(segregated_raster, col = bw_palette, legend = FALSE, axes = FALSE, main = "", box = FALSE)
plot(st_geometry(grid_500m_sf), add = TRUE, border = "black", lwd = 2)
dev.off()