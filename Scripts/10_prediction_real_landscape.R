# Load libraries
library(terra)
library(sf)
library(landscapemetrics)
library(ggplot2)

# Load prediction area shp
pred_area <- st_read("Input/Spatial Data/Prediction area/prediction_area.shp")

# Load conifer raster
conifer_LULC <- rast("Output/Spatial Data/reproj_AB_LULC/Conifer_AB_LULC_2020_10tm_reclass.tif")


####### Note: insufficient RAM to run the code below which clips the conifer_LULC raster to the extent
####### of the prediction_area. Had to run in GIS.
# --- 3. Clip conifer_LULC to the extent of pred_area ---
# First, crop the raster to the bounding box of the polygon for efficiency
#clipped_raster_cropped <- crop(conifer_LULC, vect(pred_area))

# Then, mask the cropped raster using the actual polygon shape
#clipped_raster <- mask(clipped_raster_cropped, vect(pred_area))


#print(clipped_raster_cropped)
# Verify the clipped raster
# plot(clipped_raster, main = "Clipped Conifer LULC")
# plot(st_geometry(pred_area), add = TRUE, border = "red", lwd = 2)


# Import clipped raster (clipped in GIS due to high RAM requirements that made R crash)
clipped_raster <- rast("Output/Spatial Data/reproj_LULC_clipped/clipped_LULC.tif")




########### Calculate CLUMPY for townships w/in clipped area

library(terra)
library(raster)
library(landscapemetrics)
library(dplyr)

# Load raster
clipped_raster <- rast("Output/Spatial Data/reproj_LULC_clipped/clipped_LULC.tif")
r_raster <- raster(clipped_raster)

# Ensure binary (0/1) and factor
r_raster <- round(r_raster)
r_raster <- ratify(r_raster)
levels(r_raster)[[1]] <- data.frame(ID = c(0, 1), class = c("Non-Conifer", "Conifer"))

# Define 10 km² window in cells
window_area_m2 <- 10 * 1000 * 1000
cell_size <- res(r_raster)[1]
window_side_cells <- round(sqrt(window_area_m2 / cell_size^2))
if (window_side_cells %% 2 == 0) window_side_cells <- window_side_cells + 1

# Define step (e.g., non-overlapping windows)
step <- window_side_cells

# Get number of rows/cols
nrows <- nrow(r_raster)
ncols <- ncol(r_raster)

# Initialize result list
results <- list()
i <- 1

# Loop over the raster in blocks
for (row in seq(1, nrows, by = step)) {
  for (col in seq(1, ncols, by = step)) {
    
    # Get extent of the window
    r1 <- max(1, row - floor(window_side_cells/2))
    r2 <- min(nrows, row + floor(window_side_cells/2))
    c1 <- max(1, col - floor(window_side_cells/2))
    c2 <- min(ncols, col + floor(window_side_cells/2))
    
    # Get the extent in spatial coordinates
    e <- extent(r_raster, r1, r2, c1, c2)
    r_crop <- crop(r_raster, e)
    
    # Check for sufficient data
    if (cellStats(!is.na(r_crop), sum) < 0.5 * window_side_cells^2) next
    
    # Calculate CLUMPY
    lsm <- tryCatch({
      calculate_lsm(r_crop, level = "class", metric = "clumpy", classes = 1)
    }, error = function(e) NULL)
    
    if (!is.null(lsm) && nrow(lsm) > 0) {
      center_xy <- xyFromCell(r_raster, cellFromRowCol(r_raster, row, col))
      # Count conifer cells (value == 1)
      conifer_count <- cellStats(r_crop == 1, sum)
      total_count   <- cellStats(!is.na(r_crop), sum)
      prop_conifer  <- ifelse(total_count > 0, conifer_count / total_count, NA)
      
      # Store both clumpy and proportion conifer
      results[[i]] <- data.frame(
        x = center_xy[1],
        y = center_xy[2],
        clumpy = lsm$value,
        prop_con = prop_conifer
      )
      i <- i + 1
    }
  }
}

# Combine into df
clumpy_twnshp_values <- do.call(rbind, results)

# View summary
summary(clumpy_twnshp_values$clumpy)









############ Plot townships 



library(terra)
library(landscapemetrics)
library(ggplot2)
library(sf)
library(dplyr)

# Set 10 km² window size (assuming square window)
window_side_m <- sqrt(10000000)  # 10 km² = 10,000,000 m² → ~3162.28 m

# Function to extract and plot focal window around township centroid
plot_clumpy_window <- function(r, center_xy, window_size, title) {
  # Create square extent around centroid
  half_side <- window_size / 2
  e <- ext(
    center_xy[1] - half_side, center_xy[1] + half_side,
    center_xy[2] - half_side, center_xy[2] + half_side
  )
  
  # Crop raster
  focal_r <- crop(r, e)
  
  # Convert to data.frame for ggplot
  df <- as.data.frame(focal_r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "landcover"
  
  ggplot(df, aes(x = x, y = y, fill = factor(landcover))) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_d(na.value = "white") +
    labs(title = title, fill = "Landcover") +
    theme_minimal()
}

# ---- Get centroids from selected rows in clumpy_twnshp_values ----

# If clumpy_twnshp_values is a data.frame with x and y columns:
centroid_1 <- c(clumpy_twnshp_values$x[3297], clumpy_twnshp_values$y[3297])
centroid_2 <- c(clumpy_twnshp_values$x[13594], clumpy_twnshp_values$y[13594])

# Plot both
p1 <- plot_clumpy_window(clipped_raster, centroid_1, window_side_m, "Township disp")
p2 <- plot_clumpy_window(clipped_raster, centroid_2, window_side_m, "Township agg")

# Plot 
p1 + p2

