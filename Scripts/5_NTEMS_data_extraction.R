# ---
# title: "NTEMS_data_extraction"
# author: "Leonard Patterson"
# created: "2025-04-11"
# description: "This script extracts habitat covariates from the NTEMS rasters, 
# which were downloaded from these websites:
# Tree species raster: https://gee-community-catalog.org/projects/ca_lc/#class-schema"
# Forest age: https://developers.google.com/earth-engine/datasets/catalog/CANADA_NFIS_NTEMS_CA_FOREST_AGE
# ---

# Remove objects from environment
rm(list = ls())

# Load libraries

library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(landscapemetrics)  
library(raster)


### Reproject rasters

## Import rasters
#Conifer <- rast("Input/Spatial Data/NTEMS Tree Rasters/Conifer.tif")
#Broadleaf <- rast("Input/Spatial Data/NTEMS Tree Rasters/Broadleaf.tif")
#Age <- rast("Input/Spatial Data/NTEMS Tree Rasters/ForestAge.tif")

## Load point count data
#AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/filtered_point_count_loc_2009.shp")

# Extract the WKT string from your desired CRS
#target_crs_wkt <- st_as_text(st_crs(AB_point_counts_sf))

## Create the target CRS object using terra::crs()
#target_crs_terra <- terra::crs(target_crs_wkt)

## Reproject rasters
#Conifer_10TM <- terra::project(Conifer, target_crs_terra, method = "near")
#Age_10TM <- terra::project(Age, target_crs_terra, method = "near")


## Save reprojected rasters

## Define output directory and ensure it exists
#output_dir <- "Output/Spatial Data/Reprojected Tree Rasters"

## List of rasters and corresponding filenames
#rasters <- list(
#  Conifer_reproj = Conifer_10TM,
#  Age_reproj = Age_10TM
#)

# Loop to save rasters
#for (name in names(rasters)) {
#  writeRaster(rasters[[name]], file.path(output_dir, paste0(name, ".tif")), overwrite = TRUE)
#}

## Print success message
#cat("Reprojected rasters saved to:", output_dir, "\n")










############# Calculate proportion conifer

# Load point count data
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/filtered_point_counts.shp")

# Load rasters
Conifer_10TM <- rast("Output/Spatial Data/Reprojected Tree Rasters/Conifer_reproj.tif")

# Create a 1000 meter buffer around each point
buffer1000 <- st_buffer(AB_point_counts_sf, dist = 1000)
buffer1000 <- st_transform(buffer1000, crs(Conifer_10TM))

### First, identify and remove points where 1000 m buffer extends outside of AB boundary

# Load AB boundary
AB_boundary_line <- st_read("Input/Spatial Data/Alberta/AB_boundary_line.shp")

# Ensure both are in the same coordinate reference system (CRS)
boundary_line <- st_transform(AB_boundary_line, st_crs(buffer1000))

# Check for intersections between the boundary line and the buffer polygons
intersecting_polygons <- st_intersection(buffer1000, AB_boundary_line)

# Save the result in a dataframe called points_to_remove_bndry
points_to_remove_bndry <- as.data.frame(intersecting_polygons)

# Remove intersecting polygons from AB_point_count_sf
AB_point_counts_sf <- AB_point_counts_sf[!AB_point_counts_sf$location %in% intersecting_polygons$location, ]



### Now, calculate proportion of conifer, spruce, and deciduous

# Create 150, 500, and 1000 meter buffers around each point
buffer150 <- st_buffer(AB_point_counts_sf, dist = 150)
buffer500 <- st_buffer(AB_point_counts_sf, dist = 500)
buffer1000 <- st_buffer(AB_point_counts_sf, dist = 1000)

# Ensure buffers match the raster CRS
buffer150 <- st_transform(buffer150, crs(Conifer_10TM))
buffer500 <- st_transform(buffer500, crs(Conifer_10TM))
buffer1000 <- st_transform(buffer1000, crs(Conifer_10TM))


# Function to calculate proportions
calculate_proportions <- function(buffer) {
  lapply(1:nrow(buffer), function(i) {
    buffer_area <- buffer[i, ]
    
    # Crop rasters to buffer area
    conifer_crop <- terra::crop(Conifer_10TM, buffer_area)
    
    # Calculate the number of cells occupied by each species within the buffer
    conifer_cells <- terra::global(conifer_crop, "sum", na.rm = TRUE) %>% as.numeric()
    
    # Calculate total number of cells in the cropped raster (buffer area)
    total_cells <- terra::ncell(conifer_crop)
    
    # Calculate proportions
    conifer_prop <- ifelse(total_cells > 0, conifer_cells / total_cells, 0)
    
    return(c(conifer_prop = conifer_prop))
  })
}

# Calculate proportions for each buffer size
prop_150 <- calculate_proportions(buffer150)
prop_500 <- calculate_proportions(buffer500)
prop_1000 <- calculate_proportions(buffer1000)

# Convert lists to data frames
prop_150_df <- do.call(rbind, prop_150)
prop_500_df <- do.call(rbind, prop_500)
prop_1000_df <- do.call(rbind, prop_1000)

# Add proportions to the sf object with correct names
AB_point_counts_sf <- AB_point_counts_sf %>%
  mutate(
    prop_con_1 = prop_150_df[, "conifer_prop"],
    prop_con_2 = prop_500_df[, "conifer_prop"],
    prop_con_3 = prop_1000_df[, "conifer_prop"],
  )

# Add unique ID column

AB_point_counts_sf$PointID <- seq_len(nrow(AB_point_counts_sf))

# Reorder the columns to place the ID column at the left
AB_point_counts_sf <- AB_point_counts_sf[, c("PointID", setdiff(names(AB_point_counts_sf), "PointID"))]








###### Calculate CLUMPY at class level

# Filter AB_point_counts_sf based on conifer presence in 150m buffer
AB_point_counts_sf_filtered <- AB_point_counts_sf %>%
  filter(prop_con_1 > 0)

# Load cropped conifer raster (cropped to remove NA values outside of AB boundary
# - done in ArcPro because R kept crashing when trying to mask raster)

Conifer_10TM_crop <- rast("Output/Spatial Data/Reprojected Tree Rasters/Conifer_10TM_crop.tif")

# Create filtered buffers (important to recreate buffers based on filtered points)
buffer150_filtered <- st_buffer(AB_point_counts_sf_filtered, dist = 150)
buffer500_filtered <- st_buffer(AB_point_counts_sf_filtered, dist = 500)
buffer1000_filtered <- st_buffer(AB_point_counts_sf_filtered, dist = 1000)

# Ensure buffers match the raster CRS
buffer150_filtered <- st_transform(buffer150_filtered, crs(Conifer_10TM_crop))
buffer500_filtered <- st_transform(buffer500_filtered, crs(Conifer_10TM_crop))
buffer1000_filtered <- st_transform(buffer1000_filtered, crs(Conifer_10TM_crop))

# CLUMPY function with proper categorical raster handling
calculate_clumpy <- function(buffer) {
  clumpy_values <- numeric(nrow(buffer))
  
  for (i in 1:nrow(buffer)) {
    buffer_area <- buffer[i, , drop = FALSE]  # Ensure it's still an sf object
    
    if (st_is_empty(buffer_area)) {
      message("Empty buffer for point ", i)
      clumpy_values[i] <- NA
      next
    }
    
    conifer_crop <- terra::crop(Conifer_10TM_crop, buffer_area)
    
    if (terra::ncell(conifer_crop) == 0) {
      message("No raster cells found for point ", i)
      clumpy_values[i] <- NA
      next
    }
    
    # Convert raster values to categorical (binary: 0 = non-conifer, 1 = conifer)
    reclass_matrix <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
    conifer_crop_factor <- terra::classify(conifer_crop, reclass_matrix)
    
    # Ensure only 0 and 1 values exist in raster
    unique_values <- unique(values(conifer_crop_factor))
    if (!all(unique_values %in% c(0, 1, NA))) {
      warning("Non-binary values found for point ", i)
      clumpy_values[i] <- NA
      next
    }
    
    # Ensure SpatRaster is valid and non-empty
    if (terra::ncell(conifer_crop_factor) == 0 || all(is.na(values(conifer_crop_factor)))) {
      message("Empty or invalid raster for point ", i)
      clumpy_values[i] <- NA
      next
    }
    
    # Compute CLUMPY metric
    clumpy_result <- tryCatch({
      lsm_c_clumpy(conifer_crop_factor)
    }, error = function(e) {
      warning("Error calculating CLUMPY for point ", i, ": ", e$message)
      return(data.frame(layer = 1, class = 1, value = NA))
    })
    
    # Handle results
    if (nrow(clumpy_result) == 0) {
      warning("No conifers found for point ", i)
      clumpy_values[i] <- 0
    } else {
      clumpy_values[i] <- clumpy_result$value[1]
    }
  }
  return(clumpy_values)
}

# Run CLUMPY calculations for filtered buffers
clumpy_150_vec <- calculate_clumpy(buffer150_filtered)
clumpy_500_vec <- calculate_clumpy(buffer500_filtered)
clumpy_1000_vec <- calculate_clumpy(buffer1000_filtered)

# Add results to the filtered sf object
AB_point_counts_sf_filtered <- AB_point_counts_sf_filtered %>%
  mutate(
    clumpy_1 = clumpy_150_vec,
    clumpy_2 = clumpy_500_vec,
    clumpy_3 = clumpy_1000_vec
  )










#### Extract forest age

# Load age raster
Age_10TM <- rast("Output/Spatial Data/Reprojected Tree Rasters/Age_reproj.tif")

# Function to calculate mean and median age
calculate_age_stats <- function(buffer) {
  lapply(1:nrow(buffer), function(i) {
    buffer_area <- buffer[i, ]
    age_crop <- terra::crop(Age_10TM, buffer_area)
    
    # Check for empty crop (early exit if empty)
    if (terra::ncell(age_crop) == 0) {
      return(c(mean_age = NA, median_age = NA))
    }
    
    # Calculate mean and median (using more direct approach)
    mean_age <- mean(values(age_crop), na.rm = TRUE) # Direct calculation
    median_age <- median(values(age_crop), na.rm = TRUE) # Direct calculation
    
    # Handle potential NULL or NaN (though less likely now)
    if (is.null(mean_age) || is.nan(mean_age)) mean_age <- NA
    if (is.null(median_age) || is.nan(median_age)) median_age <- NA
    
    return(c(mean_age = mean_age, median_age = median_age))
  })
}

# Calculate age statistics for each buffer size
age_stats_150 <- calculate_age_stats(buffer150_filtered)
age_stats_500 <- calculate_age_stats(buffer500_filtered)
age_stats_1000 <- calculate_age_stats(buffer1000_filtered)

# Convert lists to data frames
age_stats_150_df <- do.call(rbind, age_stats_150)
age_stats_500_df <- do.call(rbind, age_stats_500)
age_stats_1000_df <- do.call(rbind, age_stats_1000)

# Add age statistics to the sf object (UPDATED NAMING CONVENTION)
AB_point_counts_sf_filtered <- AB_point_counts_sf_filtered %>%
  mutate(
    age_mn_1 = age_stats_150_df[, "mean_age"],  # Mean
    age_md_1 = age_stats_150_df[, "median_age"], # Median
    age_mn_2 = age_stats_500_df[, "mean_age"],
    age_md_2 = age_stats_500_df[, "median_age"],
    age_mn_3 = age_stats_1000_df[, "mean_age"],
    age_md_3 = age_stats_1000_df[, "median_age"]
  )

# Save as shp file
st_write(AB_point_counts_sf_filtered, "Output/Spatial Data/AB_point_counts/point_count_locs_NTEMS.shp", append=FALSE) 

# Save as CSV
AB_point_counts_df_filtered <- st_drop_geometry(AB_point_counts_sf_filtered)
write.csv(AB_point_counts_df_filtered, "Output/Tabular Data/point_counts_NTEMS.csv", row.names = FALSE)

