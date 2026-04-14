# ---
# title: "NTEMS_data_extraction, extracting from pine and black spruce rasters"
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










############# Calculate proportion black spruce

# Load point count data
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/filtered_point_counts.shp")

# Load rasters
BS_10TM <- rast("Output/Spatial Data/Reprojected Tree Rasters/BlackSpruce_reproj.tif")

# Create a 1000 meter buffer around each point
buffer1000 <- st_buffer(AB_point_counts_sf, dist = 1000)
buffer1000 <- st_transform(buffer1000, crs(BS_10TM))

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
buffer150 <- st_transform(buffer150, crs(BS_10TM))
buffer500 <- st_transform(buffer500, crs(BS_10TM))
buffer1000 <- st_transform(buffer1000, crs(BS_10TM))


# Function to calculate proportions
calculate_proportions <- function(buffer) {
  lapply(1:nrow(buffer), function(i) {
    buffer_area <- buffer[i, ]
    
    # Crop rasters to buffer area
    bs_crop <- terra::crop(BS_10TM, buffer_area)
    
    # Calculate the number of cells occupied by each species within the buffer
    bs_cells <- terra::global(bs_crop, "sum", na.rm = TRUE) %>% as.numeric()
    
    # Calculate total number of cells in the cropped raster (buffer area)
    total_cells <- terra::ncell(bs_crop)
    
    # Calculate proportions
    bs_prop <- ifelse(total_cells > 0, bs_cells / total_cells, 0)
    
    return(c(bs_prop = bs_prop))
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
AB_point_counts_sf_bs <- AB_point_counts_sf %>%
  mutate(
    prop_bs_1 = prop_150_df[, "bs_prop"],
    prop_bs_2 = prop_500_df[, "bs_prop"],
    prop_bs_3 = prop_1000_df[, "bs_prop"],
  )

# Add unique ID column
AB_point_counts_sf_bs$PointID <- seq_len(nrow(AB_point_counts_sf_bs))

# Reorder the columns to place the ID column at the left
AB_point_counts_sf_bs <- AB_point_counts_sf_bs[, c("PointID", setdiff(names(AB_point_counts_sf_bs), "PointID"))]


# Save as CSV
AB_point_counts_df_bs <- st_drop_geometry(AB_point_counts_sf_bs)
write.csv(AB_point_counts_df_bs, "Output/Tabular Data/point_counts_NTEMS_black_spruce.csv", row.names = FALSE)
















############# Calculate proportion PINE 

# Load point count data
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/filtered_point_counts.shp")

# Load rasters
PINE_10TM <- rast("Output/Spatial Data/Reprojected Tree Rasters/BlackSpruce_reproj.tif")

# Create a 1000 meter buffer around each point
buffer1000 <- st_buffer(AB_point_counts_sf, dist = 1000)
buffer1000 <- st_transform(buffer1000, crs(PINE_10TM))

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
buffer150 <- st_transform(buffer150, crs(PINE_10TM))
buffer500 <- st_transform(buffer500, crs(PINE_10TM))
buffer1000 <- st_transform(buffer1000, crs(PINE_10TM))


# Function to calculate proportions
calculate_proportions <- function(buffer) {
  lapply(1:nrow(buffer), function(i) {
    buffer_area <- buffer[i, ]
    
    # Crop rasters to buffer area
    pine_crop <- terra::crop(PINE_10TM, buffer_area)
    
    # Calculate the number of cells occupied by each species within the buffer
    pine_cells <- terra::global(pine_crop, "sum", na.rm = TRUE) %>% as.numeric()
    
    # Calculate total number of cells in the cropped raster (buffer area)
    total_cells <- terra::ncell(pine_crop)
    
    # Calculate proportions
    pine_prop <- ifelse(total_cells > 0, pine_cells / total_cells, 0)
    
    return(c(pine_prop = pine_prop))
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
AB_point_counts_sf_pine <- AB_point_counts_sf %>%
  mutate(
    prop_pine_1 = prop_150_df[, "pine_prop"],
    prop_pine_2 = prop_500_df[, "pine_prop"],
    prop_pine_3 = prop_1000_df[, "pine_prop"],
  )

# Add unique ID column
AB_point_counts_sf_pine$PointID <- seq_len(nrow(AB_point_counts_sf_pine))

# Reorder the columns to place the ID column at the left
AB_point_counts_sf_pine <- AB_point_counts_sf_pine[, c("PointID", setdiff(names(AB_point_counts_sf_pine), "PointID"))]



# Save as CSV
AB_point_counts_df_pine <- st_drop_geometry(AB_point_counts_sf_pine)
write.csv(AB_point_counts_df_pine, "Output/Tabular Data/point_counts_NTEMS_pine.csv", row.names = FALSE)

