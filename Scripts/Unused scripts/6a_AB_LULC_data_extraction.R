# ---
# title: "LULC_data_extraction"
# author: "Leonard Patterson"
# created: "2025-04-11"
# description: "This script extracts habitat covariates from the LULC raster, 
# which were downloaded from here: https://ags.aer.ca/publications/all-publications/dig-2021-0019
# ---


# Load libraries

library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(landscapemetrics)  
library(raster)



# Load data

#LULC <-rast("Input/Spatial Data/LULC/AB_LULC_2020_10tm.tif")

### Reclassify LULC raster to conifer

# Save conifer reclassification scheme in matrix
#rcl_matrix_con <- matrix(
#  c(
#    1, 0,  # Water
#    2, 0,  # Bryoids
#    3, 0,  # Wetland – Treed
#    4, 0,  # Herbs
#    5, 0,  # Exposed/Barren Land
#    6, 0,  # Shrubland
#    7, 0,  # Wetland
#    8, 0,  # Grassland
#    9, 1,  # Coniferous
#    10, 0, # Broadleaf
#    11, 1, # Mixedwood
#    12, 0,  # Agriculture
#    13, 0,  # Developed
#    14, 0   # Burned Areas
#  ),
#  ncol = 2,
#  byrow = TRUE
#)

# Reclassify to conifer
#LULC_reclassified_con <- classify(LULC, rcl_matrix_con)

# Save the reclassified rasters
#writeRaster(LULC_reclassified_con, "Output/Spatial Data/Reprojected LULC/Conifer_AB_LULC_2020_10tm_reclass.tif", overwrite = TRUE)








############# Calculate proportion conifer 

# Load point count location from NTEMS
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/point_count_locs_NTEMS.shp")

# Remove unnecessary column
AB_point_counts_sf <- AB_point_counts_sf %>%
  dplyr::select(PointID, location, year, lat, lon, x_AEP10TM, y_AEP10TM, geometry)

# Load rasters
Conifer_10TM <- rast("Output/Spatial Data/Reprojected LULC/Conifer_AB_LULC_2020_10tm_reclass.tif")


### Calculate proportion of conifer and deciduous

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
    prop_con_3 = prop_1000_df[, "conifer_prop"]
  )



###### Calculate CLUMPY at class level

# Create filtered buffers (important to recreate buffers based on filtered points)
buffer150 <- st_buffer(AB_point_counts_sf, dist = 150)
buffer500 <- st_buffer(AB_point_counts_sf, dist = 500)
buffer1000 <- st_buffer(AB_point_counts_sf, dist = 1000)

# Ensure buffers match the raster CRS
buffer150 <- st_transform(buffer150, crs(Conifer_10TM))
buffer500 <- st_transform(buffer500, crs(Conifer_10TM))
buffer1000 <- st_transform(buffer1000, crs(Conifer_10TM))

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
    
    conifer_crop <- terra::crop(Conifer_10TM, buffer_area)
    
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
clumpy_150_vec <- calculate_clumpy(buffer150)
clumpy_500_vec <- calculate_clumpy(buffer500)
clumpy_1000_vec <- calculate_clumpy(buffer1000)

# Add results to the filtered sf object
AB_point_counts_sf <- AB_point_counts_sf %>%
  mutate(
    clumpy_1 = clumpy_150_vec,
    clumpy_2 = clumpy_500_vec,
    clumpy_3 = clumpy_1000_vec
  )

# Save as shp file
st_write(AB_point_counts_sf, "Output/Spatial Data/AB_point_counts/point_count_locs_AB_LULC.shp", append = FALSE) 

# Save as CSV
AB_point_counts_df <- st_drop_geometry(AB_point_counts_sf)
write.csv(AB_point_counts_df, "Output/Tabular Data/point_counts_AB_LULC.csv", row.names = FALSE)
