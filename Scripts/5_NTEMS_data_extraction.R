# ---
# title: "NTEMS_data_extraction"
# author: "Leonard Patterson"
# created: ""
# description: "."
# ---


# Load libraries

library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(landscapemetrics)  
library(raster)


### Reproject rasters

## Import rasters
#Broadleaf150 <- rast("Input/Spatial Data/NTEMS Tree Rasters/Broadleaf.H150.tif")
#Conifer150 <- rast("Input/Spatial Data/NTEMS Tree Rasters/Conifer.H150.tif")
#Conifer <- rast("Input/Spatial Data/NTEMS Tree Rasters/Conifer.tif")
#Picegla <- rast("Input/Spatial Data/NTEMS Tree Rasters/Conifer.tif")
#Broadleaf <- rast("Input/Spatial Data/NTEMS Tree Rasters/Broadleaf.tif")
#Age <- rast("Input/Spatial Data/NTEMS Tree Rasters/ForestAge.tif")

## Load point count data
#AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/filtered_point_count_loc_2009.shp")

## Extract the WKT string from your desired CRS
#target_crs_wkt <- st_as_text(st_crs(AB_point_counts_sf))

## Create the target CRS object using terra::crs()
#target_crs_terra <- terra::crs(target_crs_wkt)

## Reproject rasters
#Conifer_10TM <- terra::project(Conifer, target_crs_terra, method = "near")
#Broadleaf_10TM <- terra::project(Broadleaf, target_crs_terra, method = "near")
#Picegla_10TM <- terra::project(Picegla, target_crs_terra, method = "near")
#Age_10TM <- terra::project(Age, target_crs_terra, method = "near")
#Broadleaf150_10TM <- terra::project(Broadleaf150, target_crs_terra, method = "bilinear")
#Conifer150_10TM <- terra::project(Conifer150, target_crs_terra, method = "bilinear")

## Save reprojected rasters

## Define output directory and ensure it exists
output_dir <- "Output/Spatial Data/reproj_tree_rasters"

## List of rasters and corresponding filenames
#rasters <- list(
  #Broadleaf150_reproj = Broadleaf150_10TM,
  #Conifer150_reproj = Conifer150_10TM,
  #Conifer_reproj = Conifer_10TM,
  #Broadleaf_reproj = Broadleaf_10TM,
  #Age_reproj = Age_10TM
#)

# Loop to save rasters
#for (name in names(rasters)) {
#  writeRaster(rasters[[name]], file.path(output_dir, paste0(name, ".tif")), overwrite = TRUE)
#}

#writeRaster(Picegla_10TM, "Output/Spatial Data/reproj_tree_rasters/Picegla_reproj.tif", overwrite = TRUE)

## Print success message
#cat("Reprojected rasters saved to:", output_dir, "\n")











############# Calculate proportion conifer, white spruce, and broadleaf

# Load point count data
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/filtered_point_count_loc_2009.shp")

# Load rasters
Broadleaf_10TM <- rast("Output/Spatial Data/reproj_tree_rasters/Broadleaf_reproj.tif")
Conifer_10TM <- rast("Output/Spatial Data/reproj_tree_rasters/Conifer_reproj.tif")
Picegla_10TM <- rast("Output/Spatial Data/reproj_tree_rasters/Picegla_reproj.tif")


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
    broadleaf_crop <- terra::crop(Broadleaf_10TM, buffer_area)
    picegla_crop <- terra::crop(Picegla_10TM, buffer_area)
    
    
    # Calculate the number of cells occupied by each species within the buffer
    conifer_cells <- terra::global(conifer_crop, "sum", na.rm = TRUE) %>% as.numeric()
    broadleaf_cells <- terra::global(broadleaf_crop, "sum", na.rm = TRUE) %>% as.numeric()
    picegla_cells <- terra::global(picegla_crop, "sum", na.rm = TRUE) %>% as.numeric()
    
    
    # Calculate total number of cells in the cropped raster (buffer area)
    total_cells <- terra::ncell(conifer_crop)
    
    # Calculate proportions
    conifer_prop <- ifelse(total_cells > 0, conifer_cells / total_cells, 0)
    broadleaf_prop <- ifelse(total_cells > 0, broadleaf_cells / total_cells, 0)
    picegla_prop <- ifelse(total_cells > 0, picegla_cells / total_cells, 0)
    
    
    return(c(conifer_prop = conifer_prop, broadleaf_prop = broadleaf_prop, picegla_prop = picegla_prop))
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
    prop_dec_1 = prop_150_df[, "broadleaf_prop"],
    prop_con_2 = prop_500_df[, "conifer_prop"],
    prop_dec_2 = prop_500_df[, "broadleaf_prop"],
    prop_con_3 = prop_1000_df[, "conifer_prop"],
    prop_dec_3 = prop_1000_df[, "broadleaf_prop"],
    prop_Sw_1 = prop_150_df[, "picegla_prop"],
    prop_Sw_2 = prop_500_df[, "picegla_prop"],
    prop_Sw_3 = prop_1000_df[, "picegla_prop"],
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

Conifer_10TM_crop <- rast("Output/Spatial Data/reproj_tree_rasters/Conifer_10TM_crop.tif")

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
Age_10TM <- rast("Output/Spatial Data/reproj_tree_rasters/Age_reproj.tif")

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
st_write(AB_point_counts_sf_filtered, "Output/Spatial Data/AB_point_counts/point_count_locs_NTEMS.shp") 

# Save as CSV
AB_point_counts_df_filtered <- st_drop_geometry(AB_point_counts_sf_filtered)
write.csv(AB_point_counts_df_filtered, "Output/Tabular Data/point_counts_NTEMS.csv", row.names = FALSE)























########## Calculate Aggregation Index

# Load the binary conifer raster (1 = conifer, 0 = no conifer)
Conifer_10TM <- rast("Output/Spatial Data/reproj_tree_rasters/Conifer_reproj.tif")

# Function to calculate AI for each buffer size
calculate_ai <- function(buffer) {
  lapply(1:nrow(buffer), function(i) {
    buffer_raster <- terra::crop(Conifer_10TM, buffer[i, ])  # Crop raster
    
    # Ensure raster is valid and contains conifer pixels
    if (!is.null(buffer_raster) && terra::ncell(buffer_raster) > 0) {
      buffer_raster <- raster(buffer_raster)  # Convert to raster for landscapemetrics
      
      # Check if conifer pixels are present
      if (any(values(buffer_raster) == 1, na.rm = TRUE)) {
        # Calculate AI only for conifer (class 1)
        ai_values <- calculate_lsm(buffer_raster, level = "landscape", metric = "ai", class = 1)
        
        if (nrow(ai_values) > 0) {
          return(ifelse(is.na(ai_values$value), 0, ai_values$value))  # Replace NA with 0
        } else {
          return(0)  # If no AI value was calculated, return 0
        }
      } else {
        return(0)  # No conifer pixels, AI should be 0 instead of NA
      }
    } else {
      return(0)  # Return 0 for invalid buffers
    }
  })
}

# Calculate AI for each buffer size
AI_150 <- calculate_ai(buffer150)
AI_500 <- calculate_ai(buffer500)
AI_1000 <- calculate_ai(buffer1000)

# Ensure AI outputs are lists and unlist them
AI_150 <- unlist(AI_150)
AI_500 <- unlist(AI_500)
AI_1000 <- unlist(AI_1000)

# Create AI_results data frame using the ID column from AB_point_counts_sf
AI_results <- data.frame(
  ID = AB_point_counts_sf$PointID,  # Use the existing 'ID' column
  AI_150 = AI_150,             # Add AI_150 values
  AI_500 = AI_500,             # Add AI_500 values
  AI_1000 = AI_1000            # Add AI_1000 values
)

# Check if AI columns exist before appending
if (!inherits(AB_point_counts_sf, "sf")) stop("AB_point_counts_sf is not an sf object!")

# Append AI results to the shapefile, keeping the 'ID' column
AB_point_counts_sf <- st_sf(cbind(AB_point_counts_sf, AI_results[, -1]))  # Exclude the ID column from AI_results

# Rename the AI columns safely
ai_cols <- c("AI_150", "AI_500", "AI_1000")
if (all(ai_cols %in% names(AB_point_counts_sf))) {
  names(AB_point_counts_sf)[names(AB_point_counts_sf) %in% ai_cols] <- 
    paste0("NTEMS_AI_", c(1, 2, 3))
} else {
  stop("AI columns not found in AB_point_counts_sf!")
}






########## Calculate Patch Density

# Function to calculate Patch Density for each buffer size
calculate_pd <- function(buffer) {
  lapply(1:nrow(buffer), function(i) {
    buffer_raster <- terra::crop(Conifer_10TM, buffer[i, ])  # Crop raster
    
    # Ensure raster is valid and contains conifer pixels
    if (!is.null(buffer_raster) && terra::ncell(buffer_raster) > 0) {
      buffer_raster <- raster(buffer_raster)  # Convert to raster for landscapemetrics
      
      # Check if conifer pixels are present
      if (any(values(buffer_raster) == 1, na.rm = TRUE)) {
        # Calculate Patch Density only for conifer (class 1)
        pd_values <- calculate_lsm(buffer_raster, level = "landscape", metric = "pd", class = 1)
        
        if (nrow(pd_values) > 0) {
          return(ifelse(is.na(pd_values$value), 0, pd_values$value))  # Replace NA with 0
        } else {
          return(0)  # If no PD value was calculated, return 0
        }
      } else {
        return(0)  # No conifer pixels, PD should be 0 instead of NA
      }
    } else {
      return(0)  # Return 0 for invalid buffers
    }
  })
}

# Calculate Patch Density for each buffer size
PD_150 <- calculate_pd(buffer150)
PD_500 <- calculate_pd(buffer500)
PD_1000 <- calculate_pd(buffer1000)

# Ensure PD outputs are lists and unlist them
PD_150 <- unlist(PD_150)
PD_500 <- unlist(PD_500)
PD_1000 <- unlist(PD_1000)

# Create PD_results data frame using the ID column from AB_point_counts_sf
PD_results <- data.frame(
  ID = AB_point_counts_sf$PointID,  # Use the existing 'ID' column
  PD_150 = PD_150,             # Add PD_150 values
  PD_500 = PD_500,             # Add PD_500 values
  PD_1000 = PD_1000            # Add PD_1000 values
)

# Check if PD columns exist before appending
if (!inherits(AB_point_counts_sf, "sf")) stop("AB_point_counts_sf is not an sf object!")

# Append PD results to the shapefile, keeping the 'ID' column
AB_point_counts_sf <- st_sf(cbind(AB_point_counts_sf, PD_results[, -1]))  # Exclude the ID column from PD_results

# Rename the PD columns safely
pd_cols <- c("PD_150", "PD_500", "PD_1000")
if (all(pd_cols %in% names(AB_point_counts_sf))) {
  names(AB_point_counts_sf)[names(AB_point_counts_sf) %in% pd_cols] <- 
    paste0("NTEMS_PD_", c(1, 2, 3))
} else {
  stop("PD columns not found in AB_point_counts_sf!")
}








### Create index based on AI and PD

# Multiply AI and PD for each buffer size
AB_point_counts_sf$FragIndex_1 <- AB_point_counts_sf$NTEMS_AI_1 * AB_point_counts_sf$NTEMS_PD_1
AB_point_counts_sf$FragIndex_2 <- AB_point_counts_sf$NTEMS_AI_2 * AB_point_counts_sf$NTEMS_PD_2
AB_point_counts_sf$FragIndex_3 <- AB_point_counts_sf$NTEMS_AI_3 * AB_point_counts_sf$NTEMS_PD_3


# Save the updated shapefile
st_write(AB_point_counts_sf, "Output/Spatial Data/AB_point_counts/AB_point_counts_10TM_3.shp", append = FALSE)









########### Create index based on AI and PD

# Multiply AI and PD for each buffer size
AB_point_counts_sf$FragIndex_1 <- AB_point_counts_sf$NTEMS_AI_1 * AB_point_counts_sf$NTEMS_PD_1
AB_point_counts_sf$FragIndex_2 <- AB_point_counts_sf$NTEMS_AI_2 * AB_point_counts_sf$NTEMS_PD_2
AB_point_counts_sf$FragIndex_3 <- AB_point_counts_sf$NTEMS_AI_3 * AB_point_counts_sf$NTEMS_PD_3

# Check the new columns
head(AB_point_counts_sf[, c("FragIndex_1", "FragIndex_2", "FragIndex_3")])
