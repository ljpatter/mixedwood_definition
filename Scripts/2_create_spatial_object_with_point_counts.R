# ---
# title: "Create spatial object with point count locations "
# author: "Brendan Casey and Leonard Patterson"
# created: "2025-04-11"
# description: "Creating shapefiles of point count locations"
# ---


# Load packages
library(sf)
library(tidyverse)

# Load the cleaned data (ensure the object 'wildtrax_cleaned' exists in the file)
load("Output/R Data/wildtrax_cleaned_2025-09-15.rData")

# Read AB boundary and set target CRS
AB_boundary_10TM <- st_read("Input/Spatial Data/Alberta/AB_boundary.shp")
target_crs <- st_crs(AB_boundary_10TM)

# Build unique site-year table
ss_xy <- wildtrax_cleaned %>%
  ungroup() %>%
  select(location, year, lon, lat, x_AEP10TM, y_AEP10TM) %>%
  distinct() %>%
  filter(!is.na(x_AEP10TM), !is.na(y_AEP10TM))  # avoid st_as_sf errors

# Make projected points from AEP 10TM coords
ss_xy_10TM <- ss_xy %>%
  st_as_sf(coords = c("x_AEP10TM", "y_AEP10TM"), crs = target_crs, remove = FALSE)

# Keep points within Alberta boundary 
AB_point_counts <- st_filter(ss_xy_10TM, AB_boundary_10TM, .pred = st_within)

# Select final columns
AB_point_counts <- AB_point_counts %>%
  select(location, year, lon, lat, x_AEP10TM, y_AEP10TM, geometry)

# Write shapefiles
st_write(AB_point_counts, "Output/Spatial Data/AB_point_counts/AB_point_counts_10TM.shp", append = FALSE)
st_write(ss_xy_10TM,      "Output/Spatial Data/AB_point_counts/ss_xy_10TM.shp",          append = FALSE)

# Save RData objects with distinct filenames
save(AB_point_counts, file = paste0("Output/R Data/AB_point_counts_", Sys.Date(), ".rData"))
save(ss_xy_10TM,      file = paste0("Output/R Data/ss_xy_10TM_",      Sys.Date(), ".rData"))


