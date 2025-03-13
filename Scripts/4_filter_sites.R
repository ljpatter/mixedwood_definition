# ---
# title: "NTEMS_data_extraction"
# author: "Leonard Patterson"
# created: ""
# description: "."
# ---


rm(list = ls())  # Removes all objects from the environment


# Load libraries

library(dplyr)
library(sf)

########### Filter out sites that are within 50 m of a road, pipeline, wellpad, or transmission line
########### for all sites 

# Load point count data
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/AB_point_counts_10TM.shp")

# Read all vector layers
pipelines <- st_read("Input/Spatial Data/ABMI Human Footprint/o19_Pipelines_HFI_2021.shp")
wellsabnd <- st_read("Input/Spatial Data/ABMI Human Footprint/o16_WellsAbnd_HFI_2021.shp")
roads <- st_read("Input/Spatial Data/ABMI Human Footprint/o03_Roads_HFI_2021.shp")
cultivation <- st_read("Input/Spatial Data/ABMI Human Footprint/o17_Cultivation_HFI_2021.shp")
wellsactive <- st_read("Input/Spatial Data/ABMI Human Footprint/o09_WellsActive_HFI_2021.shp")
harvestareas <- st_read("Input/Spatial Data/ABMI Human Footprint/o18_HarvestAreas_HFI_2021.shp")
wildfire <- st_read("Input/Spatial Data/HistoricalWildfirePerimeters/WildfirePerimeters1931to2023.shp")
industrial <- st_read("Input/Spatial Data/ABMI Human Footprint/o08_Industrials_HFI_2021.shp")
residential <- st_read("Input/Spatial Data/ABMI Human Footprint/o15_Residentials_HFI_2021.shp")




########## Filter out sites before 2009

# Load point count data
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/AB_point_counts_10TM.shp")

# Filter points counts to only include those after 2009
AB_point_counts_sf <- AB_point_counts_sf %>% filter(year >= 2009)

# Filter out wildfire and harvest area where YEAR > 2009
harvestareas_2009 <- harvestareas %>% filter(YEAR >= 2009)
wildfire_2009 <- wildfire %>% filter(YEAR >= 2009)

# Create a unique ID column (Important!)
AB_point_counts_sf <- AB_point_counts_sf %>%
  mutate(AB_point_counts_sf_ID = 1:n())

# Apply buffer (AFTER creating the ID)
AB_point_counts_buff <- st_buffer(AB_point_counts_sf, dist = 50)

# Create points to keep (Correctly initialized)
points_to_keep_2009 <- AB_point_counts_buff # Use AB_point_counts_buff directly

# List of layers to filter against
layers <- list(pipelines, wellsabnd, roads, cultivation, wellsactive, harvestareas_2009, wildfire_2009, residential, industrial)

# Iterate and filter (optimized without st_prepare)
for (layer in layers) {
  intersected_points <- st_join(points_to_keep_2009, layer, join = st_intersects, left = FALSE)
  
  # Get the IDs of the intersecting points
  intersecting_ids <- intersected_points %>%
    st_drop_geometry() %>%
    distinct(AB_point_counts_sf_ID) %>%
    pull(AB_point_counts_sf_ID)
  
  # Filter points_to_keep by ID (Correctly using points_to_keep)
  points_to_keep_2009 <- points_to_keep_2009 %>%
    filter(!(AB_point_counts_sf_ID %in% intersecting_ids))
  
  gc()
}

# Restore original point geometries
filtered_ids <- points_to_keep_2009 %>%
  st_drop_geometry() %>%
  pull(AB_point_counts_sf_ID)

points_to_keep_2009 <- AB_point_counts_sf %>%
  filter(AB_point_counts_sf_ID %in% filtered_ids)

# Remove ID column and select final columns
points_to_keep_2009 <- points_to_keep_2009 %>%
  dplyr::select(location, year, lon, lat, x_AEP10TM, y_AEP10TM, geometry)

# Save layer
st_write(points_to_keep_2009, "Output/Spatial Data/AB_point_counts/filtered_point_count_loc_2009.shp")






































### Filter out sites before 1990

# Filter out wildfire and harvest area where YEAR > 1990
harvestareas_1990 <- harvestareas %>% filter(YEAR >= 1990)
wildfire_1990 <- wildfire %>% filter(YEAR >= 1990)

# Create a unique ID column (Important!)
AB_point_counts_sf <- AB_point_counts_sf %>%
  mutate(AB_point_counts_sf_ID = 1:n())

# Apply buffer (AFTER creating the ID)
AB_point_counts_buff <- st_buffer(AB_point_counts_sf, dist = 50)

# Create points to keep (Correctly initialized)
points_to_keep_1990 <- AB_point_counts_buff # Use AB_point_counts_buff directly

# List of layers to filter against
layers <- list(pipelines, wellsabnd, roads, cultivation, wellsactive, residential, industrial)

# Iterate and filter (optimized without st_prepare)
for (layer in layers) {
  intersected_points <- st_join(points_to_keep_1990, layer, join = st_intersects, left = FALSE)
  
  # Get the IDs of the intersecting points
  intersecting_ids <- intersected_points %>%
    st_drop_geometry() %>%
    distinct(AB_point_counts_sf_ID) %>%
    pull(AB_point_counts_sf_ID)
  
  # Filter points_to_keep by ID (Correctly using points_to_keep)
  points_to_keep_1990 <- points_to_keep_1990 %>%
    filter(!(AB_point_counts_sf_ID %in% intersecting_ids))
  
  gc()
}

st_write(points_to_keep_1990, "Output/Spatial Data/AB_point_counts/filtered_point_count_locs_1990.shp")















