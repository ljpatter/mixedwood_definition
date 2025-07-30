# ---
# title: "Filter sites"
# author: "Leonard Patterson"
# created: "2025-04-11"
# description: "This script filters point count data using disturbance shapefiles
# from ABMI"
# ---

library(dplyr)
library(sf)

# Load points
pts <- st_read("Output/Spatial Data/AB_point_counts/AB_point_counts_10TM.shp", quiet = TRUE) %>%
  filter(year >= 2009) %>%                                  
  mutate(AB_point_counts_sf_ID = dplyr::row_number())

# 50 m buffer (in projected meters)
pts_buf <- st_buffer(pts, dist = 50) %>%
  rename(point_year = year)                                  # avoid collisions in joins

# =========================
# Load disturbance layers
# =========================
pipelines    <- st_read("Input/Spatial Data/ABMI Human Footprint 2022/o19_Pipelines_HFI_2022.shp", quiet = TRUE)
wellsabnd    <- st_read("Input/Spatial Data/ABMI Human Footprint 2022/o16_WellsAbnd_HFI_2022.shp", quiet = TRUE)
roads        <- st_read("Input/Spatial Data/ABMI Human Footprint 2022/o03_Roads_HFI_2022.shp", quiet = TRUE)
cultivation  <- st_read("Input/Spatial Data/ABMI Human Footprint 2022/o17_Cultivation_HFI_2022.shp", quiet = TRUE)
wellsactive  <- st_read("Input/Spatial Data/ABMI Human Footprint 2022/o09_WellsActive_HFI2022.shp", quiet = TRUE)
harvestareas <- st_read("Input/Spatial Data/ABMI Human Footprint 2022/o18_TimberHarvestandWoodyVegetationRemoval_HFI_2022.shp", quiet = TRUE)
wildfire     <- st_read("Input/Spatial Data/Historical Wildfire Perimeters/WildfirePerimeters1931to2023.shp", quiet = TRUE)
industrial   <- st_read("Input/Spatial Data/ABMI Human Footprint 2022/o08_Industrials_HFI_2022.shp", quiet = TRUE)
transmission <- st_read("Input/Spatial Data/ABMI Human Footprint 2022/o13_TransmissionLines_HFI_2022.shp", quiet = TRUE)

# =========================
# Filtering logic
# =========================
layers <- list(
  pipelines    = pipelines,
  wellsabnd    = wellsabnd,
  roads        = roads,
  cultivation  = cultivation,
  wellsactive  = wellsactive,
  harvestareas = harvestareas,
  wildfire     = wildfire,
  transmission = transmission,
  industrial   = industrial
)

points_to_keep <- pts_buf
removed_counts <- setNames(integer(length(layers)), names(layers))

for (nm in names(layers)) {
  lyr <- layers[[nm]]
  
  # Intersections within 50 m
  inter <- st_join(points_to_keep, lyr, join = st_intersects, left = FALSE)
  if (nrow(inter) == 0) next
  
  if (nm == "harvestareas") {
    # Remove if harvest YEAR >= survey year
    if (!("YEAR" %in% names(inter))) next  # no event-year column; skip
    inter <- inter %>%
      mutate(event_year = suppressWarnings(as.integer(.data$YEAR))) %>%
      filter(!is.na(event_year) & event_year >= point_year)
  } else if (nm == "wildfire") {
    # Remove if fire YEAR >= survey year
    if (!("YEAR" %in% names(inter))) next
    inter <- inter %>%
      mutate(event_year = suppressWarnings(as.integer(.data$YEAR))) %>%
      filter(!is.na(event_year) & event_year >= point_year)
  } else {
    # Static footprint: any intersection removes the point
    # (pipelines, wells, roads, cultivation, transmission, industrial)
    # nothing additional to filter
  }
  
  if (nrow(inter) > 0) {
    bad_ids <- inter %>%
      st_drop_geometry() %>%
      distinct(AB_point_counts_sf_ID) %>%
      pull(AB_point_counts_sf_ID)
    
    removed_counts[nm] <- removed_counts[nm] + length(bad_ids)
    
    points_to_keep <- points_to_keep %>%
      filter(!(AB_point_counts_sf_ID %in% bad_ids))
  }
}

# Restore original points
kept_ids <- points_to_keep %>%
  st_drop_geometry() %>%
  pull(AB_point_counts_sf_ID)

points_final <- pts %>%
  filter(AB_point_counts_sf_ID %in% kept_ids) %>%
  select(-AB_point_counts_sf_ID)   # drop helper ID




############## Filter out sites outside of the study area

study_area <- st_read("Input/Spatial Data/Study Area/study_area.shp")

# Perform the spatial intersection
intersecting_points <- st_intersection(points_final, study_area)

# Remove unnecessary columns
intersecting_points_trim <- intersecting_points %>% select(-Shape_Leng, -Shape_Area)

# Save layer
st_write(intersecting_points_trim, "Output/Spatial Data/AB_point_counts/filtered_point_counts.shp", append=FALSE)














































# ---
# title: "NTEMS_data_extraction"
# author: "Leonard Patterson"
# created: "2025-04-11"
# description: "This script filters point count data using disturbance shapefiles
# from ABMI"
# ---
# ABMI footprint data can be found here: https://abmi.ca/data-portal/80.html
# Historical wildfire data can be found here: https://www.alberta.ca/wildfire-maps-and-data

# Load libraries

library(dplyr)
library(sf)

########### Filter out sites that are within 50 m of a road, pipeline, wellpad, harvest area, or transmission line
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
wildfire <- st_read("Input/Spatial Data/Historical Wildfire Perimeters/WildfirePerimeters1931to2023.shp")
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

















