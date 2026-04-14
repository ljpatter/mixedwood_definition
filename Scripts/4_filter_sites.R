# ---
# title: "Filter sites"
# author: "Leonard Patterson"
# created: "2025-04-11"
# description: "This script filters point count data using disturbance shapefiles
# from ABMI"
# ---

library(dplyr)
library(sf)

# =========================
# User-defined settings
# =========================
landcover_year <- 2020   # set this to the effective year of your landcover product

# =========================
# Load points
# =========================
pts <- st_read("Output/Spatial Data/AB_point_counts/AB_point_counts_10TM.shp", quiet = TRUE) %>%
  filter(year >= 2009) %>%
  mutate(AB_point_counts_sf_ID = dplyr::row_number())

# 50 m buffer (in projected meters)
pts_buf <- st_buffer(pts, dist = 50) %>%
  rename(point_year = year)

# =========================
# Load disturbance layers
# =========================
pipelines    <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/o19_Pipelines_HFI_2022.shp", quiet = TRUE)
wellsabnd    <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/o16_WellsAbnd_HFI_2022.shp", quiet = TRUE)
roads        <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/o03_Roads_HFI_2022.shp", quiet = TRUE)
cultivation  <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/o17_Cultivation_HFI_2022.shp", quiet = TRUE)
wellsactive  <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/o09_WellsActive_HFI2022.shp", quiet = TRUE)
harvestareas <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/o18_TimberHarvestandWoodyVegetationRemoval_HFI_2022.shp", quiet = TRUE)
wildfire     <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/WildfirePerimeters1931to2023.shp", quiet = TRUE)
industrial   <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/o08_Industrials_HFI_2022.shp", quiet = TRUE)
transmission <- st_read("C:/Users/leona/OneDrive/Desktop/Spatial Data/ABMI Human Footprint 2022 (extracted)/o13_TransmissionLines_HFI_2022.shp", quiet = TRUE)

# =========================
# Put layers in a list
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

# Reproject all layers to match points if needed
layers <- lapply(layers, function(x) {
  if (st_crs(x) != st_crs(pts_buf)) {
    x <- st_transform(x, st_crs(pts_buf))
  }
  x
})

# =========================
# Filtering logic
# =========================
points_to_keep <- pts_buf
removed_counts <- setNames(integer(length(layers)), names(layers))

for (nm in names(layers)) {
  
  lyr <- layers[[nm]]
  
  # Intersections within 50 m
  inter <- st_join(points_to_keep, lyr, join = st_intersects, left = FALSE)
  if (nrow(inter) == 0) next
  
  if (nm == "harvestareas") {
    
    # Remove only if harvest year falls between point year and landcover year
    # This works for both:
    # 1) point before landcover: point_year < event_year <= landcover_year
    # 2) point after landcover:  landcover_year < event_year <= point_year
    if (!("YEAR" %in% names(inter))) next
    
    inter <- inter %>%
      mutate(event_year = suppressWarnings(as.integer(.data$YEAR))) %>%
      filter(
        !is.na(event_year),
        event_year > pmin(point_year, landcover_year),
        event_year <= pmax(point_year, landcover_year)
      )
    
  } else if (nm == "wildfire") {
    
    # Same logic for wildfire
    if (!("YEAR" %in% names(inter))) next
    
    inter <- inter %>%
      mutate(event_year = suppressWarnings(as.integer(.data$YEAR))) %>%
      filter(
        !is.na(event_year),
        event_year > pmin(point_year, landcover_year),
        event_year <= pmax(point_year, landcover_year)
      )
    
  } else {
    
    # Static footprint: any intersection removes the point
    # (pipelines, wells, roads, cultivation, transmission, industrial)
    
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

# =========================
# Restore original points
# =========================
kept_ids <- points_to_keep %>%
  st_drop_geometry() %>%
  pull(AB_point_counts_sf_ID)

points_final <- pts %>%
  filter(AB_point_counts_sf_ID %in% kept_ids) %>%
  select(-AB_point_counts_sf_ID)

# =========================
# Filter out sites outside study area
# =========================
study_area <- st_read("Input/Spatial Data/Study Area/study_area.shp", quiet = TRUE)

if (st_crs(study_area) != st_crs(points_final)) {
  study_area <- st_transform(study_area, st_crs(points_final))
}

intersecting_points <- st_intersection(points_final, study_area)

# Remove unnecessary columns if they exist
drop_cols <- intersect(c("Shape_Leng", "Shape_Area"), names(intersecting_points))
intersecting_points_trim <- intersecting_points %>%
  select(-all_of(drop_cols))

# =========================
# Save layer
# =========================
st_write(
  intersecting_points_trim,
  "Output/Spatial Data/AB_point_counts/filtered_point_counts.shp",
  append = FALSE
)

# Optional: print summary
print(removed_counts)
