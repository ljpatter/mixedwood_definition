library(terra)
library(raster)

# --- Input (10 m classes, 9=conifer, 11=mixedwood) ---
LULC_10m <- rast("Input/Spatial Data/LULC/AB_LULC_2020_10tm.tif")

# --- Upsample to 5 m by replication (2×2) ---
LULC_5m <- disagg(LULC_10m, fact = 2)

# --- Output target ---
out_path <- "Output/Spatial Data/Reprojected LULC/Conifer_AB_LULC_2020_5m_reclass_mixedwood_split.tif"
out <- rast(LULC_5m)  # same grid as LULC_5m

# --- Block sizing via raster::blockSize() wrapper ---
r_compat <- raster::raster(LULC_5m)
bs <- raster::blockSize(r_compat)   # gives $n, $row, $nrows
nC <- terra::ncol(LULC_5m)

# --- Open input for reading, output for writing ---
terra::readStart(LULC_5m)
terra::writeStart(
  out, filename = out_path, overwrite = TRUE,
  wopt = list(
    datatype = "INT1U",
    gdal     = c("COMPRESS=DEFLATE","ZLEVEL=9","TILED=YES","BIGTIFF=YES")
  )
)

# --- Block loop ---
for (i in 1:bs$n) {
  # Read a block of class codes
  v <- terra::readValues(LULC_5m, row = bs$row[i], nrows = bs$nrows[i])
  
  # Build 5 m row/col indices for this block
  rows <- rep(seq(from = bs$row[i], length.out = bs$nrows[i]), each = nC)
  cols <- rep(seq_len(nC), times = bs$nrows[i])
  
  # Default output = 0 (non-conifer); preserve NAs
  outv <- rep(0L, length(v))
  na   <- is.na(v)
  outv[na] <- NA_integer_
  
  # Conifer (9) -> 1
  idx9 <- !na & (v == 9L)
  if (any(idx9)) outv[idx9] <- 1L
  
  # Mixedwood (11): exactly 2 of the 4 subcells per parent 10 m cell are 1
  idx11 <- !na & (v == 11L)
  if (any(idx11)) {
    r <- rows[idx11]; c <- cols[idx11]
    pos <- as.integer(((r - 1L) %% 2L) * 2L + ((c - 1L) %% 2L))   # {0,1,2,3}
    gr  <- ((r - 1L) %/% 2L) + 1L                                  # parent 10 m row
    gc  <- ((c - 1L) %/% 2L) + 1L                                  # parent 10 m col
    pat <- as.integer((gr * 1103515245L + gc * 12345L) %% 6L)      # 6 patterns
    sel <- ( (pat == 0L & (pos == 0L | pos == 1L)) |
               (pat == 1L & (pos == 0L | pos == 2L)) |
               (pat == 2L & (pos == 0L | pos == 3L)) |
               (pat == 3L & (pos == 1L | pos == 2L)) |
               (pat == 4L & (pos == 1L | pos == 3L)) |
               (pat == 5L & (pos == 2L | pos == 3L)) )
    outv[idx11] <- as.integer(sel)
  }
  
  # Write this block (note the nrows argument)
  terra::writeValues(out, outv, bs$row[i], nrows = bs$nrows[i])
}

# --- Close files ---
conifer_5m <- terra::writeStop(out)
terra::readStop(LULC_5m)

# # Optional quick checks:
terra::global(conifer_5m, range, na.rm = TRUE)  # expect 0..1
terra::freq(conifer_5m)                         # counts of 0/1 (excludes NA)















############# Calculate proportion conifer 

# Load libraries
library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(landscapemetrics)  
library(raster)

# Load point count location from NTEMS
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/point_count_locs_NTEMS.shp")

# Remove unnecessary column
AB_point_counts_sf <- AB_point_counts_sf %>%
  dplyr::select(PointID, location, year, x_AEP10TM, y_AEP10TM, geometry)

# Load raster
Conifer_10TM <- rast("Output/Spatial Data/Reprojected LULC/Conifer_AB_LULC_2020_5m_reclass_mixedwood_split.tif")


### Calculate proportion of conifer 

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
st_write(AB_point_counts_sf, "Output/Spatial Data/AB_point_counts/point_count_locs_AB_LULC_5m.shp", append = FALSE) 

# Save as CSV
AB_point_counts_df <- st_drop_geometry(AB_point_counts_sf)
write.csv(AB_point_counts_df, "Output/Tabular Data/point_counts_AB_LULC_5m.csv", row.names = FALSE)
