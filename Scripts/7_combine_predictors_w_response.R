# Clear environment

rm(list = ls())  # Removes all objects from the environment

# Load libraries

library(dplyr)

# Load data
NTEMS <- read.csv("Output/Tabular Data/point_counts_NTEMS.csv")
LULC <- read.csv("Output/Tabular Data/point_counts_AB_LULC.csv")
bird_data <- read.csv("Output/Tabular Data/bird_data_2.csv")



####### Prep covariate and count data for join

# Ensure that NTEMS and LULC contain the same sites
LULC_NTEMS_diff <- LULC %>%
  filter(!(location %in% NTEMS$location))
print(LULC_NTEMS_diff)

# Filter bird_data to only include species of interest
bd_filt_1 <- bird_data %>%
  dplyr::select(project, project_id, location, survey_type, x_AEP10TM,
                y_AEP10TM, survey_duration_method, survey_distance_method, 
                max_dist_band, hssr, ordinalDay, year, BTNW, BBWA, TEWA)

# Filter bird_data to only include species of interest
#bd_filt_1 <- bird_data %>%
#  dplyr::select(project, project_id, survey_type, location, x_AEP10TM,
#                y_AEP10TM, survey_duration_method, survey_distance_method, 
#                survey_effort, year, BTNW, BBWA, TEWA)



# Filter bird data to only include LULC and NTEMS sites 
bd_filt_2 <- bd_filt_1 %>%
  filter(location %in% NTEMS$location)

# Filter bird data to only include years of interest (i.e., after 2009)

bird_data_2009 <- bd_filt_2 %>% filter(year >= 2009)

# Filter NTEMS and LULC to only include unique rows based on "location"
NTEMS_unique_loc <- NTEMS %>%
  distinct(location, .keep_all = TRUE)

LULC_unique_loc <- LULC %>%
  distinct(location, .keep_all = TRUE)







### Now that bird and NTEMS have the same # of locations,
### merge NTEMS -> bird_data using location column

# Ensure "location" is a character type
bird_data_2009$location <- as.character(bird_data_2009$location)
NTEMS_unique_loc$location <- as.character(NTEMS_unique_loc$location)

# Left Join bird data and NTEMS
NTEMS_combined <- left_join(bird_data_2009, NTEMS_unique_loc, by = "location")

# Check for discrepancies between coordinates in each df
NTEMS_mismatch <- NTEMS_combined %>%
  mutate(x_coordinates_match = ifelse(x_AEP10TM.x != x_AEP10TM.y, "Mismatch", "Match"))
NTEMS_mismatch_df <- NTEMS_mismatch %>% filter(x_coordinates_match == "Mismatch")

# Clean df to only include columns of interest

NTEMS_combined <- NTEMS_combined %>%
  dplyr::select(project, location, survey_type, year.x, ordinalDay, hssr,  
         survey_duration_method, max_dist_band, x_AEP10TM.x, y_AEP10TM.x, prop_con_1, 
         prop_con_2, prop_con_3, prop_dec_1, prop_dec_2, prop_dec_3, clumpy_1, clumpy_2, 
         clumpy_3, age_mn_1, age_mn_2, age_mn_3,
         BBWA, BTNW, TEWA) %>%
  rename(year = year.x, x_AEP10TM = x_AEP10TM.x, y_AEP10TM = y_AEP10TM.x)

### Create survey effort column

NTEMS_combined$survey_duration_method_clean <- gsub("min", "", NTEMS_combined$survey_duration_method)
NTEMS_combined$survey_duration_method_clean <- gsub("\\+", "", NTEMS_combined$survey_duration_method_clean)  # Remove "+"

# Extract survey duration
NTEMS_combined$survey_effort <- sapply(strsplit(as.character(NTEMS_combined$survey_duration_method_clean), "-"), function(x) {
  x <- x[x != ""] # Remove empty strings that can occur after splitting
  if (length(x) > 0) {
    max(as.numeric(x), na.rm = TRUE)
  } else {
    NA
  }
})

# Remove cleaned column
NTEMS_combined$survey_duration_method_clean <- NULL

# Check for NA values in the survey_effort column
sum(is.na(NTEMS_combined$survey_effort))


### Adjust forest age based on the difference from 2019
NTEMS_combined_age_corr <- NTEMS_combined %>%
  mutate(
    age_mn_1 = age_mn_1 - (2019 - year),
    age_mn_2 = age_mn_2 - (2019 - year),
    age_mn_3 = age_mn_3 - (2019 - year)
  )

# Remove rows where there are NA values in any column
NTEMS_final <- NTEMS_combined_age_corr %>%
  drop_na()

length(unique(NTEMS_final$location))











### Now that bird and LULC have the same # of locations,
### merge LULC -> bird_data using location column

# Ensure "location" is a character type
bird_data_2009$location <- as.character(bird_data_2009$location)
LULC_unique_loc$location <- as.character(LULC_unique_loc$location)

# Left Join bird data and NTEMS
LULC_combined <- left_join(bird_data_2009, LULC_unique_loc, by = "location")

# Check for discrepancies between coordinates in each df
LULC_mismatch <- LULC_combined %>%
  mutate(x_coordinates_match = ifelse(x_AEP10TM.x != x_AEP10TM.y, "Mismatch", "Match"))
LULC_mismatch_df <- LULC_mismatch %>% filter(x_coordinates_match == "Mismatch")

# Clean df to only include columns of interest

LULC_combined <- LULC_combined %>%
  dplyr::select(project, location, survey_type, year.x, ordinalDay, hssr,  
                survey_duration_method, max_dist_band, x_AEP10TM.x, y_AEP10TM.x, prop_con_1, 
                prop_con_2, prop_con_3, prop_dec_1, prop_dec_2, prop_dec_3, clumpy_1, clumpy_2, 
                clumpy_3, age_mn_1, age_mn_2, age_mn_3,
                BBWA, BTNW, TEWA) %>%
  rename(year = year.x, x_AEP10TM = x_AEP10TM.x, y_AEP10TM = y_AEP10TM.x)

### Create survey effort column

LULC_combined$survey_duration_method_clean <- gsub("min", "", LULC_combined$survey_duration_method)
LULC_combined$survey_duration_method_clean <- gsub("\\+", "", LULC_combined$survey_duration_method_clean)  # Remove "+"

# Extract survey duration
LULC_combined$survey_effort <- sapply(strsplit(as.character(LULC_combined$survey_duration_method_clean), "-"), function(x) {
  x <- x[x != ""] # Remove empty strings that can occur after splitting
  if (length(x) > 0) {
    max(as.numeric(x), na.rm = TRUE)
  } else {
    NA
  }
})

# Remove cleaned column
LULC_combined$survey_duration_method_clean <- NULL


### Adjust forest age based on the difference from 2019
LULC_combined_age_corr <- LULC_combined %>%
  mutate(
    age_mn_1 = age_mn_1 - (2019 - year),
    age_mn_2 = age_mn_2 - (2019 - year),
    age_mn_3 = age_mn_3 - (2019 - year)
  )

# Remove rows where there are NA values in any column
LULC_final <- LULC_combined_age_corr %>%
  drop_na()


## NTEMS had more NA values (primarily in clumpy), due to having a 
## different pixel size. This code will filter LULC to only include
## sites that are also present in NTEMS.

# Create a unique identifier for location and year in both dataframes
NTEMS_final$location_year <- paste(NTEMS_final$location, NTEMS_final$year, sep = "_")
LULC_final$location_year <- paste(LULC_final$location, LULC_final$year, sep = "_")

# Find the intersection of location_year combinations
shared_location_years <- intersect(NTEMS_final$location_year, LULC_final$location_year)

# Filter both dataframes
NTEMS_final_2 <- NTEMS_final[NTEMS_final$location_year %in% shared_location_years, ]
LULC_final_2 <- LULC_final[LULC_final$location_year %in% shared_location_years, ]

# Remove location_year column

NTEMS_final_2 <- NTEMS_final_2 %>%
  dplyr::select(-location_year)
LULC_final_2 <- LULC_final_2 %>%
  dplyr::select(-location_year)


# Save NTEMS and LULC dfs

write.csv(NTEMS_final_2, "Output/Tabular Data/NTEMS_with_bird_data.csv")
write.csv(LULC_final_2, "Output/Tabular Data/LULC_with_bird_data.csv")




































########## Lastly, remove any point counts locations that are not located
########## in a forest.

### First, reproject NTEMS land cover raster

# Import raster
NTEMS_LC <- rast("Input/Spatial Data/NTEMS_LandCover/CA_forest_VLCE2_2019.tif")

# Extract the WKT string from your desired CRS
target_crs_wkt <- st_as_text(st_crs(LULC_final))

# Create the target CRS object using terra::crs()
target_crs_terra <- terra::crs(target_crs_wkt)

# Reproject rasters
NTEMS_LC_10TM <- terra::project(NTEMS_LC, target_crs_terra, method = "near")

## Save reprojected rasters

## Define output directory and ensure it exists
output_dir <- "Output/Spatial Data/NTEMS_LandCover"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define output file path
output_file <- file.path(output_dir, "NTEMS_LandCover_reprojected.tif")

# Save the reprojected raster
terra::writeRaster(NTEMS_LC_10TM, filename = output_file, overwrite = TRUE)


# Load reprojected NTEMS land cover raster

NTEMS_LC <- rast("Input/Spatial Data/NTEMS_LandCover/NTEMS_LandCover_reprojected.tif")












