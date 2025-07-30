# ---
# title: "Combine predictors with response"
# author: "Leonard Patterson"
# created: "2025-04-11"
# description: "This script combines covariates data with point count data
# ---

# Load libraries

library(dplyr)
library(tidyr)

# Load data
NTEMS <- read.csv("Output/Tabular Data/point_counts_NTEMS.csv")
LULC <- read.csv("Output/Tabular Data/point_counts_AB_LULC_5m.csv")
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
NTEMS_mismatch_df <- NTEMS_combined %>%
  mutate(x_ok = near(x_AEP10TM.x, x_AEP10TM.y),
         y_ok = near(y_AEP10TM.x, y_AEP10TM.y)) %>%
  filter(!(x_ok & y_ok)) %>%
  dplyr::select(location, x_AEP10TM.x, x_AEP10TM.y, y_AEP10TM.x, y_AEP10TM.y)


# Clean df to only include columns of interest

NTEMS_combined <- NTEMS_combined %>%
  dplyr::select(project, location, survey_type, year.x, ordinalDay, hssr,  
         survey_duration_method, max_dist_band, x_AEP10TM.x, y_AEP10TM.x, prop_con_1, 
         prop_con_2, prop_con_3, clumpy_1, clumpy_2, clumpy_3, age_mn_1, age_mn_2, age_mn_3,
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

# Drop any row with NA anywhere, then drop any row with negative ages
NTEMS_final <- NTEMS_combined_age_corr %>%
  drop_na() %>%                                # remove rows with any NA in any column
  filter(age_mn_1 >= 0,                        # keep only non-negative ages
         age_mn_2 >= 0,
         age_mn_3 >= 0)

length(unique(NTEMS_final$location))






### Now that bird and LULC have the same # of locations,
### merge LULC -> bird_data using location column

# Ensure "location" is a character type
bird_data_2009$location <- as.character(bird_data_2009$location)
LULC_unique_loc$location <- as.character(LULC_unique_loc$location)

# Left Join bird data and NTEMS
LULC_combined <- left_join(bird_data_2009, LULC_unique_loc, by = "location")

# Check for discrepancies between coordinates in each df
LULC_mismatch_df <- LULC_combined %>%
  mutate(x_ok = near(x_AEP10TM.x, x_AEP10TM.y),
         y_ok = near(y_AEP10TM.x, y_AEP10TM.y)) %>%
  filter(!(x_ok & y_ok)) %>%
  dplyr::select(location, x_AEP10TM.x, x_AEP10TM.y, y_AEP10TM.x, y_AEP10TM.y)

# Clean df to only include columns of interest

LULC_combined <- LULC_combined %>%
  dplyr::select(project, location, survey_type, year.x, ordinalDay, hssr,  
                survey_duration_method, max_dist_band, x_AEP10TM.x, y_AEP10TM.x, prop_con_1, 
                prop_con_2, prop_con_3, clumpy_1, clumpy_2, clumpy_3,
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

# Remove rows where there are NA values in any column
LULC_final <- LULC_combined %>%
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

write.csv(NTEMS_final_2, "Output/Tabular Data/NTEMS_with_bird_data_5m.csv")
write.csv(LULC_final_2, "Output/Tabular Data/LULC_with_bird_data_5m.csv")






############### Next, combine LULC and NTEMS dataframes

# Load data
NTEMS <- read.csv("Output/Tabular Data/NTEMS_with_bird_data_5m.csv")
LULC  <- read.csv("Output/Tabular Data/LULC_with_bird_data_5m.csv")

### Merge prop_con_, clumpy_, and age_mn_ from LULC to NTEMS

# Remove the 'X'columns from both NTEMS and LULC before renaming
NTEMS_cleaned <- NTEMS %>% dplyr::select(-c(X))
LULC_cleaned <- LULC %>% dplyr::select(-c(X))

# Rename _1, _2, and _3 to _150, _500, and _1000
scale_suffix_map <- c("_1" = "_150", "_2" = "_500", "_3" = "_1000")

rename_scale_suffixes <- function(names) {
  for (old in names(scale_suffix_map)) {
    names <- gsub(old, scale_suffix_map[[old]], names)
  }
  names
}

prefixes <- c("prop_con_", "clumpy_", "age_mn_")

NTEMS_converted <- NTEMS_cleaned %>%
  rename_with(rename_scale_suffixes, .cols = matches(paste0("^(", paste(prefixes, collapse = "|"), ")")))

LULC_converted <- LULC_cleaned %>%
  rename_with(rename_scale_suffixes, .cols = matches(paste0("^(", paste(prefixes, collapse = "|"), ")")))

NTEMS_renamed <- NTEMS_converted %>%
  rename_with(.cols = starts_with("prop_con_"), ~ paste0(., "_NTEMS")) %>%
  rename_with(.cols = starts_with("clumpy_"), ~ paste0(., "_NTEMS")) %>%
  rename_with(.cols = starts_with("age_mn_"), ~ paste0(., "_NTEMS"))

LULC_renamed <- LULC_converted %>%
  rename_with(.cols = starts_with("prop_con_"), ~ paste0(., "_LULC")) %>%
  rename_with(.cols = starts_with("clumpy_"), ~ paste0(., "_LULC"))

# Join the renamed data frames
joined_data <- left_join(NTEMS_renamed, LULC_renamed,
                         by = c("project", "location", "survey_type", "year", "ordinalDay",
                                "survey_duration_method", "max_dist_band"),
                         suffix = c(".x", ".y"))

# Check for mismatches before dropping columns
check_match <- function(df) {
  cols_to_check <- names(NTEMS_renamed)[!names(NTEMS_renamed) %in% c(
    "prop_con_1_NTEMS", "prop_con_2_NTEMS", "prop_con_3_NTEMS",
    "prop_dec_1_NTEMS", "prop_dec_2_NTEMS", "prop_dec_3_NTEMS",
    "clumpy_1_NTEMS", "clumpy_2_NTEMS", "clumpy_3_NTEMS",
    "age_mn_1_NTEMS", "age_mn_2_NTEMS", "age_mn_3_NTEMS"
  )]
  
  for (col in cols_to_check) {
    col_x <- paste0(col, ".x")
    col_y <- paste0(col, ".y")
    
    if (col_y %in% names(df)) {
      if (!all(df[[col_x]] == df[[col_y]], na.rm = TRUE)) {
        warning(paste("Mismatch found in column:", col))
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

# Add back species columns when reorganizing
species_cols <- c("BTNW", "TEWA", "BBWA") 

if (check_match(joined_data)) {
  joined_data <- joined_data %>%
    dplyr::select(-ends_with(".y"), -starts_with("prop_dec_")) %>%
    rename_with(~ gsub("\\.x$", "", .))
  
  joined_data <- joined_data %>%
    dplyr::select(
      project, location, survey_type, year, ordinalDay, hssr, survey_effort,  # <- added
      survey_duration_method, max_dist_band, x_AEP10TM, y_AEP10TM,
      all_of(species_cols),
      starts_with("prop_con_"), starts_with("clumpy_"), starts_with("age_mn_")
    )
  
  
  print("Join and replacement successful.")
} else {
  print("Mismatches found. Join and replacement aborted.")
}


# Save
write.csv(joined_data, "Output/Tabular Data/joined_data_5m.csv", row.names = FALSE)
