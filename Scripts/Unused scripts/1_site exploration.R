# ---
# title: "site_exploration"
# author: "Leonard Patterson"
# created: "2024-12-17"
# inputs: [list the required input files]
# outputs: [list the output files produced by the script]
# notes: 
#   "This script explores and selects potential site location 
#   for a chapter of my M.Sc. thesis at the University of Alberta.
#   Code in this script download main reports from all publicly available
#   ARU projects in the cloud-based bioacoustic platform WildTrax platform;
#   this code filters these project to only include unique ARU locations.
# ---



# Install and load necessary packages

#install.packages("remotes")
#remotes::install_github("ABbiodiversity/wildrtrax")
library(wildrtrax)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)
library(sf)

# Authenticate with WildTrax using environment variables for credentials
Sys.setenv(WT_USERNAME = "ljpatter", WT_PASSWORD = "Kingedwardpark13")
wt_auth()

# Create table showing all publicly available projects with ARUs in WildTrax
my_projects <- wt_get_download_summary(
  sensor_id = "ARU"
)

# Extract the project_id column from my_projects and store it in a new dataframe
project_id <- data.frame(project_id = my_projects$project_id)

aru_as_pc <- wt_download_report(project_id = 620, sensor_id = 'PC', reports = "main", weather_cols = F)

### Download all main reports 

# Define the output directory
output_dir <- "G:/Shared drives/WildTrax Main Reports"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through all project IDs and download reports
for (id in project_id$project_id) {
  # Construct the full file path
  file_name <- file.path(output_dir, paste0("main_report_", id, ".csv"))
  
  # Download the report with the specified parameters
  tryCatch({
    report <- wt_download_report(
      project_id = id,
      sensor_id = "ARU",
      report = "main",
      weather_cols = FALSE  # Exclude weather columns
    )
    write.csv(report, file_name, row.names = FALSE)
    message(paste("Successfully downloaded and saved project", id, "to", file_name))
  }, error = function(e) {
    message(paste("Failed to download project", id, ":", e$message))
    print(paste("Failed Project ID:", id))  # Print the project ID of the failed project
  })
}



### Create a list of all of the main reports

# Specify the folder path containing the CSV files
folder_path <- "G:/Shared drives/WildTrax Main Reports"

# Create a list of all the files in the folder
main_reports <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Convert main_reports to a dataframe
main_reports <- data.frame(project_id = main_reports, stringsAsFactors = FALSE)

# Extract everything after the last slash of the file_path in csv_files$project_id
main_reports <- main_reports %>%
  mutate(project_name = gsub("^.*/", "", project_id)) # Extract file name from file path

# Extract just the project number from the new project_name column
main_reports <- main_reports %>%
  mutate(project_number = as.numeric(gsub("^main_report_(\\d+)\\.csv$", "\\1", project_name)))

# Save csv_files
write.csv(main_reports, "Output/all_wildtrax_projects.csv")



### Next,  filter through all reports to identify ARUs with hazed locations. 
### Then, remove all of those ARUs from further analysis.

# Define the folder containing the main reports
folder_path <- "G:/Shared drives/WildTrax Main Reports"

# List all the CSV files in the folder
main_report_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty vector to store project IDs with hazed locations (non-NA values)
hazed_projects <- c()

# Loop through each CSV file
for (file in main_report_files) {
  # Read the CSV file
  data <- tryCatch(read.csv(file), error = function(e) NULL)
  
  # Check if the file was read successfully and if the location_buffer_m column exists
  if (!is.null(data) && "location_buffer_m" %in% colnames(data)) {
    # Check if there are any non-NA values in the location_buffer_m column
    if (any(!is.na(data$location_buffer_m))) {
      # Extract the project ID from the file name (assuming the format is 'main_report_<id>.csv')
      project_id <- sub("main_report_", "", sub("\\.csv$", "", basename(file)))
      hazed_projects <- c(hazed_projects, project_id)  # Add project ID to the list
    }
  }
}

# Convert the list of hazed projects into a data frame
hazed_projects_df <- data.frame(project_id = hazed_projects, stringsAsFactors = FALSE)

# Save the results to a CSV file
write.csv(hazed_projects_df, "Output/hazed_WT_projects.csv", row.names = FALSE)



### Remove all projects that have hazed ARU locations from filtered_WT_projects

# Read in hazed_projects_df
hazed_projects_df <- read.csv("Output/hazed_WT_projects.csv")

# Read in file with names of all wildtrax main reports
all_WT_projects <- read.csv("Output/all_wildtrax_projects.csv")

# Create a logical vector indicating which rows to remove
remove_rows <- all_WT_projects$project_number %in% hazed_projects_df$project_id

# Remove projects with buffered ARU locations from all_WT_projects
filtered_WT_projects <- all_WT_projects[!remove_rows, ]

# Save filtered file

write.csv(filtered_WT_projects, "Output/filtered_WT_projects.csv")



###  Next, parse through all of the recordings listed in filtered_WT_projects
###  and create a new dataframe that contains all ARU locations.

# Load in filtered_WT_projects
filtered_WT_projects <- read.csv("Output/filtered_WT_projects.csv")

# Define the path to your folder of CSV files
folder_path <- "G:/Shared drives/WildTrax Main Reports"

# List all CSV files in the folder
csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

# Extract the project_name (file name without path) from csv_files
csv_files <- data.frame(
  project_id = basename(csv_files),  # Get only file name (without path)
  stringsAsFactors = FALSE
)

# Get the project_name values from the filtered_WT_projects dataframe
project_names_to_process <- filtered_WT_projects$project_name

# Initialize an empty dataframe to store the unique combinations
all_unique_lat_long_year <- data.frame()

# Columns to read from the CSV files
columns_to_read <- c("longitude", "latitude", "recording_date_time", "organization", "project_id")

# Loop through each CSV file
for (file in csv_files$project_id) {
  # Read the CSV file, only the relevant columns
  df <- read_csv(file.path(folder_path, file), col_types = cols_only(
    longitude = col_double(),
    latitude = col_double(),
    recording_date_time = col_character(),
    organization = col_character(),
    project_id = col_character()
  ))
  
  # Check if there are any parsing issues
  if (nrow(problems(df)) > 0) {
    message("Parsing issues found in file: ", file)
    print(problems(df))
  }
  
  # Ensure the necessary columns exist
  if (all(c("longitude", "latitude", "recording_date_time", "organization", "project_id") %in% colnames(df))) {
    # Filter rows based on the project_name values in filtered_WT_projects dataframe
    if (file %in% project_names_to_process) {
      # Extract the first 4 characters of recording_date_time to get the year
      df$year <- substr(df$recording_date_time, 1, 4)
      
      # Create a unique combination of longitude, latitude, and year
      unique_df <- df %>%
        select(longitude, latitude, year, organization, project_id) %>%
        distinct()  # This ensures uniqueness based on longitude, latitude, and year
      
      # Combine the unique rows with the final dataframe
      all_unique_lat_long_year <- bind_rows(all_unique_lat_long_year, unique_df)
    }
  }
}



### Where ARUs were deployed in multiple years, remove duplicate lat/long rows and add
### multiple years to the 'year'column
all_unique_lat_long <- all_unique_lat_long_year %>%
  group_by(latitude, longitude) %>%
  summarise(
    data_year = paste(unique(year), collapse = ", "),
    project_id = first(project_id),
    organization = first(organization),
    .groups = "drop"
  )

# Convert project_id in all_unique_lat_long to integer for consistency
all_unique_lat_long <- all_unique_lat_long %>%
  mutate(project_id = as.integer(project_id))

# Extract project names from my_projects and add them to all_unique_lat_long
all_unique_lat_long <- all_unique_lat_long %>%
  left_join(
    my_projects %>% select(project_id, project), # Select only relevant columns
    by = "project_id" # Match on project_id
  ) 




### Remove all rows in which the lat and long values are outside of the approx.
### boundary of Alberta

# Filter to include rows within Alberta's boundaries
AB_unique_lat_long <- all_unique_lat_long %>%
  filter(latitude >= 49.0, latitude <= 60.0, 
         longitude >= -120.0, longitude <= -110.0)

# Save AB_unique_lat_long

write.csv(AB_unique_lat_long, "Output/all_ARU_locations_AB.csv")



### Create a shapefile with all ARU loocation in AB and extract and attach
### other relevant attributes (e.g., forest age)

# Read in AB_unique_lat_long
AB_unique_lat_long <- read.csv("Output/all_ARU_locations_AB.csv")

# Convert the dataframe to a spatial object in WGS84
AB_unique_lat_long_WGS84 <- st_as_sf(
  AB_unique_lat_long,
  coords = c("longitude", "latitude"), # Specify longitude and latitude columns
  crs = 4326 # Set the CRS (WGS84)
)

# Transform coordinate system to NAD83 / Alberta 10-TM (Forest) (EPSG: 3400)
AB_unique_lat_long_sf <- st_transform(AB_unique_lat_long_WGS84, crs = 3400)

# Load in ABMI Wall-to-Wall Human Footprint Data (https://abmi.ca/abmi-home/data-resources/data-portal-main/data-portal/single-data-portal-detail.html?name=46)
ABMI_human_footprint <- "G:/My Drive/Thesis/Input/Spatial Data/ABMI/ABMI_Harvest_Areas_2021.shp"
ABMI_human_footprint_sf <- st_read(ABMI_human_footprint)

# Ensure both shapefiles use the same CRS
if (st_crs(ABMI_human_footprint_sf) != st_crs(AB_unique_lat_long_sf)) {
  AB_unique_lat_long_sf <- st_transform(AB_unique_lat_long_sf, st_crs(ABMI_human_footprint_sf))
}

# Perform intersection
intersection <- st_intersection(AB_unique_lat_long_sf, ABMI_human_footprint_sf)

# Extract the "Year" column and add it as a new column to AB_unique_lat_long_sf
AB_unique_lat_long_sf$harvest_yr <- NA  # Initialize new column
intersection_data <- intersection[, c("Year", "geometry")]

# Match the intersected geometries back to the original data
for (i in seq_len(nrow(intersection_data))) {
  idx <- which(st_equals(AB_unique_lat_long_sf, intersection_data[i, ]))
  if (length(idx) > 0) {
    AB_unique_lat_long_sf$harvest_yr[idx] <- intersection_data$Year[i]
  }
}

# Define the output file path
output_path <- "G:/My Drive/Thesis/Output/Spatial Data/all_AB_ARU_locations_harvest_yr/all_AB_ARU_locations_w_harvest_yr.shp"

# Normalize the file path for proper handling
output_path_normalized <- normalizePath(output_path, mustWork = FALSE)

# Export shapefile to output file path
st_write(AB_unique_lat_long_sf, output_path_normalized)




### Import relevant Alberta Vegetation Inventory shape files

# Load AVI BAM
avi_bam_location <- "G:/My Drive/Thesis/Input/Spatial Data/AVI_BAM/avi_bam.shp"
avi_bam_sf <- st_read(avi_bam_location)

# Load AVI Crown
avi_crown_location <- "G:/My Drive/Thesis/Input/Spatial Data/AVI_Crown/avi_crown.shp"
avi_crown_sf <- st_read(avi_crown_location)

# Load ABI PVLI
avi_plvi_location <- "G:/My Drive/Thesis/Input/Spatial Data/AVI_PLVI/plvi.shp"
avi_plvi_sf <- st_read(avi_plvi_location)

# Convert shapefile into dataframe to view attribute tables

avi_bam_df <- as.data.frame(avi_bam_sf)
avi_crown_df <- as.data.frame(avi_crown_sf)
avi_plvi_df <- as.data.frame(avi_plvi_sf)

# Confirm that all of the shapefiles have the same coordinate system

st_crs(AB_unique_lat_long_sf)
st_crs(avi_bam_sf)
st_crs(avi_crown_sf)
st_crs(avi_plvi_sf)
        


### Perform intersect function between ARU locations and AVI Crown polygons

# Perform a spatial join (left join) to keep all points in AB_unique_lat_long, even if they don't intersect
AVI_Crown_ARU <- st_join(AB_unique_lat_long_sf, avi_crown_sf, join = st_intersects, left = TRUE)

# Convert result to data frame
AVI_Crown_ARU_df <- as.data.frame(AVI_Crown_ARU)

# Filter out sites that did not intersect AVI_Crown 
AVI_Crown_ARU_df <- AVI_Crown_ARU_df[!is.na(AVI_Crown_ARU_df$SP1), ]



### For all remaining rows that did not join with AVI_Crown, run intersect function with AVI_BAM

# Filter out points that did not intersect with AVI_Crown (i.e., rows where SP1 is NA)
remaining_points <- AVI_Crown_ARU[is.na(AVI_Crown_ARU$SP1), ]

# Identify common columns between AB_unique_lat_long and remaining_points
common_columns <- intersect(colnames(AB_unique_lat_long), colnames(remaining_points))

# Subset remaining_points to include only the common columns
remaining_points <- remaining_points[, common_columns]

# Perform spatial join with AVI_BAM
AVI_BAM_ARU <- st_join(remaining_points, avi_bam_sf, join = st_intersects, left = TRUE)

# Convert to a dataframe
AVI_BAM_ARU_df <- as.data.frame(AVI_BAM_ARU)

# Filter out sites that did not intersect AVI_BAM
AVI_BAM_ARU_df <- AVI_BAM_ARU_df[!is.na(AVI_BAM_ARU_df$SP1), ]



### For all remaining rows that did not join with AVI_BAM, run intersect function with AVI_PLVI

# Filter out points that did not intersect with AVI_Crown (i.e., rows where SP1 is NA)
remaining_points2 <- AVI_BAM_ARU[is.na(AVI_BAM_ARU$SP1), ]

# Identify common columns between AB_unique_lat_long and remaining_points2
common_columns2 <- intersect(colnames(AB_unique_lat_long), colnames(remaining_points2))

# Subset remaining_points to include only the common columns
remaining_points2 <- remaining_points2[, common_columns2]

# Perform spatial join with AVI_BAM
AVI_PLVI_ARU <- st_join(remaining_points2, avi_plvi_sf, join = st_intersects, left = TRUE)

# Convert to a dataframe
AVI_PLVI_ARU_df <- as.data.frame(AVI_PLVI_ARU)

# Filter out sites that did not intersect AVI_BAM
AVI_PLVI_ARU_df <- AVI_PLVI_ARU_df[!is.na(AVI_PLVI_ARU_df$SITE_TYP1), ]



### Now filter all of the sites that intersected AVI polygons to only include those that have some
### combination of Sw and Aw as their Leading and secondary species

# Filter AVI_Crown_ARU to only include rows where SP1 and SP2 are either (Aw, Sw) or (Sw, Aw)
mw_AVI_Crown_ARU <- AVI_Crown_ARU_df %>%
  filter((SP1 == "Aw" & SP2 == "Sw") | (SP1 == "Sw" & SP2 == "Aw"))

# Remove unnecessary columns
mw_AVI_Crown_ARU <- mw_AVI_Crown_ARU %>%
  select(data_year,project_id,organization,project,geometry,HEIGHT,SP1,SP1_PER,SP2,SP2_PER,SP3,SP3_PER)

# Filter AVI_BAM to only include rows where SP1 and SP2 are either (Aw, Sw) or (Sw, Aw)
mw_AVI_BAM_ARU <- AVI_BAM_ARU_df %>%
  filter((SP1 == "Aw" & SP2 == "Sw") | (SP1 == "Sw" & SP2 == "Aw"))

# Remove unnecessary columns
mw_AVI_BAM_ARU <- mw_AVI_BAM_ARU %>%
  select(data_year,project_id,organization,project,geometry,HEIGHT,SP1,SP1_PER,SP2,SP2_PER,SP3,SP3_PER)

# Rename columns to ensure compatibility with shapefile constraints
mw_AVI_Crown_ARU <- mw_AVI_Crown_ARU %>%
  rename(
    proj_id = project_id,
    org = organization,
    proj = project,
    ht = HEIGHT,
    sp1 = SP1,
    sp1_per = SP1_PER,
    sp2 = SP2,
    sp2_per = SP2_PER,
    sp3 = SP3,
    sp3_per = SP3_PER,
    year = data_year
  )

# Rename columns to ensure compatibility with shapefile constraints
mw_AVI_BAM_ARU <- mw_AVI_BAM_ARU %>%
  rename(
    proj_id = project_id,
    org = organization,
    proj = project,
    ht = HEIGHT,
    sp1 = SP1,
    sp1_per = SP1_PER,
    sp2 = SP2,
    sp2_per = SP2_PER,
    sp3 = SP3,
    sp3_per = SP3_PER,
    year = data_year
  )





### Export to shapefile for use in GIS

# Convert mw_AVI_BAM_ARU dataframe to sf object (no need for st_as_sfc if geometry is already in sfc format)
sf_mw_AVI_BAM_ARU <- st_as_sf(mw_AVI_BAM_ARU)

# Convert mw_AVI_Crown_ARU dataframe to sf object (same reasoning)
sf_mw_AVI_Crown_ARU <- st_as_sf(mw_AVI_Crown_ARU)

# Specify the output location
output_path <- "G:/My Drive/Thesis/Output/Spatial Data/mw_AVI_locations"

# Remove all files associated with the shapefile
output_path <- "G:/My Drive/Thesis/Output/Spatial Data/mw_AVI_locations"
unlink(list.files(output_path, pattern = "mw_AVI_Crown_ARU.*", full.names = TRUE))

# Save both shapefiles, overwriting if they already exist
st_write(sf_mw_AVI_BAM_ARU, file.path(output_path, "mw_AVI_BAM_ARU_2.shp"), append = FALSE)
st_write(sf_mw_AVI_Crown_ARU, file.path(output_path, "mw_AVI_Crown_ARU_2.shp"), append = FALSE)











### Create a 150 m buffer around all ARU locations to then extract forest type






### Remove all points whose 150 m intersect some type of hard human disturbance



















































































