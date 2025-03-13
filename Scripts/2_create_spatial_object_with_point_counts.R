# ---
# title: "Create spatial object with point count locations "
# author: "Brendan Casey and Leonard Patterson"
# created: "2024-02-06"
# description: "Creating shapefiles of point count locations"
# ---



##Load packages----

library(sf) # for creating spatial objects
library(tidyverse)

##Import data----
load("Output/wildtrax_cleaned_2025-03-12.rData")

# Set target crs

AB_boundary_10TM <- st_read("Input/Spatial Data/Alberta/AB_boundary.shp")
target_crs <- st_crs(AB_boundary_10TM)

# Create dataframe with x,y and site ID
ss_xy <- wildtrax_cleaned%>%
  ungroup()%>%
  dplyr::select(location, year, lon, lat, x_AEP10TM, y_AEP10TM)%>%
  distinct()

# Create spatial object in correct PCS

ss_xy_10TM<-ss_xy%>%
  st_as_sf(coords=c("x_AEP10TM", "y_AEP10TM"), crs=target_crs, remove=FALSE)

# Ensure correct GRS and PCS

st_crs(ss_xy_10TM)



### Limit point count to only those within AB

# Keep only the points that are inside the inner boundary
AB_point_counts <- st_intersection(ss_xy_10TM, AB_boundary_10TM)

# Filter out unnecessary columns

AB_point_counts <- AB_point_counts %>%
  dplyr::select(location, year, lon, lat, x_AEP10TM, y_AEP10TM, geometry)

# Save the filtered points 
st_write(AB_point_counts, "Output/Spatial Data/AB_point_counts/AB_point_counts_10TM.shp", append=FALSE)
save(AB_point_counts, file=paste0("Output/AB_point_counts", Sys.Date(), ".rData"))

# Save 
save(ss_xy_10TM, file=paste0("Output/AB_point_counts", Sys.Date(), ".rData"))
st_write(ss_xy_10TM, "Output/Spatial Data/AB_point_counts/ss_xy_10TM.shp", append=FALSE)






