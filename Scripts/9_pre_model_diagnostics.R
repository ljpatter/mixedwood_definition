# Load libraries

library(car)   
library(caret) 
library(dplyr)
library(spdep)
library(tidyverse)
library(sf)
library(fitdistrplus)
library(sp)
library(lme4)
library(MASS) 
library(gstat)
library(igraph)
library(terra)
library(ape)


# Clear environment
rm(list = ls())  # Removes all objects from the environment

# Load data
NTEMS <- read.csv("Output/Tabular Data/NTEMS_with_bird_data.csv")
LULC <- read.csv("Output/Tabular Data/LULC_with_bird_data.csv")



############# Assess correlation between covariates for NTEMS and LULC

# select covariates
df1<-NTEMS%>%
  dplyr::select(prop_con_1, prop_con_2, prop_con_3, clumpy_1, clumpy_2, clumpy_3,
                age_mn_1, age_mn_2, age_mn_3)%>%
  lapply(as.numeric)%>%
  as.data.frame()

df2<-LULC%>%
  dplyr::select(prop_con_1, prop_con_2, prop_con_3, clumpy_1, clumpy_2, clumpy_3,
                age_mn_1, age_mn_2, age_mn_3)%>%
  lapply(as.numeric)%>%
  as.data.frame()

# Compute correlation matrix
correlation_matrix_NTEMS <- cor(df1)
correlation_matrix_LULC <- cor(df2)

# Print correlation matrix
print(correlation_matrix_NTEMS)
print(correlation_matrix_LULC)





######### Calculate VIF between predictors

# Function to calculate VIF for a given dataset, response variable, and predictor suffixes
calculate_vif <- function(data, response_vars, suffixes) {
  for (response in response_vars) {
    for (suffix in suffixes) {
      formula_str <- paste(response, "~ prop_con_", suffix, "+ clumpy_", suffix, "+ age_mn_", suffix, sep = "") # Corrected line
      model <- lm(as.formula(formula_str), data = data)
      vif_result <- vif(model)
      print(paste("VIF for", response, "with suffix", suffix, ":"))
      print(vif_result)
    }
  }
}

# Define response variables and predictor suffixes
response_vars <- c("BTNW", "BBWA", "TEWA")
suffixes <- 1:3

# Calculate VIF for NTEMS and LULC
print("VIF for NTEMS:")
calculate_vif(NTEMS, response_vars, suffixes)

print("\nVIF for LULC:")
calculate_vif(LULC, response_vars, suffixes)







######### Evaluate distribution of response variable

# Response variables to check
response_vars <- c("BTNW", "TEWA", "BBWA")

for (response in response_vars) {
  print(paste("Evaluating distribution for:", response))
  
  # Histogram
  hist(NTEMS[[response]],
       main = paste("Histogram of", response),
       xlab = response,
       ylab = "Frequency",
       col = "lightblue",
       border = "black")
  
  # Distribution plots
  plotdist(NTEMS[[response]], discrete = TRUE, histo = TRUE, demp = TRUE)
  descdist(NTEMS[[response]], discrete = TRUE, boot = 1000)
  
  # Fit distributions
  fit_p <- fitdist(NTEMS[[response]], "pois")
  fit_nb <- fitdist(NTEMS[[response]], "nbinom")
  fit_norm <- fitdist(NTEMS[[response]], "norm")
  
  # Comparison plots
  par(mfrow = c(1, 1))
  plot.legend <- c("Poisson", "Negative binomial", "normal")
  denscomp(list(fit_p, fit_nb, fit_norm), legendtext = plot.legend)
  cdfcomp(list(fit_p, fit_nb, fit_norm), legendtext = plot.legend)
  qqcomp(list(fit_p, fit_nb, fit_norm), legendtext = plot.legend)
  
  # Goodness-of-fit statistics
  print("Goodness-of-fit for Normal:")
  print(gofstat(fit_norm))
  print("Goodness-of-fit for Poisson:")
  print(gofstat(fit_p))
  print("Goodness-of-fit for Negative Binomial:")
  print(gofstat(fit_nb))
  
  # Add a separator for clarity
  cat("\n----------------------------------------\n\n") #Added for better readability
}

# Calculate mean to variance ratio
BTNW_mean_var_ratio <- var(NTEMS$BTNW) / mean(NTEMS$BTNW)
TEWA_mean_var_ratio <- var(NTEMS$TEWA) / mean(NTEMS$TEWA)
BBWA_mean_var_ratio <- var(NTEMS$BBWA) / mean(NTEMS$BBWA)

## All spp. show moderate to slight overdispersion







###### Run lm with count as response and survey effort as predictor 
###### to get the residuals (i.e., the count without the effect of
###### survey effort). This is needed to calculate Moran's I

BTNW_residuals <- 
BBWA_residuals <- 
TEWA_residuals <-






# Calculate Moran's for three different thresholds

library(lme4)  # For fitting GLMM
library(dbscan) # For clustering
library(spdep)  # For spatial analysis

# Define the list of thresholds
thresholds <- c(50000)

# Define the response variables
response_vars <- c("BTNW", "TEWA", "BBWA")

# Create the output directory for saving plots
dir.create("Figures/Temp", recursive = TRUE, showWarnings = FALSE)

# Loop over each response variable (species)
for (response in response_vars) {
  
  cat("Evaluating spatial autocorrelation for:", response, "\n")
  
  # Filter out duplicate sites and calculate the mean for each site
  unique_sites <- NTEMS %>%
    group_by(x_AEP10TM, y_AEP10TM) %>%
    summarise(mean_value = mean(get(response), na.rm = TRUE), .groups = 'drop')
  
  # Extract coordinates (x_AEP10TM and y_AEP10TM)
  coords_matrix <- as.matrix(unique_sites[, c("x_AEP10TM", "y_AEP10TM")])
  
  # Loop over each threshold
  for (threshold in thresholds) {
    
    cat("Using threshold:", threshold, "meters\n")
    
    # Define a sequence of distances (from 100 to the current threshold)
    distances <- seq(100, threshold, by = 1000)
    
    # Initialize a vector to store Moran's I values
    morans_i_values <- numeric(length(distances))
    
    # Loop over the distances and calculate Moran's I
    for (i in seq_along(distances)) {
      nb <- tryCatch(dnearneigh(coords_matrix, 0, distances[i], longlat = FALSE), 
                     error = function(e) NULL)
      
      if (!is.null(nb) && length(nb) > 0) {
        # Create spatial weights matrix
        listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
        
        # Calculate Moran's I for the specified response variable
        moran_test_result <- tryCatch({
          moran.test(unique_sites$mean_value, listw, zero.policy = TRUE)
        }, error = function(e) NULL)
        
        if (!is.null(moran_test_result)) {
          morans_i_values[i] <- moran_test_result$estimate[1]
        } else {
          morans_i_values[i] <- NA
        }
      } else {
        morans_i_values[i] <- NA
      }
    }
    
    # Find the flattening point (where Moran's I approaches 0)
    flattening_index <- which.min(abs(diff(morans_i_values, na.rm = TRUE)))
    flattening_distance <- distances[flattening_index]
    flattening_moran_value <- morans_i_values[flattening_index]
    
    cat("Flattening distance:", flattening_distance, "meters\n")
    cat("Moran's I at flattening distance:", flattening_moran_value, "\n")
    
    # Save the plot for this threshold and response variable
    png(paste("Figures/Temp/spatial_autocorrelation_", response, "_threshold_", threshold, ".png", sep = ""))
    
    # Plot Moran's I values against distance
    plot(distances, morans_i_values, type = "b", pch = 16, col = "blue",
         xlab = "Distance (meters)", ylab = "Moran's I",
         main = paste("Moran's I vs. Distance Threshold\n", response, "Threshold:", threshold, "m"))
    
    abline(h = 0, col = "red", lwd = 2)
    abline(v = flattening_distance, col = "blue", lwd = 2)
    
    # Annotate the flattening point
    text_string <- paste("Flatten's at:", round(flattening_distance, 2), "\nMoran's I:", round(flattening_moran_value, 2))
    text(flattening_distance, max(morans_i_values, na.rm = TRUE) * 0.8, text_string, pos = 2, col = "blue")
    
    dev.off()
  }
}

cat("Analysis complete. Plots saved in 'Figures/Temp'.\n")


















############# Group locations based on spatial proximity with multiple distances

# Load point count data
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/filtered_point_count_loc_2009.shp")

# Extract the WKT string from your desired CRS
target_crs_wkt <- st_as_text(st_crs(AB_point_counts_sf))

# Create the target CRS object using terra::crs()
target_crs_terra <- terra::crs(target_crs_wkt)

# Remove point count file from environment
rm(list = ls()[startsWith(ls(), "AB_point_counts_sf")])

# Assuming NTEMS is already loaded as shown in your output
NTEMS_sf <- st_as_sf(NTEMS, coords = c("x_AEP10TM", "y_AEP10TM"), crs = target_crs_terra)

# Extract the coordinates
coords <- st_coordinates(NTEMS_sf$geometry)

# Convert the coordinates to a matrix
coords_matrix <- as.matrix(coords[, c("X", "Y")])

# Define the sequence of maximum distances (5000 to 15000 in increments of 1000)
max_distances <- seq(1000, 15000, by = 1000)

# Create a copy of the NTEMS dataframe
NTEMS_grouped <- NTEMS

# Loop through each maximum distance
for (max_dist in max_distances) {
  
  # Find nearest neighbors
  nb <- dnearneigh(coords_matrix, 0, max_dist)
  
  # Create a spatial weights matrix
  listw <- nb2listw(nb, style = "B", zero.policy = TRUE)
  
  # Convert the listw object to a binary matrix
  binary_matrix <- listw2mat(listw)
  
  # Compute the dissimilarity matrix
  dissimilarity_matrix <- 1 - binary_matrix
  
  # Perform hierarchical clustering
  hc <- hclust(as.dist(dissimilarity_matrix))
  
  # Cut the dendrogram to create groups
  # Adjust the h parameter to get the desired number of groups
  spatial_group <- cutree(hc, h = 0.5)
  
  # Add the spatial groups to the dataframe with unique column names
  col_name <- paste0("spatial_group_", max_dist)
  NTEMS_grouped[[col_name]] <- spatial_group
  
  # Print the number of unique groups for the current distance
  print(paste("Number of unique groups for distance", max_dist, ":", length(unique(spatial_group))))
}

# Rename output back to NTEMS
NTEMS <- NTEMS_grouped

# Save
write.csv(NTEMS, "Output/Tabular Data/NTEMS_w_clustering.csv")

cat("Clustering complete. Data saved in 'NTEMS_w_clustering.csv'.\n")













############# Calculate inverse distance weighting between sites

# Step 1: Create a data frame with coordinates only
NTEMS_coords_BTNW <- NTEMS %>%
  dplyr::select(x_AEP10TM, y_AEP10TM, BTNW)

# Step 2: Convert to SpatialPointsDataFrame
coordinates(NTEMS_coords_BTNW) <- ~x_AEP10TM+y_AEP10TM
proj4string(NTEMS_coords_BTNW) <- CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-115 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Step 3: Create a distance matrix (10km threshold)
dist_matrix <- spDists(NTEMS_coords_BTNW)
dist_threshold <- 10000 # 10 km threshold
binary_weights <- ifelse(dist_matrix <= dist_threshold & dist_matrix > 0, 1, 0)

# Step 4: Convert to spatial weights list
spatial_weights_list <- mat2listw(binary_weights, style = "W", zero.policy = TRUE)

# Step 5: Calculate the sum of spatial weights for each site
NTEMS_coords_BTNW$spatial_weight <- rowSums(binary_weights)

# Convert the SpatialPointsDataFrame back to a regular data frame
NTEMS_coords_BTNW <- as.data.frame(NTEMS_coords_BTNW)

# Step 6: One-to-many join to reattach weights
NTEMS <- NTEMS %>%
  left_join(NTEMS_coords_BTNW[, c("x_AEP10TM", "y_AEP10TM", "spatial_weight")], 
            by = c("x_AEP10TM", "y_AEP10TM"))

# Add small constant to spatial_weight

NTEMS$spatial_weight <- NTEMS$spatial_weight + 1e-5


# Save
write.csv(NTEMS, "Output/Tabular Data/NTEMS_w_clustering.csv")


