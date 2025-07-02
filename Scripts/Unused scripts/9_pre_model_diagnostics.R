# Load libraries


library(dplyr)
library(spdep)
library(ggplot2)
library(stats)
library(ggplot2)
library(dplyr)
library(fitdistrplus)


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


### Based on AIC, negative binomial is the best fit for all spp.




















###### Create bubble plot of response variables


# BTNW Plot (Size Only)
unique_sites_btnw <- NTEMS %>%
  group_by(x_AEP10TM, y_AEP10TM) %>%
  summarise(mean_BTNW = mean(BTNW, na.rm = TRUE), .groups = 'drop')

bubble_plot_btnw <- ggplot(unique_sites_btnw, aes(x = x_AEP10TM, y = y_AEP10TM, size = mean_BTNW)) +
  geom_point(alpha = 0.7, color = "blue") +
  scale_size_continuous(range = c(1, 10), name = "Mean BTNW Abundance") +
  labs(title = "Mean BTNW Abundance", x = "x_AEP10TM", y = "y_AEP10TM") +
  theme_minimal() +
  coord_equal()

print(bubble_plot_btnw)

# BBWA Plot (Size Only)
unique_sites_bbwa <- NTEMS %>%
  group_by(x_AEP10TM, y_AEP10TM) %>%
  summarise(mean_BBWA = mean(BBWA, na.rm = TRUE), .groups = 'drop')

bubble_plot_bbwa <- ggplot(unique_sites_bbwa, aes(x = x_AEP10TM, y = y_AEP10TM, size = mean_BBWA)) +
  geom_point(alpha = 0.7, color = "blue") +
  scale_size_continuous(range = c(1, 10), name = "Mean BBWA Abundance") +
  labs(title = "Mean BBWA Abundance", x = "x_AEP10TM", y = "y_AEP10TM") +
  theme_minimal() +
  coord_equal()

print(bubble_plot_bbwa)

# TEWA Plot (Size Only)
unique_sites_tewa <- NTEMS %>%
  group_by(x_AEP10TM, y_AEP10TM) %>%
  summarise(mean_TEWA = mean(TEWA, na.rm = TRUE), .groups = 'drop')

bubble_plot_tewa <- ggplot(unique_sites_tewa, aes(x = x_AEP10TM, y = y_AEP10TM, size = mean_TEWA)) +
  geom_point(alpha = 0.7, color = "blue") +
  scale_size_continuous(range = c(1, 10), name = "Mean TEWA Abundance") +
  labs(title = "Mean TEWA Abundance", x = "x_AEP10TM", y = "y_AEP10TM") +
  theme_minimal() +
  coord_equal()

print(bubble_plot_tewa)















############## kNN means approach to defining clusters

### NTEMS

# 1. Group and Summarize Data
NTEMS_grouped <- NTEMS %>%
  group_by(x_AEP10TM, y_AEP10TM) %>%
  summarise(
    mean_BTNW = mean(BTNW, na.rm = TRUE),
    mean_BBWA = mean(BBWA, na.rm = TRUE),
    mean_TEWA = mean(TEWA, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Prepare Coordinates from Grouped Data
coords <- NTEMS_grouped %>%
  dplyr::select(x_AEP10TM, y_AEP10TM)

# 3. Create Neighborhood List
nb <- dnearneigh(as.matrix(coords), d1 = 0, d2 = 30000, longlat = FALSE)

# 4. Identify Isolates
isolates <- which(lengths(nb) == 0)

# 5. Remove Isolates from Grouped Data
if (length(isolates) > 0) {
  coords_filtered <- coords[-isolates, ]
  NTEMS_grouped_filtered <- NTEMS_grouped %>%
    filter(!(paste(x_AEP10TM, y_AEP10TM) %in% paste(coords$x_AEP10TM[isolates], coords$y_AEP10TM[isolates])))
  
  cat("Removed", length(isolates), "isolates.\n")
} else {
  coords_filtered <- coords
  NTEMS_grouped_filtered <- NTEMS_grouped
  cat("No isolates found.\n")
}

# 6. Elbow Method with Filtered Grouped Data
wss <- numeric(15)
for (i in 1:15) {
  set.seed(123)
  kmeans_result <- kmeans(coords_filtered, centers = i)
  wss[i] <- kmeans_result$tot.withinss
}

# 7. Plot the Elbow Curve
elbow_plot <- ggplot(data.frame(k = 1:15, wss = wss), aes(x = k, y = wss)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow Method for Optimal k (Isolates Removed)", x = "Number of Clusters (k)", y = "Within-Cluster Sum of Squares (WSS)") +
  theme_minimal()

print(elbow_plot)

# 8. Perform k-means Clustering for k = 3 to 10
for (k in 3:10) {
  set.seed(123)
  kmeans_result <- kmeans(coords_filtered, centers = k)
  NTEMS_grouped_filtered <- NTEMS_grouped_filtered %>%
    mutate(!!paste0("cluster_k_", k) := kmeans_result$cluster[match(paste(x_AEP10TM, y_AEP10TM), paste(coords_filtered$x_AEP10TM, coords_filtered$y_AEP10TM))])
  
  # 9. Visualize Clusters (Optional)
  cluster_plot <- ggplot(NTEMS_grouped_filtered, aes(x = x_AEP10TM, y = y_AEP10TM, color = factor(.data[[paste0("cluster_k_", k)]]))) +
    geom_point() +
    labs(title = paste("k-means Clusters (k =", k, ")")) +
    theme_minimal() +
    coord_equal()
  
  print(cluster_plot)
}

# 10. Join Cluster Assignments Back to Original Data
NTEMS <- NTEMS %>%
  left_join(NTEMS_grouped_filtered %>% dplyr::select(x_AEP10TM, y_AEP10TM, starts_with("cluster_k_")), by = c("x_AEP10TM", "y_AEP10TM"))

LULC <- LULC %>%
  left_join(NTEMS_grouped_filtered %>% dplyr::select(x_AEP10TM, y_AEP10TM, starts_with("cluster_k_")), by = c("x_AEP10TM", "y_AEP10TM"))

# 11. Save

write.csv(NTEMS, "Output/Tabular Data/NTEMS_w_clustering.csv")
write.csv(LULC, "Output/Tabular Data/LULC_w_clustering.csv")














############## Moran's I approach to clustering

###### Run lm with count as response and survey effort as predictor 
###### to get the residuals (i.e., the count without the effect of
###### survey effort). This is needed to calculate Moran's I

# Load necessary libraries
library(lme4)   # For fitting GLMM
library(dbscan) # For clustering
library(spdep)  # For spatial analysis
library(dplyr)  # For data manipulation

# Run linear models for all response variables
model_BTNW <- lm(BTNW ~ log(survey_effort), data = NTEMS)
model_BBWA <- lm(BBWA ~ log(survey_effort), data = NTEMS)
model_TEMA <- lm(TEWA ~ log(survey_effort), data = NTEMS)

# Extract residuals from each model
NTEMS$BTNW_residuals <- residuals(model_BTNW)
NTEMS$BBWA_residuals <- residuals(model_BBWA)
NTEMS$TEWA_residuals <- residuals(model_TEMA)

# Create df with residuals and x,y coordinates
response_residuals <- NTEMS %>%
  dplyr::select(x_AEP10TM, y_AEP10TM, BTNW_residuals, BBWA_residuals, TEWA_residuals)

# Define the list of distance thresholds
thresholds <- c(50000)

# Define the response residual variables
response_vars <- c("BTNW_residuals", "TEWA_residuals", "BBWA_residuals")

# Create the output directory for saving plots
dir.create("Figures/Temp", recursive = TRUE, showWarnings = FALSE)

# Loop over each residual response variable
for (response in response_vars) {
  
  cat("Evaluating spatial autocorrelation for:", response, "\n")
  
  # Ensure we're working with the response_residuals dataframe
  unique_sites <- response_residuals %>%
    group_by(x_AEP10TM, y_AEP10TM) %>%
    summarise(mean_value = mean(get(response), na.rm = TRUE), .groups = 'drop')
  
  # Extract coordinates (x_AEP10TM and y_AEP10TM)
  coords_matrix <- as.matrix(unique_sites[, c("x_AEP10TM", "y_AEP10TM")])
  
  # Loop over each threshold distance
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
        
        # Calculate Moran's I for the residuals
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
    png(paste0("Figures/Temp/spatial_autocorrelation_", response, "_threshold_", threshold, ".png"))
    
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

# Define the sequence of maximum distances (5000 to 30000 in increments of 5000)
max_distances <- seq(5000, 30000, by = 5000)

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





























# Save
write.csv(NTEMS, "Output/Tabular Data/NTEMS_w_clustering.csv")
