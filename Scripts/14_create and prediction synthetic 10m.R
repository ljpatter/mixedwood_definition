######### 01_Create_Simulated_Landscapes.R #########

# Load libraries
library(sf)
library(terra)
library(landscapemetrics) # Keep this for lsm_c_clumpy as it was working

# --- DIAGNOSTIC: Check if terra is loaded correctly ---
if (!"package:terra" %in% search()) {
  stop("terra package not loaded. Please ensure library(terra) runs successfully.")
}

# Define township center coordinates (change as needed)
township_center <- c(x = 528551, y = 6000000) # Example UTM coordinates

# Define township size (10 km x 10 km)
township_size <- 10000 # in meters

# Define half_size globally or pass it to matrix_to_raster
half_size <- township_size / 2

# Set CRS
target_crs <- "EPSG:3400"

# Create a square polygon around the center
township_polygon <- st_polygon(list(rbind(
  c(township_center["x"] - half_size, township_center["y"] - half_size),
  c(township_center["x"] + half_size, township_center["y"] - half_size),
  c(township_center["x"] + half_size, township_center["y"] + half_size),
  c(township_center["x"] - half_size, township_center["y"] + half_size),
  c(township_center["x"] - half_size, township_center["y"] - half_size)
))) %>% st_sfc(crs = target_crs) # Use EPSG code here

# Convert to sf object
township_sf <- st_sf(geometry = township_polygon)

# Create output directory if it doesn't exist
if (!dir.exists("Output/Spatial Data/Township polygon")) {
  dir.create("Output/Spatial Data/Township polygon", recursive = TRUE)
}

# Export to shapefile
st_write(township_sf, "Output/Spatial Data/Township polygon/township_polygon.shp", delete_layer = TRUE)


# --- Parameters ---
cell_size <- 10 # 10m resolution for base landscapes
township_length_m <- 10000
ncol <- nrow <- township_length_m / cell_size

block_size_m <- 447
block_size_pixels <- round(block_size_m / cell_size)
blocks_per_side <- ceiling(ncol / block_size_pixels)
# --- END Parameters ---

# Helper: convert matrix to SpatRaster with CRS and correct extent
matrix_to_raster <- function(mat, crs_string, center_x, center_y, total_size) {
  half_total_size <- total_size / 2
  
  xmin_val <- center_x - half_total_size
  xmax_val <- center_x + half_total_size
  ymin_val <- center_y - half_total_size
  ymax_val <- center_y + half_total_size 
  
  rast_obj <- rast(nrows = nrow(mat), ncols = ncol(mat),
                   xmin = xmin_val, xmax = xmax_val,
                   ymin = ymin_val, ymax = ymax_val, 
                   crs = crs_string)
  values(rast_obj) <- as.vector(t(mat))
  return(rast_obj)
}

# Generate intimate landscape: fully random 0/1
generate_intimate <- function() {
  matrix(sample(c(0,1), size = nrow * ncol, replace = TRUE), nrow = nrow, ncol = ncol)
}

# Generate one block with half conifer pixels clustered on a random side
generate_half_half_block <- function(size) {
  block <- matrix(0, nrow = size, ncol = size)
  side <- sample(c("top", "bottom", "left", "right"), 1)
  
  if (side %in% c("top", "bottom")) {
    rows_to_fill <- floor(size / 2)
    if (side == "top") {
      block[1:rows_to_fill, ] <- 1
    } else {
      block[(size - rows_to_fill + 1):size, ] <- 1
    }
  } else {
    cols_to_fill <- floor(size / 2)
    if (side == "left") {
      block[, 1:cols_to_fill] <- 1
    } else {
      block[, (size - cols_to_fill + 1):size] <- 1
    }
  }
  return(block)
}

# Generate segregated landscape: stitch blocks with half-half pattern on random side
generate_segregated <- function() {
  full_mat <- matrix(0, nrow = nrow, ncol = ncol)
  for (row_block in 0:(blocks_per_side - 1)) {
    for (col_block in 0:(blocks_per_side - 1)) {
      block <- generate_half_half_block(block_size_pixels)
      
      r_start <- row_block * block_size_pixels + 1
      r_end <- min((row_block + 1) * block_size_pixels, nrow)
      c_start <- col_block * block_size_pixels + 1
      c_end <- min((col_block + 1) * block_size_pixels, ncol)
      
      block_rows <- r_end - r_start + 1
      block_cols <- c_end - c_start + 1
      
      if (block_rows != block_size_pixels || block_cols != block_size_pixels) {
        block <- block[1:block_rows, 1:block_cols]
      }
      
      full_mat[r_start:r_end, c_start:c_end] <- block
    }
  }
  return(full_mat)
}

# Generate base matrices (0=deciduous, 1=conifer)
intimate_mat <- generate_intimate()
segregated_mat <- generate_segregated()

# Convert to SpatRaster and assign CRS
intimate_raster_base <- matrix_to_raster(intimate_mat, target_crs,
                                         township_center["x"], township_center["y"],
                                         township_size)
segregated_raster_base <- matrix_to_raster(segregated_mat, target_crs,
                                           township_center["x"], township_center["y"],
                                           township_size)

# Ensure integer values for the base rasters
values(intimate_raster_base) <- as.integer(values(intimate_raster_base))
values(segregated_raster_base) <- as.integer(values(segregated_raster_base))

message("Finished generating base 10m LULC rasters.")

# Define a black and white color palette (for plotting the base layer)
bw_palette <- c("white", "black")

# Set plot margins to provide space for axis labels (removed axis-related elements)
old_par <- par(no.readonly = TRUE) # Save current par settings
par(mar = c(2, 2, 4, 2) + 0.1) # Adjusted margins for plots without axes

# Plot Intimate Mixedwood (displaying the base 0/1 layer) without axes
plot(intimate_raster_base, main = "Intimate Mixedwood (10m Base)", col = bw_palette, legend = FALSE,
     axes = FALSE) # axes set to FALSE
# Removed axis() calls
# Removed box() call

# Plot Segregated Mixedwood (displaying the base 0/1 layer) without axes
plot(segregated_raster_base, main = "Segregated Mixedwood (10m Base)", col = bw_palette, legend = FALSE,
     axes = FALSE) # axes set to FALSE
# Removed axis() calls
# Removed box() call

# Restore original par settings
par(old_par)

# --- Calculate and Print Landscape Metrics ---
message("\nCalculating landscape metrics:")

# Calculate CLUMPY (keeping this as it worked)
intimate_clumpy <- lsm_c_clumpy(intimate_raster_base)
segregated_clumpy <- lsm_c_clumpy(segregated_raster_base)

message(paste0("Intimate Landscape CLUMPY: ", round(intimate_clumpy$value[intimate_clumpy$class == 1], 3)))
message(paste0("Segregated Landscape CLUMPY: ", round(segregated_clumpy$value[segregated_clumpy$class == 1], 3)))

# Calculate Proportion of Conifer (class 1) using terra directly
# Get all values from the raster
intimate_values <- values(intimate_raster_base)
segregated_values <- values(segregated_raster_base)

# Count the number of pixels equal to 1 (conifer)
intimate_conifer_pixels <- sum(intimate_values == 1, na.rm = TRUE)
segregated_conifer_pixels <- sum(segregated_values == 1, na.rm = TRUE)

# Get the total number of pixels
total_pixels_intimate <- ncell(intimate_raster_base)
total_pixels_segregated <- ncell(segregated_raster_base)

# Calculate proportion
intimate_prop_conifer_calc <- intimate_conifer_pixels / total_pixels_intimate
segregated_prop_conifer_calc <- segregated_conifer_pixels / total_pixels_segregated

message(paste0("Intimate Landscape Proportion Conifer (class 1) [terra]: ", round(intimate_prop_conifer_calc, 3)))
message(paste0("Segregated Landscape Proportion Conifer (class 1) [terra]: ", round(segregated_prop_conifer_calc, 3)))


# --- Output and Save ---

# Create output directory if it doesn't exist
if (!dir.exists("Output/Spatial Data/Simulated landscapes")) {
  dir.create("Output/Spatial Data/Simulated landscapes", recursive = TRUE)
}

# Save as TIFF (single band)
intimate_path <- "Output/Spatial Data/Simulated landscapes/intimate_mixedwood_10m.tif"
segregated_path <- "Output/Spatial Data/Simulated landscapes/segregated_mixedwood_10m.tif"
writeRaster(intimate_raster_base, filename = intimate_path, overwrite = TRUE, datatype = "INT1U") # Use INT1U for 0/1 integer
writeRaster(segregated_raster_base, filename = segregated_path, overwrite = TRUE, datatype = "INT1U")

message(paste("Saved single-band 10m intimate raster to:", intimate_path))
message(paste("Saved single-band 10m segregated raster to:", segregated_path))

# --- Save Plots as PNG without axes ---

# Create output directory for plots if it doesn't exist
if (!dir.exists("Output/Plots")) {
  dir.create("Output/Plots", recursive = TRUE)
}

# Save Intimate Mixedwood plot as PNG
png("Output/Plots/intimate_mixedwood_10m_base.png", width = 800, height = 800, res = 100)
par(mar = c(2, 2, 4, 2) + 0.1) # Adjusted margins for plots without axes
plot(intimate_raster_base, main = "Intimate Mixedwood (10m Base)", col = bw_palette, legend = FALSE,
     axes = FALSE) # axes set to FALSE
dev.off()
par(old_par)

# Save Segregated Mixedwood plot as PNG
png("Output/Plots/segregated_mixedwood_10m_base.png", width = 800, height = 800, res = 100)
par(mar = c(2, 2, 4, 2) + 0.1) # Adjusted margins for plots without axes
plot(segregated_raster_base, main = "Segregated Mixedwood (10m Base)", col = bw_palette, legend = FALSE,
     axes = FALSE) # axes set to FALSE
dev.off()
par(old_par)
















######### 02_Predict_Abundance_Isolated_Effects.R #########

# Load necessary libraries
library(terra)
library(sf)
library(stringr)
library(ggplot2)
library(tools) # For file_path_sans_ext
library(dplyr) # For bind_rows and data manipulation
library(tidyr) # For pivot_longer
library(gbm)
library(scales) # ADDED: Required for non-scientific notation on plot axes

# --- 1. Set the folder path where your TIFF files are stored ---
# Make sure this path is correct for your system
tif_folder <- "Output/Spatial Data/Prediction landscapes" # This is for saving outputs, not inputs

# Set CRS
target_crs <- "EPSG:3400"

# Load the BASE 0/1 intimate and segregated rasters
# These are the single-band rasters generated by 01_Create_Simulated_Landscapes.R
intimate_raster_base <- rast("Output/Spatial Data/Simulated landscapes/intimate_mixedwood_10m.tif")
segregated_raster_base <- rast("Output/Spatial Data/Simulated landscapes/segregated_mixedwood_10m.tif")

# Ensure base rasters have same extent and resolution (should be from creation script)
if (!compareGeom(intimate_raster_base, segregated_raster_base)) {
  stop("Base intimate and segregated rasters do not have matching geometries.")
}

# Calculate pixel area once (assuming 10m x 10m)
pixel_area_sq_m <- prod(res(intimate_raster_base)) # e.g., 10m x 10m = 100 sq m

# --- Create a base dataframe of just X, Y coordinates for all pixels ---
# We'll use this as the spatial reference for all predictions, since LULC won't be spatial anymore.
xy_coords_df <- as.data.frame(intimate_raster_base, xy = TRUE, cell = FALSE) %>%
  select(x_AEP10TM = x, y_AEP10TM = y) # Rename columns to match model predictor names

# --- DEBUGGING: Check if TIFF files are found (for existing outputs, not for generation) ---
tif_files <- list.files(tif_folder, pattern = "\\.tif$", full.names = TRUE)
print(paste("DEBUG: Contents of tif_files (for existing prediction outputs, if any):", paste(basename(tif_files), collapse = ", ")))
if (length(tif_files) == 0) {
  message("No TIFF files found in the specified folder for *existing* prediction outputs. This is expected if you are running predictions for the first time.")
} else {
  message(paste("Found", length(tif_files), "TIFF files in the prediction output folder."))
}

# === Load model object ===
final_model_list <- readRDS("Output/Tabular Data/final_model_list_parallel_poisson.rds")

# Define spatial extents (radii for focal calculations)
extents <- c(150, 500, 1000)

# Define age classes for prediction
age_classes <- c(30, 50, 70)

# --- HELPER FUNCTION: Get means/modes of all predictors from training data ---
get_predictor_stats <- function(model) {
  train_data_full <- model$gbm.call$dataframe
  predictor_stats <- list(means = list(), modes = list())
  for (col_name in model$gbm.call$predictor.names) { # Iterate through ALL model predictors
    col_data <- train_data_full[[col_name]]
    if (is.numeric(col_data)) {
      predictor_stats$means[[col_name]] <- mean(col_data, na.rm = TRUE)
    } else if (is.factor(col_data) || is.character(col_data)) {
      mode_val <- names(sort(table(col_data), decreasing = TRUE))[1]
      predictor_stats$modes[[col_name]] <- mode_val
    }
  }
  # Also get the mean of the offset for prediction, as it's part of the model's linear predictor
  if ("offset" %in% names(train_data_full) && is.numeric(train_data_full$offset)) {
    predictor_stats$means$offset <- mean(train_data_full$offset, na.rm = TRUE)
  } else {
    warning("Offset column not found or not numeric in training data. Assuming mean offset of 0 for prediction.")
    predictor_stats$means$offset <- 0 # Default if offset is missing
  }
  return(predictor_stats)
}

# --- REVISED predict_synthetic HELPER FUNCTION for Isolated Effects ---
# This function now takes a pre-generated dataframe of X,Y coordinates
# and the layout type to hardcode the focal metrics.
predict_synthetic <- function(model, xy_coords_df, current_extent, current_age_val, predictor_stats, pixel_area_sq_m, landscape_layout) {
  all_model_predictor_names <- model$gbm.call$predictor.names
  num_pixels <- nrow(xy_coords_df)
  
  # Create a list to build columns efficiently
  cols_for_prediction <- list(
    x_AEP10TM = xy_coords_df$x_AEP10TM,
    y_AEP10TM = xy_coords_df$y_AEP10TM
  )
  
  # Determine active predictor names for current extent
  prop_con_col_name <- paste0("prop_con_", current_extent, "_LULC")
  clumpy_col_name <- paste0("clumpy_", current_extent, "_LULC")
  age_mn_col_name <- paste0("age_mn_", current_extent, "_NTEMS")
  
  # Populate columns based on model predictor names and hardcoded values
  for (v_name in all_model_predictor_names) {
    if (v_name %in% names(cols_for_prediction)) { # Skip if already added (like x, y)
      next
    }
    
    if (v_name == prop_con_col_name) {
      cols_for_prediction[[v_name]] <- 0.5 # Hardcode prop_con to 0.5 for all
    } else if (v_name == clumpy_col_name) {
      cols_for_prediction[[v_name]] <- ifelse(landscape_layout == "Intimate", 0, 1) # Hardcode clumpy based on layout
    } else if (v_name == age_mn_col_name) {
      cols_for_prediction[[v_name]] <- current_age_val # Forced age value
    } else if (v_name == "offset") {
      cols_for_prediction[[v_name]] <- predictor_stats$means$offset
    } else if (v_name == "year") {
      if ("year" %in% names(predictor_stats$modes)) {
        cols_for_prediction[[v_name]] <- factor(rep(predictor_stats$modes[["year"]], num_pixels), levels = levels(model$gbm.call$dataframe[["year"]]))
      } else {
        warning("Year is a factor in model but mode not found. Setting to first level.")
        cols_for_prediction[[v_name]] <- factor(rep(levels(model$gbm.call$dataframe[["year"]])[1], num_pixels), levels = levels(model$gbm.call$dataframe[["year"]]))
      }
    } else {
      # All other predictors: assign mean/mode from training data (vectorized)
      if (v_name %in% names(predictor_stats$means)) {
        cols_for_prediction[[v_name]] <- predictor_stats$means[[v_name]]
      } else if (v_name %in% names(predictor_stats$modes)) {
        mode_val <- predictor_stats$modes[[v_name]]
        levels_from_model <- levels(model$gbm.call$dataframe[[v_name]])
        cols_for_prediction[[v_name]] <- factor(rep(mode_val, num_pixels), levels = levels_from_model)
      } else {
        warning(paste("Predictor", v_name, "mean/mode not found in predictor_stats. Setting to 0 or first factor level."))
        if (is.numeric(model$gbm.call$dataframe[[v_name]])) {
          cols_for_prediction[[v_name]] <- 0
        } else if (is.factor(model$gbm.call$dataframe[[v_name]])) {
          cols_for_prediction[[v_name]] <- factor(rep(levels(model$gbm.call$dataframe[[v_name]])[1], num_pixels), levels = levels(model$gbm.call$dataframe[[v_name]]))
        } else {
          cols_for_prediction[[v_name]] <- NA
        }
      }
    }
  }
  
  # Convert list of columns to data frame
  df_for_prediction <- as.data.frame(cols_for_prediction)
  
  # Final check for missing columns before prediction
  missing_cols <- setdiff(all_model_predictor_names, names(df_for_prediction))
  if (length(missing_cols) > 0) {
    stop(paste("Following required columns are still missing after predict_synthetic processing:", paste(missing_cols, collapse=", ")))
  }
  
  # Ensure 'year' is a factor for prediction - this check might be redundant if handled above correctly
  if ("year" %in% all_model_predictor_names && !is.factor(df_for_prediction$year)) {
    df_for_prediction$year <- factor(df_for_prediction$year, levels = levels(model$gbm.call$dataframe$year))
  }
  
  # Get the mean offset for prediction (this is handled once and passed to predict.gbm)
  mean_offset_for_prediction <- predictor_stats$means$offset
  
  # Make initial prediction (this gives you predicted birds PER HECTARE)
  df_for_prediction$pred_raw <- predict.gbm(
    model,
    newdata = df_for_prediction,
    n.trees = model$gbm.call$best.trees,
    type = "response",
    offset = mean_offset_for_prediction
  )
  
  # --- CRITICAL FIX START ---
  # Scale prediction from "birds per hectare" to "birds per pixel (100 sq m)"
  # pred_raw is in "birds/hectare".
  # A pixel (10m x 10m) has an area of pixel_area_sq_m (which is 100 sq m).
  # 1 hectare = 10,000 sq m.
  # So, a pixel area in hectares = pixel_area_sq_m / 10000.
  df_for_prediction$pred <- df_for_prediction$pred_raw * (pixel_area_sq_m / 10000)
  # --- CRITICAL FIX END ---
  
  return(df_for_prediction)
}

# --- Run predictions ---
species_list <- c("BBWA", "BTNW", "TEWA")

all_species_preds <- list() # To collect dataframes for combined_pred_df

for (sp in species_list) {
  cat(paste0("\nProcessing predictions for species: ", sp, "\n"))
  model <- final_model_list[[sp]]
  
  current_model_predictor_stats <- get_predictor_stats(model)
  
  # Loop through all defined extents (150, 500, 1000)
  for (extent_val in extents) {
    for (age_val in age_classes) {
      cat(paste0("  - Extent: ", extent_val, "m, Age: ", age_val, " years\n"))
      
      # Predict for Intimate landscape (prop_con = 0.5, clumpy = 0)
      pred_int <- predict_synthetic(
        model, xy_coords_df, extent_val, age_val,
        current_model_predictor_stats, pixel_area_sq_m, landscape_layout = "Intimate"
      ) %>%
        mutate(layout = "Intimate", extent = extent_val, Species = sp, age_class = age_val)
      
      # Predict for Segregated landscape (prop_con = 0.5, clumpy = 1)
      pred_seg <- predict_synthetic(
        model, xy_coords_df, extent_val, age_val,
        current_model_predictor_stats, pixel_area_sq_m, landscape_layout = "Segregated"
      ) %>%
        mutate(layout = "Segregated", extent = extent_val, Species = sp, age_class = age_val)
      
      # Combine results for this species, extent, and age
      all_species_preds[[paste(sp, "int", extent_val, age_val, sep = "_")]] <- pred_int
      all_species_preds[[paste(sp, "seg", extent_val, age_val, sep = "_")]] <- pred_seg
    }
  }
}

# Combine all predictions into one dataframe
combined_pred_df <- bind_rows(all_species_preds)

# Calculate total abundance per species × layout × extent × age_class
total_abundance_df <- combined_pred_df %>%
  group_by(Species, layout, extent, age_class) %>%
  summarise(total_abundance = sum(pred, na.rm = TRUE), .groups = "drop") %>%
  mutate(total_abundance = if_else(is.na(total_abundance), 0, total_abundance))

# View summary table
print(total_abundance_df)

# --- Generate Publication-Quality Bar Chart ---

# Determine max y-value for plot limit with some padding
max_y_val <- max(total_abundance_df$total_abundance, na.rm = TRUE)
y_limit_padded <- max_y_val * 1.1 # 10% padding

bar_chart_plot <- ggplot(total_abundance_df, aes(x = factor(age_class), y = total_abundance, fill = layout)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", alpha = 0.8) +
  facet_grid(Species ~ extent, scales = "free_y") +
  scale_fill_viridis_d(name = "Landscape\nLayout", option = "plasma") +
  scale_y_continuous(labels = scales::comma, limits = c(0, y_limit_padded)) + # ADDED/FIXED: Apply non-scientific notation and set limits
  labs(
    title = "Predicted Avian Abundance by Landscape Layout, Extent, and Forest Age (Isolated Effects)",
    x = "Forest Age Class (Years)",
    y = "Total Predicted Abundance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )

print(bar_chart_plot)

# --- Save Publication-Quality Bar Chart (Optional) ---
output_dir_plots <- "Output/Plots"
if (!dir.exists(output_dir_plots)) {
  dir.create(output_dir_plots, recursive = TRUE)
}
ggsave(paste0(output_dir_plots, "/predicted_abundance_by_layout_age_isolated.png"),
       plot = bar_chart_plot, width = 12, height = 8, dpi = 300)
ggsave(paste0(output_dir_plots, "/predicted_abundance_by_layout_age_isolated.tiff"),
       plot = bar_chart_plot, width = 12, height = 8, dpi = 300)

cat("All predictions summarized and publication-quality bar chart generated and saved successfully!\n")

# This is still good for checking any new warnings generated by this run.
# lifecycle::last_lifecycle_warnings()










