# ---
# title: "Mixed regression models"
# author: "Leonard Patterson"
# created: ""
# description: ""
# ---


# Clear environment
rm(list = ls())  # Removes all objects from the environment

#Setup ----

##Load packages----

# Setup ----

## Load packages ----
library(dplyr)

# Load data
NTEMS <- read.csv("Output/Tabular Data/NTEMS_w_clustering.csv")
LULC  <- read.csv("Output/Tabular Data/LULC_w_clustering.csv")

########## Merge prop_con_, clumpy_, and age_mn_ from LULC to NTEMS

# Remove the 'X' and 'X.1' columns from both NTEMS and LULC before renaming
NTEMS_cleaned <- NTEMS %>% dplyr::select(-c(X, X.1))
LULC_cleaned <- LULC %>% dplyr::select(-c(X, X.1))

# Rename _1, _2, and _3 to _150, _500, and _1000
scale_suffix_map <- c("_1" = "_150", "_2" = "_500", "_3" = "_1000")

rename_scale_suffixes <- function(names) {
  for (old in names(scale_suffix_map)) {
    names <- gsub(old, scale_suffix_map[[old]], names)
  }
  names
}

prefixes <- c("prop_con_", "clumpy_", "age_mn_", "prop_dec_")

NTEMS_converted <- NTEMS_cleaned %>%
  rename_with(rename_scale_suffixes, .cols = matches(paste0("^(", paste(prefixes, collapse = "|"), ")")))

LULC_converted <- LULC_cleaned %>%
  rename_with(rename_scale_suffixes, .cols = matches(paste0("^(", paste(prefixes, collapse = "|"), ")")))

NTEMS_renamed <- NTEMS_converted %>%
  rename_with(.cols = starts_with("prop_con_"), ~ paste0(., "_NTEMS")) %>%
  rename_with(.cols = starts_with("prop_dec_"), ~ paste0(., "_NTEMS")) %>%
  rename_with(.cols = starts_with("clumpy_"), ~ paste0(., "_NTEMS")) %>%
  rename_with(.cols = starts_with("age_mn_"), ~ paste0(., "_NTEMS"))

LULC_renamed <- LULC_converted %>%
  rename_with(.cols = starts_with("prop_con_"), ~ paste0(., "_LULC")) %>%
  rename_with(.cols = starts_with("prop_dec_"), ~ paste0(., "_LULC")) %>%
  rename_with(.cols = starts_with("clumpy_"), ~ paste0(., "_LULC")) %>%
  rename_with(.cols = starts_with("age_mn_"), ~ paste0(., "_LULC"))

# Join the renamed data frames
joined_data <- left_join(NTEMS_renamed, LULC_renamed,
                         by = c("project", "location", "survey_type", "year", "ordinalDay", "hssr",
                                "survey_duration_method", "max_dist_band", "x_AEP10TM", "y_AEP10TM"),
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
species_cols <- c("BTNW", "TEWA", "BBWA")  # Add more if needed

if (check_match(joined_data)) {
  joined_data <- joined_data %>%
    dplyr::select(-ends_with(".y"), -starts_with("prop_dec_")) %>%
    rename_with(~ gsub("\\.x$", "", .))
  
  joined_data <- joined_data %>%
    dplyr::select(
      project, location, survey_type, year, ordinalDay, hssr, survey_effort,  # <- added
      survey_duration_method, max_dist_band,
      x_AEP10TM, y_AEP10TM,
      all_of(species_cols),
      starts_with("prop_con_"), starts_with("clumpy_"), starts_with("age_mn_"), cluster_k_10
    )

  
  print("Join and replacement successful.")
} else {
  print("Mismatches found. Join and replacement aborted.")
}



rm(list = c("LULC", "LULC_cleaned", "LULC_converted", "LULC_renamed"), "NTEMS", "NTEMS_cleaned", "NTEMS_converted", "NTEMS_renamed")














kable(top5_filtered, caption = "Top 5 Friedman H Interactions per Species (Strict Filter)")

final_model_list[["BBWA"]]$gbm.call$predictor.names
final_model_list[["TEWA"]]$gbm.call$predictor.names
final_model_list[["BTNW"]]$gbm.call$predictor.names












##### Centre and scale predictors prior to assessing functional response

# Define all predictor columns to scale
predictors_to_scale <- c(
  "prop_con_150_NTEMS", "prop_con_500_NTEMS", "prop_con_1000_NTEMS",
  "prop_con_150_LULC",  "prop_con_500_LULC",  "prop_con_1000_LULC",
  "clumpy_150_NTEMS",   "clumpy_500_NTEMS",   "clumpy_1000_NTEMS",
  "clumpy_150_LULC",    "clumpy_500_LULC",    "clumpy_1000_LULC",
  "age_mn_150_NTEMS",   "age_mn_500_NTEMS",   "age_mn_1000_NTEMS",
  "age_mn_150_LULC",    "age_mn_500_LULC",    "age_mn_1000_LULC",
  "ordinalDay"
)

# Loop and scale each one
for (predictor in predictors_to_scale) {
  scaled_name <- paste0(predictor, "_scaled")
  joined_data[[scaled_name]] <- as.vector(scale(joined_data[[predictor]]))
}

# Convert hssr to proportion of survey day (0 to 1)
joined_data$hssr_prop <- (joined_data$hssr - min(joined_data$hssr)) / 
  (max(joined_data$hssr) - min(joined_data$hssr))

# Then scale the hssr proportion
joined_data$hssr_scaled <- as.vector(scale(joined_data$hssr_prop))

# Then scale the hssr proportion
joined_data$survey_effort <- as.numeric(as.character(joined_data$survey_effort))
joined_data$log_effort <- log(joined_data$survey_effort)

# Convert survey effort to numeric
joined_data$survey_effort <- as.numeric(as.character(joined_data$survey_effort))

# Ensure year is a factor
joined_data$year <- as.factor(joined_data$year)













############ Calculate distance threshold for use in BRT

library(spdep)
library(ggplot2)
library(dplyr)

# STEP 1: Collapse to unique locations by averaging BTNW across visits
unique_sites <- joined_data %>%
  group_by(location, x_AEP10TM, y_AEP10TM) %>%
  summarise(BTNW_mean = mean(BTNW, na.rm = TRUE), .groups = "drop")

# STEP 2: Fit a simple model to generate predicted probabilities
library(xgboost)
X <- as.matrix(unique_sites[, c("x_AEP10TM", "y_AEP10TM")])  # just spatial predictors for simplicity
y <- ifelse(unique_sites$BTNW_mean > 0, 1, 0)  # convert to binary

dtrain <- xgb.DMatrix(data = X, label = y)
params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eval_metric = "logloss",
  eta = 0.01,
  max_depth = 3
)
model <- xgb.train(params = params, data = dtrain, nrounds = 100)
preds <- predict(model, dtrain)

# STEP 3: Compute residuals
residuals <- y - preds

# STEP 4: Coordinates
coords <- cbind(unique_sites$x_AEP10TM, unique_sites$y_AEP10TM)

# STEP 5: Loop through distance thresholds and compute Moran's I
dist_thresholds <- seq(500, 30000, by = 500)
moran_vals <- numeric(length(dist_thresholds))

for (i in seq_along(dist_thresholds)) {
  dist <- dist_thresholds[i]
  nb <- dnearneigh(coords, 0, dist)
  
  if (any(card(nb) > 0)) {
    lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    moran_test <- moran.test(residuals, lw, randomisation = TRUE, zero.policy = TRUE)
    moran_vals[i] <- moran_test$estimate["Moran I statistic"]
  } else {
    moran_vals[i] <- NA
  }
}

# STEP 6: Plot Moran's I vs distance
ggplot(data.frame(Distance = dist_thresholds, Moran_I = moran_vals),
       aes(x = Distance, y = Moran_I)) +
  geom_line(color = "steelblue") +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Moran's I vs Distance Threshold (BTNW, non-thinned data)",
    x = "Distance (m)", y = "Moran's I"
  )


var(joined_sf$TEWA)/mean(joined_sf$TEWA)
var(joined_sf$BBWA)/mean(joined_sf$BBWA)
var(joined_sf$BTNW)/mean(joined_sf$BTNW)






################## BRT w/ QPAD


##############################################
# QPAD-Corrected Boosted Regression Tree Models
# Species: BTNW, TEWA, BBWA
# Goal: Model species occurrence using QPAD-based detection offsets
# Author: [Your Name]
# Date: [Today’s Date]
##############################################

# Load required packages
library(dismo)       # gbm.step for BRT
library(blockCV)     # spatial blocks and folds
library(sp)          # spatial data
library(raster)      # raster prediction
library(gbm)         # base BRT functions
library(dplyr)       # data manipulation
library(ggplot2)     # plotting

# --- INPUT: joined_data must be loaded into workspace ---
# joined_data must contain: lat, lon, survey_effort (minutes), BTNW, TEWA, BBWA

# STEP 1: Create binary presence/absence columns if missing
species_list <- c("BTNW", "TEWA", "BBWA")

for (sp in species_list) {
  bin_col <- paste0(sp, "_bin")
  if (!(bin_col %in% names(joined_data))) {
    joined_data[[bin_col]] <- ifelse(joined_data[[sp]] > 0, 1, 0)
  }
}

# STEP 2: Define QPAD-like detectability offset function
# We'll use a simple linear logistic model of detectability vs survey length:
#   p = logistic(alpha + beta * survey_effort)
# Assume alpha = -1.5, beta = 0.12 as a reasonable default shape
# These would ideally be species-specific from QPAD but are approximated here

alpha <- -1.5
beta <- 0.12

logit_to_prob <- function(x) {
  exp(x) / (1 + exp(x))
}

# Compute detection probability and offset for each species
for (sp in species_list) {
  p <- logit_to_prob(alpha + beta * joined_data$survey_effort)
  offset <- -log(p)
  joined_data[[paste0(sp, "_offset")]] <- offset
}

# STEP 3: Convert to SpatialPointsDataFrame
coords <- joined_data[, c("lon", "lat")]
spdf <- SpatialPointsDataFrame(coords = coords, data = joined_data,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))
spdf_lcc <- spTransform(spdf, CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +datum=NAD83 +units=m +no_defs"))

# STEP 4: Define a function to model each species
run_species_brt <- function(sp, predictors, pred_raster, output_dir) {
  cat("\n====== Running BRT for", sp, "======\n")
  
  # Prep data
  joined_data$ABUND <- joined_data[[paste0(sp, "_bin")]]
  joined_data$logoffset <- joined_data[[paste0(sp, "_offset")]]
  
  datcombo <- joined_data
  datcombo$X <- datcombo$lon
  datcombo$Y <- datcombo$lat
  
  datcombo_sp <- SpatialPointsDataFrame(coords = datcombo[, c("lon", "lat")],
                                        data = datcombo,
                                        proj4string = CRS("+proj=longlat +datum=WGS84"))
  datcombo_sp <- spTransform(datcombo_sp, CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +datum=NAD83 +units=m +no_defs"))
  
  # Remove rows with NA in predictor columns
  datcombo <- datcombo[complete.cases(datcombo[, predictors]), ]
  datcombo_sp <- datcombo_sp[datcombo$PKEY %in% datcombo$PKEY, ]
  
  # Calculate autocorrelation range
  auto <- spatialAutoRange(pred_raster[[predictors]], doParallel = FALSE)
  
  # Spatial blocking
  set.seed(123)
  blocks <- spatialBlock(speciesData = datcombo_sp,
                         species = "ABUND",
                         rasterLayer = pred_raster[[1]],
                         theRange = auto$range,
                         k = 5,
                         selection = "random",
                         iteration = 250)
  
  # Run BRT model
  set.seed(123)
  brt <- gbm.step(data = datcombo,
                  gbm.x = predictors,
                  gbm.y = "ABUND",
                  family = "bernoulli",
                  fold.vector = blocks$foldID,
                  tree.complexity = 3,
                  learning.rate = 0.005,
                  bag.fraction = 0.5,
                  offset = datcombo$logoffset,
                  site.weights = rep(1, nrow(datcombo)), # can replace with species-specific weights
                  keep.fold.models = TRUE,
                  keep.fold.fit = TRUE,
                  verbose = TRUE)
  
  # Save outputs
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  save(brt, file = file.path(output_dir, paste0(sp, "_brt.RData")))
  save(blocks, file = file.path(output_dir, paste0(sp, "_blocks.RData")))
  
  # Save variable importance
  varimp <- summary(brt)
  write.csv(varimp, file = file.path(output_dir, paste0(sp, "_varimp.csv")), row.names = FALSE)
  
  # Predict to raster
  pred <- predict(pred_raster, brt, n.trees = brt$n.trees, type = "response")
  writeRaster(pred, filename = file.path(output_dir, paste0(sp, "_pred.tif")),
              format = "GTiff", overwrite = TRUE)
  
  # Plot partials
  pdf(file.path(output_dir, paste0(sp, "_partial.pdf")))
  gbm.plot(brt, plot.layout = c(3, 3), common.scale = FALSE, write.title = FALSE)
  dev.off()
  
  cat("✅ Finished BRT for", sp, "\n")
}

# STEP 5: Define predictors and run models
# Replace with actual predictor names from your raster
predictors <- c("Structure_Stand_Age_v1", "roads_NS", "urbag_NS")  # EXAMPLE ONLY

# Load your raster stack (must contain predictor layers)
pred_raster <- brick("0_data/1_processed/prediction dataset/abs2011_250m.grd")

# Run BRT for each species
for (sp in species_list) {
  run_species_brt(sp, predictors, pred_raster, output_dir = paste0("2_BRT_outputs/QPAD_model_", sp))
}











######################## FULL RAC BRT WORKFLOW (no year, RAC, AUC, Moran’s I, PDPs, Interactions) ########################

library(xgboost)
library(pROC)
library(dplyr)
library(spdep)
library(ggplot2)
library(pdp)
library(data.table)
library(purrr)
library(tidyr)

# Create binary columns if not present
species_list <- c("BTNW", "TEWA", "BBWA")
for (sp in species_list) {
  bin_col <- paste0(sp, "_bin")
  if (!(bin_col %in% names(joined_data))) {
    joined_data[[bin_col]] <- ifelse(joined_data[[sp]] > 0, 1, 0)
  }
}

# Predictor groups
prop_con <- grep("prop_con_.*scaled", names(joined_data), value = TRUE)
age_mn   <- grep("age_mn_.*scaled", names(joined_data), value = TRUE)
clumpy   <- grep("clumpy_.*scaled", names(joined_data), value = TRUE)

# Combine feature set (removed 'year')
features <- c(prop_con, age_mn, clumpy)
base_features <- c(features, "hssr_scaled", "ordinalDay_scaled", "x_AEP10TM", "y_AEP10TM")

# Create spatial weights matrix (0–5000 m)
coords <- cbind(joined_data$x_AEP10TM, joined_data$y_AEP10TM)
nb <- dnearneigh(coords, 0, 5000)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Containers
results_list <- list()
thresholds_list <- list()
pdp_list <- list()
interaction_list <- list()
auc_list <- list()
moran_list <- list()

# Loop ----
for (sp in species_list) {
  cat("\n====== Modeling", sp, "with RAC + offset ======\n")
  
  y <- joined_data[[paste0(sp, "_bin")]]
  offset <- log(joined_data$survey_effort)
  
  # Step 1: Initial model (no RAC)
  X_base <- model.matrix(~ . - 1, data = joined_data[, base_features])
  dtrain_base <- xgb.DMatrix(data = X_base, label = y, base_margin = offset)
  
  params <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = "logloss",
    eta = 0.005,
    max_depth = 3,
    subsample = 0.5,
    colsample_bytree = 1
  )
  
  set.seed(123)
  cv_base <- xgb.cv(params = params, data = dtrain_base, nrounds = 1000,
                    nfold = 10, early_stopping_rounds = 20, verbose = 0)
  
  best_nrounds_base <- cv_base$best_iteration
  model_base <- xgb.train(params = params, data = dtrain_base,
                          nrounds = best_nrounds_base, verbose = 0)
  
  preds_base <- predict(model_base, dtrain_base)
  residuals_base <- y - preds_base
  auc_base <- auc(y, preds_base)
  moran_base <- moran.test(residuals_base, lw, zero.policy = TRUE)$estimate[1]
  
  # Step 2: Compute RAC
  rac <- lag.listw(lw, residuals_base, zero.policy = TRUE)
  rac_name <- paste0(sp, "_RAC")
  joined_data[[rac_name]] <- rac
  
  # Step 3: Final model with RAC
  X_rac <- model.matrix(~ . - 1, data = joined_data[, c(base_features, rac_name)])
  dtrain_rac <- xgb.DMatrix(data = X_rac, label = y, base_margin = offset)
  
  cv_rac <- xgb.cv(params = params, data = dtrain_rac, nrounds = 1000,
                   nfold = 10, early_stopping_rounds = 20, verbose = 0)
  
  best_nrounds_rac <- cv_rac$best_iteration
  model_rac <- xgb.train(params = params, data = dtrain_rac,
                         nrounds = best_nrounds_rac, verbose = 0)
  
  final_preds <- predict(model_rac, dtrain_rac)
  joined_data[[paste0(sp, "_pred_prob")]] <- final_preds
  residuals_rac <- y - final_preds
  auc_rac <- auc(y, final_preds)
  moran_rac <- moran.test(residuals_rac, lw, zero.policy = TRUE)$estimate[1]
  
  auc_list[[sp]] <- c(AUC_base = auc_base, AUC_rac = auc_rac)
  moran_list[[sp]] <- c(Moran_base = moran_base, Moran_rac = moran_rac)
  
  # Variable importance
  imp <- xgb.importance(model = model_rac) %>% as_tibble()
  nuisance_vars <- c("hssr_scaled", "ordinalDay_scaled", "x_AEP10TM", "y_AEP10TM", rac_name)
  imp <- imp[!imp$Feature %in% nuisance_vars, ]
  imp$Species <- sp
  results_list[[sp]] <- imp
  
  # PDPs
  top_prop_con <- imp %>% filter(Feature %in% prop_con) %>% arrange(desc(Gain)) %>% slice_head(n = 1)
  top_age_mn   <- imp %>% filter(Feature %in% age_mn)   %>% arrange(desc(Gain)) %>% slice_head(n = 1)
  top_clumpy   <- imp %>% filter(Feature %in% clumpy)   %>% arrange(desc(Gain)) %>% slice_head(n = 1)
  top_vars <- bind_rows(top_prop_con, top_age_mn, top_clumpy)
  
  pdp_vals <- list()
  for (var in top_vars$Feature) {
    pdp_obj <- partial(model_rac, pred.var = var, train = as.data.frame(X_rac), plot = FALSE)
    pdp_vals[[var]] <- pdp_obj
  }
  pdp_list[[sp]] <- pdp_vals
  
  # Thresholds
  roc_obj <- roc(y, final_preds)
  coords_df <- coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"))
  coords_df$diff <- abs(coords_df$sensitivity - coords_df$specificity)
  best_idx <- which.min(coords_df$diff)
  thresholds_list[[sp]] <- list(
    threshold = coords_df$threshold[best_idx],
    sensitivity = coords_df$sensitivity[best_idx],
    specificity = coords_df$specificity[best_idx]
  )
  
  # Interactions
  top_vars_10 <- imp$Feature[1:10]
  interaction_scores <- expand.grid(var1 = top_vars_10, var2 = top_vars_10, stringsAsFactors = FALSE) %>%
    filter(var1 != var2) %>%
    mutate(Interaction = purrr::map2_dbl(var1, var2, ~{
      gain_int <- tryCatch({
        xgb.interaction(model_rac, data = X_rac, feature_names = colnames(X_rac)) %>%
          filter(Feature1 == .x & Feature2 == .y) %>% pull(Gain)
      }, error = function(e) NA)
      if (length(gain_int) == 0) NA else gain_int
    }))
  
  interaction_list[[sp]] <- interaction_scores %>%
    arrange(desc(Interaction)) %>%
    filter(!is.na(Interaction)) %>%
    head(5)
}

# Combine importance
imp_df <- bind_rows(results_list)

# Plot variable importance for each species
ggplot(imp_df, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  facet_wrap(~Species, scales = "free_y") +
  labs(title = "Variable Importance (RAC-adjusted)", x = "Predictor", y = "Relative Influence") +
  theme_minimal(base_size = 12)

# Print AUCs and Moran's I
print("AUCs:")
print(auc_list)

print("Residual Spatial Autocorrelation (Moran's I):")
print(moran_list)

# Print thresholds
print(thresholds_list)

# PDPs
for (sp in names(pdp_list)) {
  cat("\n=== Partial Dependence for", sp, "===\n")
  for (var in names(pdp_list[[sp]])) {
    plot(pdp_list[[sp]][[var]], main = paste(sp, "-", var), xlab = var)
  }
}

# Interactions
for (sp in names(interaction_list)) {
  cat("\n=== Top Interactions for", sp, "===\n")
  print(interaction_list[[sp]])
}

















################# Spatially thin point count data


# Load required libraries
library(sf)
library(dplyr)
library(terra)
library(xgboost)
library(spdep)

# Load ARU point data 
AB_point_counts_sf <- st_read("Output/Spatial Data/AB_point_counts/filtered_point_count_loc_2009.shp")

# Load Alberta boundary shapefile
AB_boundary_10TM <- st_read("Input/Spatial Data/Alberta/AB_boundary.shp")

# Extract and define the correct CRS
target_crs_wkt <- st_as_text(st_crs(AB_point_counts_sf))
target_crs_terra <- terra::crs(target_crs_wkt)

# Convert joined_data to an sf object
joined_data_sf <- st_as_sf(joined_data, coords = c("x_AEP10TM", "y_AEP10TM"), crs = target_crs_terra)

# Create a 1000m grid over Alberta
grid_1000m <- st_make_grid(AB_boundary_10TM, cellsize = 1000, what = "polygons")
grid_sf <- st_sf(grid_id = 1:length(grid_1000m), geometry = grid_1000m)

# Assign ARUs to grid cells
joined_data_sf <- st_join(joined_data_sf, grid_sf, left = FALSE)

# Randomly select one ARU per grid cell
set.seed(123)
selected_ARUs <- joined_data_sf %>%
  group_by(grid_id) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Create a unique ID for matching
joined_data_sf$site_year <- paste0(joined_data_sf$location, "_", joined_data_sf$year)
selected_ARUs$site_year <- paste0(selected_ARUs$location, "_", selected_ARUs$year)

# Retain remaining ARUs for validation using site_year
validation_ARUs <- joined_data_sf %>%
  filter(!paste0(location, "_", year) %in% selected_ARUs$site_year)

# Save shapefiles (optional)
# st_write(selected_ARUs, "Output/Spatial Data/Filtered_ARUs_1000m.shp", delete_layer = TRUE)
# st_write(validation_ARUs, "Output/Spatial Data/Validation_ARUs_1000m.shp", delete_layer = TRUE)

# Convert back to dataframes with coordinates restored
joined_data_thinned <- as.data.frame(selected_ARUs) %>%
  mutate(x_AEP10TM = st_coordinates(selected_ARUs)[, 1],
         y_AEP10TM = st_coordinates(selected_ARUs)[, 2]) %>%
  dplyr::select(-geometry)

joined_data_validation <- as.data.frame(validation_ARUs) %>%
  mutate(x_AEP10TM = st_coordinates(validation_ARUs)[, 1],
         y_AEP10TM = st_coordinates(validation_ARUs)[, 2]) %>%
  dplyr::select(-geometry)

cat("Total ARUs retained:", nrow(joined_data_thinned), "\n")
cat("Total ARUs in validation dataset:", nrow(joined_data_validation), "\n")

# Save

write.csv(joined_data_thinned, "Output/Tabular Data/joined_data_thinned.csv")


































######################## BRT with Residual Autocovariate (RAC) ########################

library(xgboost)
library(pROC)
library(dplyr)
library(fastDummies)
library(spdep)

# Ensure cluster and year are factors
joined_data_thinned$cluster_k_10 <- as.factor(joined_data_thinned$cluster_k_10)
joined_data_thinned$year <- as.factor(joined_data_thinned$year)

# Create dummy variables
joined_data_thinned <- fastDummies::dummy_cols(
  joined_data_thinned,
  select_columns = c("cluster_k_10", "year"),
  remove_selected_columns = TRUE
)

# Define predictors (scaled), excluding nuisance
features <- c(
  "prop_con_150_NTEMS_scaled", "prop_con_500_NTEMS_scaled", "prop_con_1000_NTEMS_scaled",
  "prop_con_150_LULC_scaled", "prop_con_500_LULC_scaled", "prop_con_1000_LULC_scaled",
  "clumpy_150_NTEMS_scaled", "clumpy_500_NTEMS_scaled", "clumpy_1000_NTEMS_scaled",
  "clumpy_150_LULC_scaled", "clumpy_500_LULC_scaled", "clumpy_1000_LULC_scaled",
  "age_mn_150_NTEMS_scaled", "age_mn_500_NTEMS_scaled", "age_mn_1000_NTEMS_scaled",
  "age_mn_150_LULC_scaled", "age_mn_500_LULC_scaled", "age_mn_1000_LULC_scaled"
)

cluster_dummies <- grep("^cluster_k_10_", names(joined_data_thinned), value = TRUE)
year_dummies <- grep("^year_", names(joined_data_thinned), value = TRUE)

# Base full feature set
base_features <- c(
  features,
  "hssr_scaled", "ordinalDay_scaled", "x_AEP10TM", "y_AEP10TM",
  cluster_dummies,
  year_dummies
)

species_list <- c("BTNW", "TEWA", "BBWA")
results_list <- list()
thresholds_list <- list()

# Build spatial weights using a 0–5000 m distance band
coords <- cbind(joined_data_thinned$x_AEP10TM, joined_data_thinned$y_AEP10TM)
nb <- dnearneigh(coords, 0, 5000)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

### Model loop with RAC
for (sp in species_list) {
  cat("\n====== Modeling", sp, "with RAC ======\n")
  
  joined_data_thinned[[paste0(sp, "_bin")]] <- ifelse(joined_data_thinned[[sp]] > 0, 1, 0)
  y <- joined_data_thinned[[paste0(sp, "_bin")]]
  
  # Fit initial model without RAC
  X_base <- as.matrix(joined_data_thinned[, base_features])
  dtrain_base <- xgb.DMatrix(data = X_base, label = y)
  
  params <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = "logloss",
    eta = 0.005,
    max_depth = 3,
    subsample = 0.5,
    colsample_bytree = 1
  )
  
  set.seed(123)
  cv_base <- xgb.cv(
    params = params,
    data = dtrain_base,
    nrounds = 1000,
    nfold = 10,
    early_stopping_rounds = 20,
    verbose = 0
  )
  
  best_nrounds_base <- cv_base$best_iteration
  model_base <- xgb.train(params = params, data = dtrain_base, nrounds = best_nrounds_base, verbose = 0)
  preds_base <- predict(model_base, dtrain_base)
  
  # Compute residuals and RAC
  residuals <- y - preds_base
  rac <- lag.listw(lw, residuals, zero.policy = TRUE)
  rac_name <- paste0(sp, "_RAC")
  joined_data_thinned[[rac_name]] <- rac
  
  # Refit with RAC
  X_rac <- as.matrix(joined_data_thinned[, c(base_features, rac_name)])
  dtrain_rac <- xgb.DMatrix(data = X_rac, label = y)
  
  cv_rac <- xgb.cv(
    params = params,
    data = dtrain_rac,
    nrounds = 1000,
    nfold = 10,
    early_stopping_rounds = 20,
    verbose = 0
  )
  
  best_nrounds_rac <- cv_rac$best_iteration
  model_rac <- xgb.train(params = params, data = dtrain_rac, nrounds = best_nrounds_rac, verbose = 0)
  final_preds <- predict(model_rac, dtrain_rac)
  
  joined_data_thinned[[paste0(sp, "_pred_prob")]] <- final_preds
  
  # Variable importance
  imp <- xgb.importance(model = model_rac)
  nuisance_vars <- c("hssr_scaled", "ordinalDay_scaled", "x_AEP10TM", "y_AEP10TM", rac_name)
  imp <- imp[!imp$Feature %in% nuisance_vars, ]
  imp$Species <- sp
  results_list[[sp]] <- imp
  
  # Optimal threshold
  roc_obj <- roc(y, final_preds)
  coords_df <- coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"))
  coords_df$diff <- abs(coords_df$sensitivity - coords_df$specificity)
  coords_df <- coords_df[!is.na(coords_df$diff), ]
  best_idx <- which.min(coords_df$diff)
  
  thresholds_list[[sp]] <- list(
    threshold = coords_df$threshold[best_idx],
    sensitivity = coords_df$sensitivity[best_idx],
    specificity = coords_df$specificity[best_idx]
  )
}

# Plot variable importance
imp_df <- do.call(rbind, results_list)
ggplot(imp_df, aes(x = reorder(Feature, Gain), y = Gain, fill = Species)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Variable Importance (RAC-adjusted)", x = "Predictor", y = "Relative Influence") +
  theme_minimal()

# Show thresholds
print(thresholds_list)





















########## Moran's I to identify distance threshold


library(spdep)
library(ggplot2)

# Example distances to try (in meters)
dist_thresholds <- seq(500, 30000, by = 500)
moran_vals <- numeric(length(dist_thresholds))

for (i in seq_along(dist_thresholds)) {
  dist <- dist_thresholds[i]
  nb <- dnearneigh(coords, 0, dist)
  
  if (any(card(nb) > 0)) {
    lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    moran_test <- moran.test(residuals, lw, randomisation = TRUE, zero.policy = TRUE)
    moran_vals[i] <- moran_test$estimate["Moran I statistic"]
  } else {
    moran_vals[i] <- NA
  }
}

# Plot
ggplot(data.frame(Distance = dist_thresholds, Moran_I = moran_vals),
       aes(x = Distance, y = Moran_I)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Moran's I vs Distance Threshold",
       x = "Distance (m)", y = "Moran's I")

















######validation
install.packages("pROC")
install.packages("PRROC")
# Load required packages
library(xgboost)
library(pROC)
library(PRROC)
library(ggplot2)
library(sf)
library(dplyr)

# Step 1: Prepare validation data
features <- c("prop_con_1_scaled", "prop_con_2_scaled", "prop_con_3_scaled",   
              "clumpy_1_scaled", "clumpy_2_scaled", "clumpy_3_scaled",         
              "age_mn_1_scaled", "age_mn_2_scaled", "age_mn_3_scaled",         
              "hssr_scaled", "ordinalDay_scaled", "x_AEP10TM", "y_AEP10TM")

X_val_bin <- as.matrix(NTEMS_validation[, features])
y_val_bin <- ifelse(NTEMS_validation$BTNW > 0, 1, 0)

# Step 2: Predict probability using trained model
dval_bin <- xgb.DMatrix(data = X_val_bin)
NTEMS_validation$pred_prob <- predict(xgb_bin_model, dval_bin)

# Step 3: Predict binary outcome with threshold (e.g. 0.3)
NTEMS_validation$pred_bin <- ifelse(NTEMS_validation$pred_prob >= 0.3, 1, 0)

# Step 4: Confusion matrix
conf_matrix_val <- table(Predicted = NTEMS_validation$pred_bin, Observed = y_val_bin)
accuracy_val <- sum(diag(conf_matrix_val)) / sum(conf_matrix_val)
cat("Validation Accuracy (0.3 threshold):", round(accuracy_val, 3), "\n")
print(conf_matrix_val)

# Step 5: ROC and AUC
roc_val <- roc(y_val_bin, NTEMS_validation$pred_prob)
auc_val <- auc(roc_val)
cat("Validation AUC (ROC):", round(auc_val, 3), "\n")
plot(roc_val, col = "darkgreen", lwd = 2, main = "ROC Curve (Validation)")
abline(a = 0, b = 1, lty = 2, col = "gray")

# Step 6: Precision-Recall Curve
pr_val <- pr.curve(scores.class0 = NTEMS_validation$pred_prob[y_val_bin == 1],
                   scores.class1 = NTEMS_validation$pred_prob[y_val_bin == 0],
                   curve = TRUE)
cat("Validation AUC (PR Curve):", round(pr_val$auc.integral, 3), "\n")
plot(pr_val, main = "Precision-Recall Curve (Validation)")

# Step 7: Spatial plot of misclassifications (optional)
NTEMS_validation$misclassified <- ifelse(NTEMS_validation$pred_bin != y_val_bin, 1, 0)

NTEMS_val_sf <- st_as_sf(NTEMS_validation, coords = c("x_AEP10TM", "y_AEP10TM"),
                         crs = 3402)  # Replace with your actual CRS if different

ggplot(NTEMS_val_sf) +
  geom_sf(aes(color = factor(misclassified)), alpha = 0.6) +
  scale_color_manual(values = c("0" = "gray60", "1" = "red"),
                     labels = c("Correct", "Misclassified")) +
  labs(title = "Misclassified Observations (Validation Set)",
       color = "Classification") +
  theme_minimal()















#### Assess diff threshold
library(pROC)
library(dplyr)
library(pROC)

# Step 1: Compute ROC object
roc_val <- roc(y_val_bin, NTEMS_validation$pred_prob)

# Step 2: Youden's J statistic
coords_j <- coords(roc_val, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
thresh_j <- as.numeric(coords_j["threshold"])
sensitivity_j <- as.numeric(coords_j["sensitivity"])
specificity_j <- as.numeric(coords_j["specificity"])

# Step 3: Sensitivity = Specificity (fixing structure)
all_coords <- as.data.frame(coords(roc_val, x = "all", ret = c("threshold", "sensitivity", "specificity")))
all_coords$diff <- abs(all_coords$sensitivity - all_coords$specificity)
best_eq_row <- all_coords[which.min(all_coords$diff), ]
thresh_eq <- best_eq_row$threshold
sensitivity_eq <- best_eq_row$sensitivity
specificity_eq <- best_eq_row$specificity

# Step 4: Maximize F1 score
compute_f1 <- function(pred, truth) {
  TP <- sum(pred == 1 & truth == 1)
  FP <- sum(pred == 1 & truth == 0)
  FN <- sum(pred == 0 & truth == 1)
  precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
  recall <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
  if ((precision + recall) == 0) return(0)
  2 * precision * recall / (precision + recall)
}

thresholds <- all_coords$threshold
f1_scores <- sapply(thresholds, function(t) {
  preds <- ifelse(NTEMS_validation$pred_prob >= t, 1, 0)
  compute_f1(preds, y_val_bin)
})

best_f1_idx <- which.max(f1_scores)
thresh_f1 <- thresholds[best_f1_idx]
best_f1 <- f1_scores[best_f1_idx]

# Step 5: Output results
cat("\n--- Threshold Comparison (Validation Set) ---\n")
cat("1. Youden’s J Statistic\n")
cat("  Threshold:", round(thresh_j, 3),
    "| Sensitivity:", round(sensitivity_j, 3),
    "| Specificity:", round(specificity_j, 3), "\n")

cat("2. Sensitivity = Specificity\n")
cat("  Threshold:", round(thresh_eq, 3),
    "| Sensitivity:", round(sensitivity_eq, 3),
    "| Specificity:", round(specificity_eq, 3), "\n")

cat("3. Max F1 Score\n")
cat("  Threshold:", round(thresh_f1, 3),
    "| F1 Score:", round(best_f1, 3), "\n")


















 ###############Boosted regression tree (GMB)
install.packages("gbm")
install.packages("dismo")
install.packages("xgboost")

# Load required libraries
library(gbm)
library(dismo)
library(ggplot2)
#-------------------------------------------------
# 1. RUN BRT MODEL WITH ALL SPATIAL SCALES
#-------------------------------------------------
# Load required libraries
library(gbm)
library(dismo)
library(ggplot2)

BRT_tweedie <- gbm.step(
  data = NTEMS,
  gbm.x = c("prop_con_1_scaled", "prop_con_2_scaled", "prop_con_3_scaled",   
            "clumpy_1_scaled", "clumpy_2_scaled", "clumpy_3_scaled",         
            "age_mn_1_scaled", "age_mn_2_scaled", "age_mn_3_scaled",         
            "hssr_scaled", "ordinalDay_scaled"),  # Use transformed predictors
  gbm.y = "BTNW",
  offset = log(NTEMS$survey_effort),  # Keep offset in log scale (no standardization)
  family = "poisson",
  tree.complexity = 3,   # Recommended value for modeling interactions
  learning.rate = 0.001,  # Slow learning for better model generalization
  bag.fraction = 0.75,     # Recommended value to introduce stochasticity
  n.trees = 1000,         # At least 1000, will be adjusted by cross-validation
  step.size = 50,         # Standard step size
  cv.folds = 10           # Cross-validation for model selection
)


#-------------------------------------------------
# 2. IDENTIFY MOST IMPORTANT SCALE FOR EACH PREDICTOR TYPE
#-------------------------------------------------
# Extract variable importance from the BRT model
var_importance <- summary(BRT_tweedie)

# Find the most important spatial scale for each predictor type
best_prop_con <- var_importance$var[which.max(var_importance$rel.inf[var_importance$var %in% c("prop_con_1", "prop_con_2", "prop_con_3")])]
best_clumpy <- var_importance$var[which.max(var_importance$rel.inf[var_importance$var %in% c("clumpy_1", "clumpy_2", "clumpy_3")])]
best_age_mn <- var_importance$var[which.max(var_importance$rel.inf[var_importance$var %in% c("age_mn_1", "age_mn_2", "age_mn_3")])]

# Print selected best scales
print(paste("Best scale for prop_con:", best_prop_con))
print(paste("Best scale for clumpy:", best_clumpy))
print(paste("Best scale for age_mn:", best_age_mn))

#-------------------------------------------------
# 3. PLOT PARTIAL DEPENDENCE FOR BEST PREDICTORS
#-------------------------------------------------
par(mfrow = c(1,3))  # Arrange plots in a single row
plot(BRT_tweedie, i.var = best_prop_con, main = paste("Partial Dependence of", best_prop_con))
plot(BRT_tweedie, i.var = best_clumpy, main = paste("Partial Dependence of", best_clumpy))
plot(BRT_tweedie, i.var = best_age_mn, main = paste("Partial Dependence of", best_age_mn))

#-------------------------------------------------
# 4. IDENTIFY THRESHOLD VALUES FOR BEST PREDICTORS
#-------------------------------------------------
find_brt_threshold <- function(model, predictor, data) {
  x_vals <- seq(min(data[[predictor]], na.rm = TRUE), max(data[[predictor]], na.rm = TRUE), length.out = 100)
  grid <- data.frame(x = x_vals)
  names(grid) <- predictor
  
  # Keep other predictors at mean values
  for (var in setdiff(names(data), c(predictor, "BTNW"))) {
    grid[[var]] <- mean(data[[var]], na.rm = TRUE)
  }
  
  # Predict across the range of the predictor
  preds <- predict.gbm(model, grid, n.trees = model$gbm.call$best.trees, type = "response")
  
  # Manually apply offset correction
  preds_corrected <- preds + log(data$survey_effort[1])  # Use mean effort if necessary
  preds_final <- exp(preds_corrected)  # Convert back to count scale
  
  # Find threshold where predicted abundance is maximized
  threshold_val <- x_vals[which.max(preds_final)]
  return(threshold_val)
}

# Compute threshold values for best predictors
threshold_best_prop_con <- find_brt_threshold(BRT_tweedie, best_prop_con, NTEMS)
threshold_best_clumpy <- find_brt_threshold(BRT_tweedie, best_clumpy, NTEMS)
threshold_best_age_mn <- find_brt_threshold(BRT_tweedie, best_age_mn, NTEMS)

# Print threshold values
print(paste("Optimal", best_prop_con, ":", round(threshold_best_prop_con, 3)))
print(paste("Optimal", best_clumpy, ":", round(threshold_best_clumpy, 3)))
print(paste("Optimal", best_age_mn, ":", round(threshold_best_age_mn, 3)))

#-------------------------------------------------
# 5. CHECK INTERACTIONS BETWEEN PREDICTORS
#-------------------------------------------------
interaction_results <- gbm.interactions(BRT_tweedie)

# View ranked interactions
print(interaction_results$rank.list)

# Plot strongest interaction (e.g., best coniferous cover scale × best spatial config scale)
gbm.perspec(BRT_tweedie, best_prop_con, best_clumpy, z.range = c(0, max(NTEMS$BTNW)))

#-------------------------------------------------
# 6. MODEL PERFORMANCE & RESIDUAL CHECKS
#-------------------------------------------------
# Check deviance explained
deviance_explained <- 100 * (BRT_tweedie$self.statistics$mean.null - BRT_tweedie$self.statistics$mean.resid) / BRT_tweedie$self.statistics$mean.null
paste("Deviance explained:", round(deviance_explained, 2), "%")

# Predict on training data
NTEMS$predicted <- predict.gbm(BRT_tweedie, NTEMS, n.trees = BRT_tweedie$gbm.call$best.trees, type = "response")

# Manually apply offset correction to predictions
NTEMS$predicted_corrected <- NTEMS$predicted + log(NTEMS$survey_effort)
NTEMS$predicted_final <- exp(NTEMS$predicted_corrected)  # Convert back to count scale

# Compute RMSE
RMSE <- sqrt(mean((NTEMS$BTNW - NTEMS$predicted_final)^2))
paste("Root Mean Squared Error (RMSE):", round(RMSE, 4))

# Residual diagnostics
ggplot(NTEMS, aes(x = predicted_final, y = BTNW - predicted_final)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residual Plot", x = "Predicted Abundance", y = "Residuals")


















library(mgcv)
library(spdep)

# Fit the GAM with both random effect and spatial smooth
gam_model <- gam(
  BTNW ~ 
    s(prop_con_1_scaled) + s(prop_con_2_scaled) + s(prop_con_3_scaled) +  
    s(clumpy_1_scaled) + s(clumpy_2_scaled) + s(clumpy_3_scaled) +         
    s(age_mn_1_scaled) + s(age_mn_2_scaled) + s(age_mn_3_scaled) +         
    s(hssr_scaled) + s(ordinalDay_scaled) +
    s(spatial_group_5000, bs = "re") +  # Cluster random effect
    s(x_AEP10TM, y_AEP10TM, bs = "tp", k = 20),  # Flexible spatial smoother
  family = nb(),
  data = NTEMS,
  method = "REML"
)

# Extract residuals
NTEMS$gam_residuals_spatial <- resid(gam_model)

# Generate nearest neighbors (you only need to do this once)
coords <- NTEMS[, c("x_AEP10TM", "y_AEP10TM")]
neighbors <- knearneigh(coords, k = 4)
neighbor_list <- knn2nb(neighbors)

# Run Moran’s I test on the GAM residuals
moran_test_spatial <- moran.test(NTEMS$gam_residuals_spatial, nb2listw(neighbor_list))

# Print Moran's I results
print(moran_test_spatial)


































##Import data----
# load("0_data/manual/bd_cov_scaled_2024-05-03.rData")
load("0_data/manual/bd_cov_scaled_grouped_2024-05-07.rData")
d<-bd_cov_scaled_grouped%>%
  mutate(spatial_group=as.factor(spatial_group))


##////////////////////////////////////////////////////////////////

#Shannon diversity ----

##Evaluate nonlinear effects----
### Fractional mixedwood----
nl1<-lme4::lmer(shan ~prcB_1000_mean   + (1|location),  
                data=d,  na.action = na.exclude) 
nl2<-lme4::lmer(shan ~poly(prcB_1000_mean,2)   + (1|location),  
                data=d,  na.action = na.exclude) 
nl3<-lme4::lmer(shan ~poly(prcB_1000_mean,3)   + (1|location), 
                data=d,  na.action = na.exclude) 
nl4<-lme4::lmer(shan ~poly(prcB_1000_mean,4)   + (1|location), 
                data=d,  na.action = na.exclude) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$prcB_1000_mean, 'scaled:scale')
ce_1 <- attr(d$prcB_1000_mean, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl2, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

# Find the x-coordinate corresponding to the peak of the trend line
peak_x <- dd2$x[which.max(dd2$predicted)]

# Create the plot
plot <- ggplot(dd2) +
  aes(x, predicted, ymin = conf.low, ymax = conf.high) +
  # Add ribbon for confidence interval
  geom_ribbon(alpha = 0.3) +  
  # Add line for the predicted values
  geom_line() +  
  # Add axis labels
  labs(x = "Fractional Broadleaf", y = "Shannon diversity" ) +  
  # Apply minimal theme
  theme_minimal() +  
  theme(# Remove major gridlines
    panel.grid.major = element_blank(),  
    # Remove minor gridlines
    panel.grid.minor = element_blank(), 
    # Set panel background color to white
    panel.background = element_rect(fill = "white"), 
    # Set axis line color to black
    axis.line = element_line(colour = "black")) + 
  # Add vertical line
  geom_vline(xintercept = peak_x, linetype = "dashed", color = "red")  

# Save the plot
ggsave("3_output/figures/shan_fractional_broadleaf.png", plot, 
       width = 6, height = 4, dpi = 300)


#^2

### Forest age----
nl1<-lme4::lmer(shan ~age_mean_500   + (1|location),  
                data=d,  na.action = na.exclude) 
nl2<-lme4::lmer(shan ~poly(age_mean_500,2)   + (1|location),  
                data=d,  na.action = na.exclude) 
nl3<-lme4::lmer(shan ~poly(age_mean_500,3)   + (1|location), 
                data=d,  na.action = na.exclude) 
nl4<-lme4::lmer(shan ~poly(age_mean_500,4)   + (1|location), 
                data=d,  na.action = na.exclude) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$age_mean_500, 'scaled:scale')
ce_1 <- attr(d$age_mean_500, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl2, terms = c("age_mean_500 [all]"))%>%
  as_tibble()%>%
  mutate(x = x* sc_1 + ce_1)

ggplot(dd2)+
  aes(x, predicted, ymin = conf.low, ymax = conf.high)+
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              se = TRUE)

#^2

### Height upper----
nl1<-lme4::lmer(shan ~height_500_mean   + (1|location),  
                data=d,  na.action = na.exclude) 
nl2<-lme4::lmer(shan ~poly(height_500_mean,2)   + (1|location),  
                data=d,  na.action = na.exclude) 
nl3<-lme4::lmer(shan ~poly(height_500_mean,3)   + (1|location), 
                data=d,  na.action = na.exclude) 
nl4<-lme4::lmer(shan ~poly(height_500_mean,4)   + (1|location), 
                data=d,  na.action = na.exclude) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$height_500_mean, 'scaled:scale')
ce_1 <- attr(d$height_500_mean, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl2, terms = c("height_500_mean [all]"))%>%
  as_tibble()%>%
  mutate(x = x* sc_1 + ce_1)

ggplot(dd2)+
  aes(x, predicted, ymin = conf.low, ymax = conf.high)+
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              se = TRUE)

#^2
##////////////////////////////////////////////////////////////////

## 150----
shan_m_1_150<-lme4::lmer(shan ~  poly(age_mean_150,2, raw=TRUE) 
                         * poly(prcB_150_mean,2, raw=TRUE) 
                         + poly(height_150_mean,2, raw=TRUE)
                         + (1|spatial_group/location)
                         + offset(log(survey_effort)),
                         data=d,  na.action = na.fail)


shan_dredge_150<-MuMIn::dredge(shan_m_1_150,
                               fixed = c("offset(log(survey_effort))",
                                         "(1|spatial_group/location)"),                    
                               extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))


shan_m_2_150 <- lme4::lmer(shan ~ poly(prcB_150_mean,2, raw=TRUE) 
                           + poly(height_150_mean,2, raw=TRUE)
                           + (1|spatial_group/location)
                           + offset(log(survey_effort)),
                           data=d,  na.action = na.fail)

shan_m_2_150_boot <- bootMer(shan_m_2_150, function(x) fixef(x), nsim = 100)


## 500----
shan_m_1_500<-lme4::lmer(shan ~  poly(age_mean_500,2, raw=TRUE) 
                         * poly(prcB_500_mean,2, raw=TRUE) 
                         + poly(height_500_mean,2, raw=TRUE)
                         + (1|spatial_group/location)
                         + offset(log(survey_effort)),
                         data=d,  na.action = na.fail)


shan_dredge_500<-MuMIn::dredge(shan_m_1_500,
                               fixed = c("offset(log(survey_effort))",
                                         "(1|spatial_group/location)"),                    
                               extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

shan_m_2_500 <- lme4::lmer(shan ~  poly(age_mean_500,2, raw=TRUE) 
                           * poly(prcB_500_mean,2, raw=TRUE) 
                           + poly(height_500_mean,2, raw=TRUE)
                           + (1|spatial_group/location)
                           + offset(log(survey_effort)),
                           data=d,  na.action = na.fail)

shan_m_2_500_boot <- bootMer(shan_m_2_500, function(x) fixef(x), nsim = 100)

## 1000----
shan_m_1_1000<-lme4::lmer(shan ~  poly(age_mean_1000,2, raw=TRUE) 
                          * poly(prcB_1000_mean,2, raw=TRUE) 
                          + poly(height_1000_mean,2, raw=TRUE)
                          + (1|spatial_group/location)
                          + offset(log(survey_effort)),
                          data=d,  na.action = na.fail)


shan_dredge_1000<-MuMIn::dredge(shan_m_1_500,
                                fixed = c("offset(log(survey_effort))",
                                          "(1|spatial_group/location)"),                    
                                extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

shan_m_2_1000 <- lme4::lmer(shan ~  poly(age_mean_1000,2, raw=TRUE) 
                            * poly(prcB_1000_mean,2, raw=TRUE) 
                            + poly(height_1000_mean,2, raw=TRUE)
                            + (1|spatial_group/location)
                            + offset(log(survey_effort)),
                            data=d,  na.action = na.fail)

shan_m_2_1000_boot <- bootMer(shan_m_2_1000, function(x) fixef(x), nsim = 100)

## Save----
save(shan_dredge_1000, shan_dredge_150, shan_dredge_500, 
     shan_m_1_1000, shan_m_1_150, shan_m_1_500, shan_m_2_1000, 
     shan_m_2_1000_boot, shan_m_2_150, shan_m_2_150_boot, 
     shan_m_2_500, shan_m_2_500_boot, 
     file = "3_output/models/shan_models.RData")

##////////////////////////////////////////////////////////////////
#Richness ----

##Evaluate nonlinear effects----
### Fractional mixedwood
nl1<-lme4::glmer(richness ~prcB_1000_mean + (1|location),
                 data=d, family=poisson,  na.action = na.exclude) 

nl2<-lme4::glmer(richness ~poly(prcB_500_mean,2) + (1|location), 
                 data=d, family=poisson, na.action = na.exclude, 
                 control=glmerControl(optCtrl=list(maxfun=1e9))) 

nl3<-lme4::glmer(richness ~poly(prcB_1000_mean,3) + (1|location), 
                 data=d, family=poisson,  na.action = na.exclude, 
                 control=glmerControl(optCtrl=list(maxfun=1e9))) 

nl4<-lme4::glmer(richness ~poly(prcB_1000_mean,4) + (1|location), 
                 data=d,family=poisson,  na.action = na.exclude, 
                 control=glmerControl(optCtrl=list(maxfun=1e9))) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$prcB_500_mean, 'scaled:scale')
ce_1 <- attr(d$prcB_500_mean, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl2, terms = c("prcB_500_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

dd3<-ggpredict(nl3, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

dd4<-ggpredict(nl4, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)
# Find the x-coordinate corresponding to the peak of the trend line
peak_x <- dd2$x[which.max(dd2$predicted)]

# Create the plot
plot <- ggplot(dd2) +
  aes(x, predicted, ymin = conf.low, ymax = conf.high) +
  # Add ribbon for confidence interval
  geom_ribbon(alpha = 0.3) +  
  # Add line for the predicted values
  geom_line() +  
  # Add axis labels
  labs(x = "Fractional Broadleaf", y = "Richness") +  
  # Apply minimal theme
  theme_minimal() +  
  theme(# Remove major gridlines
    panel.grid.major = element_blank(),  
    # Remove minor gridlines
    panel.grid.minor = element_blank(), 
    # Set panel background color to white
    panel.background = element_rect(fill = "white"), 
    # Set axis line color to black
    axis.line = element_line(colour = "black")) + 
  # Add vertical line
  geom_vline(xintercept = peak_x, linetype = "dashed", color = "red")  

# Save the plot
ggsave("3_output/figures/richness_fractional_broadleaf.png", plot, width = 6, 
       height = 4, dpi = 300)

##////////////////////////////////////////////////////////////////

## 150----
richness_m_1_150<-lme4::glmer(richness ~  poly(age_mean_150,2, raw=TRUE) 
                              * poly(prcB_150_mean,2, raw=TRUE) 
                              + poly(height_150_mean,2, raw=TRUE)
                              + (1|spatial_group/location)
                              + offset(log(survey_effort)),
                              data=d, family=poisson, na.action = na.fail,
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


richness_dredge_150<-MuMIn::dredge(richness_m_1_150,
                                   fixed = c("offset(log(survey_effort))",
                                             "(1|spatial_group/location)"),                    
                                   extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

richness_m_2_150 <-lme4::glmer(richness ~  poly(age_mean_150,2, raw=TRUE) 
                               * poly(prcB_150_mean,2, raw=TRUE) 
                               + poly(height_150_mean,2, raw=TRUE)
                               + (1|spatial_group/location)
                               + offset(log(survey_effort)),
                               data=d, family=poisson, na.action = na.fail)

richness_m_2_150_boot <- bootMer(richness_m_2_150, function(x) fixef(x), nsim = 100)

## 500----
richness_m_1_500<-lme4::glmer(richness ~  poly(age_mean_500,2, raw=TRUE) 
                              * poly(prcB_500_mean,2, raw=TRUE) 
                              + poly(height_500_mean,2, raw=TRUE)
                              + (1|spatial_group/location)
                              + offset(log(survey_effort)),
                              data=d, family=poisson, na.action = na.fail,
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


richness_dredge_500<-MuMIn::dredge(richness_m_1_500,
                                   fixed = c("offset(log(survey_effort))",
                                             "(1|spatial_group/location)"),                    
                                   extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

richness_m_2_500 <- lme4::glmer(richness ~  poly(age_mean_500,2, raw=TRUE) 
                                + poly(prcB_500_mean,2, raw=TRUE) 
                                + poly(height_500_mean,2, raw=TRUE)
                                + (1|spatial_group/location)
                                + offset(log(survey_effort)),
                                data=d, family=poisson, na.action = na.fail,
                                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


richness_m_2_500_boot <- bootMer(richness_m_2_500, function(x) fixef(x), nsim = 100)

## 1000----
richness_m_1_1000<-lme4::glmer(richness ~  poly(age_mean_1000,2, raw=TRUE) 
                               * poly(prcB_1000_mean,2, raw=TRUE) 
                               + poly(height_1000_mean,2, raw=TRUE)
                               + (1|spatial_group/location)
                               + offset(log(survey_effort)),
                               data=d, family=poisson, na.action = na.fail,
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


richness_dredge_1000<-MuMIn::dredge(richness_m_1_500,
                                    fixed = c("offset(log(survey_effort))",
                                              "(1|spatial_group/location)"),                    
                                    extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

richness_m_2_1000 <- lme4::glmer(richness ~  poly(age_mean_1000,2, raw=TRUE) 
                                 * poly(prcB_1000_mean,2, raw=TRUE) 
                                 + poly(height_1000_mean,2, raw=TRUE)
                                 + (1|spatial_group/location)
                                 + offset(log(survey_effort)),
                                 data=d, family=poisson, na.action = na.fail)	
richness_m_2_1000_boot <- bootMer(richness_m_2_1000, function(x) fixef(x), nsim = 100)

## Save----
save(richness_dredge_1000, richness_dredge_150, richness_dredge_500, 
     richness_m_1_1000, richness_m_1_150, richness_m_1_500, 
     richness_m_2_1000, richness_m_2_1000_boot, richness_m_2_150, 
     richness_m_2_150_boot, richness_m_2_500, richness_m_2_500_boot,
     file = "3_output/models/richness_models.RData")

##////////////////////////////////////////////////////////////////
# Evenness ----

##Evaluate nonlinear effects----
### Fractional mixedwood
nl1<-lme4::lmer(pilou_even ~prcB_1000_mean  + (1|spatial_group)+(1|location),  
                data=d,  na.action = na.exclude) 
nl2<-lme4::lmer(pilou_even ~poly(prcB_1000_mean,2)+ (1|spatial_group)+(1|location),  
                data=d,  na.action = na.exclude) 
nl3<-lme4::lmer(pilou_even ~poly(prcB_1000_mean,3)+ (1|spatial_group)+(1|location), 
                data=d,  na.action = na.exclude) 
nl4<-lme4::lmer(pilou_even ~poly(prcB_1000_mean,4)+ (1|spatial_group)+(1|location), 
                data=d,  na.action = na.exclude) 

print(AIC(nl1, nl2, nl3, nl4))

sc_1 <- attr(d$prcB_1000_mean, 'scaled:scale')
ce_1 <- attr(d$prcB_1000_mean, 'scaled:center')

# make predictions from model using ggpredict
dd2<-ggpredict(nl1, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

dd3<-ggpredict(nl3, terms = c("prcB_1000_mean[all]"))%>%
  as_tibble()%>%
  # unscale the predictor
  mutate(x = x* sc_1 + ce_1)

# Find the x-coordinate corresponding to the peak of the trend line
peak_x <- dd2$x[which.min(dd2$predicted)]

# Create the plot
plot <- ggplot(dd2) +
  aes(x, predicted, ymin = conf.low, ymax = conf.high) +
  # Add ribbon for confidence interval
  geom_ribbon(alpha = 0.3) +  
  # Add line for the predicted values
  geom_line() +  
  # Add axis labels
  labs(x = "Fractional Broadleaf", y = "Pielou's evenness") +  
  # Apply minimal theme
  theme_minimal() +  
  theme(# Remove major gridlines
    panel.grid.major = element_blank(),  
    # Remove minor gridlines
    panel.grid.minor = element_blank(), 
    # Set panel background color to white
    panel.background = element_rect(fill = "white"), 
    # Set axis line color to black
    axis.line = element_line(colour = "black")) + 
  # Add vertical line
  geom_vline(xintercept = peak_x, linetype = "dashed", color = "red")  

# Save the plot
ggsave("3_output/figures/pilou_even_fractional_broadleaf.png", plot, width = 6, 
       height = 4, dpi = 300)

##////////////////////////////////////////////////////////////////

## 150----
pilou_even_m_1_150<-lme4::lmer(pilou_even ~  poly(age_mean_150,2, raw=TRUE) 
                               * poly(prcB_150_mean,2, raw=TRUE) 
                               + poly(height_150_mean,2, raw=TRUE)
                               + (1|spatial_group/location)
                               + offset(log(survey_effort)),
                               data=d,  na.action = na.fail,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


pilou_even_dredge_150<-MuMIn::dredge(pilou_even_m_1_150,
                                     fixed = c("offset(log(survey_effort))",
                                               "(1|spatial_group/location)"),                    
                                     extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

pilou_even_m_2_150 <- lme4::lmer(pilou_even ~  poly(age_mean_150,2, raw=TRUE) 
                                 + prcB_150_mean
                                 + poly(height_150_mean,2, raw=TRUE)
                                 + (1|location)
                                 + offset(log(survey_effort)),
                                 data=d,  na.action = na.fail,
                                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

pilou_even_m_2_150_boot <- bootMer(pilou_even_m_2_150, function(x) fixef(x), nsim = 100)

## 500----
pilou_even_m_1_500<-lme4::lmer(pilou_even ~  poly(age_mean_500,2, raw=TRUE) 
                               + poly(prcB_500_mean,2, raw=TRUE) 
                               + poly(height_500_mean,2, raw=TRUE)
                               + (1|spatial_group)
                               + offset(log(survey_effort)),
                               data=d,  na.action = na.fail,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


pilou_even_dredge_500<-MuMIn::dredge(pilou_even_m_1_500,
                                     fixed = c("offset(log(survey_effort))",
                                               "(1|spatial_group/location)"),                    
                                     extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

pilou_even_m_2_500 <- lme4::lmer(pilou_even ~  poly(age_mean_500,2, raw=TRUE) 
                                 + prcB_500_mean
                                 + poly(height_500_mean,2, raw=TRUE)
                                 + (1|location)
                                 + offset(log(survey_effort)),
                                 data=d,  na.action = na.fail,
                                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


pilou_even_m_2_500_boot <- bootMer(pilou_even_m_2_500, function(x) fixef(x), nsim = 100)

## 1000----
pilou_even_m_1_1000<-lme4::lmer(pilou_even ~  poly(age_mean_1000,2, raw=TRUE) 
                                * poly(prcB_1000_mean,2, raw=TRUE) 
                                + poly(height_1000_mean,2, raw=TRUE)
                                + (1|spatial_group)
                                + offset(log(survey_effort)),
                                data=d,  na.action = na.fail,
                                control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))


pilou_even_dredge_1000<-MuMIn::dredge(pilou_even_m_1_500,
                                      fixed = c("offset(log(survey_effort))",
                                                "(1|spatial_group/location)"),                    
                                      extra = c(r2=function(x) round(r.squaredGLMM(x)[1,c(1,2)],3)))

pilou_even_m_2_1000 <- lme4::lmer(pilou_even ~  poly(age_mean_1000,2, raw=TRUE) 
                                  + prcB_1000_mean
                                  + poly(height_1000_mean,2, raw=TRUE)
                                  + (1|location)
                                  + offset(log(survey_effort)),
                                  data=d,  na.action = na.fail,
                                  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))



pilou_even_m_2_1000_boot <- bootMer(pilou_even_m_2_1000, function(x) fixef(x), nsim = 100)


## Save----
# save(pilou_even_dredge_1000, pilou_even_dredge_150, 
#      pilou_even_dredge_500, pilou_even_m_1_1000, 
#      pilou_even_m_1_150, pilou_even_m_1_500, 
#      pilou_even_m_2_1000, pilou_even_m_2_1000_boot, 
#      pilou_even_m_2_150, pilou_even_m_2_150_boot, 
#      pilou_even_m_2_500, pilou_even_m_2_500_boot,
#      file = "3_output/models/pilou_even_models.RData")

save(   pilou_even_m_2_1000,  
        pilou_even_m_2_150, 
        pilou_even_m_2_500, 
        # pilou_even_m_2_150_boot, 
        # pilou_even_m_2_500_boot,
        # pilou_even_m_2_1000_boot,
        # pilou_even_dredge_1000, pilou_even_dredge_150, 
        # pilou_even_dredge_500, pilou_even_m_1_1000, 
        # pilou_even_m_1_150, pilou_even_m_1_500, 
        
        file = "3_output/models/pilou_even_models.RData")

##////////////////////////////////////////////////////////////////



library(caret)
# Define the control parameters for cross-validation
ctrl <- trainControl(method = "cv", number = 3)

# Perform cross-validation
cv_results <- train(nl2, method = "lmer", trControl = ctrl)

cv_results <- bootMer(nl2, FUN = fixef, nsim = 3)  # Use appropriate number of nsim

boot_gm1 <- bootMer(gm1, FUN = function(x) fixef(x), nsim = 100)

pilou_even_m_2_150_boot <- bootMer(pilou_even_m_2_150, function(x) fixef(x), nsim = 100)




































library(MASS) # For glm.nb
library(lme4) # For glmer.nb


summary(BTNW_comp_grad)
summary(BTNW_spatial_config)
summary(BTNW_forest_age)
#-------------------------------------------------------------------------------------------------------------------
# BTNW Models
#-------------------------------------------------------------------------------------------------------------------

# Baseline Model (Random Effects Only)
BTNW_baseline <- glmer.nb(BTNW ~ offset(log(survey_effort)) + hssr + ordinalDay + (1 | cluster_k_10), data = NTEMS)

# Compositional Gradient
BTNW_comp_grad <- glmer.nb(BTNW ~ poly(prop_con_1_scaled, 3, raw = TRUE) +
                             hssr + ordinalDay + offset(log(survey_effort)) + (1 | cluster_k_10),
                           data = NTEMS,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
# Spatial Configuration 
BTNW_spatial_config <- glmer.nb(BTNW ~ I(clumpy_1_scaled^(1/3)) + hssr + ordinalDay +
                                  offset(log(survey_effort)) + (1 | cluster_k_10), 
                                data = NTEMS,
                                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
# Forest Age 
BTNW_forest_age <- glmer.nb(BTNW ~ poly(age_mn_1_scaled,3, raw=TRUE) + hssr + ordinalDay +
                              offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS,
                            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Species Composition + Spatial Configuration
BTNW_comp_spatial <- glmer.nb(BTNW ~ poly(prop_con_1_scaled,3, raw=TRUE) + hssr + ordinalDay + I(clumpy_1_scaled^(1/3)) + 
                                offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS,
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Species Composition x Spatial Configuration
BTNW_comp_x_spatial <- glmer.nb(BTNW ~ poly(prop_con_1_scaled,3,raw=TRUE)*I(clumpy_1_scaled^(1/3)) + hssr + ordinalDay +
                                  offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS,
                                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Species Composition x Age
BTNW_comp_x_age <- glmer.nb(BTNW ~ poly(prop_con_1_scaled,3,raw=TRUE)*poly(age_mn_1_scaled,3,raw=TRUE) + hssr + ordinalDay +
                              offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS,
                            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Species Composition x Spatial Configuration x Age
BTNW_comp_x_spatial_x_age <- glmer.nb(BTNW ~ poly(prop_con_1_scaled,2,raw=TRUE)*poly(age_mn_1_scaled,2,raw=TRUE)*I(clumpy_1_scaled^(1/3)) + hssr + ordinalDay +
                                        offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS,
                                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Global Model (All terms included)
BTNW_global <- glmer.nb(BTNW ~ poly(prop_con_1_scaled,3,raw=TRUE) + poly(age_mn_1_scaled,3,raw=TRUE) + I(clumpy_1_scaled^(1/3)) + hssr + ordinalDay +
                          offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Extract and print AIC values for BTNW
BTNW_aic <- data.frame(
  Model = c("Baseline", "Comp_Grad", "Spatial_Config", "Forest_Age", "Comp_Spatial", "Comp_x_Spatial", 
            "Comp_x_Age", "Comp_x_Spatial_x_Age", "Global"),
  AIC = c(AIC(BTNW_baseline), AIC(BTNW_comp_grad), AIC(BTNW_spatial_config), AIC(BTNW_forest_age), 
          AIC(BTNW_comp_spatial), AIC(BTNW_comp_x_spatial), AIC(BTNW_comp_x_age), 
          AIC(BTNW_comp_x_spatial_x_age), AIC(BTNW_global))
)

print(BTNW_aic[order(BTNW_aic$AIC), ])  # Sort models by AIC (lower is better)












summary(BTNW_global)


#-------------------------------------------------------------------------------------------------------------------
# TEWA Models
#-------------------------------------------------------------------------------------------------------------------

# Null Model
TEWA_null <- glmer.nb(TEWA ~ 1 + (1 | cluster_k_10), data = NTEMS)

# Baseline Model (Random Effects Only)
TEWA_baseline <- glmer.nb(TEWA ~ offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Compositional Gradient
TEWA_comp_grad <- glmer.nb(TEWA ~ prop_con_1_scaled + offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Spatial Configuration
TEWA_spatial_config <- glmer.nb(TEWA ~ I(clumpy_1_scaled^(1/3)) + offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Forest Age
TEWA_forest_age <- glmer.nb(TEWA ~ poly(age_mn_1_scaled, 4, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Species Composition + Spatial Configuration
TEWA_comp_spatial <- glmer.nb(TEWA ~ prop_con_1_scaled + I(clumpy_1_scaled^(1/3)) + offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Species Composition x Spatial Configuration
TEWA_comp_x_spatial <- glmer.nb(TEWA ~ prop_con_1_scaled * I(clumpy_1_scaled^(1/3)) + offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Species Composition x Age
TEWA_comp_x_age <- glmer.nb(TEWA ~ prop_con_1_scaled * poly(age_mn_1_scaled, 4, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Species Composition x Spatial Configuration x Age
TEWA_comp_x_spatial_x_age <- glmer.nb(TEWA ~ prop_con_1_scaled * I(clumpy_1_scaled^(1/3)) * poly(age_mn_1_scaled, 4, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Global Model
TEWA_global <- glmer.nb(TEWA ~ prop_con_1_scaled + I(clumpy_1_scaled^(1/3)) + poly(age_mn_1_scaled, 4, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Extract and print AIC values for TEWA
TEWA_aic <- data.frame(
  Model = c("Null", "Baseline", "Comp_Grad", "Spatial_Config", "Forest_Age", "Comp_Spatial", "Comp_x_Spatial", "Comp_x_Age", "Comp_x_Spatial_x_Age", "Global"),
  AIC = c(AIC(TEWA_null), AIC(TEWA_baseline), AIC(TEWA_comp_grad), AIC(TEWA_spatial_config), AIC(TEWA_forest_age), AIC(TEWA_comp_spatial), AIC(TEWA_comp_x_spatial), AIC(TEWA_comp_x_age), AIC(TEWA_comp_x_spatial_x_age), AIC(TEWA_global))
)
print(TEWA_aic[order(TEWA_aic$AIC), ])

#-------------------------------------------------------------------------------------------------------------------
# BBWA Models
#-------------------------------------------------------------------------------------------------------------------

# Null Model
BBWA_null <- glmer.nb(BBWA ~ 1 + (1 | cluster_k_10), data = NTEMS)

# Baseline Model (Random Effects Only)
BBWA_baseline <- glmer.nb(BBWA ~ offset(log(survey_effort)) + (1 | cluster_k_10), data = NTEMS)

# Compositional Gradient
BBWA_comp_grad <- glmer.nb(BBWA ~ poly(prop_con_1
