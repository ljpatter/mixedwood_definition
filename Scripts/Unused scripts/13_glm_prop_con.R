# ---
# title: "GLM prop con"
# author: "Leonard Patterson"
# created: "2025-08-11"
# description: This code runs a GLMM To compare the outputs with those created from the BRT.
# ---


# Load libraries
library(dplyr)
library(spdep)
library(ggplot2)
library(stats)
library(ggplot2)
library(dplyr)
library(fitdistrplus)

# Load data
joined_data <- read.csv("Output/Tabular Data/joined_data.csv")


######### Evaluate distribution of response variable

# Response variables to check
response_vars <- c("BTNW")

for (response in response_vars) {
  print(paste("Evaluating distribution for:", response))
  
  # Histogram
  hist(joined_data[[response]],
       main = paste("Histogram of", response),
       xlab = response,
       ylab = "Frequency",
       col = "lightblue",
       border = "black")
  
  # Distribution plots
  plotdist(joined_data[[response]], discrete = TRUE, histo = TRUE, demp = TRUE)
  descdist(joined_data[[response]], discrete = TRUE, boot = 1000)
  
  # Fit distributions
  fit_p <- fitdist(joined_data[[response]], "pois")
  fit_nb <- fitdist(joined_data[[response]], "nbinom")
  fit_norm <- fitdist(joined_data[[response]], "norm")
  
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


### Negative binomial is the best fit 







############## Account for SAC in response by using a kNN means approach to defining clusters
############## which can then be used as a random effect in a mixed effect model


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







######## Combine NTEMS and LULC w/ clustering


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
      survey_duration_method, max_dist_band, cluster_k_10,
      x_AEP10TM, y_AEP10TM,
      all_of(species_cols),
      starts_with("prop_con_"), starts_with("clumpy_"), starts_with("age_mn_")
    )
  
  
  print("Join and replacement successful.")
} else {
  print("Mismatches found. Join and replacement aborted.")
}

write.csv(joined_data, "Output/Tabular Data/joined_data_w_clustering.csv")











################ MODELS


library(glmmTMB)
library(dplyr)
library(ggplot2)

# Load data
joined_data <- read.csv("Output/Tabular Data/joined_data_w_clustering.csv")

# Ensure factor for random effect
joined_data$cluster_k_10 <- as.factor(joined_data$cluster_k_10)

### --- Models with prop_con_150_NTEMS ---

model_lin_ntems <- glmmTMB(
  BTNW ~ prop_con_150_NTEMS + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_quad_ntems <- glmmTMB(
  BTNW ~ poly(prop_con_150_NTEMS, 2, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_cubic_ntems <- glmmTMB(
  BTNW ~ poly(prop_con_150_NTEMS, 3, raw = TRUE) + x_AEP10TM + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_log_ntems <- glmmTMB(
  BTNW ~ log(prop_con_150_NTEMS + 0.001) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_exp_ntems <- glmmTMB(
  BTNW ~ poly(prop_con_150_NTEMS, 2, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),  # same as quadratic
  family = nbinom2,
  data = joined_data
)

threshold_ntems <- median(joined_data$prop_con_150_NTEMS, na.rm = TRUE)
joined_data$above_thresh_ntems <- as.numeric(joined_data$prop_con_150_NTEMS > threshold_ntems)

model_piecewise_ntems <- glmmTMB(
  BTNW ~ prop_con_150_NTEMS * above_thresh_ntems + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)


### --- Models with prop_con_150_LULC ---

model_lin_lulc <- glmmTMB(
  BTNW ~ prop_con_150_LULC + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_quad_lulc <- glmmTMB(
  BTNW ~ poly(prop_con_150_LULC, 2, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_cubic_lulc <- glmmTMB(
  BTNW ~ poly(prop_con_150_LULC, 3, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_log_lulc <- glmmTMB(
  BTNW ~ log(prop_con_150_LULC + 0.001) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_exp_lulc <- glmmTMB(
  BTNW ~ poly(prop_con_150_LULC, 2, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),  # same as quadratic
  family = nbinom2,
  data = joined_data
)

threshold_lulc <- median(joined_data$prop_con_150_LULC, na.rm = TRUE)
joined_data$above_thresh_lulc <- as.numeric(joined_data$prop_con_150_LULC > threshold_lulc)

model_piecewise_lulc <- glmmTMB(
  BTNW ~ prop_con_150_LULC * above_thresh_lulc + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

### --- Compare All Models by AIC ---
aic_values <- data.frame(
  Model = c(
    "NTEMS_Linear", "NTEMS_Quadratic", "NTEMS_Cubic", "NTEMS_Log", "NTEMS_Quadratic(Exp)", "NTEMS_Piecewise",
    "LULC_Linear", "LULC_Quadratic", "LULC_Cubic", "LULC_Log", "LULC_Quadratic(Exp)", "LULC_Piecewise"
  ),
  AIC = c(
    AIC(model_lin_ntems), AIC(model_quad_ntems), AIC(model_cubic_ntems), AIC(model_log_ntems), AIC(model_exp_ntems), AIC(model_piecewise_ntems),
    AIC(model_lin_lulc), AIC(model_quad_lulc), AIC(model_cubic_lulc), AIC(model_log_lulc), AIC(model_exp_lulc), AIC(model_piecewise_lulc)
  )
)







# Load required library
library(ggplot2)

# Build prediction data
newdata <- data.frame(
  prop_con_150_NTEMS = seq(
    min(joined_data$prop_con_150_NTEMS, na.rm = TRUE),
    max(joined_data$prop_con_150_NTEMS, na.rm = TRUE),
    length.out = 200
  )
)

# Add polynomial terms for raw = TRUE
newdata$`poly(prop_con_150_NTEMS, 3, raw = TRUE)1` <- newdata$prop_con_150_NTEMS
newdata$`poly(prop_con_150_NTEMS, 3, raw = TRUE)2` <- newdata$prop_con_150_NTEMS^2
newdata$`poly(prop_con_150_NTEMS, 3, raw = TRUE)3` <- newdata$prop_con_150_NTEMS^3

# Add offset term (mean survey effort)
newdata$survey_effort <- mean(joined_data$survey_effort, na.rm = TRUE)

# Predict values
newdata$predicted <- predict(
  model_cubic_ntems,
  newdata = newdata,
  type = "response",
  re.form = NA
)

# Plot with raw data + predicted line
ggplot() +
  geom_point(data = joined_data, aes(x = prop_con_150_NTEMS, y = BTNW),
             alpha = 0.3, color = "black") +
  geom_line(data = newdata, aes(x = prop_con_150_NTEMS, y = predicted),
            color = "blue", linewidth = 1.2) +
  labs(
    title = "BTNW Abundance vs. Proportion Conifer (150m, NTEMS)",
    x = "Proportion Conifer (150m, NTEMS)",
    y = "BTNW Relative Abundance"
  ) +
  theme_minimal(base_size = 14)











########### Assess model residuals

# Load necessary packages
library(glmmTMB)
library(DHARMa) # For residual analysis

# Ensure 'joined_data' is available in your R environment
# (Assuming 'joined_data' and the variables within it are already loaded)

# 1. Fit the glmmTMB model
# This is the model provided by the user
model_cubic_ntems <- glmmTMB(
  BTNW ~ poly(prop_con_150_NTEMS, 3, raw = TRUE) + x_AEP10TM + y_AEP10TM + offset(log(survey_effort)),
  family = nbinom2,
  data = joined_data
)

# Print model summary (optional, but good for checking model fit)
summary(model_cubic_ntems)

# 2. Simulate residuals using DHARMa
# The simulateResiduals function takes the fitted model as input.
# It simulates new data from the fitted model and compares the observed data
# to these simulations to create standardized residuals.
# 'n = 250' is a good starting number of simulations; increase for more precision if needed.
simulated_residuals <- simulateResiduals(model_cubic_ntems, n = 250)

# 3. Plot the DHARMa residuals
# This generates a set of diagnostic plots:
# - Quantile-quantile (QQ) plot: Checks for overall fit and distribution of residuals.
# - Residuals vs. predicted values: Checks for homogeneity of variance and patterns.
# - Residuals vs. each predictor: Checks for patterns with individual predictors.
plot(simulated_residuals)

# 4. Run specific DHARMa tests (recommended for formal assessment)

# Test for dispersion (overdispersion or underdispersion)
# For negative binomial models (nbinom2), this checks for residual dispersion
# beyond what the negative binomial distribution already accounts for.
# A p-value < 0.05 suggests significant over/underdispersion.
testDispersion(simulated_residuals)

# Test for zero-inflation
# This checks if there are more zeros in the observed data than predicted by the model.
# A p-value < 0.05 suggests significant zero-inflation not captured by the model.
testZeroInflation(simulated_residuals)

# You can also test for specific predictors if you suspect issues
# For example, residuals vs. 'prop_con_150_NTEMS'
plotResiduals(simulated_residuals, form = joined_data$prop_con_150_NTEMS)

# Or residuals vs. the fitted values
plotResiduals(simulated_residuals, form = fitted(model_cubic_ntems))



























################ BBWA MODELS


library(glmmTMB)
library(dplyr)
library(ggplot2)

# Ensure factor for random effect
joined_data$cluster_k_10 <- as.factor(joined_data$cluster_k_10)

### --- Models with prop_con_150_NTEMS ---

model_lin_ntems <- glmmTMB(
  BBWA ~ prop_con_150_NTEMS + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_quad_ntems <- glmmTMB(
  BBWA ~ poly(prop_con_150_NTEMS, 2, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_cubic_ntems <- glmmTMB(
  BBWA ~ poly(prop_con_150_NTEMS, 3, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_log_ntems <- glmmTMB(
  BBWA ~ log(prop_con_150_NTEMS + 0.001) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_exp_ntems <- glmmTMB(
  BBWA ~ poly(prop_con_150_NTEMS, 2, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),  # same as quadratic
  family = nbinom2,
  data = joined_data
)

threshold_ntems <- median(joined_data$prop_con_150_NTEMS, na.rm = TRUE)
joined_data$above_thresh_ntems <- as.numeric(joined_data$prop_con_150_NTEMS > threshold_ntems)

model_piecewise_ntems <- glmmTMB(
  BBWA ~ prop_con_150_NTEMS * above_thresh_ntems + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)


### --- Models with prop_con_150_LULC ---

model_lin_lulc <- glmmTMB(
  BBWA ~ prop_con_150_LULC + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_quad_lulc <- glmmTMB(
  BBWA ~ poly(prop_con_150_LULC, 2, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_cubic_lulc <- glmmTMB(
  BBWA ~ poly(prop_con_150_LULC, 3, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_log_lulc <- glmmTMB(
  BBWA ~ log(prop_con_150_LULC + 0.001) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_exp_lulc <- glmmTMB(
  BBWA ~ poly(prop_con_150_LULC, 2, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),  # same as quadratic
  family = nbinom2,
  data = joined_data
)

threshold_lulc <- median(joined_data$prop_con_150_LULC, na.rm = TRUE)
joined_data$above_thresh_lulc <- as.numeric(joined_data$prop_con_150_LULC > threshold_lulc)

model_piecewise_lulc <- glmmTMB(
  BBWA ~ prop_con_150_LULC * above_thresh_lulc + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

### --- Compare All Models by AIC ---
aic_values <- data.frame(
  Model = c(
    "NTEMS_Linear", "NTEMS_Quadratic", "NTEMS_Cubic", "NTEMS_Log", "NTEMS_Quadratic(Exp)", "NTEMS_Piecewise",
    "LULC_Linear", "LULC_Quadratic", "LULC_Cubic", "LULC_Log", "LULC_Quadratic(Exp)", "LULC_Piecewise"
  ),
  AIC = c(
    AIC(model_lin_ntems), AIC(model_quad_ntems), AIC(model_cubic_ntems), AIC(model_log_ntems), AIC(model_exp_ntems), AIC(model_piecewise_ntems),
    AIC(model_lin_lulc), AIC(model_quad_lulc), AIC(model_cubic_lulc), AIC(model_log_lulc), AIC(model_exp_lulc), AIC(model_piecewise_lulc)
  )
)

aic_values <- aic_values[order(aic_values$AIC), ]
print(aic_values)





## Plot best model (LULC_cubic)

# 2. Create 'newdata' for predictions
# Define the range for the predictor variable (prop_con_150_LULC)
min_prop <- min(joined_data$prop_con_150_LULC, na.rm = TRUE)
max_prop <- max(joined_data$prop_con_150_LULC, na.rm = TRUE)

# Create a sequence of values for the predictor to get a smooth prediction line
plot_prop_con_150_LULC <- seq(min_prop, max_prop, length.out = 100)

# Calculate the mean of 'survey_effort' to use for the offset in predictions.
# This represents a typical survey effort across which to predict.
mean_survey_effort <- mean(joined_data$survey_effort, na.rm = TRUE)

# Create the newdata dataframe for prediction.
# We only need the fixed effect predictor and the offset variable.
# The random effect 'cluster_k_10' is handled by 're.form = NA' in predict().
newdata <- data.frame(
  prop_con_150_LULC = plot_prop_con_150_LULC,
  survey_effort = mean_survey_effort
)

# 3. Predict 'BBWA' values using the model
# type = "response" ensures predictions are on the original scale of BBWA (counts).
# re.form = NA ensures that only the fixed effects are used for prediction,
# effectively giving the population-level trend, ignoring cluster-specific variations.
newdata$predicted <- predict(model_cubic_lulc, newdata = newdata, type = "response", re.form = NA)

# 4. Plotting the raw data and the predicted line
# We apply the log10(x+1) transformation to the y-axis, consistent with previous plots,
# to handle potential zero values in 'BBWA' and visualize the data on a log scale.
ggplot() +
  # Plot raw data points. Add 1 to BBWA before logging to handle zeros.
  geom_point(data = joined_data, aes(x = prop_con_150_LULC, y = BBWA + 1),
             alpha = 0.3, color = "black") +
  # Plot the predicted line. Add 1 to predicted values for consistency with the y-axis scale.
  geom_line(data = newdata, aes(x = prop_con_150_LULC, y = predicted + 1),
            color = "blue", linewidth = 1.2) +
  # Apply log10 transformation to the y-axis.
  # Use trans_breaks and trans_format from the 'scales' package for cleaner log-scale labels
  # (e.g., 10^0, 10^1, etc.).
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(
    title = "BBWA Abundance vs. Proportion Conifer (150m, LULC)",
    x = "Proportion Conifer (150m, LULC)",
    y = "BBWA Relative Abundance (Log10(x+1) Scale)" # Reflect the transformation in the label
  ) +
  theme_minimal(base_size = 14)


















########## Convert easting and northing to lat and long so model can solve

# Install the 'sf' package if you haven't already
# install.packages("sf")

# Load the 'sf' package
library(sf)

utm_crs <- target_crs_terra  # NAD83 / UTM zone 11N
latlon_crs <- st_crs(4326) # WGS84 (Latitude/Longitude)

# 2. Create an 'sf' object from your UTM coordinates
#    st_as_sf takes your data frame and specifies which columns are X and Y,
#    and what their current CRS is.
sf_points <- st_as_sf(joined_data,
                      coords = c("x_AEP10TM", "y_AEP10TM"),
                      crs = utm_crs)

# 3. Transform the 'sf' object to the desired Latitude/Longitude CRS
transformed_points <- st_transform(sf_points, crs = latlon_crs)

# 4. Extract the new latitude and longitude coordinates
#    st_coordinates extracts the X and Y coordinates from the 'sf' object.
#    The order will be Longitude (X) then Latitude (Y).
coords_matrix <- st_coordinates(transformed_points)

# 5. Add the new latitude and longitude columns to your original data frame
joined_data$longitude <- coords_matrix[, "X"]
joined_data$latitude <- coords_matrix[, "Y"]

# 6. View the updated data frame with the new columns
print(head(joined_data))

# You can also check the class of the new columns
# print(class(joined_data$latitude))
# print(class(joined_data$longitude))









########### Models including lat and long

model_1 <- glmmTMB(
  BBWA ~ poly(prop_con_150_NTEMS, 3, raw = TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_2 <- glmmTMB(
  BBWA ~ poly(prop_con_150_NTEMS, 3, raw = TRUE) + latitude + longitude + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

model_3 <- glmmTMB(
  BBWA ~ poly(prop_con_150_NTEMS, 3, raw = TRUE) + age_mn_150_NTEMS + latitude + longitude + latitude*longitude + offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)


### --- Compare All Models by AIC ---
aic_values <- data.frame(
  Model = c(
    "model_1", "model_2", "model_3"
  ),
  AIC = c(
    AIC(model_1), AIC(model_2), AIC(model_3)
  )
)

aic_values <- aic_values[order(aic_values$AIC), ]
print(aic_values)
















########### Assess model residuals

# Load necessary packages
library(glmmTMB)
library(DHARMa) # For residual analysis

# 2. Simulate residuals using DHARMa
# The simulateResiduals function takes the fitted model as input.
# It simulates new data from the fitted model and compares the observed data
# to these simulations to create standardized residuals.
# 'n = 250' is a good starting number of simulations; increase for more precision if needed.
simulated_residuals <- simulateResiduals(model_3, n = 250)

# 3. Plot the DHARMa residuals
# This generates a set of diagnostic plots:
# - Quantile-quantile (QQ) plot: Checks for overall fit and distribution of residuals.
# - Residuals vs. predicted values: Checks for homogeneity of variance and patterns.
# - Residuals vs. each predictor: Checks for patterns with individual predictors.
plot(simulated_residuals)

# 4. Run specific DHARMa tests (recommended for formal assessment)

# Test for dispersion (overdispersion or underdispersion)
# For negative binomial models (nbinom2), this checks for residual dispersion
# beyond what the negative binomial distribution already accounts for.
# A p-value < 0.05 suggests significant over/underdispersion.
testDispersion(simulated_residuals)

# Test for zero-inflation
# This checks if there are more zeros in the observed data than predicted by the model.
# A p-value < 0.05 suggests significant zero-inflation not captured by the model.
testZeroInflation(simulated_residuals)

# You can also test for specific predictors if you suspect issues
# For example, residuals vs. 'prop_con_150_NTEMS'
plotResiduals(simulated_residuals, form = joined_data$prop_con_150_NTEMS)

# Or residuals vs. the fitted values
plotResiduals(simulated_residuals, form = fitted(model_cubic_ntems))




######## Plot model_3

# Create a sequence of prop_con_150_NTEMS values for prediction
new_prop_con_values <- seq(min(joined_data$prop_con_150_NTEMS, na.rm = TRUE),
                           max(joined_data$prop_con_150_NTEMS, na.rm = TRUE),
                           length.out = 100)

# Create a data frame for prediction
# Set other covariates to their mean, median, or a sensible reference value.
# For random effects, set to NA if you want to predict only fixed effects.
predict_data <- data.frame(
  prop_con_150_NTEMS = new_prop_con_values,
  age_mn_150_NTEMS = mean(joined_data$age_mn_150_NTEMS, na.rm = TRUE),
  latitude = mean(joined_data$latitude, na.rm = TRUE),
  longitude = mean(joined_data$longitude, na.rm = TRUE),
  # For offset, use the mean of survey_effort
  survey_effort = mean(joined_data$survey_effort, na.rm = TRUE),
  # For random effects, set to NA to ignore them in prediction (predict fixed effects only)
  cluster_k_10 = NA
)

# 2. Predict from the model
# type = "response" gives predictions on the original scale (e.g., counts for nbinom2)
# re.form = NA ensures predictions are only for the fixed effects (ignoring random effects)
predictions <- predict(model_3,
                       newdata = predict_data,
                       type = "response",
                       re.form = NA)

# Add predictions to the predict_data data frame
predict_data$predicted_BBWA <- predictions

# 3. Plot the model's predictions
plot_model <- ggplot(predict_data, aes(x = prop_con_150_NTEMS, y = predicted_BBWA)) +
  geom_line(color = "blue", size = 1) + # Plot the predicted line
  geom_point(data = joined_data, aes(x = prop_con_150_NTEMS, y = BBWA),
             alpha = 0.5, color = "darkgreen") + # Overlay actual data points
  labs(title = "Predicted BBWA Abundance vs. Prop_con_150_NTEMS",
       x = "Proportion of Coniferous Forest within 150m (NTEMS)",
       y = "Predicted BBWA Abundance") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) # Center the title

print(plot_model)














# BTNW models in function

library(glmmTMB)
library(dplyr)

# Vector of all prop_con_ variables you want to test
prop_vars <- c(
  "prop_con_150_NTEMS", "prop_con_500_NTEMS", "prop_con_1000_NTEMS",
  "prop_con_150_LULC", "prop_con_500_LULC", "prop_con_1000_LULC"
)

run_models_for_var <- function(var_name, data) {
  # Make a local copy
  data_mod <- data
  
  # Make sure cluster_k_10 is factor
  if(!is.factor(data_mod$cluster_k_10)) {
    data_mod$cluster_k_10 <- as.factor(data_mod$cluster_k_10)
  }
  
  # Create above_thresh column based on median of var_name
  data_mod$above_thresh <- as.numeric(data_mod[[var_name]] > median(data_mod[[var_name]], na.rm = TRUE))
  
  formulas <- list(
    linear = as.formula(paste0("BTNW ~ ", var_name, " + offset(log(survey_effort)) + (1 | cluster_k_10)")),
    quadratic = as.formula(paste0("BTNW ~ poly(", var_name, ", 2, raw=TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10)")),
    cubic = as.formula(paste0("BTNW ~ poly(", var_name, ", 3, raw=TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10)")),
    log = as.formula(paste0("BTNW ~ log(", var_name, " + 0.001) + offset(log(survey_effort)) + (1 | cluster_k_10)")),
    quadratic_exp = as.formula(paste0("BTNW ~ poly(", var_name, ", 2, raw=TRUE) + offset(log(survey_effort)) + (1 | cluster_k_10)")),
    piecewise = as.formula(paste0("BTNW ~ ", var_name, " * above_thresh + offset(log(survey_effort)) + (1 | cluster_k_10)"))
  )
  
  fits <- lapply(formulas, function(fml) {
    glmmTMB(fml, family = nbinom2, data = data_mod)
  })
  
  aic_vals <- sapply(fits, AIC)
  
  data.frame(
    variable = var_name,
    model_type = names(aic_vals),
    AIC = aic_vals,
    row.names = NULL
  )
}

# Run for all variables and combine results
all_results <- do.call(rbind, lapply(prop_vars, run_models_for_var, data = joined_data))

# View sorted by AIC to find best models
all_results <- all_results %>% arrange(AIC)
print(all_results)


summary(model_top)

### Plot model

library(glmmTMB)

# Fit cubic model for prop_con_1000_NTEMS
model_top <- glmmTMB(
  BTNW ~ poly(prop_con_1000_NTEMS, 3, raw = TRUE) + latitude + longitude + latitude*longitude +
    offset(log(survey_effort)) + (1 | cluster_k_10),
  family = nbinom2,
  data = joined_data
)

# Create prediction data frame across observed range
newdata <- data.frame(
  prop_con_1000_NTEMS = seq(
    min(joined_data$prop_con_1000_NTEMS, na.rm = TRUE),
    max(joined_data$prop_con_1000_NTEMS, na.rm = TRUE),
    length.out = 200
  ),
  latitude = mean(joined_data$latitude, na.rm = TRUE),
  longitude = mean(joined_data$longitude, na.rm = TRUE),
  survey_effort = 1,
  cluster_k_10 = NA
)

# Predict fitted values
newdata$pred <- predict(model_top, newdata = newdata, type = "response")


library(ggplot2)

library(ggplot2)

ggplot(newdata, aes(x = prop_con_1000_NTEMS, y = pred)) +
  geom_line(color = "blue", linewidth = 1.2) +
  labs(
    x = "Proportion Conifer (1000m, NTEMS)",
    y = "Predicted BTNW Abundance",
    title = "Cubic Model Fit: BTNW ~ prop_con_1000_NTEMS + Spatial Terms"
  ) +
  theme_minimal()



