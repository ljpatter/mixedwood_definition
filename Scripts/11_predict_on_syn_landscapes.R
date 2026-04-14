# ---
# title: "Predictions on synthetic landscapes"
# author: "Leonard Patterson"
# created: "2025-08-11"
# description: This script output from final multiscale BRT model to make abundance predictions for focal 
#              species on synthetic intimate and segregated mixed wood forests.
# ---

# -----------------
# 1. Load package
# -----------------
library(dplyr)
library(tidyr)
library(readr)
library(sf)
library(gbm)
library(stringr)
library(ggplot2)

# -----------------
# 2. Settings
# -----------------
TOWNSHIP_HA <- 1000   # Area of the township for scaling
CI_LO <- 0.025        # 95% Confidence Interval (Lower)
CI_HI <- 0.975        # 95% Confidence Interval (Upper)

species_order  <- c("BTNW", "TEWA", "BBWA")
layout_levels  <- c("Intimate", "Segregated")
age_vals       <- c(30, 60, 90)


# -----------------
# 3. Load point count locations and bootstrap model data
# -----------------
points_sf <- st_read("Output/Spatial Data/Simulated landscapes/point_grid_2km_inset1km.shp")   # Point counts
boot_model_list <- readRDS("Output/Tabular Data/bootstrapped_multiscale_model_list_5m.rds")      # Bootstrap models

# -----------------
# 4. Preconditions & Data Loading
# -----------------
# Check for your specific model object and spatial points
if (!exists("boot_model_list")) stop("Object 'boot_model_list' not found.")
if (!exists("points_sf")) stop("Object 'points_sf' not found.")

# Load spatial metrics if not in memory
if (!exists("metrics_points_df")) {
  metrics_path <- "Output/Tabular Data/propcon_clumpy_by_point_extent.rds"
  if (file.exists(metrics_path)) {
    metrics_points_df <- readRDS(metrics_path)
  } else {
    stop("Metric data file not found at: ", metrics_path)
  }
}

# Importance table for SoE predictor selection
imp_table_path <- "Output/Tabular Data/imp_table_5m.rds"
if (!file.exists(imp_table_path)) stop("Importance table missing.")
imp_use <- readRDS(imp_table_path)

# -----------------
# 5. Utilities
# -----------------
`%||%` <- function(a,b) if (is.null(a)) b else a

# Function to predict density (Birds per Hectare)
# This ignores offsets because the model scale (link) + exp = Density
predict_density <- function(model, newdata) {
  ntrees <- model$gbm.call$best.trees %||% model$n.trees
  if (is.null(ntrees) || is.na(ntrees) || ntrees <= 0) return(rep(NA_real_, nrow(newdata)))
  
  # QPAD logic: link scale is log(Density). Exp(link) is Density.
  eta <- gbm::predict.gbm(model, newdata = newdata, n.trees = ntrees, type = "link")
  exp(eta)
}

# Variable name builder
varname_for_soe <- function(prefix, extent, suffix) sprintf("%s_%d_%s", prefix, extent, suffix)

# -----------------
# 6. Prepare Predictors (SoE and Coordinates)
# -----------------
# Extract the best Scale of Effect (SoE) for each species/metric
soe_tbl <- imp_use %>%
  mutate(var_type = case_when(
    str_starts(Feature, "prop_con") ~ "prop_con",
    str_starts(Feature, "clumpy")   ~ "clumpy",
    str_starts(Feature, "age_mn")   ~ "age_mn",
    TRUE ~ "other"
  )) %>%
  filter(var_type %in% c("prop_con", "clumpy", "age_mn")) %>%
  group_by(Species, var_type) %>%
  slice_max(MeanGain, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    soe_extent = as.integer(sub(".*_(\\d+)_.*$", "\\1", Feature)),
    soe_suffix = sub("^.*_\\d+_([A-Za-z0-9]+)$", "\\1", Feature)
  )

# Point Coordinate Setup
coords <- st_coordinates(points_sf)
pts_key <- tibble(x = as.numeric(coords[,1]), y = as.numeric(coords[,2])) %>%
  mutate(xr = round(x, 3), yr = round(y, 3))

metrics_points_df_join <- metrics_points_df %>%
  mutate(xr = round(as.numeric(x), 3), yr = round(as.numeric(y), 3))

# Metric extractor function
get_metric_vec <- function(layout_name, target_extent, metric_name) {
  filt <- metrics_points_df_join %>%
    filter(layout == layout_name, extent == target_extent)
  joined <- pts_key %>% left_join(filt, by = c("xr", "yr"))
  out <- joined[[metric_name]]
  if (anyNA(out)) stop("Alignment error in spatial metrics.")
  as.numeric(out)
}

# -----------------
# 7. Prediction Loop
# -----------------
results <- list()

for (sp in intersect(names(boot_model_list), species_order)) {
  message("Processing Species: ", sp)
  models_sp <- boot_model_list[[sp]]
  
  # Determine Species SoE extents
  s_info <- soe_tbl %>% filter(Species == sp)
  v_prop_name <- s_info %>% filter(var_type == "prop_con") %>% slice(1) %>% pull(Feature)
  v_clmp_name <- s_info %>% filter(var_type == "clumpy")   %>% slice(1) %>% pull(Feature)
  v_age_name  <- s_info %>% filter(var_type == "age_mn")   %>% slice(1) %>% pull(Feature)
  
  ext_prop <- s_info %>% filter(var_type == "prop_con") %>% slice(1) %>% pull(soe_extent)
  ext_clmp <- s_info %>% filter(var_type == "clumpy")   %>% slice(1) %>% pull(soe_extent)
  
  for (layout in layout_levels) {
    prop_vals <- get_metric_vec(layout, ext_prop, "prop_con")
    clmp_vals <- get_metric_vec(layout, ext_clmp, "clumpy")
    
    for (age_val in age_vals) {
      boot_means <- vapply(models_sp, function(m) {
        if (is.null(m)) return(NA_real_)
        
        need <- m$gbm.call$predictor.names
        
        # Skip bootstrap if it does not contain the exact SoE predictors
        if (!(v_prop_name %in% need && v_clmp_name %in% need && v_age_name %in% need)) {
          return(NA_real_)
        }
        
        nd <- pts_key %>% transmute(x_AEP10TM = x, y_AEP10TM = y)
        nd[[v_prop_name]] <- prop_vals
        nd[[v_clmp_name]] <- clmp_vals
        nd[[v_age_name]]  <- age_val
        
        if ("year" %in% need) {
          train_df <- m$gbm.call$dataframe
          top_year <- names(sort(table(train_df$year), decreasing = TRUE))[1]
          nd$year <- factor(top_year, levels = levels(train_df$year))
        }
        
        miss <- setdiff(need, names(nd))
        for (v in miss) {
          mdf <- m$gbm.call$dataframe
          if (is.factor(mdf[[v]])) {
            nd[[v]] <- factor(levels(mdf[[v]])[1], levels = levels(mdf[[v]]))
          } else {
            nd[[v]] <- 0
          }
        }
        
        dens_vec <- predict_density(m, nd[, need, drop = FALSE])
        mean(dens_vec, na.rm = TRUE)
      }, numeric(1))
      
      valid_boots <- boot_means[is.finite(boot_means)]
      if (length(valid_boots) == 0) next
      
      m_dens <- mean(valid_boots)
      lwr_dens <- as.numeric(quantile(valid_boots, CI_LO))
      upr_dens <- as.numeric(quantile(valid_boots, CI_HI))
      
      results[[paste(sp, layout, age_val, sep = "_")]] <- tibble(
        Species = sp,
        layout = layout,
        age = age_val,
        mean_density_ha = m_dens,
        lwr_density_ha = lwr_dens,
        upr_density_ha = upr_dens,
        total_abundance = m_dens * TOWNSHIP_HA,
        lwr_abundance = lwr_dens * TOWNSHIP_HA,
        upr_abundance = upr_dens * TOWNSHIP_HA
      )
    }
  }
}

# -----------------
# 8. Output & Visualization
# -----------------

# 1. Combine list into dataframe
pred_summary_clean <- bind_rows(results)

# 2. Check if dataframe is empty before proceeding
if (nrow(pred_summary_clean) == 0) {
  stop("The 'results' list is empty. Check the prediction loop logic.")
}

# 3. Ensure Factors and Column Names match ggplot exactly
# Note: We use 'age' here because that is what was assigned in the results list above
pred_summary_clean <- pred_summary_clean %>%
  mutate(
    Species = factor(Species, levels = species_order),
    layout  = factor(layout, levels = layout_levels),
    age_f   = factor(age, levels = age_vals) # Creating a factor for discrete x-axis
  )

# 4. Plotting using the verified column names from the results list
p <- ggplot(pred_summary_clean, aes(x = age_f, y = total_abundance, fill = layout)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = lwr_abundance, ymax = upr_abundance),
    position = position_dodge(width = 0.7),
    width = 0.25, linewidth = 0.5
  ) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_fill_manual(values = c("Intimate" = "forestgreen", "Segregated" = "cornflowerblue")) +
  labs(
    x = "Forest age (years)",
    y = "Predicted mean abundance",
    fill = "Spatial configuration"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.text       = element_text(size = 13, face = "bold"),
    # Add these lines to force black text:
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    axis.ticks       = element_line(color = "black") # Optional: adds tick marks
  )

print(p)
out_fig <- "Output/Figures/synth_preds_SoE_only_5m.tiff"
ggsave(
  filename    = out_fig,
  plot        = p,
  device      = "tiff",
  width       = 42,
  height      = 18,
  units       = "cm",
  dpi         = 600,
  compression = "lzw"
)
