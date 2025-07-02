# ---
# title: "BRT model"
# author: "Leonard Patterson"
# created: "2025-04-11"
# description: This script runs boosted regression model with boostrapping and associated plots and diagnostics
# ---


# === LIBRARIES ===

library(cowplot)
library(sf)
library(sp)
library(dplyr)
library(gbm)
library(dismo)
library(QPAD)
library(spdep)
library(pROC)
library(purrr)
library(ggplot2)
library(stringr)
library(readr)
library(doParallel)
library(foreach)
library(tidyr)
library(viridis)
library(knitr)
library(tibble)
library(terra)
library(dplyr)
library(landscapemetrics)
library(patchwork)        



# Load study area
study_area <- st_read("Input/Spatial Data/Study Area/study_area.shp")

# Generate a 1 km grid over the study area
grid_1km <- st_make_grid(study_area, cellsize = 1000, square = TRUE) %>%
  st_sf(grid_id = 1:length(.), geometry = .)

# Convert point data to sf and assign CRS
coordinates(joined_data) <- ~x_AEP10TM + y_AEP10TM
proj4string(joined_data) <- CRS(st_crs(study_area)$proj4string)
joined_sf <- st_as_sf(joined_data)

# Filter points inside the grid extent
joined_sf <- joined_sf[st_within(joined_sf, st_union(grid_1km), sparse = FALSE) %>% apply(1, any), ]

# Spatial join: assign grid cell (block) to each point
joined_with_blocks <- st_join(joined_sf, grid_1km, join = st_within)
joined_with_blocks <- joined_with_blocks %>% filter(!is.na(grid_id))
joined_with_blocks$block_id <- as.factor(joined_with_blocks$grid_id)

# Convert back to SpatialPointsDataFrame
joined_sp <- as(joined_with_blocks, "Spatial")

# Restore coordinate column names to expected format
colnames(joined_sp@coords) <- c("x_AEP10TM", "y_AEP10TM")

















###################### BRT POISSON W/ PARALLELL BOOTSTRAPPING


# === SETUP ===
species_list <- c("BTNW", "TEWA", "BBWA")
results_list <- list()
friedmanH_list <- list()
final_model_list <- list()

# === Start Parallel Backend ===
cores_to_use <- 10
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

# === Offset Function ===
load_BAM_QPAD(version = 2)
calculate_offsets <- function(data, species_code) {
  localBAMcorrections(
    species = species_code,
    t = data$survey_effort,
    r = Inf,
    jday = data$ordinalDay,
    tssr = data$hssr
  ) %>% corrections2offset()
}

# === Friedman H Wrapper ===
.gbm.interactions_full <- function(model, model_data) {
  var_names <- model$gbm.call$predictor.names
  var_pairs <- combn(var_names, 2, simplify = FALSE)
  
  interaction_list <- lapply(var_pairs, function(pair) {
    score <- tryCatch({
      gbm::interact.gbm(
        x = model,
        data = model_data,
        i.var = pair,
        n.trees = model$n.trees
      )
    }, error = function(e) NA_real_)
    
    tibble(var1 = pair[1], var2 = pair[2], H = score)
  })
  
  bind_rows(interaction_list)
}

# === BOOTSTRAPPED BRT LOOP ===
set.seed(123)
n_boot <- 25

for (sp in species_list) {
  cat("\n\n=== Species:", sp, "===\n")
  
  joined_sp@data$offset <- calculate_offsets(joined_sp@data, sp)
  
  # 🔒 Force 'year' to be a factor globally
  joined_sp@data$year <- as.factor(joined_sp@data$year)
  
  boot_results <- foreach(i = 1:n_boot, .packages = c("dismo", "gbm", "pROC", "dplyr", "tidyr", "tibble")) %dopar% {
    sampled_blocks <- sample(unique(joined_sp@data$block_id), replace = TRUE)
    sp_subset <- joined_sp[joined_sp@data$block_id %in% sampled_blocks, ]
    block_subset <- cbind(as.data.frame(sp_subset@coords), sp_subset@data)
    
    block_subset$year <- as.factor(block_subset$year)  # 👈 enforce again inside parallel loop
    
    sampled_data <- data.frame(
      response = block_subset[[sp]],
      offset = block_subset$offset,
      year = block_subset$year,
      x_AEP10TM = block_subset$x_AEP10TM,
      y_AEP10TM = block_subset$y_AEP10TM,
      block_subset %>%
        dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
        dplyr::select(where(is.numeric))
    )
    
    if (sum(sampled_data$response, na.rm = TRUE) == 0) return(NULL)
    stopifnot(is.factor(sampled_data$year))  # 🔥 crash hard if not a factor
    
    folds <- sample(rep(1:10, length.out = nrow(sampled_data)))
    
    brt_model <- tryCatch({
      gbm.step(
        data = sampled_data,
        gbm.y = 1,
        gbm.x = 3:ncol(sampled_data),
        family = "poisson",
        tree.complexity = 3,
        learning.rate = 0.005,
        bag.fraction = 0.5,
        offset = sampled_data$offset,
        fold.vector = folds,
        n.folds = 10,
        keep.fold.models = FALSE,
        keep.fold.fit = FALSE,
        verbose = FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(brt_model)) return(NULL)
    
    imp_df <- data.frame(
      Feature = summary(brt_model)$var,
      Gain = summary(brt_model)$rel.inf
    )
    
    model_data <- sampled_data %>%
      dplyr::select(all_of(brt_model$gbm.call$predictor.names)) %>%
      mutate(across(everything(), ~tidyr::replace_na(., 0)))
    
    interaction_df <- .gbm.interactions_full(brt_model, model_data)
    interaction_df$Boot <- i
    
    list(imp = imp_df, interaction = interaction_df)
  }
  
  boot_results <- boot_results[!sapply(boot_results, is.null)]
  
  imp_all <- lapply(boot_results, function(x) x$imp)
  int_all <- lapply(boot_results, function(x) x$interaction)
  
  results_list[[sp]] <- bind_rows(imp_all) %>%
    group_by(Feature) %>%
    summarise(MeanGain = mean(Gain, na.rm = TRUE), .groups = "drop") %>%
    mutate(Species = sp)
  
  friedmanH_list[[sp]] <- bind_rows(int_all) %>%
    group_by(var1, var2) %>%
    summarise(H = mean(H, na.rm = TRUE), .groups = "drop") %>%
    mutate(Species = sp)
  
  # === Final model on full data ===
  joined_sp@data$year <- as.factor(joined_sp@data$year)  # 🔁 ensure factor again
  
  full_data <- cbind(as.data.frame(joined_sp@coords), joined_sp@data)
  model_data <- data.frame(
    response = full_data[[sp]],
    offset = full_data$offset,
    year = full_data$year,
    x_AEP10TM = full_data$x_AEP10TM,
    y_AEP10TM = full_data$y_AEP10TM,
    full_data %>%
      dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
      dplyr::select(where(is.numeric))
  )
  
  stopifnot(is.factor(model_data$year))  # 🔥 hard fail if not a factor
  
  final_model <- gbm.step(
    data = model_data,
    gbm.y = 1,
    gbm.x = 3:(ncol(model_data)-1),
    family = "poisson",
    tree.complexity = 3,
    learning.rate = 0.005,
    bag.fraction = 0.5,
    offset = model_data$offset,
    n.folds = 10,
    keep.fold.models = FALSE,
    keep.fold.fit = FALSE,
    verbose = FALSE
  )
  
  final_model_list[[sp]] <- final_model
}

# === Stop Parallel Cluster ===
stopCluster(cl)

# === Save Results ===
imp_df <- bind_rows(results_list)
friedman_df <- bind_rows(friedmanH_list)

saveRDS(results_list, "Output/Tabular Data/results_list_final_parallel_poisson.rds")
saveRDS(imp_df, "Output/Tabular Data/imp_df_final_parallel_poisson.rds")
saveRDS(friedman_df, "Output/Tabular Data/friedman_interactions_bootstrapped_parallel_final_poisson.rds")
saveRDS(final_model_list, "Output/Tabular Data/final_model_list_parallel_poisson.rds")

cat("\n✅ Saved all model outputs for Poisson models with year as FACTOR!\n")











##################### Poisson BRT Variable Importance Analysis #####################


# === Load model object to verify 'year' as factor ===
final_model_list <- readRDS("Output/Tabular Data/final_model_list_parallel_poisson.rds")

for (sp in names(final_model_list)) {
  cat(paste0("\n🔍 Checking model for: ", sp, "\n"))
  model_data <- final_model_list[[sp]]$gbm.call$dataframe
  if (!"year" %in% names(model_data)) {
    stop(paste0("❌ ERROR: 'year' not found in model data for ", sp))
  }
  if (!is.factor(model_data$year)) {
    stop(paste0("❌ ERROR: 'year' is NOT a factor for ", sp))
  }
  cat("✅ 'year' is included and treated as a FACTOR.\n")
}

# === Load variable importance results ===
imp_df <- readRDS("imp_df_final_parallel_poisson.rds")

# === Nuisance variables to exclude ===
nuisance_vars <- c("year", "x_AEP10TM", "y_AEP10TM")

# === Color map and ordering ===
fill_map <- c(
  "prop_con" = "forestgreen",
  "clumpy"   = "cornflowerblue",
  "age_mn"   = "lightcoral"
)
var_order <- c("prop_con", "clumpy", "age_mn")

# === Clean, classify, and round ===
imp_df_clean <- imp_df %>%
  filter(!Feature %in% nuisance_vars) %>%
  filter(!grepl("^age_mn_\\d+_LULC$", Feature)) %>%
  mutate(
    var_type = case_when(
      str_detect(Feature, "prop_con") ~ "prop_con",
      str_detect(Feature, "clumpy")   ~ "clumpy",
      str_detect(Feature, "age_mn")   ~ "age_mn",
      TRUE                            ~ "other"
    ),
    var_type = factor(var_type, levels = var_order)
  )

# === Print Variable Importance Tables ===
imp_table <- imp_df_clean %>%
  arrange(Species, desc(MeanGain)) %>%
  mutate(MeanGain = round(MeanGain, 2))

cat("\n🔎 Variable Importance Table by Species:\n")
print(imp_table, n = Inf)

species_list <- unique(imp_table$Species)

for (sp in species_list) {
  cat(paste0("\n=== Variable Importance for ", sp, " ===\n"))
  print(imp_table %>% filter(Species == sp), n = Inf)
}

# === Plot VI per species ===
max_gain <- ceiling(max(imp_df_clean$MeanGain, na.rm = TRUE) / 5) * 5

for (sp in species_list) {
  df_sp <- imp_df_clean %>%
    filter(Species == sp) %>%
    group_by(var_type) %>%
    arrange(desc(MeanGain), .by_group = TRUE) %>%
    ungroup() %>%
    mutate(Feature = factor(Feature, levels = rev(Feature)))
  
  if (nrow(df_sp) == 0) {
    cat(paste0("⚠️ No variables found for ", sp, ". Skipping plot.\n"))
    next
  }
  
  p <- ggplot(df_sp, aes(x = Feature, y = MeanGain, fill = var_type)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_manual(
      values = fill_map,
      drop = FALSE,
      name = "Predictor Type",
      labels = c(
        "prop_con" = "Proportion Conifer",
        "clumpy" = "Clumpiness",
        "age_mn" = "Forest Age"
      )
    ) +
    scale_y_continuous(limits = c(0, max_gain)) +
    coord_flip() +
    labs(
      title = paste0(sp, " (Poisson BRT)"),
      x = "Predictor Variable",
      y = "Mean Relative Influence"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(hjust = 1, size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.position = "right"  # ✅ Show legend here
    )
  
  print(p)  # ✅ This line ensures the plot actually renders
}












# === FULL PDPs FOR ALL VARIABLES USED IN FINAL MODELS ===


# Load saved objects
results_list <- readRDS("Output/Tabular Data/results_list_final_parallel_poisson.rds")
final_model_list <- readRDS("Output/Tabular Data/final_model_list_parallel_poisson.rds")

# Define nuisance variables to exclude
nuisance_vars <- c("x_AEP10TM", "y_AEP10TM", "year")

# Loop through each species
for (sp in names(final_model_list)) {
  cat("\n=== PDPs for", sp, "===\n")
  
  model <- final_model_list[[sp]]
  var_importance <- results_list[[sp]]
  
  # 🔒 Check that year is in model and is a factor
  if ("year" %in% names(model$gbm.call$dataframe)) {
    is_factor <- is.factor(model$gbm.call$dataframe$year)
    cat("✅ 'year' is present and is.factor(year) =", is_factor, "\n")
    stopifnot(is_factor)
  } else {
    warning("'year' is not in the model data — skipping check")
  }
  
  # Get all predictor names used in the final model (excluding nuisance)
  valid_vars <- model$gbm.call$predictor.names
  valid_vars <- valid_vars[!valid_vars %in% nuisance_vars]
  
  for (var in valid_vars) {
    # Safely generate PDP
    pd <- try(gbm::plot.gbm(model, i.var = var, return.grid = TRUE), silent = TRUE)
    
    if (inherits(pd, "try-error") || !("y" %in% names(pd)) || all(is.na(pd$y))) {
      message("Skipping variable: ", var, " (could not generate valid PDP)")
      next
    }
    
    # Center y values on log scale
    pd$y_centered <- pd$y - mean(pd$y, na.rm = TRUE)
    
    # Extract relative contribution
    contrib <- var_importance %>%
      filter(Feature == var) %>%
      pull(MeanGain) %>%
      round(1)
    
    title_text <- paste0(sp, " - ", var, " (", contrib, "%)")
    
    # Plot
    p <- ggplot(pd, aes_string(x = var, y = "y_centered")) +
      geom_line(color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(
        title = title_text,
        x = var,
        y = "Centered Partial Dependence (log scale)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
    
    print(p)
  }
}

















################## AUC and model deviance explained


# If not already loaded:
# dev_list <- readRDS("Output/Tabular Data/dev_list_parallel.rds")  # You saved this one

# Summarize deviance explained
dev_summary <- lapply(names(dev_list), function(sp) {
  tibble(
    Species = sp,
    MeanDeviance = round(mean(dev_list[[sp]], na.rm = TRUE), 2),
    SDDeviance = round(sd(dev_list[[sp]], na.rm = TRUE), 2)
  )
}) %>% bind_rows()

# Summarize AUC (if auc_list is still in memory)
if (exists("auc_list")) {
  auc_summary <- lapply(names(auc_list), function(sp) {
    tibble(
      Species = sp,
      MeanAUC = round(mean(auc_list[[sp]], na.rm = TRUE), 3),
      SDAUC = round(sd(auc_list[[sp]], na.rm = TRUE), 3)
    )
  }) %>% bind_rows()
  
  # Join both summaries
  model_summary <- left_join(dev_summary, auc_summary, by = "Species")
  print(model_summary)
  
} else {
  message("⚠️ 'auc_list' not found in memory. You must rerun the model loop or restore auc_list from memory if saved.")
  print(dev_summary)
}














################# H-statistics top interactions


# Load interactions
all_friedman_interactions <- readRDS("Output/Tabular Data/friedman_interactions_bootstrapped_parallel_final_poisson.rds")

# Updated whitelist WITHOUT "predictors."
allowed_pairs <- list(
  c("prop_con_150_NTEMS", "clumpy_150_NTEMS"),
  c("prop_con_150_NTEMS", "age_mn_150_NTEMS"),
  c("prop_con_500_NTEMS", "clumpy_500_NTEMS"),
  c("prop_con_500_NTEMS", "age_mn_500_NTEMS"),
  c("prop_con_1000_NTEMS", "clumpy_1000_NTEMS"),
  c("prop_con_1000_NTEMS", "age_mn_1000_NTEMS"),
  c("prop_con_150_LULC", "clumpy_150_LULC"),
  c("prop_con_150_LULC", "age_mn_150_NTEMS"),
  c("prop_con_500_LULC", "clumpy_500_LULC"),
  c("prop_con_500_LULC", "age_mn_500_NTEMS"),
  c("prop_con_1000_LULC", "clumpy_1000_LULC"),
  c("prop_con_1000_LULC", "age_mn_1000_NTEMS")
)

# Filter to whitelist only
filtered_df <- all_friedman_interactions %>%
  rowwise() %>%
  filter(any(sapply(allowed_pairs, function(pair) {
    (var1 == pair[1] & var2 == pair[2]) |
      (var1 == pair[2] & var2 == pair[1])
  }))) %>%
  ungroup()

# Top 5 per species
top5_filtered <- filtered_df %>%
  group_by(Species) %>%
  arrange(desc(H)) %>%
  slice_head(n = 5) %>%
  ungroup()

# Show the pretty table
kable(top5_filtered, caption = "Top 5 Friedman H Interactions per Species (Strict Filter)")


write.csv(top5_filtered, "Output/Tabular Data/interactions_table.csv")

















##################### 3d Interation Plots


final_model_list <- readRDS("Output/Tabular Data/final_model_list_parallel.rds")


# Helper: extract range of a predictor from the model
get_var_range <- function(model, var_index) {
  var_name <- model$gbm.call$predictor.names[var_index]
  data <- model$gbm.call$dataframe
  range(data[[var_name]], na.rm = TRUE)
}

# BBWA age × prop_con
gbm.perspec(
  gbm.object = final_model_list[["BBWA"]],
  x = 15,  # prop_con_150_LULC
  y = 6, # age_mn_150_NTEMS
  y.range = get_var_range(final_model_list[["BBWA"]], 6),
  z.range = c(0, 0.6),
  main = "BBWA: prop_con_150_LULC × age_mn_150_NTEMS"
)



# TEWA clumpy × prop_con
gbm.perspec(
  gbm.object = final_model_list[["TEWA"]],
  x = 8,  # prop_con_150_LULC
  y = 14, # age_mn_150_NTEMS
  y.range = get_var_range(final_model_list[["TEWA"]],14),
  z.range = c(0, 1),
  main = "TEWA: prop_con_1000_LULC × clumpy_1000_NTEMS"
)

# Or this version
final_model_list[["TEWA"]]$gbm.call$predictor.names



#BTNW

# Extract variable names from the model
x_var <- final_model_list[["BTNW"]]$gbm.call$predictor.names[11]
y_var <- final_model_list[["BTNW"]]$gbm.call$predictor.names[5]

# Use your original model data (the one used to train the model)
model_data <- final_model_list[["BTNW"]]$gbm.call$dataframe

# Automatically get the range for y.axis
y_range_auto <- range(model_data[[y_var]], na.rm = TRUE)

# Generate the 3D interaction plot
gbm.perspec(
  gbm.object = final_model_list[["BTNW"]],
  x = 11,
  y = 5,
  y.range = y_range_auto,
  z.range = c(0, 1),
  main = paste("BTNW:", x_var, "×", y_var)
)













################ Moran's I for Poisson Final Models ################

# === Coordinates setup ===
coords_sp <- coordinates(joined_sp)
rownames(coords_sp) <- rownames(joined_sp@data)

# Build neighborhood
nb <- dnearneigh(coords_sp, 0, 30000)
nb <- make.sym.nb(nb)
has_neighbors <- which(card(nb) > 0)

coords_use <- coords_sp[has_neighbors, , drop = FALSE]
data_use <- joined_sp@data[has_neighbors, ]

# Ensure 'year' is factor
data_use$year <- as.factor(data_use$year)

# Rebuild neighborhood and weights
nb_use <- dnearneigh(coords_use, 0, 30000)
nb_use <- make.sym.nb(nb_use)
lw_use <- nb2listw(nb_use, style = "W", zero.policy = TRUE)

# === Moran's I Results ===
moran_results <- list()

for (sp in species_list) {
  cat(paste0("\n=== Moran's I for ", sp, " ===\n"))
  
  model <- final_model_list[[sp]]
  if (is.null(model)) {
    warning(paste("Model for", sp, "is NULL. Skipping."))
    moran_results[[sp]] <- NA
    next
  }
  
  # Get raw response and required predictors
  data_use$response <- data_use[[sp]]
  data_use$x_AEP10TM <- coords_use[, 1]
  data_use$y_AEP10TM <- coords_use[, 2]
  
  # Ensure predictors match the model
  predictors_needed <- model$gbm.call$predictor.names
  pred_data <- data_use %>%
    dplyr::select(all_of(predictors_needed)) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  # Add offset and predict on link scale
  preds_link <- predict(model, pred_data, n.trees = model$n.trees, type = "link")
  preds_mu <- exp(preds_link + data_use$offset)
  
  resids <- data_use$response - preds_mu
  
  # Quick check on residuals
  print(summary(resids))
  print(sum(is.na(resids)))
  
  # Skip invalid residuals
  if (all(is.na(resids)) || var(resids, na.rm = TRUE) == 0) {
    warning(paste("Residuals for", sp, "are all NA or zero variance."))
    moran_results[[sp]] <- NA
    next
  }
  
  # Moran’s I
  moran_val <- tryCatch({
    unname(moran.test(resids, lw_use, zero.policy = TRUE)$estimate[1])
  }, error = function(e) {
    warning(paste("Moran's I failed for", sp, ":", e$message))
    NA
  })
  
  moran_results[[sp]] <- moran_val
  print(moran_val)
}

cat("\n✅ === Moran's I (Final Poisson Model Residuals) ===\n")
print(moran_results)


















# === FULL PDPs FOR ALL VARIABLES USED IN FINAL MODELS ===


# Load saved model objects
results_list <- readRDS("Output/Tabular Data/results_list_final_parallel_poisson.rds")
final_model_list <- readRDS("Output/Tabular Data/final_model_list_parallel_poisson.rds")

# Define nuisance variables to exclude
nuisance_vars <- c("x_AEP10TM", "y_AEP10TM", "year")

# Loop through each species
for (sp in names(final_model_list)) {
  cat("\n=== PDPs for", sp, "===\n")
  
  model <- final_model_list[[sp]]
  var_importance <- results_list[[sp]]
  
  # Get all predictor names used in the final model (excluding nuisance)
  valid_vars <- model$gbm.call$predictor.names
  valid_vars <- valid_vars[!valid_vars %in% nuisance_vars]
  
  for (var in valid_vars) {
    # Safely generate PDP
    pd <- try(gbm::plot.gbm(model, i.var = var, return.grid = TRUE), silent = TRUE)
    
    if (inherits(pd, "try-error") || !("y" %in% names(pd)) || all(is.na(pd$y))) {
      message("Skipping variable: ", var, " (could not generate valid PDP)")
      next
    }
    
    # === Back-transform from log scale ===
    pd$expected_count <- exp(pd$y)
    
    # === Clean and center/scale ===
    pd <- pd %>%
      filter(is.finite(expected_count)) %>%
      mutate(expected_scaled = as.numeric(scale(expected_count))) %>%
      filter(!is.na(expected_scaled))
    
    # Extract relative contribution from variable importance
    contrib <- var_importance %>%
      filter(Feature == var) %>%
      pull(MeanGain) %>%
      round(1)
    
    title_text <- paste0(sp, " - ", var, " (", contrib, "% RI)")
    
    # Plot
    p <- ggplot(pd, aes_string(x = var, y = "expected_scaled")) +
      geom_line(color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(
        title = title_text,
        x = var,
        y = "Centered & Scaled Expected Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
    
    print(p)
  }
}











############## Synthetic landscape


# === Load model ===
final_model_list <- readRDS("Output/Tabular Data/final_model_list_parallel_poisson.rds")

# === Scale-of-effect predictors per species ===
scale_effect_predictors <- list(
  BTNW = c("prop_con_150_LULC", "clumpy_1000_NTEMS", "age_mn_150_NTEMS"),
  TEWA = c("prop_con_1000_LULC", "clumpy_1000_LULC", "age_mn_150_NTEMS"),
  BBWA = c("prop_con_150_LULC", "clumpy_150_LULC", "age_mn_150_NTEMS")
)

# === Function to generate prediction surface ===
predict_on_simulated_surface <- function(species, model, predictor_set, grid_steps = 50) {
  
  prop_vals <- seq(0, 1, length.out = grid_steps)
  clumpy_vals <- seq(0, 1, length.out = grid_steps)
  age_vals <- c(30, 60, 90, 120)
  
  grid <- expand.grid(
    prop_con = prop_vals,
    clumpy = clumpy_vals,
    age = age_vals
  )
  
  # Rename for prediction
  names(grid)[names(grid) == "prop_con"] <- predictor_set[1]
  names(grid)[names(grid) == "clumpy"] <- predictor_set[2]
  names(grid)[names(grid) == "age"] <- predictor_set[3]
  
  # Fill missing vars
  all_model_vars <- model$gbm.call$predictor.names
  missing_vars <- setdiff(all_model_vars, names(grid))
  for (v in missing_vars) {
    grid[[v]] <- 0
  }
  
  # Predict
  grid$prediction <- predict.gbm(
    object = model,
    newdata = grid,
    n.trees = model$gbm.call$best.trees,
    type = "response"
  )
  
  # Rename back for plotting
  names(grid)[names(grid) == predictor_set[1]] <- "prop_con"
  names(grid)[names(grid) == predictor_set[2]] <- "clumpy"
  names(grid)[names(grid) == predictor_set[3]] <- "age"
  
  grid$Species <- factor(species, levels = c("BTNW", "TEWA", "BBWA")) # Ensure correct order
  grid$AgeGroup <- factor(grid$age)
  
  return(grid)
}

# === Generate predictions for all species ===
simulated_predictions <- map_dfr(
  names(final_model_list),
  ~predict_on_simulated_surface(
    species = .x,
    model = final_model_list[[.x]],
    predictor_set = scale_effect_predictors[[.x]]
  )
)

# === Rescale within species
simulated_predictions <- simulated_predictions %>%
  group_by(Species) %>%
  mutate(scaled_prediction = (prediction - min(prediction)) / (max(prediction) - min(prediction))) %>%
  ungroup()

# === Final Plot with Magma Color Scale ===
# Create a new factor for the facet labels
simulated_predictions$Species_Label <- factor(simulated_predictions$Species, 
                                              levels = c("BTNW", "TEWA", "BBWA"),
                                              labels = c("a)", "b)", "c)"))

ggplot(simulated_predictions, aes(x = prop_con, y = clumpy, fill = scaled_prediction)) +
  geom_tile() +
  facet_grid(Species_Label ~ AgeGroup) + # Use the new factor here
  scale_fill_viridis_c(
    option = "magma",
    name = "Scaled\nPrediction",
    limits = c(0, 1)
  ) +
  labs(
    x = "Proportion Conifer",
    y = "Spatial Configuration",
    title = "                                                                                                Forest Age"
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.text.y = element_text(angle = 0, hjust = 1)) # Keep labels horizontal















library(terra)

# --- Create 10x10 raster matrix for each layout ---

# INTIMATE: Random mix of 0s and 1s with ~equal frequency for CLUMPY ≈ 0
set.seed(123)  # For reproducibility
intimate_mat <- matrix(sample(rep(c(0, 1), each = 50)), nrow = 10)

# SEGREGATED: Left half conifer (1), right half deciduous (0)
segregated_mat <- matrix(0, nrow = 10, ncol = 10)
segregated_mat[, 1:5] <- 1  # Left side = conifer

# --- Convert to SpatRaster objects ---

make_raster <- function(mat) {
  r <- rast(mat)
  ext(r) <- ext(0, 10, 0, 10)  # Set extent (arbitrary units)
  crs(r) <- "EPSG:32612"       # Projected CRS
  names(r) <- "conifer"
  return(r)
}

r_intimate <- make_raster(intimate_mat)
r_segregated <- make_raster(segregated_mat)

# --- Optional: Save for visual checks ---
writeRaster(r_intimate, "intimate_layout.tif", overwrite = TRUE)
writeRaster(r_segregated, "segregated_layout.tif", overwrite = TRUE)

















### Visualize raster


library(ggplot2)
library(terra)

# Load and convert raster to data.frame
r <- rast("segregated_layout.tif")
df <- as.data.frame(r, xy = TRUE)
names(df)[3] <- "cover"  # Rename layer to 'cover'

# Plot using ggplot
ggplot(df, aes(x = x, y = y, fill = factor(cover))) +
  geom_raster() +
  scale_fill_manual(
    values = c("lightgreen", "darkgreen"),
    labels = c("Deciduous", "Conifer"),
    name = "Cover Type"
  ) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(title = "Synthetic Segregated Forest Layout", x = "X", y = "Y")




library(ggplot2)
library(terra)

# Load and convert intimate raster to data.frame
r <- rast("intimate_layout.tif")
df <- as.data.frame(r, xy = TRUE)
names(df)[3] <- "cover"  # Rename raster layer to 'cover'

# Plot using ggplot
ggplot(df, aes(x = x, y = y, fill = factor(cover))) +
  geom_raster() +
  scale_fill_manual(
    values = c("lightgreen", "darkgreen"),
    labels = c("Deciduous", "Conifer"),
    name = "Cover Type"
  ) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(title = "Synthetic Intimate Forest Layout (CLUMPY ≈ 0)", x = "X", y = "Y")









### Calculate clumpy for intimate mixedwood

library(landscapemetrics)

# Ensure values are categorical (0 = deciduous, 1 = conifer)
r_int_cat <- classify(r_intimate, matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))
r_seg_cat <- classify(r_segregated, matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))

# Calculate CLUMPY for all classes
clumpy_intimate_all <- lsm_c_clumpy(r_int_cat)
clumpy_segregated_all <- lsm_c_clumpy(r_seg_cat)

# ✅ Filter to conifer pixels only (class == 1)
clumpy_val_intimate <- clumpy_intimate_all %>%
  dplyr::filter(class == 1) %>%
  dplyr::pull(value)

clumpy_val_segregated <- clumpy_segregated_all %>%
  dplyr::filter(class == 1) %>%
  dplyr::pull(value)

# Print results
print(clumpy_val_intimate)
print(clumpy_val_segregated)












library(terra)
library(dplyr)
library(gbm)

# Load your synthetic rasters
r_int <- rast("intimate_layout.tif")
r_seg <- rast("segregated_layout.tif")

# Assign names for consistency
names(r_int) <- "prop_con"
names(r_seg) <- "prop_con"

# Define conifer CLUMPY values from earlier calculations
clumpy_int_val <- 0  # randomized layout ~ 0
clumpy_seg_val <- 1  # perfectly segregated ~ 1

# Define age values
ages <- c(30, 60, 90)

# Expand raster into prediction dataframe
make_prediction_df <- function(r, clumpy_val, age_val) {
  df <- as.data.frame(r, xy = TRUE)
  df$clumpy <- clumpy_val
  df$age <- age_val
  return(df)
}

# Predict for one species on one raster at one age
predict_synthetic <- function(species, model, predictors, landscape_df) {
  df <- landscape_df %>%
    rename(
      !!predictors[1] := prop_con,
      !!predictors[2] := clumpy,
      !!predictors[3] := age
    )
  
  # Fill in zeros for other variables
  for (v in setdiff(model$gbm.call$predictor.names, names(df))) {
    df[[v]] <- 0
  }
  
  df$pred <- predict.gbm(model,
                         newdata = df,
                         n.trees = model$gbm.call$best.trees,
                         type = "response")
  
  df$Species <- species
  return(df)
}

# Loop over all combinations: species × layout × age
synthetic_preds <- list()

for (sp in names(final_model_list)) {
  model <- final_model_list[[sp]]
  predictors <- scale_effect_predictors[[sp]]
  
  for (age in ages) {
    
    # Create two data.frames: intimate and segregated
    df_int <- make_prediction_df(r_int, clumpy_val = clumpy_int_val, age_val = age)
    df_seg <- make_prediction_df(r_seg, clumpy_val = clumpy_seg_val, age_val = age)
    
    # Predict
    pred_int <- predict_synthetic(sp, model, predictors, df_int) %>% mutate(layout = "Intimate", age_class = age)
    pred_seg <- predict_synthetic(sp, model, predictors, df_seg) %>% mutate(layout = "Segregated", age_class = age)
    
    # Combine
    synthetic_preds[[paste(sp, age, "int", sep = "_")]] <- pred_int
    synthetic_preds[[paste(sp, age, "seg", sep = "_")]] <- pred_seg
  }
}

# Combine all into one data frame
synthetic_pred_df <- bind_rows(synthetic_preds)








### Plot

library(ggplot2)

plot_synthetic_predictions <- function(pred_df, species_code) {
  
  # Filter for one species
  df_plot <- pred_df %>%
    filter(Species == species_code)
  
  # Rescale within species for color contrast
  df_plot <- df_plot %>%
    mutate(scaled_pred = (pred - min(pred)) / (max(pred) - min(pred)))
  
  # Plot
  ggplot(df_plot, aes(x = x, y = y, fill = scaled_pred)) +
    geom_raster() +
    facet_grid(age_class ~ layout) +
    scale_fill_viridis_c(
      option = "magma",
      limits = c(0, 1),
      name = "Relative\nAbundance"
    ) +
    coord_equal() +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Predicted Abundance for", species_code),
      x = "X Coordinate", y = "Y Coordinate"
    )
}

plot_synthetic_predictions(synthetic_pred_df, "BBWA")
plot_synthetic_predictions(synthetic_pred_df, "TEWA")
plot_synthetic_predictions(synthetic_pred_df, "BTNW")














############### Synthetic aggregate and dispersed prediction plots



# Set font for graphics
theme_set(theme_minimal(base_size = 14) +
            theme(text = element_text(family = "")))

### Create 10x10 raster matrix for each configuration

# INTIMATE: Random mix of 0s and 1s with ~equal frequency for CLUMPY ≈ 0
set.seed(123) # For reproducibility
intimate_mat <- matrix(sample(c(0, 1), size = 100, replace = TRUE), nrow = 10)

# SEGREGATED: Left half conifer (1), right half deciduous (0)
segregated_mat <- matrix(0, nrow = 10, ncol = 10)
segregated_mat[, 1:5] <- 1 # Left side = conifer

# --- Convert to SpatRaster objects ---

make_raster <- function(mat) {
  r <- rast(mat)
  ext(r) <- ext(0, 10, 0, 10) # Set extent (arbitrary units)
  crs(r) <- "EPSG:32612"       # Projected CRS
  names(r) <- "conifer"        # Name the layer
  return(r)
}

r_intimate <- make_raster(intimate_mat)
r_segregated <- make_raster(segregated_mat)

# --- Optional: Save for visual checks (already in your original code) ---
# writeRaster(r_intimate, "intimate_layout.tif", overwrite = TRUE)
# writeRaster(r_segregated, "segregated_layout.tif", overwrite = TRUE)


# --- Calculate clumpy for intimate mixedwood ---

# Ensure values are categorical (0 = deciduous, 1 = conifer)
r_int_cat <- classify(r_intimate, matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))
r_seg_cat <- classify(r_segregated, matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))

# Calculate CLUMPY for all classes
clumpy_intimate_all <- lsm_c_clumpy(r_int_cat)
clumpy_segregated_all <- lsm_c_clumpy(r_seg_cat)

# Filter to conifer pixels only (class == 1)
clumpy_val_intimate <- clumpy_intimate_all %>%
  dplyr::filter(class == 1) %>%
  dplyr::pull(value)

clumpy_val_segregated <- clumpy_segregated_all %>%
  dplyr::filter(class == 1) %>%
  dplyr::pull(value)

# Print results (for verification, will be close to your assumed 0 and 1)
print(paste("Calculated CLUMPY for Intimate (conifer):", clumpy_val_intimate))
print(paste("Calculated CLUMPY for Segregated (conifer):", clumpy_val_segregated))

# Define conifer CLUMPY values from earlier calculations (or derived from above)
clumpy_int_val <- 0  # randomized layout ~ 0
clumpy_seg_val <- 1  # perfectly segregated ~ 1

# --- Define dummy model and predictors for demonstration ---
# IMPORTANT: Since final_model_list and scale_effect_predictors are not provided in your code,
# I'm creating dummy structures so the prediction loop can run.
# YOU MUST REPLACE THESE WITH YOUR ACTUAL LOADED MODELS AND PREDICTORS for real results!
if (!exists("final_model_list")) {
  message("WARNING: 'final_model_list' not found. Creating dummy models for demonstration.")
  message("         Please replace these with your actual trained gbm models for accurate results.")
  
  # Dummy training data (simple, illustrative)
  dummy_train_data <- data.frame(
    prop_con = sample(c(0, 1), 100, replace = TRUE),
    clumpy = sample(c(0, 1), 100, replace = TRUE, prob = c(0.5, 0.5)),
    age = sample(c(30, 60, 90), 100, replace = TRUE),
    other_var1 = runif(100),
    other_var2 = runif(100),
    abundance = runif(100, min = 0, max = 10) # Dummy abundance values
  )
  
  # A simple relationship: abundance increases with conifer, decreases slightly with clumpy
  dummy_train_data$abundance <- 5 + 3 * dummy_train_data$prop_con - 2 * dummy_train_data$clumpy + 0.1 * dummy_train_data$age
  dummy_train_data$abundance <- pmax(0, dummy_train_data$abundance) # Ensure non-negative predictions
  
  # Create a dummy gbm model for each species
  dummy_model_base <- gbm(abundance ~ prop_con + clumpy + age + other_var1 + other_var2,
                          data = dummy_train_data,
                          distribution = "gaussian",
                          n.trees = 10, # Very few trees for dummy, your actual models will have many more
                          interaction.depth = 1,
                          verbose = FALSE)
  # Set best.trees (this would be determined by cross-validation in a real scenario)
  dummy_model_base$gbm.call$best.trees <- 10 
  
  final_model_list <- list(
    "BBWA" = dummy_model_base,
    "TEWA" = dummy_model_base, # Using the same dummy for all for demonstration
    "BTNW" = dummy_model_base
  )
  
  scale_effect_predictors <- list(
    "BBWA" = c("prop_con", "clumpy", "age"),
    "TEWA" = c("prop_con", "clumpy", "age"),
    "BTNW" = c("prop_con", "clumpy", "age")
  )
}
# END DUMMY MODEL DEFINITION

# --- Define age values ---
ages <- c(30, 60, 90)

# Expand raster into prediction dataframe
make_prediction_df <- function(r, clumpy_val, age_val) {
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "prop_con" # Ensure the column name is 'prop_con' for the model
  df$clumpy <- clumpy_val
  df$age <- age_val
  return(df)
}

# Predict for one species on one raster at one age
predict_synthetic <- function(species, model, predictors, landscape_df) {
  df <- landscape_df %>%
    rename(
      !!predictors[1] := prop_con,
      !!predictors[2] := clumpy,
      !!predictors[3] := age
    )
  
  # Fill in zeros for other variables expected by the model but not in landscape_df
  for (v in setdiff(model$gbm.call$predictor.names, names(df))) {
    df[[v]] <- 0
  }
  
  df$pred <- predict.gbm(model,
                         newdata = df,
                         n.trees = model$gbm.call$best.trees,
                         type = "response")
  
  df$Species <- species
  return(df)
}

# Loop over all combinations: species × layout × age
synthetic_preds <- list()

for (sp in names(final_model_list)) {
  model <- final_model_list[[sp]]
  predictors <- scale_effect_predictors[[sp]]
  
  for (age in ages) {
    
    # Create two data.frames: intimate and segregated
    df_int <- make_prediction_df(r_intimate, clumpy_val = clumpy_int_val, age_val = age)
    df_seg <- make_prediction_df(r_segregated, clumpy_val = clumpy_seg_val, age_val = age)
    
    # Predict
    pred_int <- predict_synthetic(sp, model, predictors, df_int) %>% mutate(layout = "Intimate", age_class = age)
    pred_seg <- predict_synthetic(sp, model, predictors, df_seg) %>% mutate(layout = "Segregated", age_class = age)
    
    # Combine
    synthetic_preds[[paste(sp, age, "int", sep = "_")]] <- pred_int
    synthetic_preds[[paste(sp, age, "seg", sep = "_")]] <- pred_seg
  }
}

# Combine all into one data frame
synthetic_pred_df <- bind_rows(synthetic_preds)


# --- Main plotting function for all species ---
plot_all_species_predictions <- function(pred_df, r_intimate_orig, r_segregated_orig) {
  
  # --- 1. Create reference plots for original layouts (Pure Black & White) ---
  
  # Dataframe for Intimate reference
  df_intimate_ref <- as.data.frame(r_intimate_orig, xy = TRUE)
  names(df_intimate_ref)[3] <- "cover" 
  p_intimate_ref <- ggplot(df_intimate_ref, aes(x = x, y = y, fill = factor(cover))) +
    geom_raster() +
    scale_fill_manual(
      values = c("white", "black"), # Pure white for deciduous, pure black for conifer
      labels = c("Deciduous", "Conifer"),
      name = "Cover Type"
    ) +
    coord_equal() + # Ensures square pixels
    # Apply theme settings. base_size is inherited from theme_set.
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_blank(), # Remove x axis text
      axis.text.y = element_blank(), # Remove y axis text
      axis.ticks = element_blank(),  # Ensure axis ticks are removed
      axis.title.x = element_blank(), # Remove x axis title
      axis.title.y = element_blank(), # Remove y axis title
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), # Adjust margins
      legend.position = "right", # Allow legend to be on this plot for collection
      legend.justification = "center",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) + # ENSURE no ticks/labels and no padding
    scale_y_continuous(expand = c(0, 0), breaks = NULL)   # ENSURE no ticks/labels and no padding
  
  # Dataframe for Segregated reference
  df_segregated_ref <- as.data.frame(r_segregated_orig, xy = TRUE)
  names(df_segregated_ref)[3] <- "cover" 
  p_segregated_ref <- ggplot(df_segregated_ref, aes(x = x, y = y, fill = factor(cover))) +
    geom_raster() +
    scale_fill_manual(
      values = c("white", "black"), # Pure white for deciduous, pure black for conifer
      labels = c("Deciduous", "Conifer"),
      name = "Cover Type" # Keep name for collection
    ) +
    coord_equal() + # Ensures square pixels
    # Apply theme settings. base_size is inherited from theme_set.
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_blank(), # Remove x axis text
      axis.text.y = element_blank(), # Remove y axis text
      axis.ticks = element_blank(), # Ensure axis ticks are removed
      axis.title.x = element_blank(), # Remove x axis title
      axis.title.y = element_blank(), # Remove y axis title
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), # Adjust margins
      legend.position = "right" # Allow legend to be on this plot for collection
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) + # ENSURE no ticks/labels and no padding
    scale_y_continuous(expand = c(0, 0), breaks = NULL)   # ENSURE no ticks/labels and no padding
  
  # Combine reference plots side-by-side with spacers for centering
  # Collect guides and place them to the right of the entire top row
  p_top_refs_combined <- wrap_plots(
    plot_spacer(),
    p_intimate_ref,
    p_segregated_ref, # Both now contribute their legends
    plot_spacer(),
    widths = c(0.5, 1, 1, 0.5), # Relative widths for (spacer, plot, plot, spacer) - Total 3 units
    nrow = 1
  ) + plot_layout(guides = "collect") & theme(legend.position = "right") # Collect all guides and place them on the right
  
  
  # --- 2. Create prediction plots for each species ---
  
  species_plots <- list()
  all_species_codes <- unique(pred_df$Species)
  
  for (i in seq_along(all_species_codes)) {
    sp <- all_species_codes[i]
    df_plot_sp <- pred_df %>%
      filter(Species == sp) %>%
      mutate(scaled_pred = (pred - min(pred)) / (max(pred) - min(pred))) # Rescale for consistent color map
    
    p_sp <- ggplot(df_plot_sp, aes(x = x, y = y, fill = scaled_pred)) +
      geom_raster() +
      # Use facet_grid for layout and age class, add "years" to age labels
      facet_grid(age_class ~ layout, labeller = labeller(age_class = function(x) paste0(x, " years"))) + 
      scale_fill_viridis_c(
        option = "magma",
        limits = c(0, 1),
        name = "Relative\nAbundance"
      ) +
      coord_equal() + # Ensures square pixels for all panels
      # Apply theme settings. base_size is inherited from theme_set.
      theme(
        # Remove all axis text and ticks for prediction plots
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(), # Ensure axis ticks are removed here too
        axis.title.x = element_blank(), # Also remove titles
        axis.title.y = element_blank(), # Also remove titles
        # Increase space around facets to prevent label overlap
        panel.spacing = unit(0.5, "cm"),
        # Adjust facet label text size by 2 points (from base_size 14 + default 0, so 14+2=16 is too large).
        # Original `strip.text` size was 10. So +2 means 12.
        strip.text = element_text(size = 12, face = "bold"), # INCREASED SIZE HERE
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm") # Reduce margins around each plot
      ) +
      scale_x_continuous(expand = c(0, 0), breaks = NULL) + # ENSURE no ticks/labels and no padding
      scale_y_continuous(expand = c(0, 0), breaks = NULL)   # ENSURE no ticks/labels and no padding
    
    # Only keep one legend for relative abundance, placed on the rightmost plot
    if (i != length(all_species_codes)) { # If not the last species (rightmost)
      p_sp <- p_sp + theme(legend.position = "none")
    }
    
    species_plots[[sp]] <- p_sp
  }
  
  # Combine species plots side by side in a single patchwork object
  # Apply species names as tags here, using tag_levels (which internally uses tag_positions logic)
  p_combined_species_with_tags <- wrap_plots(species_plots, nrow = 1) +
    plot_annotation(tag_levels = list(all_species_codes)) & # Use list() for custom levels
    theme(
      plot.tag = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0), # Style tags
      plot.tag.position = "top" # Tags will be positioned at the top of each sub-plot column
    )
  
  # --- 3. Final combination of all plots ---
  
  # Final combination: top row (centered reference) over bottom row (combined species with tags)
  final_plot_object <- p_top_refs_combined / p_combined_species_with_tags
  
  # Apply layout settings for heights to the entire combined plot.
  final_plot_object <- final_plot_object + plot_layout(
    heights = c(1, 3) # Height ratio: 1 for reference row, 3 for prediction row
  )
  
  # Apply overall main title to the *entire* figure.
  final_plot_object <- final_plot_object + plot_annotation(
    title = "Forest Layouts and Predicted Species Abundance"
  ) & theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
  
  return(final_plot_object) # Return the final plot object
}

# --- Call the new plotting function to generate the complete figure ---
final_plot <- plot_all_species_predictions(synthetic_pred_df, r_intimate, r_segregated)
print(final_plot) # Print the final plot




