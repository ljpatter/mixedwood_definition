# ---
# title: "BRT model clipped with no wetland (but including swamp)"
# author: "Leonard Patterson"
# created: "2025-09-01"
# description: 
# ---

# Install QPAD
#remotes::install_github("borealbirds/QPAD")

# === LIBRARIES ===

library(cowplot)
library(sf)
library(sp)
library(dplyr)
library(gbm)
library(dismo)
library(remotes) # Needed for installing QPAD first time
library(QPAD)
library(spdep)
library(pROC)
library(purrr)
library(ggplot2)
library(stringr)
library(readr)
library(doParallel)
library(doRNG)
library(foreach)
library(tidyr)
library(viridis)
library(knitr)
library(tibble)
library(terra)
library(landscapemetrics)
library(patchwork)
library(ggh4x)
library(scales)
library(grid)





########## 1. Run model for BTNW, clipping PC to BTNW range (including points located in wetlands)

# Load data
BTNW_data <- read.csv("Input/Tabular Data/BTNW_clipped.csv")

# Load study area
study_area <- st_read("Input/Spatial Data/Study Area/study_area.shp")

# Generate a 1 km grid over the study area
grid_1km <- st_make_grid(study_area, cellsize = 1000, square = TRUE) %>%
  st_sf(grid_id = 1:length(.), geometry = .)

# Convert point data to sf and assign CRS
coordinates(BTNW_data) <- ~x_AEP10TM + y_AEP10TM
proj4string(BTNW_data) <- CRS(st_crs(study_area)$proj4string)
joined_sf <- st_as_sf(BTNW_data)

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

# >>> NEW: precompute row indices by block for TRUE block bootstrap ----------
idx_by_block <- split(seq_len(nrow(joined_sp)), joined_sp@data$block_id)
block_ids    <- names(idx_by_block)           # only blocks that actually have data
K_blocks     <- length(block_ids)
# ---------------------------------------------------------------------------


###################### BRT POISSON W/ PARALLEL BOOTSTRAPPING

# === SETUP ===
species_list <- c("BTNW")
results_list <- list()
friedmanH_list <- list()
final_model_list <- list()

# === Start Parallel Backend ===
cores_to_use <- 6
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)
registerDoRNG(123)

# === Offsets (QPAD/BAM) ===
load_BAM_QPAD(version = 3)
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

# PDP helpers
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")

fit_gam_pdp <- function(df, n = 300) {
  # df columns: x, y (link), var, boot, Species
  df <- df %>% mutate(y_resp = exp(y))  # inverse-link for Poisson
  xr <- range(df$x, na.rm = TRUE)
  xnew <- data.frame(x = seq(xr[1], xr[2], length.out = n))
  g <- mgcv::gam(y_resp ~ s(x), data = df, family = gaussian())
  pr <- as.data.frame(predict(g, xnew, se.fit = TRUE))
  tibble(
    x   = xnew$x,
    fit = pr$fit,
    lwr = pr$fit - 1.96 * pr$se.fit,
    upr = pr$fit + 1.96 * pr$se.fit
  )
}

pdp_dir <- "Output/Tabular Data"
if (!dir.exists(pdp_dir)) dir.create(pdp_dir, recursive = TRUE)

# === BOOTSTRAPPED BRT LOOP ===
set.seed(123)
n_boot <- 25

for (sp in species_list) {
  cat("\n\n=== Species:", sp, "===\n")
  
  joined_sp@data$offset <- calculate_offsets(joined_sp@data, sp)
  joined_sp@data$year <- as.factor(joined_sp@data$year)
  
  boot_results <- foreach(
    i = 1:n_boot,
    .packages = c("dismo", "gbm", "pROC", "dplyr", "tidyr", "tibble", "purrr")
  ) %dopar% {
    
    # ---- TRUE BLOCK BOOTSTRAP (WITH REPLACEMENT) --------------------------
    boot_blocks <- sample(block_ids, size = K_blocks, replace = TRUE)
    boot_idx    <- unlist(idx_by_block[boot_blocks], use.names = FALSE)  # keep duplicates
    sp_subset   <- joined_sp[boot_idx, ]
    block_subset <- cbind(as.data.frame(sp_subset@coords), sp_subset@data)
    block_subset$year <- as.factor(block_subset$year)
    
    sampled_data <- data.frame(
      response  = block_subset[[sp]],
      offset    = block_subset$offset,
      year      = block_subset$year,
      x_AEP10TM = block_subset$x_AEP10TM,
      y_AEP10TM = block_subset$y_AEP10TM,
      block_subset %>%
        dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
        dplyr::select(where(is.numeric))
    )
    
    if (sum(sampled_data$response, na.rm = TRUE) == 0) return(NULL)
    stopifnot(is.factor(sampled_data$year))
    
    # ---- BLOCKED CV FOLDS (NO LEAKAGE) ------------------------------------
    uniq_blocks <- sort(unique(as.character(sp_subset@data$block_id)))
    k_folds     <- min(10L, length(uniq_blocks))
    fold_map    <- setNames(((seq_along(uniq_blocks)-1L) %% k_folds) + 1L, uniq_blocks)
    folds       <- as.integer(fold_map[as.character(sp_subset@data$block_id)])
    
    # ---- FIT ---------------------------------------------------------------
    brt_model <- tryCatch({
      gbm.step(
        data             = sampled_data,
        gbm.y            = 1,
        gbm.x            = 3:ncol(sampled_data),
        family           = "poisson",
        tree.complexity  = 3,
        learning.rate    = 0.005,
        bag.fraction     = 0.5,
        offset           = sampled_data$offset,
        fold.vector      = folds,     # BLOCKED folds
        n.folds          = k_folds,   # match actual number of folds
        keep.fold.models = FALSE,
        keep.fold.fit    = FALSE,
        verbose          = FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(brt_model)) return(NULL)
    
    imp_df <- data.frame(
      Feature = summary(brt_model)$var,
      Gain    = summary(brt_model)$rel.inf
    )
    
    model_data <- sampled_data %>%
      dplyr::select(all_of(brt_model$gbm.call$predictor.names)) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ tidyr::replace_na(., 0)))
    
    interaction_df <- .gbm.interactions_full(brt_model, model_data)
    interaction_df$Boot <- i
    
    # PDPs for ALL predictors
    var_names <- brt_model$gbm.call$predictor.names
    pdp_df <- purrr::map_dfr(seq_along(var_names), function(k) {
      g <- try(gbm::plot.gbm(brt_model, i.var = k, return.grid = TRUE), silent = TRUE)
      if (inherits(g, "try-error") || !"y" %in% names(g)) return(NULL)
      colnames(g)[1] <- "x"
      tibble(
        x    = as.numeric(g$x),
        y    = as.numeric(g$y),  # link scale (log for Poisson)
        var  = var_names[k],
        boot = i
      )
    })
    
    list(imp = imp_df, interaction = interaction_df, pdp = pdp_df)
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
  
  # Knight-style PDP summarization for this species
  pdp_all <- lapply(boot_results, function(x) x$pdp)
  pdp_all <- pdp_all[!sapply(pdp_all, is.null)]
  if (length(pdp_all) > 0) {
    pdp_raw_sp <- bind_rows(pdp_all)
    if (nrow(pdp_raw_sp) > 0) {
      pdp_raw_sp <- pdp_raw_sp %>%
        filter(!var %in% nuisance_vars) %>%
        mutate(Species = sp)
      saveRDS(pdp_raw_sp, file.path(pdp_dir, paste0("BTNW_pdp_bootstrap_raw_", sp, ".rds")))
      
      pdp_gam_sp <- pdp_raw_sp %>%
        group_by(Species, var) %>%
        group_modify(~ fit_gam_pdp(.x)) %>%
        ungroup()
      
      saveRDS(pdp_gam_sp, file.path(pdp_dir, paste0("BTNW_pdp_gam_summary_", sp, ".rds")))
    }
  }
  
  # === Final model on full data (with BLOCKED CV) ===
  joined_sp@data$year <- as.factor(joined_sp@data$year)
  
  full_data <- cbind(as.data.frame(joined_sp@coords), joined_sp@data)
  model_data <- data.frame(
    response  = full_data[[sp]],
    offset    = full_data$offset,
    year      = full_data$year,
    x_AEP10TM = full_data$x_AEP10TM,
    y_AEP10TM = full_data$y_AEP10TM,
    full_data %>%
      dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
      dplyr::select(where(is.numeric))
  )
  stopifnot(is.factor(model_data$year))
  
  # ---- BLOCKED CV ON FULL DATA (DETERMINISTIC) ----
  # For full-data CV
  uniq_blocks_full <- sort(unique(as.character(joined_sp@data$block_id)))
  k_folds_full     <- min(10L, length(uniq_blocks_full))
  fold_map_full    <- setNames(((seq_along(uniq_blocks_full)-1L) %% k_folds_full) + 1L, uniq_blocks_full)
  folds_full       <- as.integer(fold_map_full[as.character(joined_sp@data$block_id)])
  
  
  final_model <- gbm.step(
    data             = model_data,
    gbm.y            = 1,
    gbm.x            = 3:ncol(model_data),   # <<< DO NOT drop last column
    family           = "poisson",
    tree.complexity  = 3,
    learning.rate    = 0.005,
    bag.fraction     = 0.5,
    offset           = model_data$offset,
    fold.vector      = folds_full,           # <<< BLOCKED folds on full data
    n.folds          = k_folds_full,
    keep.fold.models = FALSE,
    keep.fold.fit    = FALSE,
    verbose          = FALSE
  )
  
  final_model_list[[sp]] <- final_model
}

# === Stop Parallel Cluster ===
stopCluster(cl)

# === Save Results ===
imp_df <- bind_rows(results_list)
friedman_df <- bind_rows(friedmanH_list)

imp_table <- imp_df %>%
  dplyr::mutate(var_type = dplyr::case_when(
    stringr::str_starts(Feature, "prop_con") ~ "prop_con",
    stringr::str_starts(Feature, "clumpy")   ~ "clumpy",
    stringr::str_starts(Feature, "age_mn")   ~ "age_mn",
    TRUE ~ "other"
  ))
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)

saveRDS(results_list, "Output/Tabular Data/BTNW_results_list_final_parallel_poisson_5m.rds")
saveRDS(imp_df, "Output/Tabular Data/BTNW_imp_df_final_parallel_poisson_5m.rds")
saveRDS(friedman_df, "Output/Tabular Data/BTNW_friedman_interactions_bootstrapped_parallel_final_poisson_5m.rds")
saveRDS(final_model_list, "Output/Tabular Data/BTNW_final_model_list_parallel_poisson_5m.rds")
saveRDS(imp_table, file.path(pdp_dir, "BTNW_imp_table_5m.rds"))

cat("\n✅ Saved all model outputs for Poisson models with year as FACTOR!\n")






# ============================
# GAM-PDPs at Scale of Effect
# ============================


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(mgcv)
  library(readr)
})

# --- Locations ---
pdp_dir           <- "Output/Tabular Data"
imp_table_path    <- file.path("Output/Tabular Data", "BTNW_imp_table_5m.rds")
imp_fallback_path <- file.path("Output/Tabular Data", "BTNW_imp_df_final_parallel_poisson_5m.rds")

# --- Options ---
n_points      <- 300    # x-grid resolution for smooth curves
trim_prop     <- 0.01   # trim 1% tails to avoid edge wiggles
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")

# --- Outputs ---
out_tab_dir <- "Output/Tables/BTNW_GAM_PDP_5m"
out_fig_dir <- "Output/Figures/BTNW_GAM_PDP_5m"
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----
nice_lab <- function(v) {
  v %>%
    str_replace("_NTEMS$", " (NTEMS)") %>%
    str_replace("_LULC$",  " (LULC)")  %>%
    str_replace("^prop_con", "Proportion conifer") %>%
    str_replace("^clumpy",   "Clumpiness") %>%
    str_replace("^age_mn",   "Forest age")
}

var_type_from_name <- function(v) case_when(
  str_starts(v, "prop_con") ~ "prop_con",
  str_starts(v, "clumpy")   ~ "clumpy",
  str_starts(v, "age_mn")   ~ "age_mn",
  TRUE ~ "other"
)

slugify <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace("^_|_$", "")
}

`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a

# ---------- Load or rebuild imp_table ----------
if (exists("imp_table") && is.data.frame(imp_table)) {
  message("Using in-memory imp_table")
} else if (file.exists(imp_table_path)) {
  imp_table <- readRDS(imp_table_path)
  message("Loaded imp_table from ", imp_table_path)
} else if (file.exists(imp_fallback_path)) {
  message("imp_table.rds not found; deriving from ", imp_fallback_path)
  imp_table <- readRDS(imp_fallback_path) %>%
    mutate(var_type = var_type_from_name(Feature)) %>%
    filter(var_type %in% c("prop_con","clumpy","age_mn"))
} else {
  stop("Couldn't find an importance table. Expected one of:\n",
       " - ", imp_table_path, "\n",
       " - ", imp_fallback_path, "\n",
       "Or define imp_table in memory before running.")
}

imp_table <- imp_table %>%
  mutate(var_type = factor(var_type, levels = c("prop_con","clumpy","age_mn")))

# -------------------------------------------
# 1) Choose Scale-of-Effect (SoE) per species
#    = argmax MeanGain within each var_type
# -------------------------------------------
soe_tbl <- imp_table %>%
  filter(var_type %in% c("prop_con","clumpy","age_mn")) %>%
  group_by(Species, var_type) %>%
  arrange(desc(MeanGain), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(Species, var_type, Feature_soe = Feature)

print(soe_tbl, n = Inf)

# Species list (use what’s actually present in SoE table)
species_list <- unique(soe_tbl$Species) %||% c("BTNW")

# ---- Fit one population-level GAM with random intercept per bootstrap ----
fit_gam_pdp_one <- function(df, n = n_points, trim = trim_prop) {
  # df columns: x, y (link), var, boot, Species
  df <- df %>%
    mutate(
      boot   = factor(boot),
      y_resp = exp(y)               # Poisson inverse link
    ) %>%
    filter(is.finite(x), is.finite(y_resp))
  
  if (nrow(df) < 10) return(NULL)
  
  # Trim extremes for stability
  q   <- quantile(df$x, probs = c(trim, 1 - trim), na.rm = TRUE)
  xlo <- q[1]; xhi <- q[2]
  if (!is.finite(xlo) || !is.finite(xhi) || xlo >= xhi) return(NULL)
  
  # Include a dummy boot factor with correct levels for predict()
  xnew <- tibble(
    x    = seq(xlo, xhi, length.out = n),
    boot = factor(df$boot[1], levels = levels(df$boot))
  )
  
  g <- mgcv::gam(
    y_resp ~ s(x) + s(boot, bs = "re"),
    data   = df,
    family = gaussian(),
    method = "REML",
    select = TRUE
  )
  
  pr <- as.data.frame(predict(g, newdata = xnew, se.fit = TRUE, exclude = "s(boot)"))
  
  tibble(
    x   = xnew$x,
    fit = pr$fit,
    lwr = pr$fit - 1.96 * pr$se.fit,
    upr = pr$fit + 1.96 * pr$se.fit,
    edf = summary(g)$s.table[rownames(summary(g)$s.table) == "s(x)", "edf"] %||% NA_real_,
    r2  = summary(g)$r.sq %||% NA_real_
  )
}

# -------------------------------------------
# 2) Build GAM-PDPs only for those SoE features
#    and SAVE ONE FILE PER PREDICTOR
# -------------------------------------------
all_plots <- list()

for (sp in species_list) {
  pdp_file <- file.path(pdp_dir, paste0("BTNW_pdp_bootstrap_raw_", sp, ".rds"))
  if (!file.exists(pdp_file)) {
    message("Missing PDP file for ", sp, " — skipping.")
    next
  }
  
  pdp_raw <- readRDS(pdp_file) %>%
    filter(!var %in% nuisance_vars) %>%
    drop_na(x, y)
  
  # Features selected as SoE for this species
  soe_vars_sp <- soe_tbl %>% filter(Species == sp) %>% pull(Feature_soe)
  if (length(soe_vars_sp) == 0) {
    message("No SoE features for ", sp, " — skipping.")
    next
  }
  
  pdp_soe <- pdp_raw %>% filter(var %in% soe_vars_sp)
  if (nrow(pdp_soe) == 0) {
    message("No PDP rows for SoE features in ", sp, " — skipping.")
    next
  }
  
  # Diagnostics: counts per SoE variable
  message("Counts for ", sp, " (SoE):")
  print(pdp_soe %>% count(var), n = Inf)
  
  # Fit one GAM per selected feature and SAVE EACH SEPARATELY
  unique_vars <- unique(pdp_soe$var)
  for (v in unique_vars) {
    df_v <- pdp_soe %>% filter(var == v)
    sm   <- try(fit_gam_pdp_one(df_v), silent = TRUE)
    if (inherits(sm, "try-error") || is.null(sm)) {
      message("Skipping ", sp, " / ", v, " (insufficient data or fit error).")
      next
    }
    
    # Tidy + save table per variable
    out_tab <- sm %>%
      mutate(Species = sp, var = v, var_clean = nice_lab(v)) %>%
      relocate(Species, var, var_clean)
    csv_name <- file.path(out_tab_dir, paste0("GAM_PDP_SoE_", sp, "_", slugify(v), ".csv"))
    write_csv(out_tab, csv_name)
    
    # Subtitle shows chosen SoE for this var_type
    vt <- var_type_from_name(v)
    subtitle_txt <- paste0("Selected SoE (", vt, "): ", v)
    
    p <- ggplot(out_tab, aes(x = x, y = fit)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey30", colour = NA, alpha = 0.25) +
      geom_line() +
      labs(
        title    = paste0("GAM-smoothed PDP at Scale of Effect (", sp, ")"),
        subtitle = subtitle_txt,
        x = "Predictor value",
        y = "Expected count"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank())
    
    png_name <- file.path(out_fig_dir, paste0("GAM_PDP_SoE_", sp, "_", slugify(v), ".png"))
    ggsave(png_name, p, width = 7.5, height = 5.0, dpi = 300)
    
    all_plots[[paste(sp, v, sep = "_")]] <- p
  }
}

invisible(all_plots)

print(all_plots)




























rm(list = ls())




########## 2. Run model for BTNW, clipping PC to BTNW range AND removing PC in wetlands (except swamps)

# Load data
BTNW_data <- read.csv("Input/Tabular Data/BTNW_clipped_swamp.csv")

# Load study area
study_area <- st_read("Input/Spatial Data/Study Area/study_area.shp")

# Generate a 1 km grid over the study area
grid_1km <- st_make_grid(study_area, cellsize = 1000, square = TRUE) %>%
  st_sf(grid_id = 1:length(.), geometry = .)

# Convert point data to sf and assign CRS
coordinates(BTNW_data) <- ~x_AEP10TM + y_AEP10TM
proj4string(BTNW_data) <- CRS(st_crs(study_area)$proj4string)
joined_sf <- st_as_sf(BTNW_data)

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

# >>> NEW: precompute row indices by block for TRUE block bootstrap ----------
idx_by_block <- split(seq_len(nrow(joined_sp)), joined_sp@data$block_id)
block_ids    <- names(idx_by_block)           # only blocks that actually have data
K_blocks     <- length(block_ids)
# ---------------------------------------------------------------------------


###################### BRT POISSON W/ PARALLEL BOOTSTRAPPING

# === SETUP ===
species_list <- c("BTNW")
results_list <- list()
friedmanH_list <- list()
final_model_list <- list()

# === Start Parallel Backend ===
cores_to_use <- 6
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)
registerDoRNG(123)

# === Offsets (QPAD/BAM) ===
load_BAM_QPAD(version = 3)
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

# PDP helpers
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")

fit_gam_pdp <- function(df, n = 300) {
  # df columns: x, y (link), var, boot, Species
  df <- df %>% mutate(y_resp = exp(y))  # inverse-link for Poisson
  xr <- range(df$x, na.rm = TRUE)
  xnew <- data.frame(x = seq(xr[1], xr[2], length.out = n))
  g <- mgcv::gam(y_resp ~ s(x), data = df, family = gaussian())
  pr <- as.data.frame(predict(g, xnew, se.fit = TRUE))
  tibble(
    x   = xnew$x,
    fit = pr$fit,
    lwr = pr$fit - 1.96 * pr$se.fit,
    upr = pr$fit + 1.96 * pr$se.fit
  )
}

pdp_dir <- "Output/Tabular Data"
if (!dir.exists(pdp_dir)) dir.create(pdp_dir, recursive = TRUE)

# === BOOTSTRAPPED BRT LOOP ===
set.seed(123)
n_boot <- 50

for (sp in species_list) {
  cat("\n\n=== Species:", sp, "===\n")
  
  joined_sp@data$offset <- calculate_offsets(joined_sp@data, sp)
  joined_sp@data$year <- as.factor(joined_sp@data$year)
  
  boot_results <- foreach(
    i = 1:n_boot,
    .packages = c("dismo", "gbm", "pROC", "dplyr", "tidyr", "tibble", "purrr")
  ) %dopar% {
    
    # ---- TRUE BLOCK BOOTSTRAP (WITH REPLACEMENT) --------------------------
    boot_blocks <- sample(block_ids, size = K_blocks, replace = TRUE)
    boot_idx    <- unlist(idx_by_block[boot_blocks], use.names = FALSE)  # keep duplicates
    sp_subset   <- joined_sp[boot_idx, ]
    block_subset <- cbind(as.data.frame(sp_subset@coords), sp_subset@data)
    block_subset$year <- as.factor(block_subset$year)
    
    sampled_data <- data.frame(
      response  = block_subset[[sp]],
      offset    = block_subset$offset,
      year      = block_subset$year,
      x_AEP10TM = block_subset$x_AEP10TM,
      y_AEP10TM = block_subset$y_AEP10TM,
      block_subset %>%
        dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
        dplyr::select(where(is.numeric))
    )
    
    if (sum(sampled_data$response, na.rm = TRUE) == 0) return(NULL)
    stopifnot(is.factor(sampled_data$year))
    
    # ---- BLOCKED CV FOLDS (NO LEAKAGE) ------------------------------------
    uniq_blocks <- sort(unique(as.character(sp_subset@data$block_id)))
    k_folds     <- min(10L, length(uniq_blocks))
    fold_map    <- setNames(((seq_along(uniq_blocks)-1L) %% k_folds) + 1L, uniq_blocks)
    folds       <- as.integer(fold_map[as.character(sp_subset@data$block_id)])
    
    # ---- FIT ---------------------------------------------------------------
    brt_model <- tryCatch({
      gbm.step(
        data             = sampled_data,
        gbm.y            = 1,
        gbm.x            = 3:ncol(sampled_data),
        family           = "poisson",
        tree.complexity  = 3,
        learning.rate    = 0.005,
        bag.fraction     = 0.5,
        offset           = sampled_data$offset,
        fold.vector      = folds,     # BLOCKED folds
        n.folds          = k_folds,   # match actual number of folds
        keep.fold.models = FALSE,
        keep.fold.fit    = FALSE,
        verbose          = FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(brt_model)) return(NULL)
    
    imp_df <- data.frame(
      Feature = summary(brt_model)$var,
      Gain    = summary(brt_model)$rel.inf
    )
    
    model_data <- sampled_data %>%
      dplyr::select(all_of(brt_model$gbm.call$predictor.names)) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ tidyr::replace_na(., 0)))
    
    interaction_df <- .gbm.interactions_full(brt_model, model_data)
    interaction_df$Boot <- i
    
    # PDPs for ALL predictors
    var_names <- brt_model$gbm.call$predictor.names
    pdp_df <- purrr::map_dfr(seq_along(var_names), function(k) {
      g <- try(gbm::plot.gbm(brt_model, i.var = k, return.grid = TRUE), silent = TRUE)
      if (inherits(g, "try-error") || !"y" %in% names(g)) return(NULL)
      colnames(g)[1] <- "x"
      tibble(
        x    = as.numeric(g$x),
        y    = as.numeric(g$y),  # link scale (log for Poisson)
        var  = var_names[k],
        boot = i
      )
    })
    
    list(imp = imp_df, interaction = interaction_df, pdp = pdp_df)
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
  
  # Knight-style PDP summarization for this species
  pdp_all <- lapply(boot_results, function(x) x$pdp)
  pdp_all <- pdp_all[!sapply(pdp_all, is.null)]
  if (length(pdp_all) > 0) {
    pdp_raw_sp <- bind_rows(pdp_all)
    if (nrow(pdp_raw_sp) > 0) {
      pdp_raw_sp <- pdp_raw_sp %>%
        filter(!var %in% nuisance_vars) %>%
        mutate(Species = sp)
      saveRDS(pdp_raw_sp, file.path(pdp_dir, paste0("BTNW_pdp_bootstrap_raw_", sp, ".rds")))
      
      pdp_gam_sp <- pdp_raw_sp %>%
        group_by(Species, var) %>%
        group_modify(~ fit_gam_pdp(.x)) %>%
        ungroup()
      
      saveRDS(pdp_gam_sp, file.path(pdp_dir, paste0("BTNW_pdp_gam_summary_", sp, ".rds")))
    }
  }
  
  # === Final model on full data (with BLOCKED CV) ===
  joined_sp@data$year <- as.factor(joined_sp@data$year)
  
  full_data <- cbind(as.data.frame(joined_sp@coords), joined_sp@data)
  model_data <- data.frame(
    response  = full_data[[sp]],
    offset    = full_data$offset,
    year      = full_data$year,
    x_AEP10TM = full_data$x_AEP10TM,
    y_AEP10TM = full_data$y_AEP10TM,
    full_data %>%
      dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
      dplyr::select(where(is.numeric))
  )
  stopifnot(is.factor(model_data$year))
  
  # ---- BLOCKED CV ON FULL DATA (DETERMINISTIC) ----
  # For full-data CV
  uniq_blocks_full <- sort(unique(as.character(joined_sp@data$block_id)))
  k_folds_full     <- min(10L, length(uniq_blocks_full))
  fold_map_full    <- setNames(((seq_along(uniq_blocks_full)-1L) %% k_folds_full) + 1L, uniq_blocks_full)
  folds_full       <- as.integer(fold_map_full[as.character(joined_sp@data$block_id)])
  
  
  final_model <- gbm.step(
    data             = model_data,
    gbm.y            = 1,
    gbm.x            = 3:ncol(model_data),   # <<< DO NOT drop last column
    family           = "poisson",
    tree.complexity  = 3,
    learning.rate    = 0.005,
    bag.fraction     = 0.5,
    offset           = model_data$offset,
    fold.vector      = folds_full,           # <<< BLOCKED folds on full data
    n.folds          = k_folds_full,
    keep.fold.models = FALSE,
    keep.fold.fit    = FALSE,
    verbose          = FALSE
  )
  
  final_model_list[[sp]] <- final_model
}

# === Stop Parallel Cluster ===
stopCluster(cl)

# === Save Results ===
imp_df <- bind_rows(results_list)
friedman_df <- bind_rows(friedmanH_list)

imp_table <- imp_df %>%
  dplyr::mutate(var_type = dplyr::case_when(
    stringr::str_starts(Feature, "prop_con") ~ "prop_con",
    stringr::str_starts(Feature, "clumpy")   ~ "clumpy",
    stringr::str_starts(Feature, "age_mn")   ~ "age_mn",
    TRUE ~ "other"
  ))
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)

saveRDS(results_list, "Output/Tabular Data/BTNW_no_wetlands_results_list_final_parallel_poisson_5m.rds")
saveRDS(imp_df, "Output/Tabular Data/BTNW_no_wetlands_imp_df_final_parallel_poisson_5m.rds")
saveRDS(friedman_df, "Output/Tabular Data/BTNW_no_wetlands_friedman_interactions_bootstrapped_parallel_final_poisson_5m.rds")
saveRDS(final_model_list, "Output/Tabular Data/BTNW_no_wetlands_final_model_list_parallel_poisson_5m.rds")
saveRDS(imp_table, file.path(pdp_dir, "BTNW_no_wetlands_imp_table_5m.rds"))

cat("\n✅ Saved all model outputs for Poisson models with year as FACTOR!\n")








### ============================
### BTNW GAM-PDPs at Scale of Effect (no wetlands except swamp)
### ============================


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(mgcv)
  library(readr)
})

# --- Locations ---
pdp_dir           <- "Output/Tabular Data"
imp_table_path    <- file.path("Output/Tabular Data", "BTNW_no_wetlands_imp_table_5m.rds")
imp_fallback_path <- file.path("Output/Tabular Data", "BTNW_no_wetlands_imp_df_final_parallel_poisson_5m.rds")

# --- Options ---
n_points      <- 300    # x-grid resolution for smooth curves
trim_prop     <- 0.01   # trim 1% tails to avoid edge wiggles
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")

# --- Outputs ---
out_tab_dir <- "Output/Tables/BTNW_no_wetlands_GAM_PDP_5m"
out_fig_dir <- "Output/Figures/BTNW_no_wetlands_GAM_PDP_5m"
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----
nice_lab <- function(v) {
  v %>%
    str_replace("_NTEMS$", " (NTEMS)") %>%
    str_replace("_LULC$",  " (LULC)")  %>%
    str_replace("^prop_con", "Proportion conifer") %>%
    str_replace("^clumpy",   "Clumpiness") %>%
    str_replace("^age_mn",   "Forest age")
}

var_type_from_name <- function(v) case_when(
  str_starts(v, "prop_con") ~ "prop_con",
  str_starts(v, "clumpy")   ~ "clumpy",
  str_starts(v, "age_mn")   ~ "age_mn",
  TRUE ~ "other"
)

slugify <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace("^_|_$", "")
}

`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a

# ---------- Load or rebuild imp_table ----------
if (exists("imp_table") && is.data.frame(imp_table)) {
  message("Using in-memory imp_table")
} else if (file.exists(imp_table_path)) {
  imp_table <- readRDS(imp_table_path)
  message("Loaded imp_table from ", imp_table_path)
} else if (file.exists(imp_fallback_path)) {
  message("imp_table.rds not found; deriving from ", imp_fallback_path)
  imp_table <- readRDS(imp_fallback_path) %>%
    mutate(var_type = var_type_from_name(Feature)) %>%
    filter(var_type %in% c("prop_con","clumpy","age_mn"))
} else {
  stop("Couldn't find an importance table. Expected one of:\n",
       " - ", imp_table_path, "\n",
       " - ", imp_fallback_path, "\n",
       "Or define imp_table in memory before running.")
}

imp_table <- imp_table %>%
  mutate(var_type = factor(var_type, levels = c("prop_con","clumpy","age_mn")))

# -------------------------------------------
# 1) Choose Scale-of-Effect (SoE) per species
#    = argmax MeanGain within each var_type
# -------------------------------------------
soe_tbl <- imp_table %>%
  filter(var_type %in% c("prop_con","clumpy","age_mn")) %>%
  group_by(Species, var_type) %>%
  arrange(desc(MeanGain), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(Species, var_type, Feature_soe = Feature)

print(soe_tbl, n = Inf)

# Species list (use what’s actually present in SoE table)
species_list <- unique(soe_tbl$Species) %||% c("BTNW")

# ---- Fit one population-level GAM with random intercept per bootstrap ----
fit_gam_pdp_one <- function(df, n = n_points, trim = trim_prop) {
  # df columns: x, y (link), var, boot, Species
  df <- df %>%
    mutate(
      boot   = factor(boot),
      y_resp = exp(y)               # Poisson inverse link
    ) %>%
    filter(is.finite(x), is.finite(y_resp))
  
  if (nrow(df) < 10) return(NULL)
  
  # Trim extremes for stability
  q   <- quantile(df$x, probs = c(trim, 1 - trim), na.rm = TRUE)
  xlo <- q[1]; xhi <- q[2]
  if (!is.finite(xlo) || !is.finite(xhi) || xlo >= xhi) return(NULL)
  
  # Include a dummy boot factor with correct levels for predict()
  xnew <- tibble(
    x    = seq(xlo, xhi, length.out = n),
    boot = factor(df$boot[1], levels = levels(df$boot))
  )
  
  g <- mgcv::gam(
    y_resp ~ s(x) + s(boot, bs = "re"),
    data   = df,
    family = gaussian(),
    method = "REML",
    select = TRUE
  )
  
  pr <- as.data.frame(predict(g, newdata = xnew, se.fit = TRUE, exclude = "s(boot)"))
  
  tibble(
    x   = xnew$x,
    fit = pr$fit,
    lwr = pr$fit - 1.96 * pr$se.fit,
    upr = pr$fit + 1.96 * pr$se.fit,
    edf = summary(g)$s.table[rownames(summary(g)$s.table) == "s(x)", "edf"] %||% NA_real_,
    r2  = summary(g)$r.sq %||% NA_real_
  )
}

# -------------------------------------------
# 2) Build GAM-PDPs only for those SoE features
#    and SAVE ONE FILE PER PREDICTOR
# -------------------------------------------
all_plots <- list()

for (sp in species_list) {
  pdp_file <- file.path(pdp_dir, paste0("BTNW_pdp_bootstrap_raw_", sp, ".rds"))
  if (!file.exists(pdp_file)) {
    message("Missing PDP file for ", sp, " — skipping.")
    next
  }
  
  pdp_raw <- readRDS(pdp_file) %>%
    filter(!var %in% nuisance_vars) %>%
    drop_na(x, y)
  
  # Features selected as SoE for this species
  soe_vars_sp <- soe_tbl %>% filter(Species == sp) %>% pull(Feature_soe)
  if (length(soe_vars_sp) == 0) {
    message("No SoE features for ", sp, " — skipping.")
    next
  }
  
  pdp_soe <- pdp_raw %>% filter(var %in% soe_vars_sp)
  if (nrow(pdp_soe) == 0) {
    message("No PDP rows for SoE features in ", sp, " — skipping.")
    next
  }
  
  # Diagnostics: counts per SoE variable
  message("Counts for ", sp, " (SoE):")
  print(pdp_soe %>% count(var), n = Inf)
  
  # Fit one GAM per selected feature and SAVE EACH SEPARATELY
  unique_vars <- unique(pdp_soe$var)
  for (v in unique_vars) {
    df_v <- pdp_soe %>% filter(var == v)
    sm   <- try(fit_gam_pdp_one(df_v), silent = TRUE)
    if (inherits(sm, "try-error") || is.null(sm)) {
      message("Skipping ", sp, " / ", v, " (insufficient data or fit error).")
      next
    }
    
    # Tidy + save table per variable
    out_tab <- sm %>%
      mutate(Species = sp, var = v, var_clean = nice_lab(v)) %>%
      relocate(Species, var, var_clean)
    csv_name <- file.path(out_tab_dir, paste0("GAM_PDP_SoE_", sp, "_", slugify(v), ".csv"))
    write_csv(out_tab, csv_name)
    
    # Subtitle shows chosen SoE for this var_type
    vt <- var_type_from_name(v)
    subtitle_txt <- paste0("Selected SoE (", vt, "): ", v)
    
    p <- ggplot(out_tab, aes(x = x, y = fit)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey30", colour = NA, alpha = 0.25) +
      geom_line() +
      labs(
        title    = paste0("GAM-smoothed PDP at Scale of Effect (", sp, ")"),
        subtitle = subtitle_txt,
        x = "Predictor value",
        y = "Expected count"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank())
    
    png_name <- file.path(out_fig_dir, paste0("GAM_PDP_SoE_", sp, "_", slugify(v), ".png"))
    ggsave(png_name, p, width = 7.5, height = 5.0, dpi = 300)
    
    all_plots[[paste(sp, v, sep = "_")]] <- p
  }
}

# Print
print(all_plots)








### BTNW VI PLOT - No wetlands (but including swamp)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggtext)
  library(cowplot)
})

# --- Load variable importance results ---
imp_df <- readRDS("Output/Tabular Data/BTNW_no_wetlands_imp_df_final_parallel_poisson_5m.rds")

# --- Config ---
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")
fill_map <- c("prop_con" = "forestgreen",
              "clumpy"   = "cornflowerblue",
              "age_mn"   = "lightcoral")
var_order <- c("prop_con","clumpy","age_mn")

# --- Clean & classify (BTNW only) ---
imp_df_clean <- imp_df %>%
  filter(Species == "BTNW") %>%
  filter(!Feature %in% nuisance_vars) %>%
  # drop LULC age terms if undesired; comment out next line to keep them
  filter(!grepl("^age_mn_\\d+_5$", Feature)) %>%
  mutate(
    var_type = case_when(
      str_detect(Feature, "^prop_con") ~ "prop_con",
      str_detect(Feature, "^clumpy")   ~ "clumpy",
      str_detect(Feature, "^age_mn")   ~ "age_mn",
      TRUE                             ~ NA_character_
    ),
    var_type = factor(var_type, levels = var_order)
  ) %>%
  filter(!is.na(var_type))

# --- Auto-scale axis to data (this is the plotted "y", but becomes x after coord_flip) ---
y_max_data <- max(imp_df_clean$MeanGain, na.rm = TRUE)
y_max <- max(1, ceiling(y_max_data * 1.08))  # ~8% headroom, at least 1
y_breaks <- pretty(c(0, y_max), n = 6)
y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_max]

# --- One-species plot (grouped vertically by var_type; left group labels hidden) ---
make_vi_plot <- function(df_sp) {
  df_sp <- df_sp %>%
    group_by(var_type) %>%
    arrange(desc(MeanGain), .by_group = TRUE) %>%
    mutate(is_top = row_number() == 1L) %>%
    ungroup() %>%
    mutate(
      Feature       = fct_rev(fct_inorder(Feature)),  # order for coord_flip
      Feature_label = ifelse(is_top,
                             paste0("**", as.character(Feature), "**"),
                             as.character(Feature))
    )
  
  lab_map <- setNames(df_sp$Feature_label, df_sp$Feature)
  
  ggplot(df_sp, aes(x = Feature, y = MeanGain, fill = var_type)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_manual(
      values = fill_map, drop = FALSE,
      name = "Predictor type",
      labels = c("prop_con" = "Proportion conifer",
                 "clumpy"   = "Clumpiness",
                 "age_mn"   = "Forest age")
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.06))
    ) +
    scale_x_discrete(labels = lab_map) +
    coord_flip() +
    facet_grid(var_type ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = "Predictor variables",
      y = "Mean relative influence"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      axis.text          = element_text(color = "black"),
      axis.text.y        = ggtext::element_markdown(hjust = 1, size = 11),
      
      strip.text.y       = element_blank(),
      strip.background   = element_blank(),
      
      legend.position    = "right",
      legend.title       = element_text(face = "bold")
    )
}

final_plot <- make_vi_plot(imp_df_clean)

# --- Save ---
out_fig <- "Output/Figures/BTNW_NO_WETLAND_5m_VI.png"
dir.create(dirname(out_fig), showWarnings = FALSE, recursive = TRUE)
ggsave(out_fig, final_plot, width = 9, height = 10, dpi = 300)
message("✅ Saved: ", out_fig)

print(out_fig)

















































rm(list = ls())



########## 3. Run model for BBWA, clipping PC to BBWA range AND removing PC in wetlands (except swamps)

# Load data
BBWA_data <- read.csv("Input/Tabular Data/BBWA_clipped_swamp.csv")

# Load study area
study_area <- st_read("Input/Spatial Data/Study Area/study_area.shp")

# Generate a 1 km grid over the study area
grid_1km <- st_make_grid(study_area, cellsize = 1000, square = TRUE) %>%
  st_sf(grid_id = 1:length(.), geometry = .)

# Convert point data to sf and assign CRS
coordinates(BBWA_data) <- ~x_AEP10TM + y_AEP10TM
proj4string(BBWA_data) <- CRS(st_crs(study_area)$proj4string)
joined_sf <- st_as_sf(BBWA_data)

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

# >>> NEW: precompute row indices by block for TRUE block bootstrap ----------
idx_by_block <- split(seq_len(nrow(joined_sp)), joined_sp@data$block_id)
block_ids    <- names(idx_by_block)           # only blocks that actually have data
K_blocks     <- length(block_ids)
# ---------------------------------------------------------------------------


###################### BRT POISSON W/ PARALLEL BOOTSTRAPPING

# === SETUP ===
species_list <- c("BBWA")
results_list <- list()
friedmanH_list <- list()
final_model_list <- list()

# === Start Parallel Backend ===
cores_to_use <- 6
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)
registerDoRNG(123)

# === Offsets (QPAD/BAM) ===
load_BAM_QPAD(version = 3)
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

# PDP helpers
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")

fit_gam_pdp <- function(df, n = 300) {
  # df columns: x, y (link), var, boot, Species
  df <- df %>% mutate(y_resp = exp(y))  # inverse-link for Poisson
  xr <- range(df$x, na.rm = TRUE)
  xnew <- data.frame(x = seq(xr[1], xr[2], length.out = n))
  g <- mgcv::gam(y_resp ~ s(x), data = df, family = gaussian())
  pr <- as.data.frame(predict(g, xnew, se.fit = TRUE))
  tibble(
    x   = xnew$x,
    fit = pr$fit,
    lwr = pr$fit - 1.96 * pr$se.fit,
    upr = pr$fit + 1.96 * pr$se.fit
  )
}

pdp_dir <- "Output/Tabular Data"
if (!dir.exists(pdp_dir)) dir.create(pdp_dir, recursive = TRUE)

# === BOOTSTRAPPED BRT LOOP ===
set.seed(123)
n_boot <- 50

for (sp in species_list) {
  cat("\n\n=== Species:", sp, "===\n")
  
  joined_sp@data$offset <- calculate_offsets(joined_sp@data, sp)
  joined_sp@data$year <- as.factor(joined_sp@data$year)
  
  boot_results <- foreach(
    i = 1:n_boot,
    .packages = c("dismo", "gbm", "pROC", "dplyr", "tidyr", "tibble", "purrr")
  ) %dopar% {
    
    # ---- TRUE BLOCK BOOTSTRAP (WITH REPLACEMENT) --------------------------
    boot_blocks <- sample(block_ids, size = K_blocks, replace = TRUE)
    boot_idx    <- unlist(idx_by_block[boot_blocks], use.names = FALSE)  # keep duplicates
    sp_subset   <- joined_sp[boot_idx, ]
    block_subset <- cbind(as.data.frame(sp_subset@coords), sp_subset@data)
    block_subset$year <- as.factor(block_subset$year)
    
    sampled_data <- data.frame(
      response  = block_subset[[sp]],
      offset    = block_subset$offset,
      year      = block_subset$year,
      x_AEP10TM = block_subset$x_AEP10TM,
      y_AEP10TM = block_subset$y_AEP10TM,
      block_subset %>%
        dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
        dplyr::select(where(is.numeric))
    )
    
    if (sum(sampled_data$response, na.rm = TRUE) == 0) return(NULL)
    stopifnot(is.factor(sampled_data$year))
    
    # ---- BLOCKED CV FOLDS (NO LEAKAGE) ------------------------------------
    uniq_blocks <- sort(unique(as.character(sp_subset@data$block_id)))
    k_folds     <- min(10L, length(uniq_blocks))
    fold_map    <- setNames(((seq_along(uniq_blocks)-1L) %% k_folds) + 1L, uniq_blocks)
    folds       <- as.integer(fold_map[as.character(sp_subset@data$block_id)])
    
    # ---- FIT ---------------------------------------------------------------
    brt_model <- tryCatch({
      gbm.step(
        data             = sampled_data,
        gbm.y            = 1,
        gbm.x            = 3:ncol(sampled_data),
        family           = "poisson",
        tree.complexity  = 3,
        learning.rate    = 0.005,
        bag.fraction     = 0.5,
        offset           = sampled_data$offset,
        fold.vector      = folds,     # BLOCKED folds
        n.folds          = k_folds,   # match actual number of folds
        keep.fold.models = FALSE,
        keep.fold.fit    = FALSE,
        verbose          = FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(brt_model)) return(NULL)
    
    imp_df <- data.frame(
      Feature = summary(brt_model)$var,
      Gain    = summary(brt_model)$rel.inf
    )
    
    model_data <- sampled_data %>%
      dplyr::select(all_of(brt_model$gbm.call$predictor.names)) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ tidyr::replace_na(., 0)))
    
    interaction_df <- .gbm.interactions_full(brt_model, model_data)
    interaction_df$Boot <- i
    
    # PDPs for ALL predictors
    var_names <- brt_model$gbm.call$predictor.names
    pdp_df <- purrr::map_dfr(seq_along(var_names), function(k) {
      g <- try(gbm::plot.gbm(brt_model, i.var = k, return.grid = TRUE), silent = TRUE)
      if (inherits(g, "try-error") || !"y" %in% names(g)) return(NULL)
      colnames(g)[1] <- "x"
      tibble(
        x    = as.numeric(g$x),
        y    = as.numeric(g$y),  # link scale (log for Poisson)
        var  = var_names[k],
        boot = i
      )
    })
    
    list(imp = imp_df, interaction = interaction_df, pdp = pdp_df)
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
  
  # Knight-style PDP summarization for this species
  pdp_all <- lapply(boot_results, function(x) x$pdp)
  pdp_all <- pdp_all[!sapply(pdp_all, is.null)]
  if (length(pdp_all) > 0) {
    pdp_raw_sp <- bind_rows(pdp_all)
    if (nrow(pdp_raw_sp) > 0) {
      pdp_raw_sp <- pdp_raw_sp %>%
        filter(!var %in% nuisance_vars) %>%
        mutate(Species = sp)
      saveRDS(pdp_raw_sp, file.path(pdp_dir, paste0("BBWA_pdp_bootstrap_raw_", sp, ".rds")))
      
      pdp_gam_sp <- pdp_raw_sp %>%
        group_by(Species, var) %>%
        group_modify(~ fit_gam_pdp(.x)) %>%
        ungroup()
      
      saveRDS(pdp_gam_sp, file.path(pdp_dir, paste0("BBWA_pdp_gam_summary_", sp, ".rds")))
    }
  }
  
  # === Final model on full data (with BLOCKED CV) ===
  joined_sp@data$year <- as.factor(joined_sp@data$year)
  
  full_data <- cbind(as.data.frame(joined_sp@coords), joined_sp@data)
  model_data <- data.frame(
    response  = full_data[[sp]],
    offset    = full_data$offset,
    year      = full_data$year,
    x_AEP10TM = full_data$x_AEP10TM,
    y_AEP10TM = full_data$y_AEP10TM,
    full_data %>%
      dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
      dplyr::select(where(is.numeric))
  )
  stopifnot(is.factor(model_data$year))
  
  # ---- BLOCKED CV ON FULL DATA (DETERMINISTIC) ----
  # For full-data CV
  uniq_blocks_full <- sort(unique(as.character(joined_sp@data$block_id)))
  k_folds_full     <- min(10L, length(uniq_blocks_full))
  fold_map_full    <- setNames(((seq_along(uniq_blocks_full)-1L) %% k_folds_full) + 1L, uniq_blocks_full)
  folds_full       <- as.integer(fold_map_full[as.character(joined_sp@data$block_id)])
  
  
  final_model <- gbm.step(
    data             = model_data,
    gbm.y            = 1,
    gbm.x            = 3:ncol(model_data),   # <<< DO NOT drop last column
    family           = "poisson",
    tree.complexity  = 3,
    learning.rate    = 0.005,
    bag.fraction     = 0.5,
    offset           = model_data$offset,
    fold.vector      = folds_full,           # <<< BLOCKED folds on full data
    n.folds          = k_folds_full,
    keep.fold.models = FALSE,
    keep.fold.fit    = FALSE,
    verbose          = FALSE
  )
  
  final_model_list[[sp]] <- final_model
}

# === Stop Parallel Cluster ===
stopCluster(cl)

# === Save Results ===
imp_df <- bind_rows(results_list)
friedman_df <- bind_rows(friedmanH_list)

imp_table <- imp_df %>%
  dplyr::mutate(var_type = dplyr::case_when(
    stringr::str_starts(Feature, "prop_con") ~ "prop_con",
    stringr::str_starts(Feature, "clumpy")   ~ "clumpy",
    stringr::str_starts(Feature, "age_mn")   ~ "age_mn",
    TRUE ~ "other"
  ))
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)

saveRDS(results_list, "Output/Tabular Data/BBWA_no_wetlands_results_list_final_parallel_poisson_5m.rds")
saveRDS(imp_df, "Output/Tabular Data/BBWA_no_wetlands_imp_df_final_parallel_poisson_5m.rds")
saveRDS(friedman_df, "Output/Tabular Data/BBWA_no_wetlands_friedman_interactions_bootstrapped_parallel_final_poisson_5m.rds")
saveRDS(final_model_list, "Output/Tabular Data/BBWA_no_wetlands_final_model_list_parallel_poisson_5m.rds")
saveRDS(imp_table, file.path(pdp_dir, "BBWA_no_wetlands_imp_table_5m.rds"))

cat("\n✅ Saved all model outputs for Poisson models with year as FACTOR!\n")








### ============================
### BBWA GAM-PDPs at Scale of Effect (no wetlands except swamp)
### ============================


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(mgcv)
  library(readr)
})

# --- Locations ---
pdp_dir           <- "Output/Tabular Data"
imp_table_path    <- file.path("Output/Tabular Data", "BBWA_no_wetlands_imp_table_5m.rds")
imp_fallback_path <- file.path("Output/Tabular Data", "BBWA_no_wetlands_imp_df_final_parallel_poisson_5m.rds")

# --- Options ---
n_points      <- 300    # x-grid resolution for smooth curves
trim_prop     <- 0.01   # trim 1% tails to avoid edge wiggles
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")

# --- Outputs ---
out_tab_dir <- "Output/Tables/BBWA_no_wetlands_GAM_PDP_5m"
out_fig_dir <- "Output/Figures/BBWA_no_wetlands_GAM_PDP_5m"
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----
nice_lab <- function(v) {
  v %>%
    str_replace("_NTEMS$", " (NTEMS)") %>%
    str_replace("_LULC$",  " (LULC)")  %>%
    str_replace("^prop_con", "Proportion conifer") %>%
    str_replace("^clumpy",   "Clumpiness") %>%
    str_replace("^age_mn",   "Forest age")
}

var_type_from_name <- function(v) case_when(
  str_starts(v, "prop_con") ~ "prop_con",
  str_starts(v, "clumpy")   ~ "clumpy",
  str_starts(v, "age_mn")   ~ "age_mn",
  TRUE ~ "other"
)

slugify <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace("^_|_$", "")
}

`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a

# ---------- Load or rebuild imp_table ----------
if (exists("imp_table") && is.data.frame(imp_table)) {
  message("Using in-memory imp_table")
} else if (file.exists(imp_table_path)) {
  imp_table <- readRDS(imp_table_path)
  message("Loaded imp_table from ", imp_table_path)
} else if (file.exists(imp_fallback_path)) {
  message("imp_table.rds not found; deriving from ", imp_fallback_path)
  imp_table <- readRDS(imp_fallback_path) %>%
    mutate(var_type = var_type_from_name(Feature)) %>%
    filter(var_type %in% c("prop_con","clumpy","age_mn"))
} else {
  stop("Couldn't find an importance table. Expected one of:\n",
       " - ", imp_table_path, "\n",
       " - ", imp_fallback_path, "\n",
       "Or define imp_table in memory before running.")
}

imp_table <- imp_table %>%
  mutate(var_type = factor(var_type, levels = c("prop_con","clumpy","age_mn")))

# -------------------------------------------
# 1) Choose Scale-of-Effect (SoE) per species
#    = argmax MeanGain within each var_type
# -------------------------------------------
soe_tbl <- imp_table %>%
  filter(var_type %in% c("prop_con","clumpy","age_mn")) %>%
  group_by(Species, var_type) %>%
  arrange(desc(MeanGain), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(Species, var_type, Feature_soe = Feature)

print(soe_tbl, n = Inf)

# Species list (use what’s actually present in SoE table)
species_list <- unique(soe_tbl$Species) %||% c("BBWA")

# ---- Fit one population-level GAM with random intercept per bootstrap ----
fit_gam_pdp_one <- function(df, n = n_points, trim = trim_prop) {
  # df columns: x, y (link), var, boot, Species
  df <- df %>%
    mutate(
      boot   = factor(boot),
      y_resp = exp(y)               # Poisson inverse link
    ) %>%
    filter(is.finite(x), is.finite(y_resp))
  
  if (nrow(df) < 10) return(NULL)
  
  # Trim extremes for stability
  q   <- quantile(df$x, probs = c(trim, 1 - trim), na.rm = TRUE)
  xlo <- q[1]; xhi <- q[2]
  if (!is.finite(xlo) || !is.finite(xhi) || xlo >= xhi) return(NULL)
  
  # Include a dummy boot factor with correct levels for predict()
  xnew <- tibble(
    x    = seq(xlo, xhi, length.out = n),
    boot = factor(df$boot[1], levels = levels(df$boot))
  )
  
  g <- mgcv::gam(
    y_resp ~ s(x) + s(boot, bs = "re"),
    data   = df,
    family = gaussian(),
    method = "REML",
    select = TRUE
  )
  
  pr <- as.data.frame(predict(g, newdata = xnew, se.fit = TRUE, exclude = "s(boot)"))
  
  tibble(
    x   = xnew$x,
    fit = pr$fit,
    lwr = pr$fit - 1.96 * pr$se.fit,
    upr = pr$fit + 1.96 * pr$se.fit,
    edf = summary(g)$s.table[rownames(summary(g)$s.table) == "s(x)", "edf"] %||% NA_real_,
    r2  = summary(g)$r.sq %||% NA_real_
  )
}

# -------------------------------------------
# 2) Build GAM-PDPs only for those SoE features
#    and SAVE ONE FILE PER PREDICTOR
# -------------------------------------------
all_plots <- list()

for (sp in species_list) {
  pdp_file <- file.path(pdp_dir, paste0("BBWA_pdp_bootstrap_raw_", sp, ".rds"))
  if (!file.exists(pdp_file)) {
    message("Missing PDP file for ", sp, " — skipping.")
    next
  }
  
  pdp_raw <- readRDS(pdp_file) %>%
    filter(!var %in% nuisance_vars) %>%
    drop_na(x, y)
  
  # Features selected as SoE for this species
  soe_vars_sp <- soe_tbl %>% filter(Species == sp) %>% pull(Feature_soe)
  if (length(soe_vars_sp) == 0) {
    message("No SoE features for ", sp, " — skipping.")
    next
  }
  
  pdp_soe <- pdp_raw %>% filter(var %in% soe_vars_sp)
  if (nrow(pdp_soe) == 0) {
    message("No PDP rows for SoE features in ", sp, " — skipping.")
    next
  }
  
  # Diagnostics: counts per SoE variable
  message("Counts for ", sp, " (SoE):")
  print(pdp_soe %>% count(var), n = Inf)
  
  # Fit one GAM per selected feature and SAVE EACH SEPARATELY
  unique_vars <- unique(pdp_soe$var)
  for (v in unique_vars) {
    df_v <- pdp_soe %>% filter(var == v)
    sm   <- try(fit_gam_pdp_one(df_v), silent = TRUE)
    if (inherits(sm, "try-error") || is.null(sm)) {
      message("Skipping ", sp, " / ", v, " (insufficient data or fit error).")
      next
    }
    
    # Tidy + save table per variable
    out_tab <- sm %>%
      mutate(Species = sp, var = v, var_clean = nice_lab(v)) %>%
      relocate(Species, var, var_clean)
    csv_name <- file.path(out_tab_dir, paste0("GAM_PDP_SoE_", sp, "_", slugify(v), ".csv"))
    write_csv(out_tab, csv_name)
    
    # Subtitle shows chosen SoE for this var_type
    vt <- var_type_from_name(v)
    subtitle_txt <- paste0("Selected SoE (", vt, "): ", v)
    
    p <- ggplot(out_tab, aes(x = x, y = fit)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey30", colour = NA, alpha = 0.25) +
      geom_line() +
      labs(
        title    = paste0("GAM-smoothed PDP at Scale of Effect (", sp, ")"),
        subtitle = subtitle_txt,
        x = "Predictor value",
        y = "Expected count"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank())
    
    png_name <- file.path(out_fig_dir, paste0("GAM_PDP_SoE_", sp, "_", slugify(v), ".png"))
    ggsave(png_name, p, width = 7.5, height = 5.0, dpi = 300)
    
    all_plots[[paste(sp, v, sep = "_")]] <- p
  }
}

# Print
print(all_plots)








### BBWA VI PLOT - No wetlands (but including swamp)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggtext)
  library(cowplot)
})

# --- Load variable importance results ---
imp_df <- readRDS("Output/Tabular Data/BBWA_no_wetlands_imp_df_final_parallel_poisson_5m.rds")

# --- Config ---
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")
fill_map <- c("prop_con" = "forestgreen",
              "clumpy"   = "cornflowerblue",
              "age_mn"   = "lightcoral")
var_order <- c("prop_con","clumpy","age_mn")

# --- Clean & classify (BBWA only) ---
imp_df_clean <- imp_df %>%
  filter(Species == "BBWA") %>%
  filter(!Feature %in% nuisance_vars) %>%
  # drop LULC age terms if undesired; comment out next line to keep them
  filter(!grepl("^age_mn_\\d+_5$", Feature)) %>%
  mutate(
    var_type = case_when(
      str_detect(Feature, "^prop_con") ~ "prop_con",
      str_detect(Feature, "^clumpy")   ~ "clumpy",
      str_detect(Feature, "^age_mn")   ~ "age_mn",
      TRUE                             ~ NA_character_
    ),
    var_type = factor(var_type, levels = var_order)
  ) %>%
  filter(!is.na(var_type))

# --- Auto-scale axis to data (this is the plotted "y", but becomes x after coord_flip) ---
y_max_data <- max(imp_df_clean$MeanGain, na.rm = TRUE)
y_max <- max(1, ceiling(y_max_data * 1.08))  # ~8% headroom, at least 1
y_breaks <- pretty(c(0, y_max), n = 6)
y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_max]

# --- One-species plot (grouped vertically by var_type; left group labels hidden) ---
make_vi_plot <- function(df_sp) {
  df_sp <- df_sp %>%
    group_by(var_type) %>%
    arrange(desc(MeanGain), .by_group = TRUE) %>%
    mutate(is_top = row_number() == 1L) %>%
    ungroup() %>%
    mutate(
      Feature       = fct_rev(fct_inorder(Feature)),  # order for coord_flip
      Feature_label = ifelse(is_top,
                             paste0("**", as.character(Feature), "**"),
                             as.character(Feature))
    )
  
  lab_map <- setNames(df_sp$Feature_label, df_sp$Feature)
  
  ggplot(df_sp, aes(x = Feature, y = MeanGain, fill = var_type)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_manual(
      values = fill_map, drop = FALSE,
      name = "Predictor type",
      labels = c("prop_con" = "Proportion conifer",
                 "clumpy"   = "Clumpiness",
                 "age_mn"   = "Forest age")
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.06))
    ) +
    scale_x_discrete(labels = lab_map) +
    coord_flip() +
    facet_grid(var_type ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = "Predictor variables",
      y = "Mean relative influence"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      axis.text          = element_text(color = "black"),
      axis.text.y        = ggtext::element_markdown(hjust = 1, size = 11),
      
      strip.text.y       = element_blank(),
      strip.background   = element_blank(),
      
      legend.position    = "right",
      legend.title       = element_text(face = "bold")
    )
}

final_plot <- make_vi_plot(imp_df_clean)

# --- Save ---
out_fig <- "Output/Figures/BBWA_NO_WETLAND_5m_VI.png"
dir.create(dirname(out_fig), showWarnings = FALSE, recursive = TRUE)
ggsave(out_fig, final_plot, width = 9, height = 10, dpi = 300)
message("✅ Saved: ", out_fig)

print(out_fig)


