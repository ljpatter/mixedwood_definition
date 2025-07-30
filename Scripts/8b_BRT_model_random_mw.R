# ---
# title: "BRT model w/ random mixedwood"
# author: "Leonard Patterson"
# created: "2025-09-01"
# description: This script runs boosted regression model with boostrapping and associated plots and diagnostics
#              For prop con, this code uses the LULC layer where mixedwood pixels are randomly assigned to conifer and deciduous.
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

# Load data
joined_data <- read.csv("Output/Tabular Data/joined_data_NEW.csv")

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

# >>> NEW: precompute row indices by block for TRUE block bootstrap ----------
idx_by_block <- split(seq_len(nrow(joined_sp)), joined_sp@data$block_id)
block_ids    <- names(idx_by_block)           # only blocks that actually have data
K_blocks     <- length(block_ids)
# ---------------------------------------------------------------------------


###################### BRT POISSON W/ PARALLEL BOOTSTRAPPING

# === SETUP ===
species_list <- c("BTNW", "TEWA", "BBWA")
results_list <- list()
friedmanH_list <- list()
final_model_list <- list()

# === Start Parallel Backend ===
cores_to_use <- 10
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

# === Offsets (QPAD/BAM) ===
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

# Knight-style PDP helpers
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
    uniq_blocks <- unique(sp_subset@data$block_id)
    k_folds     <- min(10L, length(uniq_blocks))
    fold_labels <- rep(seq_len(k_folds), length.out = length(uniq_blocks))
    fold_map    <- setNames(sample(fold_labels), as.character(uniq_blocks))
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
      mutate(across(everything(), ~tidyr::replace_na(., 0)))
    
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
      saveRDS(pdp_raw_sp, file.path(pdp_dir, paste0("pdp_bootstrap_raw_", sp, ".rds")))
      
      pdp_gam_sp <- pdp_raw_sp %>%
        group_by(Species, var) %>%
        group_modify(~ fit_gam_pdp(.x)) %>%
        ungroup()
      
      saveRDS(pdp_gam_sp, file.path(pdp_dir, paste0("pdp_gam_summary_", sp, ".rds")))
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
  
  # Blocked folds for full-data CV
  uniq_blocks_full <- unique(joined_sp@data$block_id)
  k_folds_full     <- min(10L, length(uniq_blocks_full))
  fold_labels_full <- rep(seq_len(k_folds_full), length.out = length(uniq_blocks_full))
  fold_map_full    <- setNames(sample(fold_labels_full), as.character(uniq_blocks_full))
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

saveRDS(results_list, "Output/Tabular Data/results_list_final_parallel_poisson.rds")
saveRDS(imp_df, "Output/Tabular Data/imp_df_final_parallel_poisson.rds")
saveRDS(friedman_df, "Output/Tabular Data/friedman_interactions_bootstrapped_parallel_final_poisson.rds")
saveRDS(final_model_list, "Output/Tabular Data/final_model_list_parallel_poisson.rds")

cat("\n✅ Saved all model outputs for Poisson models with year as FACTOR!\n")















########### GAM PDPs


# ============================
# GAM-PDPs at Scale of Effect
# ============================

# Run this AFTER your modeling script has saved:
#   - "Output/Tabular Data/pdp_bootstrap_raw_<Species>.rds"
#   - "Output/Tabular Data/imp_table.rds"  (optional)
#   - "Output/Tabular Data/imp_df_final_parallel_poisson.rds" (fallback)
# ============================
# GAM-PDPs at Scale of Effect
# ============================

# Run this AFTER your modeling script has saved:
#   - "Output/Tabular Data/pdp_bootstrap_raw_<Species>.rds"
#   - "Output/Tabular Data/imp_table.rds"  (optional)
#   - "Output/Tabular Data/imp_df_final_parallel_poisson.rds" (fallback)

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
imp_table_path    <- file.path("Output/Tabular Data", "imp_table.rds")
imp_fallback_path <- file.path("Output/Tabular Data", "imp_df_final_parallel_poisson.rds")

# --- Options ---
n_points      <- 300    # x-grid resolution for smooth curves
trim_prop     <- 0.01   # trim 1% tails to avoid edge wiggles
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")

# --- Outputs ---
out_tab_dir <- "Output/Tables/GAM_PDP_SoE"
out_fig_dir <- "Output/Figures/GAM_PDP_SoE"
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
species_list <- unique(soe_tbl$Species) %||% c("BBWA","BTNW","TEWA")

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
  pdp_file <- file.path(pdp_dir, paste0("pdp_bootstrap_raw_", sp, ".rds"))
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
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) +
      geom_line() +
      labs(
        title    = paste0("GAM-smoothed PDP at Scale of Effect (", sp, ")"),
        subtitle = subtitle_txt,
        x = "Predictor value",
        y = "Expected count (Poisson, inverse link)"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank())
    
    png_name <- file.path(out_fig_dir, paste0("GAM_PDP_SoE_", sp, "_", slugify(v), ".png"))
    ggsave(png_name, p, width = 7.5, height = 5.0, dpi = 300)
    
    all_plots[[paste(sp, v, sep = "_")]] <- p
  }
}

invisible(all_plots)




















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
imp_df <- readRDS("Output/Tabular Data/imp_df_final_parallel_poisson.rds")

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










######### PDP ON RAW SCALE

# === FULL PDPs ON RAW SCALE, CENTERED ===

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
    
    # Back-transform to raw expected count and center
    pd$y_raw <- exp(pd$y)
    pd$y_raw_centered <- pd$y_raw / mean(pd$y_raw, na.rm = TRUE)
    
    # Extract relative contribution
    contrib <- var_importance %>%
      filter(Feature == var) %>%
      pull(MeanGain) %>%
      round(1)
    
    title_text <- paste0(sp, " - ", var, " (", contrib, "%)")
    
    # Plot
    p <- ggplot(pd, aes_string(x = var, y = "y_raw_centered")) +
      geom_line(color = "black", linewidth = 1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      labs(
        title = title_text,
        x = var,
        y = "Centered Expected Count (raw scale)"
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

# Top 10 per species
top5_filtered <- filtered_df %>%
  group_by(Species) %>%
  arrange(desc(H)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Show the pretty table
kable(top5_filtered, caption = "Top 10 Friedman H Interactions per Species (Strict Filter)")


write.csv(top5_filtered, "Output/Tabular Data/interactions_table.csv")















# ===========================
# 3-D Interaction Perspec Plots (by NAME)
# ===========================

library(dismo)
library(dplyr)
library(stringr)

# 1) Load the final models from your latest run
fm_path <- "Output/Tabular Data/final_model_list_parallel_poisson.rds"
if (!file.exists(fm_path)) stop("Can't find: ", fm_path)
final_model_list <- readRDS(fm_path)

# 2) Helper: list predictors with indices 
list_predictors <- function(model) {
  tibble(
    idx  = seq_along(model$gbm.call$predictor.names),
    name = model$gbm.call$predictor.names,
    cls  = vapply(model$gbm.call$dataframe[model$gbm.call$predictor.names], function(x) class(x)[1], character(1))
  )
}

# Example: see predictors for each species
cat("\nBBWA predictors:\n"); print(list_predictors(final_model_list[["BBWA"]]), n=Inf)
cat("\nTEWA predictors:\n"); print(list_predictors(final_model_list[["TEWA"]]), n=Inf)
cat("\nBTNW predictors:\n"); print(list_predictors(final_model_list[["BTNW"]]), n=Inf)

# 3) Helper: safe perspec by variable NAMES
gbm_perspec_by_name <- function(model, x_var, y_var,
                                x_range = NULL, y_range = NULL, z_range = NULL,
                                main = NULL, ...) {
  # sanity
  stopifnot(!is.null(model), is.character(x_var), is.character(y_var))
  pn <- model$gbm.call$predictor.names
  if (!(x_var %in% pn)) stop("x_var not in model predictors: ", x_var)
  if (!(y_var %in% pn)) stop("y_var not in model predictors: ", y_var)
  
  # reject factors (perspec assumes continuous surfaces)
  df <- model$gbm.call$dataframe
  if (is.factor(df[[x_var]]) || is.factor(df[[y_var]]))
    stop("One of the variables is a factor. Choose two continuous predictors.")
  
  # indices
  xi <- match(x_var, pn)
  yi <- match(y_var, pn)
  
  # auto ranges if not supplied
  if (is.null(x_range)) x_range <- range(df[[x_var]], na.rm = TRUE)
  if (is.null(y_range)) y_range <- range(df[[y_var]], na.rm = TRUE)
  
  # sensible main
  if (is.null(main)) main <- paste(x_var, "×", y_var)
  
  # Draw the surface (type='response' is what dismo::gbm.perspec uses for gbm)
  dismo::gbm.perspec(
    gbm.object = model,
    x = xi,
    y = yi,
    y.range = y_range,
    z.range = z_range,
    main = main,
    ...
  )
}

# 4) Examples (use names you actually have!)
# Tip: from your importance table and earlier code, names look like:
#   prop_con_1_LULC / prop_con_2_LULC / prop_con_3_LULC
#   clumpy_150_NTEMS / clumpy_500_NTEMS / clumpy_1000_NTEMS
#   age_mn_150_NTEMS / age_mn_500_NTEMS / age_mn_1000_NTEMS

## BBWA: prop_con_1_LULC × clumpy_1_LULC
gbm_perspec_by_name(
  model = final_model_list[["BBWA"]],
  x_var = "clumpy_1_LULC",
  y_var = "prop_con_1_LULC",
  z_range = NULL,                         # let it auto-scale, or set e.g. c(0, 0.6)
  main = "BBWA: prop_con_1_LULC × clumpy_1_LULC"
)

## TEWA: clumpy_1000_NTEMS × prop_con_1000_NTEMS (adjust if names differ)
gbm_perspec_by_name(
  model = final_model_list[["TEWA"]],
  x_var = "prop_con_1000_NTEMS",
  y_var = "clumpy_1000_NTEMS",
  z_range = NULL,
  main = "TEWA: prop_con_1000_NTEMS × clumpy_1000_NTEMS"
)

## BTNW: pick any two continuous predictors you care about
# For example, choose by looking at list_predictors(final_model_list[["BTNW"]])
gbm_perspec_by_name(
  model = final_model_list[["BTNW"]],
  x_var = "prop_con_500_NTEMS",
  y_var = "clumpy_1_LULC",
  z_range = NULL,
  main = "BTNW: prop_con_500_NTEMS × clumpy_1_LULC"
)












################ Moran's I for Poisson Final Models ################

################ Moran's I for Poisson BRTs — Average Across Bootstraps ################

suppressPackageStartupMessages({
  library(sp)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(foreach)
  library(doParallel)
  library(gbm)
  library(dismo)
})

# ---------------------------
# 0) Neighbor graph (once)
# ---------------------------
coords_sp <- coordinates(joined_sp)
rownames(coords_sp) <- rownames(joined_sp@data)

# 30 km neighbor radius (adjust if desired)
nb_all <- dnearneigh(coords_sp, 0, 30000)
nb_all <- make.sym.nb(nb_all)
has_neighbors <- which(card(nb_all) > 0)

coords_use <- coords_sp[has_neighbors, , drop = FALSE]
data_use   <- joined_sp@data[has_neighbors, , drop = FALSE]
data_use$year <- as.factor(data_use$year)

# neighbor weights for the fixed evaluation set
nb_use <- dnearneigh(coords_use, 0, 30000)
nb_use <- make.sym.nb(nb_use)
lw_use <- nb2listw(nb_use, style = "W", zero.policy = TRUE)

# ---------------------------
# 1) TRUE block bootstrap helper
# ---------------------------
if (!exists("idx_by_block") || !exists("block_ids") || !exists("K_blocks")) {
  idx_by_block <- split(seq_len(nrow(joined_sp)), joined_sp@data$block_id)
  block_ids    <- names(idx_by_block)
  K_blocks     <- length(block_ids)
}

# ---------------------------
# 2) Parallel backend (start if needed)
# ---------------------------
stop_after <- FALSE
if (!foreach::getDoParRegistered()) {
  cores_to_use <- max(1L, parallel::detectCores() - 1L)
  cl_mi <- makeCluster(cores_to_use)
  registerDoParallel(cl_mi)
  stop_after <- TRUE
}

# ---------------------------
# 3) Settings (reuse your n_boot / hyperparams)
# ---------------------------
n_boot <- 25L  # set to your bootstrap count

# ---------------------------
# 4) Core worker: one bootstrap → Moran's I
# ---------------------------
compute_moran_one_boot <- function(sp, boot_id) {
  # ---- TRUE BLOCK BOOTSTRAP (train set) ----
  boot_blocks <- sample(block_ids, size = K_blocks, replace = TRUE)
  boot_idx    <- unlist(idx_by_block[boot_blocks], use.names = FALSE)
  sp_subset   <- joined_sp[boot_idx, ]
  df_train    <- cbind(as.data.frame(sp_subset@coords), sp_subset@data)
  df_train$year <- as.factor(df_train$year)
  
  # QPAD offset for training subset (species-specific)
  # (joined_sp@data$offset may be set elsewhere; recompute to be safe)
  df_train$offset <- calculate_offsets(df_train, sp)
  
  # build modeling frame
  sampled_data <- data.frame(
    response  = df_train[[sp]],
    offset    = df_train$offset,
    year      = df_train$year,
    x_AEP10TM = df_train$x_AEP10TM,
    y_AEP10TM = df_train$y_AEP10TM,
    df_train %>%
      dplyr::select(matches("(_LULC|_NTEMS)$")) %>%
      dplyr::select(where(is.numeric))
  )
  
  # skip empty-positive samples
  if (sum(sampled_data$response, na.rm = TRUE) == 0) return(NA_real_)
  
  # ---- BLOCKED CV mapping (no leakage) ----
  uniq_blocks <- unique(sp_subset@data$block_id)
  k_folds     <- min(10L, length(uniq_blocks))
  fold_labels <- rep(seq_len(k_folds), length.out = length(uniq_blocks))
  fold_map    <- setNames(sample(fold_labels), as.character(uniq_blocks))
  folds       <- as.integer(fold_map[as.character(sp_subset@data$block_id)])
  
  # ---- Fit Poisson BRT ----
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
      fold.vector      = folds,
      n.folds          = k_folds,
      keep.fold.models = FALSE,
      keep.fold.fit    = FALSE,
      verbose          = FALSE
    )
  }, error = function(e) NULL)
  
  if (is.null(brt_model)) return(NA_real_)
  
  # ---- Predict on a FIXED evaluation set (all sites with neighbors) ----
  # species response & offset for eval set
  y_eval      <- data_use[[sp]]
  offset_eval <- calculate_offsets(data_use, sp)
  
  # predictors required by this bootstrap model
  pred_names <- brt_model$gbm.call$predictor.names
  eval_pred  <- data_use %>%
    dplyr::select(all_of(pred_names)) %>%
    mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))
  
  # link-scale prediction and mean on response scale (Poisson with offset)
  eta  <- predict(brt_model, eval_pred, n.trees = brt_model$n.trees, type = "link")
  mu   <- exp(eta + offset_eval)
  
  # residuals (Pearson-like)
  resids <- y_eval - mu
  
  # guardrails
  if (!all(is.finite(resids))) {
    ok <- which(is.finite(resids))
    if (length(ok) < 3) return(NA_real_)
    # rebuild weights for the OK subset (keeps Moran computation valid)
    nb_ok <- dnearneigh(coords_use[ok, , drop = FALSE], 0, 30000)
    nb_ok <- make.sym.nb(nb_ok)
    lw_ok <- nb2listw(nb_ok, style = "W", zero.policy = TRUE)
    v <- var(resids[ok], na.rm = TRUE)
    if (is.na(v) || v == 0) return(NA_real_)
    return(tryCatch(
      unname(moran.test(resids[ok], lw_ok, zero.policy = TRUE)$estimate[1]),
      error = function(e) NA_real_
    ))
  } else {
    # use precomputed weights on full eval set
    v <- var(resids, na.rm = TRUE)
    if (is.na(v) || v == 0) return(NA_real_)
    return(tryCatch(
      unname(moran.test(resids, lw_use, zero.policy = TRUE)$estimate[1]),
      error = function(e) NA_real_
    ))
  }
}

# ---------------------------
# 5) Run across bootstraps per species (in parallel)
# ---------------------------
moran_boot_summary <- list()

for (sp in species_list) {
  cat("\n=== Moran's I across bootstraps — ", sp, " ===\n")
  
  moran_vals <- foreach(i = seq_len(n_boot), .combine = c,
                        .packages = c("sp", "spdep", "dplyr", "gbm", "dismo", "tidyr")) %dopar% {
                          compute_moran_one_boot(sp, i)
                        }
  
  # summarize
  mi <- moran_vals[is.finite(moran_vals)]
  n_ok <- length(mi)
  mean_mi <- if (n_ok > 0) mean(mi) else NA_real_
  sd_mi   <- if (n_ok > 1)  sd(mi)   else NA_real_
  se_mi   <- if (n_ok > 1)  sd_mi / sqrt(n_ok) else NA_real_
  ci95_lo <- if (n_ok > 1)  mean_mi - 1.96 * se_mi else NA_real_
  ci95_hi <- if (n_ok > 1)  mean_mi + 1.96 * se_mi else NA_real_
  
  moran_boot_summary[[sp]] <- list(
    moran_vals = moran_vals,
    summary = tibble(
      Species = sp,
      n_boot_total = n_boot,
      n_boot_ok    = n_ok,
      mean_MI      = mean_mi,
      sd_MI        = sd_mi,
      se_MI        = se_mi,
      CI95_low     = ci95_lo,
      CI95_high    = ci95_hi
    )
  )
  
  print(moran_boot_summary[[sp]]$summary)
  # save per-species vectors
  saveRDS(moran_vals, file = file.path("Output/Tabular Data",
                                       paste0("moran_bootstrap_vals_", sp, ".rds")))
}

# combined summary table
moran_summary_tbl <- bind_rows(lapply(moran_boot_summary, `[[`, "summary"))
readr::write_csv(moran_summary_tbl,
                 file.path("Output/Tabular Data", "moran_bootstrap_summary.csv"))
print(moran_summary_tbl)

# ---------------------------
# 6) Stop cluster if we started it
# ---------------------------
if (isTRUE(stop_after)) {
  stopCluster(cl_mi)
  registerDoSEQ()
}

cat("\n✅ Moran's I across bootstraps computed and saved.\n")


















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











############## OG Synthetic landscape (smoothed surface)


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





