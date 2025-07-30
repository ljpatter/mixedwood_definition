# ---
# title: "BRT model w/ random mixedwood 5m"
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
library(doRNG)
library(foreach)
library(tidyr)
library(viridis)
library(knitr)
library(tibble)
library(terra)
library(dplyr)
library(landscapemetrics)
library(patchwork)
library(ggh4x)
library(scales)
library(grid)


# Load data
joined_data <- read.csv("Output/Tabular Data/joined_data_5m.csv")

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
registerDoRNG(123)

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
n_boot <- 100

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

saveRDS(results_list, "Output/Tabular Data/results_list_final_parallel_poisson_5m.rds")
saveRDS(imp_df, "Output/Tabular Data/imp_df_final_parallel_poisson_5m.rds")
saveRDS(friedman_df, "Output/Tabular Data/friedman_interactions_bootstrapped_parallel_final_poisson_5m.rds")
saveRDS(final_model_list, "Output/Tabular Data/final_model_list_parallel_poisson_5m.rds")
saveRDS(imp_table, file.path(pdp_dir, "imp_table_5m.rds"))

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
imp_table_path    <- file.path("Output/Tabular Data", "imp_table_5m.rds")
imp_fallback_path <- file.path("Output/Tabular Data", "imp_df_final_parallel_poisson_5m.rds")

# --- Options ---
n_points      <- 300    # x-grid resolution for smooth curves
trim_prop     <- 0.01   # trim 1% tails to avoid edge wiggles
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")

# --- Outputs ---
out_tab_dir <- "Output/Tables/GAM_PDP_SoE_5m"
out_fig_dir <- "Output/Figures/GAM_PDP_SoE_5m"
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





























######## Plot PDP in single window

# ---- One-frame overview with per-panel y-scaling: rows = Species; cols = prop_con | clumpy | age ----

# Rebuild the long table if needed (from the CSVs you just wrote)
if (!exists("all_tabs") || length(all_tabs) == 0) {
  csvs <- list.files(out_tab_dir, pattern = "^GAM_PDP_SoE_.*\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) stop("No per-predictor CSVs found in: ", out_tab_dir)
  all_tabs_df <- purrr::map_dfr(csvs, readr::read_csv, show_col_types = FALSE)
} else {
  all_tabs_df <- dplyr::bind_rows(all_tabs)
}

order_cols <- c("prop_con","clumpy","age")
order_species <- c("BTNW","TEWA","BBWA")   # desired top-to-bottom order

all_tabs_df <- all_tabs_df %>%
  dplyr::mutate(
    var_col = dplyr::recode(var_type_from_name(var),
                            "prop_con" = "prop_con",
                            "clumpy"   = "clumpy",
                            "age_mn"   = "age",
                            .default   = NA_character_),
    var_col = factor(var_col, levels = order_cols)
  ) %>%
  dplyr::filter(!is.na(var_col), var_col %in% order_cols) %>%
  # Apply factor order only to species that actually exist in the data
  dplyr::mutate(Species = factor(Species,
                                 levels = intersect(order_species,
                                                    unique(Species))))

# --- Build per-panel labels from SoE + MeanGain (as before) ---
imp_table_path    <- file.path("Output/Tabular Data", "imp_table_5.rds")
imp_fallback_path <- file.path("Output/Tabular Data", "imp_df_final_parallel_poisson_5m.rds")

var_type_from_name <- function(v) dplyr::case_when(
  stringr::str_starts(v, "prop_con") ~ "prop_con",
  stringr::str_starts(v, "clumpy")   ~ "clumpy",
  stringr::str_starts(v, "age_mn")   ~ "age_mn",
  TRUE ~ "other"
)

if (exists("imp_table") && is.data.frame(imp_table)) {
  imp_use <- imp_table
} else if (file.exists(imp_table_path)) {
  imp_use <- readRDS(imp_table_path)
} else {
  imp_use <- readRDS(imp_fallback_path) |>
    dplyr::mutate(var_type = var_type_from_name(Feature)) |>
    dplyr::filter(var_type %in% c("prop_con","clumpy","age_mn"))
}

soe_tbl <- imp_use |>
  dplyr::filter(var_type %in% c("prop_con","clumpy","age_mn")) |>
  dplyr::group_by(Species, var_type) |>
  dplyr::slice_max(MeanGain, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    Species,
    var_col = dplyr::recode(var_type,
                            "prop_con" = "prop_con",
                            "clumpy"   = "clumpy",
                            "age_mn"   = "age"),
    panel_label = sprintf("%s (%.1f%%)", Feature, MeanGain)
  )

# --- X midpoints for centering text above each panel ---
x_mids <- all_tabs_df |>
  dplyr::group_by(Species, var_col) |>
  dplyr::summarise(x_mid = mean(range(x, na.rm = TRUE)), .groups = "drop")

label_df <- dplyr::inner_join(soe_tbl, x_mids, by = c("Species","var_col"))

# --- Factor orders (same as before) ---
order_cols    <- c("prop_con","clumpy","age")
order_species <- c("BTNW","TEWA","BBWA")

all_tabs_df <- all_tabs_df |>
  dplyr::mutate(
    var_col  = factor(var_col, levels = order_cols),
    Species  = factor(Species, levels = intersect(order_species, unique(Species)))
  )

label_df <- label_df |>
  dplyr::mutate(
    var_col = factor(var_col, levels = order_cols),
    Species = factor(Species, levels = levels(all_tabs_df$Species))
  )

p_all <- ggplot(all_tabs_df, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) +
  geom_line() +
  # ensure ~5 y-axis ticks in each facet (even with free scales)
  scale_y_continuous(breaks = breaks_extended(n = 5)) +  # <-- added '+' here
  geom_text(
    data = label_df,
    aes(x = x_mid, y = Inf, label = panel_label),
    inherit.aes = FALSE,
    vjust = -0.55,
    fontface = "bold",
    size = 3.6
  ) +
  ggh4x::facet_grid2(
    rows = vars(Species),
    cols = vars(var_col),
    scales = "free",
    independent = "all"
  ) +
  labs(x = "Predictor value", y = "Expected count") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.6),
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_rect(colour = "grey30", fill = "grey95"),
    panel.spacing.y = unit(24, "pt"),
    plot.margin = margin(t = 24, r = 10, b = 10, l = 10),
    plot.title = element_blank()
  )

print(p_all)
ggsave(file.path(out_fig_dir, "GAM_PDP_SoE__combined_grid_freeY_5m.png"),
       p_all, width = 12,
       height = 4 * length(unique(all_tabs_df$Species)), dpi = 300)


















##################### Poisson BRT Variable Importance Analysis #####################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggtext)
  library(cowplot)  
})

# --- Load variable importance results ---
imp_df <- readRDS("Output/Tabular Data/imp_df_final_parallel_poisson_5m.rds")

# --- Config ---
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")
fill_map <- c("prop_con" = "forestgreen",
              "clumpy"   = "cornflowerblue",
              "age_mn"   = "lightcoral")
var_order     <- c("prop_con","clumpy","age_mn")
species_order <- c("BTNW","TEWA","BBWA")   # set order as desired

# --- Clean & classify ---
imp_df_clean <- imp_df %>%
  filter(!Feature %in% nuisance_vars) %>%
  # drop LULC age terms if undesired; comment out next line to keep them
  filter(!grepl("^age_mn_\\d+_LULC$", Feature)) %>%
  mutate(
    var_type = case_when(
      str_detect(Feature, "^prop_con") ~ "prop_con",
      str_detect(Feature, "^clumpy")   ~ "clumpy",
      str_detect(Feature, "^age_mn")   ~ "age_mn",
      TRUE                             ~ "other"
    ),
    var_type = factor(var_type, levels = var_order),
    Species  = factor(Species, levels = intersect(species_order, unique(Species)))
  ) %>% 
  filter(!is.na(var_type), !is.na(Species))

# --- Fixed x axis across all panels ---
x_max    <- 25
x_breaks <- seq(0, x_max, by = 5)

# --- One-species plot (grouped vertically by var_type; left group labels hidden) ---
make_vi_plot <- function(df_sp, show_bottom_title = FALSE) {
  df_sp <- df_sp %>%
    group_by(var_type) %>%
    arrange(desc(MeanGain), .by_group = TRUE) %>%
    mutate(is_top = row_number() == 1L) %>%           # top per var_type
    ungroup() %>%
    mutate(
      Feature       = fct_rev(fct_inorder(Feature)),  # keep order for coord_flip
      Feature_label = ifelse(is_top, paste0("**", as.character(Feature), "**"),
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
    scale_y_continuous(limits = c(0, x_max),
                       breaks = x_breaks,
                       expand = expansion(mult = c(0, 0.06))) +
    scale_x_discrete(labels = lab_map) +
    coord_flip() +
    facet_grid(var_type ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = NULL,                                            # global left label added later
      y = if (show_bottom_title) "Mean relative influence" else NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      # remove vertical grid lines (after coord_flip these are x-grid)
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      # remove horizontal grid lines too (clean look)
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      axis.text          = element_text(color = "black"),
      axis.text.y        = ggtext::element_markdown(hjust = 1, size = 11),
      axis.title.y       = element_text(margin = margin(r = 10)),
      # hide left facet (group) labels
      strip.text.y       = element_blank(),
      strip.background   = element_blank(),
      # legend kept only on one plot (extracted below)
      legend.position    = "right",
      legend.title       = element_text(face = "bold"),
      plot.title         = element_blank()
    )
}

# --- Build species panels (only bottom shows horizontal axis title) ---
species_levels <- levels(imp_df_clean$Species)
plots_with_legend <- lapply(seq_along(species_levels), function(i) {
  sp <- species_levels[i]
  make_vi_plot(filter(imp_df_clean, Species == sp),
               show_bottom_title = (i == length(species_levels)))
})
names(plots_with_legend) <- species_levels

# --- Extract legend from the first plot, then remove from panels ---
legend_g <- cowplot::get_legend(plots_with_legend[[1]])
plots <- lapply(plots_with_legend, function(p) p + theme(legend.position = "none"))

# --- Add a)/b)/c) as a small tag above each species panel ---
tag_strip <- function(lbl) {
  ggplot() +
    annotate("text", x = 0.01, y = 0.5, label = lbl,
             hjust = 0, vjust = 0.5, size = 6, fontface = "bold") +
    theme_void() + theme(plot.margin = margin(b = 0))
}

panels_with_tags <- Map(function(p, taglbl) {
  cowplot::plot_grid(tag_strip(taglbl), p, ncol = 1,
                     rel_heights = c(0.08, 0.92), align = "v")
}, plots, paste0(letters[seq_along(plots)], ")"))

# --- Stack the tagged species panels vertically (right side) ---
right_stack <- cowplot::plot_grid(plotlist = panels_with_tags, ncol = 1,
                                  rel_heights = rep(1, length(panels_with_tags)),
                                  align = "v")

# --- Global left label: "Predictor variables" centered on stack ---
left_ylabel <- cowplot::ggdraw() +
  cowplot::draw_label("Predictor variables", angle = 90,
                      x = 0.5, y = 0.5, size = 14)

main_body <- cowplot::plot_grid(left_ylabel, right_stack, ncol = 2,
                                rel_widths = c(0.05, 0.95), align = "h")

# --- Top row: put legend on the right; spacer on the left so legend is top-right ---
top_row <- cowplot::plot_grid(NULL, cowplot::ggdraw(legend_g),
                              ncol = 2, rel_widths = c(0.70, 0.30), align = "h")

# --- Final assembly ---
final_plot <- cowplot::plot_grid(top_row, main_body, ncol = 1,
                                 rel_heights = c(0.12, 0.88), align = "v")

# --- Save ---
out_fig <- "Output/Figures/VI_all_species_5m_combined.png"
dir.create(dirname(out_fig), showWarnings = FALSE, recursive = TRUE)
ggsave(out_fig, final_plot, width = 9, height = 13, dpi = 300)
message("✅ Saved: ", out_fig)

print(out_fig)




















################# H-statistics top interactions


# Load interactions
all_friedman_interactions <- readRDS("Output/Tabular Data/friedman_interactions_bootstrapped_parallel_final_poisson_5m.rds")

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

write.csv(top5_filtered, "Output/Tabular Data/interactions_table_5m.csv")















# ===========================
# 3-D Interaction Perspec Plots (by NAME)
# ===========================

library(dismo)
library(dplyr)
library(stringr)

# 1) Load the final models from your latest run
fm_path <- "Output/Tabular Data/final_model_list_parallel_poisson_5m.rds"
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
  x_var = "clumpy_150_LULC",
  y_var = "prop_con_150_LULC",
  z_range = NULL,                         # let it auto-scale, or set e.g. c(0, 0.6)
  main = "BBWA: prop_con_1_LULC × clumpy_1_LULC"
)

## TEWA: clumpy_1000_LULC × prop_con_1000_LULC (adjust if names differ)
gbm_perspec_by_name(
  model = final_model_list[["TEWA"]],
  x_var = "clumpy_1000_LULC",
  y_var = "prop_con_1000_LULC",
  z_range = NULL,
  main = "TEWA: prop_con_1000_NTEMS × clumpy_1000_NTEMS"
)

## BTNW: pick any two continuous predictors you care about
# For example, choose by looking at list_predictors(final_model_list[["BTNW"]])
gbm_perspec_by_name(
  model = final_model_list[["BTNW"]],
  x_var = "clumpy_150_LULC",
  y_var = "prop_con_150_LULC",
  z_range = NULL,
  main = "BTNW: prop_con_500_NTEMS × clumpy_1_LULC"
)









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
  library(QPAD)
  library(readr)
})

# ---------------------------
# -1) QPAD setup + offset helpers
# ---------------------------
if (is.null(getOption("BAMversion"))) {
  load_BAM_QPAD(version = 2)
}

species_list <- c("BTNW", "TEWA", "BBWA")

species_code_map <- c(
  "BTNW" = "BTNW",
  "TEWA" = "TEWA",
  "BBWA" = "BBWA"
)

calculate_offsets <- function(df, species_label) {
  sp_code <- species_code_map[[species_label]]
  if (is.na(sp_code)) stop("Unknown species label: ", species_label)
  
  t_col    <- dplyr::coalesce(df$survey_effort, df$t, df$duration)
  jday_col <- dplyr::coalesce(df$ordinalDay, df$jday, df$yday)
  tssr_col <- dplyr::coalesce(df$hssr, df$tssr, rep(NA_real_, nrow(df)))
  
  if (is.null(t_col))    stop("No survey effort column found (expected one of: survey_effort, t, duration).")
  if (is.null(jday_col)) stop("No ordinal day column found (expected one of: ordinalDay, jday, yday).")
  
  t_col[!is.finite(t_col)]       <- median(t_col[is.finite(t_col)], na.rm = TRUE)
  jday_col[!is.finite(jday_col)] <- median(jday_col[is.finite(jday_col)], na.rm = TRUE)
  if (all(!is.finite(tssr_col))) tssr_col <- rep(0, length(t_col))
  tssr_col[!is.finite(tssr_col)] <- 0
  
  corr <- localBAMcorrections(
    species = sp_code,
    t       = t_col,
    r       = Inf,
    jday    = jday_col,
    tssr    = tssr_col
  )
  corrections2offset(corr)
}

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

# >>> inject coordinates as predictors <<<
data_use$x_AEP10TM <- coords_use[, 1]
data_use$y_AEP10TM <- coords_use[, 2]

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
  
  # Load needed packages and QPAD on each worker
  clusterEvalQ(cl_mi, {
    suppressPackageStartupMessages({
      library(sp); library(spdep); library(dplyr); library(tidyr)
      library(gbm); library(dismo); library(QPAD); library(tibble)
    })
    if (is.null(getOption("BAMversion"))) load_BAM_QPAD(version = 2)
    TRUE
  })
  
  # Export helper objects/functions and data to workers
  clusterExport(cl_mi, c(
    "joined_sp", "idx_by_block", "block_ids", "K_blocks",
    "coords_use", "data_use", "lw_use",
    "species_list", "species_code_map", "calculate_offsets"
  ), envir = environment())
}

# ---------------------------
# 3) Settings
# ---------------------------
n_boot <- 100L  # set to your bootstrap count

# ---------------------------
# 4) Core worker: one bootstrap → Moran's I
# ---------------------------
compute_moran_one_boot <- function(sp, boot_id) {
  # Safety: ensure QPAD is loaded on the worker
  if (is.null(getOption("BAMversion"))) {
    load_BAM_QPAD(version = 2)
  }
  
  # ---- TRUE BLOCK BOOTSTRAP (train set) ----
  boot_blocks <- sample(block_ids, size = K_blocks, replace = TRUE)
  boot_idx    <- unlist(idx_by_block[boot_blocks], use.names = FALSE)
  sp_subset   <- joined_sp[boot_idx, ]
  
  # Build training data: take attributes, append coords
  coords_train <- sp_subset@coords
  df_train <- sp_subset@data %>%
    mutate(
      year      = as.factor(year),
      x_AEP10TM = coords_train[, 1],
      y_AEP10TM = coords_train[, 2]
    )
  
  # QPAD offset for training subset (species-specific)
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
  y_eval      <- data_use[[sp]]
  offset_eval <- calculate_offsets(data_use, sp)
  
  # predictors required by this bootstrap model
  pred_names <- brt_model$gbm.call$predictor.names
  
  # Build eval predictors and ensure any missing columns exist
  eval_pred <- data_use %>%
    dplyr::select(any_of(pred_names)) %>%
    mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))
  
  # If some predictors are still missing (e.g., rare splits), add zero columns
  missing_preds <- setdiff(pred_names, names(eval_pred))
  if (length(missing_preds)) {
    for (mp in missing_preds) eval_pred[[mp]] <- 0
    eval_pred <- eval_pred[, pred_names, drop = FALSE]
  }
  
  # link-scale prediction and mean on response scale (Poisson with offset)
  eta  <- predict(brt_model, eval_pred, n.trees = brt_model$n.trees, type = "link")
  mu   <- exp(eta + offset_eval)
  
  # residuals (Pearson-like simple)
  resids <- y_eval - mu
  
  # guardrails
  if (!all(is.finite(resids))) {
    ok <- which(is.finite(resids))
    if (length(ok) < 3) return(NA_real_)
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
                        .packages = c("sp", "spdep", "dplyr", "gbm", "dismo", "tidyr", "QPAD")) %dopar% {
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

