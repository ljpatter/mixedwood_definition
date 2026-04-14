# ---
# title: "Multiscale model"
# author: "Leonard Patterson"
# created: "2025-08-11"
# description: This script runs a final multi scale model using scale of effect predictors from the global BRT model
# ---

# === FINAL MULTISCALE BRT MODEL USING ONLY SCALE-OF-EFFECT PREDICTORS (BLOCK-BOOT + BLOCK-CV) ===

# --- CLEAN ENV ---
rm(list = ls())

# --- LIBRARIES ---
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
library(landscapemetrics)
library(patchwork)
library(dplyr)
library(gbm)


# --- Load variable importance (used to pick scale-of-effect predictors) ---
imp_df <- readRDS("Output/Tabular Data/imp_df_final_parallel_poisson_5m.rds") %>%
  mutate(
    var_type = case_when(
      str_detect(Feature, "prop_con") ~ "composition",
      str_detect(Feature, "clumpy")   ~ "configuration",
      str_detect(Feature, "age_mn")   ~ "structure",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(var_type))

# --- Data & Study Area ---
study_area  <- st_read("Input/Spatial Data/Study Area/study_area.shp")

# Create sf points directly (no proj4string / sp conversion needed)
joined_sf <- st_as_sf(
  readr::read_csv("Output/Tabular Data/joined_data_5m.csv"),
  coords = c("x_AEP10TM", "y_AEP10TM"),
  crs = st_crs(study_area)
)

# --- 1 km block grid ---
grid_1km <- st_make_grid(study_area, cellsize = 1000, square = TRUE) %>%
  st_sf(grid_id = 1:length(.), geometry = .)

# --- Clip to grid, attach blocks ---
joined_sf <- joined_sf[
  st_within(joined_sf, st_union(grid_1km), sparse = FALSE) %>% apply(1, any),
]

joined_with_blocks <- st_join(joined_sf, grid_1km, join = st_within) %>%
  dplyr::filter(!is.na(grid_id))
joined_with_blocks$block_id <- as.factor(joined_with_blocks$grid_id)

# --- Back to Spatial and restore coord names (needed for downstream code) ---
joined_sp <- as(joined_with_blocks, "Spatial")
colnames(joined_sp@coords) <- c("x_AEP10TM", "y_AEP10TM")

# >>> Precompute row indices by block (for TRUE block bootstrap with replacement) ---
idx_by_block <- split(seq_len(nrow(joined_sp)), joined_sp@data$block_id)
block_ids    <- names(idx_by_block)              # only blocks that have data
K_blocks     <- length(block_ids)
# -----------------------------------------------------------------------------------

# ====================== BOOTSTRAPPED MULTISCALE BRT ======================

# --- SETUP ---
species_list <- c("BTNW", "TEWA", "BBWA")
n_boot       <- 100
cores_to_use <- 6

# --- Parallel Backend + deterministic RNG ---
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)
registerDoRNG(123)  # reproducible parallel RNG

# --- QPAD Offsets ---
load_BAM_QPAD(version = 3)
calculate_offsets <- function(data, species_code) {
  localBAMcorrections(
    species = species_code,
    t   = data$survey_effort,
    r   = Inf,
    jday= data$ordinalDay,
    tssr= data$hssr
  ) %>% corrections2offset()
}

# --- Pick top (scale-of-effect) predictor per var_type per species ---
best_by_type <- imp_df %>%
  group_by(Species, var_type) %>%
  slice_max(MeanGain, n = 1, with_ties = FALSE) %>%
  ungroup()

# --- Output containers ---
boot_model_list <- list()
boot_model_data <- list()

# --- Species loop ---
for (sp in species_list) {
  cat("\n\n=== Bootstrapped Multiscale (Scale-of-Effect) Model for:", sp, "===\n")
  
  # Top predictors for this species (1 per var_type)
  top_vars <- best_by_type %>%
    filter(Species == sp) %>%
    pull(Feature)
  
  if (length(top_vars) == 0) {
    warning(sprintf("No top vars found in imp_df for %s; skipping.", sp))
    next
  }
  
  # Offsets + factor year
  joined_sp@data$offset <- calculate_offsets(joined_sp@data, sp)
  joined_sp@data$year   <- as.factor(joined_sp@data$year)
  
  # Parallel block-bootstraps
  boot_models <- foreach(i = 1:n_boot,
                         .packages = c("dismo", "gbm", "dplyr", "tibble", "tidyr")) %dopar% {
                           # TRUE Block Bootstrap (WITH replacement): draw K_blocks blocks
                           boot_blocks <- sample(block_ids, size = K_blocks, replace = TRUE)
                           boot_idx    <- unlist(idx_by_block[boot_blocks], use.names = FALSE)  # keep duplicates
                           sp_subset   <- joined_sp[boot_idx, ]
                           
                           # Build modeling frame in row order of sp_subset
                           block_subset <- cbind(as.data.frame(sp_subset@coords), sp_subset@data)
                           block_subset$year <- as.factor(block_subset$year)
                           
                           # Keep only selected (scale-of-effect) predictors; fail cleanly if any missing
                           missing_vars <- setdiff(top_vars, names(block_subset))
                           if (length(missing_vars) > 0) {
                             warning(sprintf("Bootstrap %d %s: missing predictors: %s", i, sp, paste(missing_vars, collapse = ", ")))
                             return(NULL)
                           }
                           
                           sampled_data <- data.frame(
                             response  = block_subset[[sp]],
                             offset    = block_subset$offset,
                             year      = block_subset$year,
                             x_AEP10TM = block_subset$x_AEP10TM,
                             y_AEP10TM = block_subset$y_AEP10TM,
                             block_subset[, top_vars, drop = FALSE]
                           )
                           
                           # --- NA & degenerate checks (avoid gbm.step failures) ---
                           keep <- stats::complete.cases(sampled_data[, c("response", "offset", top_vars), drop = FALSE])
                           if (!all(keep)) sampled_data <- sampled_data[keep, , drop = FALSE]
                           if (nrow(sampled_data) == 0L || sum(sampled_data$response, na.rm = TRUE) == 0) return(NULL)
                           stopifnot(is.factor(sampled_data$year))
                           
                           # ---- BLOCKED CV (no leakage across blocks) ----
                           uniq_blocks <- sort(unique(as.character(sp_subset@data$block_id)))
                           k_folds     <- min(10L, length(uniq_blocks))
                           if (k_folds < 2L) return(NULL)   # not enough blocks for CV
                           fold_map    <- setNames(((seq_along(uniq_blocks) - 1L) %% k_folds) + 1L, uniq_blocks)
                           folds       <- as.integer(fold_map[as.character(sp_subset@data$block_id)])
                           
                           # Fit BRT (Poisson)
                           model_fit <- tryCatch({
                             gbm.step(
                               data             = sampled_data,
                               gbm.y            = 1,
                               gbm.x            = 3:ncol(sampled_data),   # year + coords + top_vars
                               family           = "poisson",
                               tree.complexity  = 3,
                               learning.rate    = 0.005,
                               bag.fraction     = 0.5,
                               offset           = sampled_data$offset,
                               fold.vector      = folds,                   # BLOCKED folds
                               n.folds          = k_folds,
                               keep.fold.models = FALSE,
                               keep.fold.fit    = FALSE,
                               verbose          = FALSE
                             )
                           }, error = function(e) NULL)
                           
                           if (is.null(model_fit)) return(NULL)
                           list(model = model_fit, data = sampled_data)
                         }
  
  # Keep successful fits
  boot_models <- boot_models[!sapply(boot_models, is.null)]
  boot_model_list[[sp]] <- lapply(boot_models, `[[`, "model")
  boot_model_data[[sp]] <- lapply(boot_models, `[[`, "data")
  
  cat(sprintf("%d successful bootstraps for %s\n", length(boot_models), sp))
}

# --- Stop cluster ---
stopCluster(cl)

# --- Save outputs ---
saveRDS(boot_model_list, "Output/Tabular Data/bootstrapped_multiscale_model_list_5m.rds")
saveRDS(boot_model_data, "Output/Tabular Data/bootstrapped_multiscale_data_list_5m.rds")
cat("\n Saved all bootstrapped multiscale models (scale-of-effect predictors + true block bootstrap + blocked CV).\n")











########## Extract variable importance values

# === Re-load saved bootstrapped models ===
boot_model_list <- readRDS("Output/Tabular Data/bootstrapped_multiscale_model_list_5m.rds")

# === Extract variable importance from each model ===
imp_df <- lapply(names(boot_model_list), function(sp) {
  bind_rows(lapply(boot_model_list[[sp]], function(mod) {
    data.frame(
      Feature = summary(mod)$var,
      Gain = summary(mod)$rel.inf,
      Species = sp
    )
  }))
}) %>% bind_rows()

# === Tag variable types ===
imp_df <- imp_df %>%
  mutate(
    var_type = case_when(
      str_detect(Feature, "prop_con") ~ "composition",
      str_detect(Feature, "clumpy") ~ "configuration",
      str_detect(Feature, "age_mn") ~ "structure",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(var_type))

# === Select top predictor per type/species ===
best_by_type <- imp_df %>%
  group_by(Species, var_type, Feature) %>%
  summarise(MeanGain = mean(Gain), .groups = "drop") %>%
  group_by(Species, var_type) %>%
  slice_max(MeanGain, n = 1, with_ties = FALSE) %>%
  ungroup()

# === Rebuild per-species predictor list ===
top_vars_per_species <- best_by_type %>%
  group_by(Species) %>%
  summarise(vars = list(Feature), .groups = "drop") %>%
  deframe()









################### Run final model 

# === Define offset function ===
calculate_offsets <- function(data, species_code) {
  localBAMcorrections(
    species = species_code,
    t = data$survey_effort,
    r = Inf,
    jday = data$ordinalDay,
    tssr = data$hssr
  ) %>% corrections2offset()
}

# === Load top variables ===
top_vars_per_species <- readRDS("Output/Tabular Data/top_vars_per_species.rds")

# === Create model list to store final BRTs ===
final_model_list <- list()

# === Loop over species ===
for (sp in names(top_vars_per_species)) {
  cat("\nTraining final model for:", sp, "\n")
  
  top_vars <- top_vars_per_species[[sp]]
  
  # Generate offsets per species
  joined_sp@data$offset <- calculate_offsets(joined_sp@data, sp)
  joined_sp@data$year <- as.factor(joined_sp@data$year)
  
  # Build model data
  model_data <- data.frame(
    response = joined_sp@data[[sp]],
    offset = joined_sp@data$offset,
    year = joined_sp@data$year,
    joined_sp@coords,
    joined_sp@data[, top_vars, drop = FALSE]
  )
  
  # Fit final BRT model with 10-fold CV
  set.seed(123)
  folds <- sample(rep(1:10, length.out = nrow(model_data)))
  
  model_fit <- gbm.step(
    data = model_data,
    gbm.y = 1,
    gbm.x = 4:ncol(model_data),
    family = "poisson",
    tree.complexity = 3,
    learning.rate = 0.005,
    bag.fraction = 0.5,
    offset = model_data$offset,
    fold.vector = folds,
    n.folds = 10,
    keep.fold.models = FALSE,
    keep.fold.fit = FALSE,
    verbose = FALSE
  )
  
  final_model_list[[sp]] <- model_fit
}


# Save
saveRDS(final_model_list, "Output/Tabular Data/final_model_list.rds")
saveRDS(imp_df, "Output/Tabular Data/imp_df_multiscale_bootstrapped.rds")
saveRDS(top_vars_per_species, "Output/Tabular Data/top_vars_per_species.rds")











################ PDPs for SOE  MODEL

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr)
  library(stringr); library(tibble)
  library(ggplot2); library(gbm); library(mgcv); library(readr)
})

# -----------------------------
# Paths & output locations
# -----------------------------
boot_path    <- "Output/Tabular Data/bootstrapped_multiscale_model_list_5m.rds"
out_tab_dir  <- "Output/Tables/GAM_PDP_SoE"
out_fig_dir  <- "Output/Figures/GAM_PDP_SoE"
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load bootstrapped SoE models
# -----------------------------
if (!file.exists(boot_path)) stop("Can't find: ", boot_path)
boot_model_list <- readRDS(boot_path)  # list: Species -> list of gbm.step objects

# -----------------------------
# Settings
# -----------------------------
nuisance_vars <- c("year","x_AEP10TM","y_AEP10TM")
# Optionally limit boots for speed while testing; set to Inf to use all
max_boot_per_sp <- Inf
# trim tails when smoothing GAM (avoids edge wiggles)
trim_prop <- 0.01
# x grid points for GAM smooth
n_points <- 300

nice_lab <- function(v) {
  v %>%
    str_replace("_NTEMS$", " (NTEMS)") %>%
    str_replace("_LULC$",  " (LULC)")  %>%
    str_replace("^prop_con", "Proportion conifer") %>%
    str_replace("^clumpy",   "Clumpiness") %>%
    str_replace("^age_mn",   "Forest age")
}

slugify <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace("^_|_$", "")
}

`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a

# -----------------------------
# Helper: fit one GAM to PDPs
# df columns: x, y (link), var, boot, Species
# Returns x, fit, lwr, upr on RESPONSE scale (exp(link))
# -----------------------------
fit_gam_pdp_one <- function(df, n = n_points, trim = trim_prop) {
  df <- df %>%
    mutate(
      boot   = factor(boot),
      y_resp = exp(y)               # back-transform from link (Poisson)
    ) %>%
    filter(is.finite(x), is.finite(y_resp))
  if (nrow(df) < 10) return(NULL)
  
  q   <- quantile(df$x, probs = c(trim, 1 - trim), na.rm = TRUE)
  xlo <- q[1]; xhi <- q[2]
  if (!is.finite(xlo) || !is.finite(xhi) || xlo >= xhi) return(NULL)
  
  # Dummy boot (levels retained) for prediction excluding boot RE
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

# -----------------------------
# Extract PDPs across bootstraps
# -----------------------------
all_outputs <- list()

for (sp in names(boot_model_list)) {
  models_sp <- boot_model_list[[sp]]
  if (!length(models_sp)) next
  
  # Keep only non-null models and limit boots if desired
  ok_idx <- which(!vapply(models_sp, is.null, logical(1)))
  if (!length(ok_idx)) next
  if (is.finite(max_boot_per_sp)) ok_idx <- head(ok_idx, max_boot_per_sp)
  models_sp <- models_sp[ok_idx]
  
  # Grab predictor set from the first successful model
  m_ref <- models_sp[[1]]
  pred_names <- setdiff(m_ref$gbm.call$predictor.names, nuisance_vars)
  pred_names <- pred_names[nchar(pred_names) > 0]
  if (!length(pred_names)) next
  
  message("== ", sp, " : predictors = ", paste(pred_names, collapse=", "))
  
  # ---- RAW PDP aggregation across boots ----
  pdp_raw <- purrr::imap_dfr(models_sp, function(m, bidx) {
    # For each predictor, call plot.gbm (returns link-scale y)
    purrr::map_dfr(seq_along(pred_names), function(k) {
      v <- pred_names[k]
      gr <- try(
        gbm::plot.gbm(m, i.var = v, n.trees = m$gbm.call$best.trees,
                      return.grid = TRUE),
        silent = TRUE
      )
      if (inherits(gr, "try-error") || !("y" %in% names(gr))) return(NULL)
      # Standardize to x / y columns
      if (!"x" %in% names(gr)) {
        # gbm::plot.gbm names first column with predictor's name
        colnames(gr)[1] <- "x"
      } else {
        # ensure numeric x
        gr$x <- as.numeric(gr$x)
      }
      tibble(
        x      = as.numeric(gr$x),
        y      = as.numeric(gr$y),  # link scale
        var    = v,
        boot   = as.integer(bidx),
        Species= sp
      ) %>% filter(is.finite(x), is.finite(y))
    })
  })
  
  if (nrow(pdp_raw) == 0) {
    warning("No PDP rows for ", sp, " — skipping.")
    next
  }
  
  # Save raw PDPs (optional audit)
  saveRDS(pdp_raw, file.path(out_tab_dir, paste0("pdp_bootstrap_raw_SOE_", sp, ".rds")))
  
  # ---- Fit one GAM per predictor using boot RE ----
  pdp_summaries <- pdp_raw %>%
    group_by(Species, var) %>%
    group_modify(~ fit_gam_pdp_one(.x, n = n_points, trim = trim_prop) %||% tibble()) %>%
    ungroup()
  
  if (nrow(pdp_summaries) == 0) {
    warning("No GAM summaries for ", sp, " — skipping.")
    next
  }
  
  # Tidy labels and save per predictor CSV + PNG
  vars_sp <- unique(pdp_summaries$var)
  for (v in vars_sp) {
    df_v <- pdp_summaries %>% filter(Species == sp, var == v)
    if (nrow(df_v) == 0) next
    
    df_v <- df_v %>% mutate(var_clean = nice_lab(v)) %>%
      relocate(Species, var, var_clean)
    
    csv_name <- file.path(out_tab_dir, paste0("GAM_PDP_SOE_boot_", sp, "_", slugify(v), ".csv"))
    readr::write_csv(df_v, csv_name)
    
    p <- ggplot(df_v, aes(x = x, y = fit)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) +
      geom_line() +
      labs(
        title = paste0("GAM-smoothed bootstrap PDP (", sp, ")"),
        subtitle = v,
        x = "Predictor value",
        y = "Expected count (Poisson, inverse link)"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank())
    
    png_name <- file.path(out_fig_dir, paste0("GAM_PDP_SOE_boot_", sp, "_", slugify(v), ".png"))
    ggsave(png_name, p, width = 7.5, height = 5.0, dpi = 300)
  }
  
  # Combined multi-panel per species
  p_all <- pdp_summaries %>%
    mutate(var_clean = nice_lab(var)) %>%
    ggplot(aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) +
    geom_line() +
    facet_wrap(~ var_clean, scales = "free_x") +
    labs(
      title = paste0("GAM-smoothed bootstrap PDPs — ", sp),
      x = "Predictor value", y = "Expected count"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())
  
  ggsave(file.path(out_fig_dir, paste0("GAM_PDP_SOE_boot__combined_", sp, ".png")),
         p_all, width = 11, height = 7, dpi = 300)
  
  all_outputs[[sp]] <- list(raw = pdp_raw, summary = pdp_summaries)
}

invisible(all_outputs)

