# ---
# title: "Predictions on synthetic landscapes"
# author: "Leonard Patterson"
# created: "2025-08-11"
# description: This script output from final multiscale BRT model to make abundance predictions for focal 
#              species on synthetic intimate and segregated mixed wood forests.
# ---



# =========================================================
# SYNTHETIC PREDICTIONS (THESIS-SAFE)
# Uses: boot_model_list, points_sf, metrics_points_df
# Only predicts with species-specific SoE predictors
# Applies Poisson offset correctly: exp(eta + offset)
# Averages across point-counts, then summarizes across bootstraps
# Scales to 10 km² township (1000 ha)
# =========================================================
# =========================================================
# SYNTHETIC PREDICTIONS (THESIS-SAFE)
# Uses: boot_model_list, points_sf, metrics_points_df
# Only predicts with species-specific SoE predictors
# Applies Poisson offset correctly: exp(eta + offset)
# Averages across point-counts, then summarizes across bootstraps
# Scales to 10 km² township (1000 ha)
# =========================================================


library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(sf)
library(gbm)
library(ggplot2)
library(readr)


# -----------------
# Settings
# -----------------
TOWNSHIP_HA <- 1000
CI_LO <- 0.025
CI_HI <- 0.975
species_order <- c("BBWA","BTNW","TEWA")
layout_levels <- c("Intimate","Segregated")
age_vals      <- c(30, 60, 90)

# -----------------
# Preconditions
# -----------------
stopifnot(exists("boot_model_list"), is.list(boot_model_list), length(boot_model_list) > 0)
stopifnot(exists("points_sf"), inherits(points_sf, "sf"))

if (!exists("metrics_points_df")) {
  metrics_points_df <- readRDS("Output/Tabular Data/propcon_clumpy_by_point_extent.rds")
}
stopifnot(all(c("x","y","layout","extent","prop_con","clumpy") %in% names(metrics_points_df)))

# -----------------
# Load importance table to derive species-specific SoE predictors
# -----------------
imp_table_path    <- "Output/Tabular Data/imp_table_5m.rds"
imp_fallback_path <- "Output/Tabular Data/imp_df_final_parallel_poisson_5m.rds"

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
} else if (file.exists(imp_fallback_path)) {
  imp_use <- readRDS(imp_fallback_path)
} else {
  stop("No importance table found (imp_table_5m.rds or imp_df_final_parallel_poisson_5m.rds).")
}

stopifnot(all(c("Species","Feature") %in% names(imp_use)))

imp_use <- imp_use %>%
  mutate(var_type = var_type_from_name(Feature)) %>%
  filter(var_type %in% c("prop_con","clumpy","age_mn")) %>%
  group_by(Species, var_type) %>%
  slice_max(MeanGain, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    soe_extent = as.integer(sub(".*_(\\d+)_.*$", "\\1", Feature)),
    soe_suffix = sub("^.*_\\d+_([A-Za-z0-9]+)$", "\\1", Feature)  # LULC / NTEMS
  )

# Safety check: every focal species must have all three SoE rows
needed <- tidyr::expand_grid(Species = species_order,
                             var_type = c("prop_con","clumpy","age_mn"))
chk <- needed %>% left_join(imp_use %>% select(Species,var_type,Feature), by = c("Species","var_type"))
if (any(is.na(chk$Feature))) {
  print(chk %>% filter(is.na(Feature)))
  stop("Missing SoE Feature for some Species×var_type. See rows above.")
}

soe_tbl <- imp_use %>% select(Species, var_type, Feature, soe_extent, soe_suffix)

# Helper to grab a species' SoE triplet
get_soe_for_species <- function(sp) {
  s <- soe_tbl %>% filter(Species == sp)
  list(
    prop = s %>% filter(var_type=="prop_con") %>% slice(1) %>% select(soe_extent, soe_suffix) %>% as.list(),
    clmp = s %>% filter(var_type=="clumpy")   %>% slice(1) %>% select(soe_extent, soe_suffix) %>% as.list(),
    age  = s %>% filter(var_type=="age_mn")   %>% slice(1) %>% select(soe_extent, soe_suffix) %>% as.list()
  )
}

# -----------------
# Utilities
# -----------------
`%||%` <- function(a,b) if (is.null(a)) b else a

first_non_null <- function(x) {
  if (is.null(x)) return(NULL)
  idx <- which(!vapply(x, is.null, logical(1)))
  if (length(idx) == 0) NULL else x[[idx[1]]]
}

get_training_df <- function(sp) {
  m <- first_non_null(boot_model_list[[sp]])
  if (!is.null(m) && !is.null(m$gbm.call$dataframe)) return(m$gbm.call$dataframe)
  if (exists("boot_model_data") && !is.null(boot_model_data) && !is.null(boot_model_data[[sp]]))
    return(first_non_null(boot_model_data[[sp]]))
  stop(sprintf("No training dataframe available for %s", sp))
}

# year factor: pick most common level + carry levels for predictions
pick_global_year <- function(df) {
  if (!"year" %in% names(df) || !is.factor(df$year)) return(NULL)
  tab <- sort(table(df$year), decreasing = TRUE)
  list(level = names(tab)[1], levels = levels(df$year))
}

# offset to use for synthetic predictions (no per-point offsets exist here)
mean_training_offset <- function(df) {
  if ("offset" %in% names(df) && is.numeric(df$offset)) {
    mean(df$offset, na.rm = TRUE)
  } else 0
}

# predict on link, add offset, then exp to response
predict_mu_with_offset <- function(model, newdata, offset_vec) {
  ntrees <- model$gbm.call$best.trees %||% model$n.trees
  if (is.null(ntrees) || is.na(ntrees) || ntrees <= 0) return(rep(NA_real_, nrow(newdata)))
  eta <- gbm::predict.gbm(model, newdata = newdata, n.trees = ntrees, type = "link")
  # NOTE: predict.gbm() does NOT add offsets; we add them on the link scale:
  exp(eta + offset_vec)
}

# Build model var name for a specific SoE feature (exact extent + suffix)
varname_for_soe <- function(prefix, extent, suffix) sprintf("%s_%d_%s", prefix, extent, suffix)

# -----------------
# Align per-point metrics to points_sf order (by coordinates)
# -----------------
coords <- sf::st_coordinates(points_sf)
pts_key <- tibble(
  x  = as.numeric(coords[,1]),
  y  = as.numeric(coords[,2])
) %>% mutate(xr = round(x, 3), yr = round(y, 3))

metrics_points_df_join <- metrics_points_df %>%
  mutate(xr = round(as.numeric(x), 3), yr = round(as.numeric(y), 3))

get_metric_vec <- function(layout_name, target_extent, metric = c("prop_con","clumpy")) {
  metric <- match.arg(metric)
  filt <- metrics_points_df_join %>%
    filter(layout == layout_name, extent == target_extent) %>%
    select(xr, yr, prop_con, clumpy)
  joined <- pts_key %>% left_join(filt, by = c("xr","yr"))
  out <- joined[[metric]]
  # Hard fail if any NA or length mismatch
  if (length(out) != nrow(points_sf)) {
    stop("Length mismatch for ", layout_name, " @ ", target_extent, "m (", metric, "). ",
         "Got ", length(out), " but expected ", nrow(points_sf), ".")
  }
  if (anyNA(out)) {
    stop("NA values in ", layout_name, " @ ", target_extent, "m (", metric, 
         "). Check points vs metrics alignment.")
  }
  as.numeric(out)
}

# -----------------
# MAIN LOOP
# -----------------
results <- list()

for (sp in intersect(names(boot_model_list), species_order)) {
  message("\nSpecies: ", sp)
  models_sp <- boot_model_list[[sp]]
  model_ref <- first_non_null(models_sp)
  if (is.null(model_ref)) { warning("No usable models for ", sp); next }
  
  df_train <- get_training_df(sp)
  year_fix <- pick_global_year(df_train)
  off_fix  <- mean_training_offset(df_train)
  
  soe <- get_soe_for_species(sp)
  stopifnot(!is.null(soe), length(soe$prop$soe_extent)==1, length(soe$clmp$soe_extent)==1, length(soe$age$soe_extent)==1)
  
  # Exact model predictor names that must exist for a bootstrap to be used
  v_prop_exact <- varname_for_soe("prop_con", soe$prop$soe_extent, soe$prop$soe_suffix)
  v_clmp_exact <- varname_for_soe("clumpy",   soe$clmp$soe_extent, soe$clmp$soe_suffix)
  v_age_exact  <- varname_for_soe("age_mn",   soe$age$soe_extent,  soe$age$soe_suffix)
  
  # sanity: ensure these predictors appear somewhere among the species' bootstraps
  if (!any(vapply(models_sp, function(m) !is.null(m) && all(c(v_prop_exact,v_clmp_exact,v_age_exact) %in% m$gbm.call$predictor.names), logical(1)))) {
    stop("For ", sp, ", no bootstrap contains all required SoE predictors: ",
         paste(c(v_prop_exact, v_clmp_exact, v_age_exact), collapse = ", "))
  }
  
  for (layout in layout_levels) {
    # per-point vectors at the SoE extents (these are the *actual* synthetic metrics)
    prop_vec <- get_metric_vec(layout, soe$prop$soe_extent, "prop_con")
    clmp_vec <- get_metric_vec(layout, soe$clmp$soe_extent, "clumpy")
    
    # base newdata (coordinates always included; many models include them)
    base_nd <- pts_key %>% transmute(x_AEP10TM = x, y_AEP10TM = y)
    
    for (age_val in age_vals) {
      per_boot_means <- vapply(models_sp, function(m) {
        if (is.null(m)) return(NA_real_)
        need <- m$gbm.call$predictor.names
        
        # Strict rule: skip bootstrap if it lacks any SoE predictor for this species
        if (!(v_prop_exact %in% need && v_clmp_exact %in% need && v_age_exact %in% need)) {
          return(NA_real_)
        }
        
        nd <- base_nd
        
        # year factor (use most frequent level from training; preserve levels)
        if ("year" %in% need && !is.null(year_fix)) {
          nd$year <- factor(year_fix$level, levels = year_fix$levels)
        }
        
        # attach SoE predictor columns with actual per-point values / selected age
        nd[[v_prop_exact]] <- prop_vec
        nd[[v_clmp_exact]] <- clmp_vec
        nd[[v_age_exact]]  <- age_val
        
        # add any other required predictors with safe defaults (0 numeric, first level for factors)
        miss <- setdiff(need, names(nd))
        if (length(miss)) {
          mdf <- m$gbm.call$dataframe
          for (v in miss) {
            if (!is.null(mdf) && v %in% names(mdf) && is.factor(mdf[[v]])) {
              nd[[v]] <- factor(rep(levels(mdf[[v]])[1], nrow(nd)), levels = levels(mdf[[v]]))
            } else {
              nd[[v]] <- 0
            }
          }
        }
        
        # column order to match model
        nd <- nd[, need, drop = FALSE]
        
        # Predict: link → add offset → exp
        mu <- predict_mu_with_offset(m, nd, rep(off_fix, nrow(nd)))
        
        # per-ha mean over all point-counts for this bootstrap
        mean(mu, na.rm = TRUE)
      }, numeric(1))
      
      # Keep only finite bootstraps
      valid <- is.finite(per_boot_means)
      n_contrib <- sum(valid)
      if (n_contrib == 0) {
        message(sprintf("… %s / %s / age=%d → 0 contributing bootstraps; skipping.",
                        sp, layout, age_val))
        next
      }
      per_boot_means <- per_boot_means[valid]
      
      mean_per_ha <- mean(per_boot_means, na.rm = TRUE)
      lwr_per_ha  <- as.numeric(quantile(per_boot_means, probs = CI_LO, na.rm = TRUE))
      upr_per_ha  <- as.numeric(quantile(per_boot_means, probs = CI_HI, na.rm = TRUE))
      
      results[[paste(sp, layout, age_val, sep = "_")]] <- tibble(
        Species   = sp,
        layout    = layout,
        age_class = age_val,
        mean_abundance_per_ha = mean_per_ha,
        lwr_abundance_per_ha  = lwr_per_ha,
        upr_abundance_per_ha  = upr_per_ha,
        total_abundance_mean  = mean_per_ha * TOWNSHIP_HA,
        total_abundance_lwr   = lwr_per_ha  * TOWNSHIP_HA,
        total_abundance_upr   = upr_per_ha  * TOWNSHIP_HA,
        n_boot_contrib        = n_contrib
      )
    }
  }
}

# -----------------
# Summarize & Plot
# -----------------
pred_summary <- bind_rows(results) %>%
  mutate(
    Species   = factor(Species, levels = species_order),
    layout    = factor(layout, levels = layout_levels),
    age_class = factor(age_class, levels = age_vals)
  )

if (!nrow(pred_summary)) stop("No scenarios produced any valid predictions — nothing to plot.")

message("\nContributing bootstraps by scenario:")
print(pred_summary %>% arrange(Species, layout, age_class) %>%
        select(Species, layout, age_class, n_boot_contrib), n = Inf)

pred_summary_clean <- pred_summary %>% filter(is.finite(total_abundance_mean))

dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)
write_csv(pred_summary_clean, "Output/Tabular Data/synthetic_predictions_SoE_point_based_5m.csv")

p <- ggplot(pred_summary_clean, aes(x = age_class, y = total_abundance_mean, fill = layout)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = total_abundance_lwr, ymax = total_abundance_upr),
    position = position_dodge(width = 0.7),
    width = 0.25, linewidth = 0.5
  ) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_fill_manual(values = c("Intimate" = "forestgreen", "Segregated" = "cornflowerblue")) +
  labs(
    x = "Forest age",
    y = "Predicted mean abundance",
    fill = "Spatial configuration"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.text       = element_text(size = 13, face = "bold")
  )

print(p)
dir.create("Output/Figures", recursive = TRUE, showWarnings = FALSE)
ggsave("Output/Figures/synth_preds_SoE_only_5m.png", p, width = 9, height = 7, dpi = 300)

# -----------------
# Sanity checks (optional, but good to keep)
# -----------------
# 1) Offset handling identity on one model_ref & scenario (uses average offset)
#    mean(exp(eta+off)) == mean(exp(eta))*exp(off) only if off is constant
do_offset_check <- TRUE
if (do_offset_check) {
  sp_ex <- intersect(names(boot_model_list), species_order)[1]
  m_ex  <- first_non_null(boot_model_list[[sp_ex]])
  if (!is.null(m_ex)) {
    df_train_ex <- get_training_df(sp_ex)
    off_fix_ex  <- mean_training_offset(df_train_ex)
    year_fix_ex <- pick_global_year(df_train_ex)
    soe_ex <- get_soe_for_species(sp_ex)
    v_prop_ex <- varname_for_soe("prop_con", soe_ex$prop$soe_extent, soe_ex$prop$soe_suffix)
    v_clmp_ex <- varname_for_soe("clumpy",   soe_ex$clmp$soe_extent, soe_ex$clmp$soe_suffix)
    v_age_ex  <- varname_for_soe("age_mn",   soe_ex$age$soe_extent,  soe_ex$age$soe_suffix)
    # choose a layout/age
    layout_ex <- layout_levels[1]; age_ex <- age_vals[1]
    prop_vec_ex <- get_metric_vec(layout_ex, soe_ex$prop$soe_extent, "prop_con")
    clmp_vec_ex <- get_metric_vec(layout_ex, soe_ex$clmp$soe_extent, "clumpy")
    nd <- pts_key %>% transmute(x_AEP10TM = x, y_AEP10TM = y)
    need <- m_ex$gbm.call$predictor.names
    if ("year" %in% need && !is.null(year_fix_ex)) {
      nd$year <- factor(year_fix_ex$level, levels = year_fix_ex$levels)
    }
    nd[[v_prop_ex]] <- prop_vec_ex
    nd[[v_clmp_ex]] <- clmp_vec_ex
    nd[[v_age_ex]]  <- age_ex
    miss <- setdiff(need, names(nd))
    if (length(miss)) {
      mdf <- m_ex$gbm.call$dataframe
      for (v in miss) {
        if (!is.null(mdf) && v %in% names(mdf) && is.factor(mdf[[v]])) {
          nd[[v]] <- factor(rep(levels(mdf[[v]])[1], nrow(nd)), levels = levels(mdf[[v]]))
        } else {
          nd[[v]] <- 0
        }
      }
    }
    nd <- nd[, need, drop = FALSE]
    eta <- gbm::predict.gbm(m_ex, newdata = nd, n.trees = (m_ex$gbm.call$best.trees %||% m_ex$n.trees), type = "link")
    mu_no_off <- exp(eta)
    mu_with   <- exp(eta + off_fix_ex)
    stopifnot(all.equal(mean(mu_with), mean(mu_no_off) * exp(off_fix_ex), tol = 1e-6))
    message("Offset sanity check passed for ", sp_ex, " (constant mean offset applied).")
  }
}



suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# --- 1) Build a tidy table with both layouts side-by-side and ratios ---
ratio_tbl <- pred_summary_clean %>%
  select(Species, age_class, layout, total_abundance_mean) %>%
  tidyr::pivot_wider(
    names_from = layout,
    values_from = total_abundance_mean
  ) %>%
  # keep only rows where both layouts exist
  filter(!is.na(Intimate), !is.na(Segregated)) %>%
  # avoid divide-by-zero; drop scenarios with zero denominators
  mutate(
    ratio_int_over_seg = ifelse(Segregated > 0, Intimate / Segregated, NA_real_),
    ratio_seg_over_int = ifelse(Intimate  > 0, Segregated / Intimate, NA_real_),
    higher_layout = case_when(
      is.finite(Intimate)   & is.finite(Segregated) & Intimate > Segregated ~ "Intimate higher",
      is.finite(Intimate)   & is.finite(Segregated) & Segregated > Intimate ~ "Segregated higher",
      is.finite(Intimate)   & is.finite(Segregated) & Intimate == Segregated ~ "Tie",
      TRUE ~ NA_character_
    )
  )

# Peek at scenario-level ratios (optional)
ratio_tbl_print <- ratio_tbl %>%
  arrange(Species, age_class) %>%
  select(Species, age_class, Intimate, Segregated,
         ratio_int_over_seg, ratio_seg_over_int, higher_layout)
print(ratio_tbl_print, n = Inf)

# --- 2) Per-species ranges of "times greater" when a given layout wins ---
# When Intimate is higher, report Intimate/Segregated range
range_int_wins <- ratio_tbl %>%
  filter(higher_layout == "Intimate higher", is.finite(ratio_int_over_seg)) %>%
  group_by(Species) %>%
  summarise(
    n_scenarios = dplyr::n(),
    min_times_greater = min(ratio_int_over_seg, na.rm = TRUE),
    max_times_greater = max(ratio_int_over_seg, na.rm = TRUE),
    .groups = "drop"
  )

# When Segregated is higher, report Segregated/Intimate range
range_seg_wins <- ratio_tbl %>%
  filter(higher_layout == "Segregated higher", is.finite(ratio_seg_over_int)) %>%
  group_by(Species) %>%
  summarise(
    n_scenarios = dplyr::n(),
    min_times_greater = min(ratio_seg_over_int, na.rm = TRUE),
    max_times_greater = max(ratio_seg_over_int, na.rm = TRUE),
    .groups = "drop"
  )

# --- 3) Overall ranges across all species/ages (optional) ---
overall_int_wins <- ratio_tbl %>%
  filter(higher_layout == "Intimate higher", is.finite(ratio_int_over_seg)) %>%
  summarise(
    n_scenarios = dplyr::n(),
    min_times_greater = min(ratio_int_over_seg, na.rm = TRUE),
    max_times_greater = max(ratio_int_over_seg, na.rm = TRUE)
  ) %>% mutate(layout = "Intimate higher")

overall_seg_wins <- ratio_tbl %>%
  filter(higher_layout == "Segregated higher", is.finite(ratio_seg_over_int)) %>%
  summarise(
    n_scenarios = dplyr::n(),
    min_times_greater = min(ratio_seg_over_int, na.rm = TRUE),
    max_times_greater = max(ratio_seg_over_int, na.rm = TRUE)
  ) %>% mutate(layout = "Segregated higher")

overall_ranges <- bind_rows(overall_int_wins, overall_seg_wins) %>%
  select(layout, n_scenarios, min_times_greater, max_times_greater)

# Print concise summaries
cat("\n=== Range when INTIMATE is higher (Intimate / Segregated) ===\n")
print(range_int_wins, n = Inf)

cat("\n=== Range when SEGREGATED is higher (Segregated / Intimate) ===\n")
print(range_seg_wins, n = Inf)

cat("\n=== Overall ranges across all species/ages ===\n")
print(overall_ranges, n = Inf)

# --- 4) (Optional) Save tables ---
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)
write_csv(ratio_tbl_print,  "Output/Tabular Data/intimate_segregated_ratios_by_scenario.csv")
write_csv(range_int_wins,   "Output/Tabular Data/range_when_intimate_higher.csv")
write_csv(range_seg_wins,   "Output/Tabular Data/range_when_segregated_higher.csv")
write_csv(overall_ranges,   "Output/Tabular Data/overall_ratio_ranges.csv")











############# SANITY CHECK

# --- Requirements: points_sf, metrics_points_df, soe_tbl already in memory ---

library(dplyr)
library(sf)
library(ggplot2)

# Stable point order (matches points_sf)
# ---------- Stable point key (ordered like points_sf) ----------
coords  <- sf::st_coordinates(points_sf)
pts_key <- tibble::tibble(
  x  = as.numeric(coords[,1]),
  y  = as.numeric(coords[,2])
) %>%
  mutate(xr = round(x, 3), yr = round(y, 3))  # small rounding for stable joins

# Round coords in the metrics table as well (non-destructive copy)
metrics_points_df_join <- metrics_points_df %>%
  mutate(
    xr = round(as.numeric(x), 3),
    yr = round(as.numeric(y), 3)
  )

# ---------- Pull a vector of per-point values, preserving points_sf order ----------
# metric_col must be one of: "prop_con" or "clumpy"
get_metric <- function(metric_col, target_extent, layout_name) {
  # Filter to the required layout + extent, then left-join to pts_key
  filt <- metrics_points_df_join %>%
    dplyr::filter(layout == layout_name, extent == target_extent) %>%
    dplyr::select(xr, yr, prop_con, clumpy)
  
  joined <- pts_key %>%
    dplyr::left_join(filt, by = c("xr", "yr"))
  
  # Return the requested metric as a numeric vector aligned to points_sf order
  out <- joined[[metric_col]]
  as.numeric(out)
}


# Helper: pull a vector of per-point metrics aligned to points_sf order
get_metric_vec <- function(layout_name = "Intimate", extent_m, metric = c("prop_con","clumpy")) {
  metric <- match.arg(metric)
  vals <- metrics_points_df %>%
    filter(layout == layout_name, extent == extent_m) %>%
    select(x, y, prop_con, clumpy) %>%
    right_join(pts_key, by = c("x","y")) %>%  # preserves points_sf order; may introduce NA
    pull({{ metric }})
  as.numeric(vals)
}

# --------------------------
# A) View for a chosen extent
# --------------------------
EXT <- 150  # change to 500 or 1000 as needed

prop_vec_I <- get_metric_vec("Intimate", extent_m = EXT, metric = "prop_con")
clmp_vec_I <- get_metric_vec("Intimate", extent_m = EXT, metric = "clumpy")

cat("\nINTIMATE @", EXT, "m — head(prop_con):\n"); print(head(prop_vec_I))
cat("\nINTIMATE @", EXT, "m — head(clumpy):\n");    print(head(clmp_vec_I))
cat("\nLengths (should equal nrow(points_sf)):\n",
    length(prop_vec_I), length(clmp_vec_I), "\n")

# Quick summaries & histos
summary_df <- tibble::tibble(prop_con = prop_vec_I, clumpy = clmp_vec_I)
print(summary_df %>% summarise(across(everything(),
                                      list(mean = ~mean(.x, na.rm=TRUE),
                                           sd   = ~sd(.x,   na.rm=TRUE),
                                           min  = ~min(.x,  na.rm=TRUE),
                                           max  = ~max(.x,  na.rm=TRUE)))))

g_prop <- ggplot(summary_df, aes(prop_con)) + geom_histogram(bins = 20, color="black") +
  labs(title = paste0("Intimate prop_con @ ", EXT, " m"))
g_clmp <- ggplot(summary_df, aes(clumpy))   + geom_histogram(bins = 20, color="black") +
  labs(title = paste0("Intimate CLUMPY @ ", EXT, " m"))
print(g_prop); print(g_clmp)

# -----------------------------------------------------------------
# B) View the *Intimate* vectors specifically at a species' SoE
# -----------------------------------------------------------------
# soe_tbl must contain: Species, var_type (prop_con/clumpy/age_mn), Feature (e.g., "prop_con_150_LULC")
# and parsed columns: soe_extent, soe_suffix (if not present, add them quickly)
if (!all(c("soe_extent","soe_suffix") %in% names(soe_tbl))) {
  soe_tbl <- soe_tbl %>%
    mutate(
      soe_extent = as.integer(sub(".*_(\\d+)_.*$", "\\1", Feature)),
      soe_suffix = sub("^.*_\\d+_([A-Za-z0-9]+)$", "\\1", Feature)
    )
}

view_intimate_for_species <- function(sp) {
  s <- soe_tbl %>% filter(Species == sp)
  if (nrow(s) == 0) stop("No SoE rows for species: ", sp)
  ext_prop <- s %>% filter(var_type == "prop_con") %>% slice(1) %>% pull(soe_extent)
  ext_clmp <- s %>% filter(var_type == "clumpy")   %>% slice(1) %>% pull(soe_extent)
  
  p_vec <- get_metric_vec("Intimate", extent_m = ext_prop, metric = "prop_con")
  c_vec <- get_metric_vec("Intimate", extent_m = ext_clmp, metric = "clumpy")
  
  cat("\n=== ", sp, " — INTIMATE @ SoE ===\n", sep = "")
  cat("prop_con extent:", ext_prop, "m; head:\n"); print(head(p_vec))
  cat("clumpy   extent:", ext_clmp, "m; head:\n"); print(head(c_vec))
  cat("Lengths:", length(p_vec), length(c_vec), " (nrow(points_sf) =", nrow(points_sf), ")\n")
  
  tibble::tibble(prop_con = p_vec, clumpy = c_vec)
}

# Example usage:
df_BTNW <- view_intimate_for_species("BTNW")
df_TEWA <- view_intimate_for_species("TEWA")
df_BBWA <- view_intimate_for_species("BBWA")

# ------------------------------------------------------
# C) Compare Intimate vs Segregated quickly at one extent
# ------------------------------------------------------
quick_compare <- function(ext) {
  tibble::tibble(
    layout = rep(c("Intimate","Segregated"), each = nrow(points_sf)),
    prop_con = c(get_metric_vec("Intimate", ext, "prop_con"),
                 get_metric_vec("Segregated", ext, "prop_con")),
    clumpy   = c(get_metric_vec("Intimate", ext, "clumpy"),
                 get_metric_vec("Segregated", ext, "clumpy"))
  ) %>%
    group_by(layout) %>%
    summarise(
      prop_con_mean = mean(prop_con, na.rm=TRUE),
      prop_con_sd   = sd(prop_con,   na.rm=TRUE),
      clumpy_mean   = mean(clumpy,   na.rm=TRUE),
      clumpy_sd     = sd(clumpy,     na.rm=TRUE),
      .groups = "drop"
    )
}
# quick_compare(150); quick_compare(500); quick_compare(1000)



















# Confirm CRS of points and that they sit in your township
sf::st_crs(points_sf)    # should be EPSG:3400
print(nrow(points_sf))   # number of point-count locations

# The exact X/Y used as predictors:
coords <- sf::st_coordinates(points_sf)
summary(coords[,1])  # x_AEP10TM range
summary(coords[,2])  # y_AEP10TM range
head(coords)         # first few XY pairs

# Quick sanity: do they fall within the township you defined?
all(sf::st_within(points_sf, township_sf, sparse = FALSE)[,1])




