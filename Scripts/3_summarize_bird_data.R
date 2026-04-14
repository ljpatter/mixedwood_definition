# ---
# title: "Summarize bird data"
# author: "Brendan Casey"
# created: "2024-01-13"
# description: "Summarize bird data from WildTrax."
# ---

## Load packages ----
library(tidyverse)
library(vegan)
library(dplyr)
library(wildrtrax) # for species table

# Authenticate with WildTrax using environment variables for credentials
config <- "Scripts/login.R"
source(config)
wt_auth()

## Import data ----
load("Output/R Data/wildtrax_cleaned_2025-09-15.rData")

## ////////////////////////////////////////////////////////////////
## Filter species ---
## Filter flyover species and water birds
aves <- wt_get_species() %>%
  dplyr::filter(
    species_class == "AVES",
    species_order != "Podicipediformes",
    species_order != "Pelecaniformes",
    species_order != "Anseriformes",
    species_order != "Gaviiformes"
  )

bd1 <- wildtrax_cleaned %>% semi_join(aves)

# view species list
species_list <- bd1 %>%
  dplyr::select(species_code) %>%
  distinct() %>%
  left_join(aves)

## ////////////////////////////////////////////////////////////////
## Prepare data ---
## Convert table to wide format
bd2 <- bd1 %>%
  dplyr::select(-c("detection_time", "detection_distance", "vocalization")) %>% # control for individuals counted multiple times in same survey
  distinct() %>%
  pivot_wider(
    names_from  = species_code,
    values_from = individual_count,
    values_fn   = list(individual_count = sum)
  )

## Subset columns
colnms <- c(
  "organization", "location_buffer_m",
  "lon", "lat",
  "tz", "sunrise", "task_method",
  "max_duration", "n_surveys", "individual_order"
)

bd3 <- bd2 %>%
  dplyr::select(-c(all_of(colnms)))

## Select one point count per year per site (reproducible)
set.seed(1234)

select_random_point_count <- function(df) {
  # Check if the required columns exist
  if (!all(c("location", "year") %in% names(df))) {
    stop("DataFrame must contain 'location' and 'year' columns.")
  }
  df %>%
    dplyr::group_by(location, year) %>%
    dplyr::slice_sample(n = 1) %>%
    dplyr::ungroup()
}

# Filter
bd4 <- select_random_point_count(bd3)

# Replace NAs in spp columns w/ 0s
replace_na_with_zero <- function(df) {
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], function(x) { x[is.na(x)] <- 0; x })
  return(df)
}

bd5 <- replace_na_with_zero(bd4)

## FINAL CHECK: ensure no NA values in species columns
species_cols <- setdiff(names(bd5), c("location", "year"))
if (anyNA(bd5[, species_cols, drop = FALSE])) {
  stop("NA values detected in species columns.")
}

## Save
write.csv(bd5, "Output/Tabular Data/bird_data_2.csv", row.names = FALSE)

