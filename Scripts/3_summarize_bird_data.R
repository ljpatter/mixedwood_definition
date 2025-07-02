# ---
# title: "Calculate taxonomic diversity"
# author: "Brendan Casey"
# created: "2024-01-13"
# description: "Calculate taxonomic diversity from WildTrax detection data."
# ---



##Load packages----
library(tidyverse)
library(vegan)
library(dplyr)
library(wildrtrax)#for species table

# Authenticate with WildTrax using environment variables for credentials
Sys.setenv(WT_USERNAME = "ljpatter", WT_PASSWORD = "Kingedwardpark13")
wt_auth()

##Import data----
load("Output/wildtrax_cleaned_2025-03-12.rData")

##////////////////////////////////////////////////////////////////
# Filter species---
## Filter flyover species and water birds
aves<-wt_get_species()%>%
  dplyr::filter(species_class=="AVES",
                species_order!="Podicipediformes",
                species_order!="Pelecaniformes",
                species_order!="Anseriformes",
                species_order!="Gaviiformes")  

bd1<-wildtrax_cleaned%>%semi_join(aves)

#view species list
species_list<-bd1 %>%
  dplyr::select(species_code) %>%
  distinct() %>%
  left_join(aves)

##////////////////////////////////////////////////////////////////
# Prepare data----
## Convert table to wide format----
bd2<-bd1%>%
  dplyr::select(-c("detection_time", "detection_distance", "vocalization"))%>% #control for individuals that were counted multiple times in the same survey
  distinct()%>%
  pivot_wider(names_from = species_code, values_from = individual_count, values_fn = list(individual_count = sum))

## Subset columns----
colnms<- c("organization", "location_buffer_m", 
           "lon", "lat", 
           "tz", "sunrise", "task_method", 
           "max_duration", "n_surveys", "individual_order") 

bd3<-bd2%>%
  dplyr::select(-c(all_of(colnms)))

### Select one point count per year per sites

select_random_point_count <- function(df) {
  # Check if the required columns exist
  if (!all(c("location", "year") %in% names(df))) {
    stop("DataFrame must contain 'location' and 'year' columns.")
  }
  
  # Group by location and year, then sample one row from each group
  selected_rows <- df %>%
    group_by(location, year) %>%
    slice_sample(n = 1) %>%
    ungroup() # Remove grouping
  
  return(selected_rows)
}

# Filter
bd4 <- select_random_point_count(bd3)

### Replace any NA values with 0s

# Assuming your data frame is named 'bd3'

replace_na_with_zero <- function(df) {
  # Replace NA with 0 in all columns
  df[is.na(df)] <- 0
  return(df)
}

# Example Usage (replace 'bd3' with your actual data frame name)
bd5 <- replace_na_with_zero(bd4)

#Save
write.csv(bd5, "Output/Tabular Data/bird_data_2.csv")
