# ---
# title: "Get WT reports"
# author: "Brendan Casey and Leonard Patterson"
# created: "2025-04-11"
# description: This code downloads point count data from WildTrax and combines 
# data into a single dataframe
# ---

##Load packages----
# install.packages("remotes")
#remotes::install_github("ABbiodiversity/wildRtrax")
# remotes::install_github("ABbiodiversity/wildRtrax@development")
library(wildrtrax)
library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Authenticate with WildTrax using environment variables for credentials
config <- "Scripts/login.R"
source(config)
wt_auth()

# Initialize a vector to store project IDs that cause errors
error_projects <- c()

# ARU only
my_report_aru <- wt_get_download_summary(sensor_id = "PC") %>%
  as_tibble() %>%
  filter(sensor == "ARU") %>%
  mutate(data = purrr::map(
    .x = project_id,
    .f = ~tryCatch(
      {
        # Attempt to download the report
        result <- wt_download_report(
          project_id = .x,
          sensor_id = "ARU",
          weather_cols = FALSE,
          reports = "main"
        )
        # Ensure the result is a tibble
        if (!inherits(result, "tbl_df")) {
          result <- tibble::tibble()
        }
        
        # Explicitly convert 'location' column to character if it exists
        if ("location" %in% colnames(result)) {
          result$location <- as.character(result$location)
        }
        
        result
      },
      error = function(e) {
        error_projects <<- c(error_projects, .x)
        NULL
      }
    )
  )) %>%
  filter(!sapply(data, is.null)) %>%
  select(project, data) %>%
  unnest(cols = data, keep_empty = TRUE) %>%
  mutate(survey_type = "ARU")


# Print the IDs of projects that caused errors
if (length(error_projects) > 0) {
  message("The following project IDs caused errors:")
  print(error_projects)
} else {
  message("No errors encountered.")
}





###PC only---

my_report_pc <- wt_get_download_summary(sensor_id = "PC") %>%
  tibble::as_tibble() %>%
  filter(sensor == "PC") %>%
  dplyr::mutate(data = purrr::map(.x = project_id, 
                                  .f = ~wt_download_report(project_id = .x, 
                                                           sensor_id = "PC", 
                                                           weather_cols = FALSE, 
                                                           reports = "main"))) %>%
  dplyr::select(c(data)) %>%
  mutate(class = sapply(data, function(df) class(df$individual_count))) %>%
  filter(class != "logical") %>% 
  filter(class != "NULL") %>%
  mutate(data = lapply(data, function(df) {
    df %>%
      mutate(individual_count = as.character(individual_count),
             observer = as.character(observer)) # Convert observer to character
  })) %>%
  dplyr::select(-class) %>%
  unnest(col = data) %>%
  mutate(survey_type = "PC")

 # Save----
wildtrax_raw_pc<-my_report_pc
save(wildtrax_raw_pc, file=paste0("Output/R Data/wildtrax_raw_pc_", Sys.Date(), ".rData"))

wildtrax_raw_aru<-my_report_aru
save(wildtrax_raw_aru, file=paste0("Output/R Data/wildtrax_raw_aru_", Sys.Date(), ".rData"))
