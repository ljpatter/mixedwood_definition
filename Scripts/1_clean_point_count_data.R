# ---
# title: "Clean bird data"
# author: "Brendan Casey and Leonard Patterson"
# created: "2025-04-11"
# description: "Code to clean data brought in from wildtrax. 
# Adapted some of the script contained in https://github.com/
# borealbirds/QPADestimation/blob/main/script/02.CleanData.R"
# ---

#Setup ----

##Load packages----
library(tidyverse)
library(wildrtrax)
library(lubridate)
library(suncalc) #sunrise time retrieval
library(lutz) #get timezone
library(expss) #for calulating TMTT values
library(sf) #spatial projections and getting new coordinates


##Import data----
load("Output/R Data/wildtrax_raw_aru_2025-01-18.rData")
load("Output/R Data/wildtrax_raw_pc_2025-01-18.rData")

##////////////////////////////////////////////////////////////////
# Login to WildTrax----
# NOTE: Edit the 'loginexample.R' script to include your WildTrax 
# login details and rename to 'login.R'. 
# DO NOT PUSH YOUR LOGIN TO GITHUB
#config <- "1_code/r_scripts/login.R"
#source(config)
#wt_auth()


# Authenticate with WildTrax using environment variables for credentials
Sys.setenv(WT_USERNAME = "ljpatter", WT_PASSWORD = "Kingedwardpark13")
wt_auth()

##////////////////////////////////////////////////////////////////
# Clean methods columns----

## define columns of interest
colnms<- c("organization","project","project_id","survey_type", 
           "location", "location_buffer_m", "longitude","latitude",
           "date","survey_duration_method","max_duration",
           "survey_distance_method","detection_time",
           "detection_distance","task_method","species_code", 
           "individual_order","individual_count", "vocalization")


## ARU----
wildtrax_raw_aru_1<-wildtrax_raw_aru%>%
  #filter by task_duration values that are multiples of 60 %>%
  filter(task_duration %% 60 == 0 & task_duration!=0)%>% 
  rename(date=recording_date_time, max_duration=task_duration)%>%
  mutate(minutes=max_duration/60)%>%
  mutate(survey_duration_method = case_when(minutes==1 ~ "0-1min",
                                            minutes==2 ~ "0-1-2min",
                                            minutes==3 ~ "0-1-2-3min",
                                            minutes==4 ~ "0-1-2-3-4min",
                                            minutes==5 ~ "0-1-2-3-4-5min",
                                            minutes==6 ~ "0-1-2-3-4-5-6min",
                                            minutes==7 ~ "0-1-2-3-4-5-6-7min",
                                            minutes==8 ~ "0-1-2-3-4-5-6-7-8min",
                                            minutes==9 ~ "0-1-2-3-4-5-6-7-8-9min",
                                            minutes==10 ~ "0-1-2-3-4-5-6-7-8-9-10min"))%>%
  mutate(survey_distance_method="0m-INF-ARU")%>%
  mutate(detection_distance="0m-INF-ARU")%>%
  mutate(detection_time=as.character(detection_time))%>%
  dplyr::select(any_of(colnms))

## PC----
wildtrax_raw_pc_1<-wildtrax_raw_pc%>%
  dplyr::select(-c(observer))%>%
  mutate(max_duration=as.numeric(
    str_extract(survey_duration_method, 
                "(?<=-)(\\d+)(?=min)"))*60)%>%
  rename(date=survey_date)%>%
  dplyr::select(any_of(colnms))

## Combine PC and ARU dataframes----
wildtrax_all<-dplyr::bind_rows(wildtrax_raw_pc_1, wildtrax_raw_aru_1)



##////////////////////////////////////////////////////////////////
# Filter data---- 
w<-wildtrax_all

## Remove non bird species----
aves<-wt_get_species()%>%
  filter(species_class=="AVES")

w2<-w%>%semi_join(aves)

## Remove single species projects----
w3 <- w2%>%group_by(project_id) %>%
  #single species projects only have one species or none.
  filter(n_distinct(species_code) > 3)%>% 
  ungroup()%>%
  distinct()

## Exclude other projects that shouldn't be used----
# Exclude BU training, ABMI, & all 
# "DO NOT USE" projects listed in the below .csv
instructions <- read.csv(file.path(
  paste0("Input/Tabular Data/projectInstructions.csv")))

w4 <- w3 %>% 
  anti_join(instructions) %>%
  dplyr::filter(organization!="BU-TRAINING")%>%
  dplyr::filter(organization!="ABMI")


## Remove missing data, etc.----
#Remove surveys with no location
#Standardize sig figs for location to remove duplicates
#Remove outliers for day of year (use 99% quantile)
#Take out BBS because isn't useful for removal or distance sampling
#Remove none tag methods
#Remove surveys with unknown survey time - specific coding
#Remove midnight PC surveys
#Remove locations that are within a buffer
w5 <- w4 %>%
  dplyr::filter(!is.na(date),
                project!="BAM-BBS") %>% 
  mutate(location_buffer_m = ifelse(is.na(location_buffer_m), 0, 
                                    location_buffer_m),
         ordinalDay = yday(date)) %>% 
  unique() %>% 
  dplyr::filter(!is.na(latitude),
                latitude > 0, 
                longitude < 0,
                !is.na(date),
                year(date) > 1900,
                ordinalDay > quantile(ordinalDay, 0.005),
                ordinalDay < quantile(ordinalDay, 0.995),
                !str_sub(date, -8, -1) %in% c("00:00:01", "00:01:01"),
                !(hour(date)==0 & survey_type=="PC"),
                location_buffer_m==0)




##////////////////////////////////////////////////////////////////

#Convert TMTT counts to numeric----
## Some of the Wildtrax abundance data is labeled TMTT (too many to count). 
## For each species, I'll convert this to a number by taking the average number of 
## individuals in a group of 3 or more and using that average in place of TMTT.


## Seperate numerical and "too many to count" abundance values----
w6<-w5%>%
  filter(individual_count!="CI 3"&individual_count!="CI 1"&individual_count!="CI 2"&individual_count!="N/A"&individual_count!="NA")
#add new columns for abundance data (convert words to numbers)
w6$ABUND_TMTT[w6$individual_count=="TMTT"]<-"TMTT"

#Seperate dataframe by species with a numeric ABUND vs too many to count
w_TMTT<-filter(w6, individual_count=="TMTT")
# w_TMTC<- distinct(w_TMTT)

#filter to only those where there are numeric ABUND estimates
w_ABUND<-w6%>%
  filter(individual_count!="TMTT")%>%
  mutate(individual_count=as.numeric(individual_count))%>%
  dplyr::select(-c(ABUND_TMTT))

## Convert bird data to wide format----
# ## using the helper function in the wildrtrax package,
# wt_wide<- wt_make_wide(data = wildtrax_raw,sensor = "PC")
# wt_wide_2<-wt_wide%>%
#   dplyr::select(-c(1:14))

wt_wide<-w_ABUND%>%
  pivot_wider(names_from = species_code, values_from = individual_count, values_fn = list(individual_count = sum))
wt_wide_2<-wt_wide%>%
  dplyr::select(-c(1:17))
#1159504

## Calculate mean TMTT----

bmean<-round(mean_col_if(gt(3), wt_wide_2), digits = 0)
bm<-as.data.frame(bmean)
bm <- cbind(rownames(bm), data.frame(bm, row.names=NULL))
colnames(bm)<- c("species_code","ABUND_TMTT_2")

#IF a species hasn't been detected in groups of 4 or more, convert TMTT to 4
bm[is.na(bm)]<-4

## Create a new dataframe with estimated TMTT values----
w7<-w_TMTT%>%
  #join TMTT abundance values with the subsetted TMTT Wiltrax data
  left_join(bm)%>%
  #add TMTT abundance estimates to individual count
  mutate(individual_count=(ABUND_TMTT_2))%>%
  dplyr::select(-c(ABUND_TMTT, ABUND_TMTT_2))%>%
  ###### Bind back together the TMTT and numerical abundance subsets of the Wildtrax data {-}
  bind_rows(w_ABUND)

##////////////////////////////////////////////////////////////////
# Temporal variables----
## Add month and year columns----
w8<-w7%>%
  mutate(year = year(date),
         month=month(date),
         ordinalDay = yday(date),
         start_time = hour(date) + minute(date)/60)

## Get local timezone----
w9 <- w8 %>%
  mutate(date = ymd(str_sub(date, 1, 10))) %>%
  mutate(tz=tz_lookup_coords(latitude, longitude, method="accurate"))%>%
  rename(lat=latitude, lon=longitude)

## Calculate time since sunrise----
#Loop through timezones
tzs <- unique(w9$tz)
tzs<- Filter(function(x) !is.na(x), tzs) #remove NA values

sun.list <- list()
for(i in 1:length(tzs)){
  
  all.i <- w9 %>% 
    dplyr::filter(tz==tzs[i])
  
  all.i$sunrise <- getSunlightTimes(data=all.i, keep="sunrise", tz=tzs[i])$sunrise
  all.i$hssr <- all.i$start_time-(hour(all.i$sunrise) + minute(all.i$sunrise)/60)
  
  sun.list[[i]] <- all.i
}

sun <- do.call(rbind, sun.list)

ggplot(sun) +
  geom_histogram(aes(x=hssr, fill=survey_type))

w10<-sun

##////////////////////////////////////////////////////////////////
# Calculate survey effort----  
w11<-w10%>% 
  dplyr::select(project_id, survey_type, location, survey_duration_method, survey_distance_method, year, date, max_duration)%>%
  distinct()%>%
  group_by(project_id, survey_type, location,survey_duration_method, survey_distance_method, year)%>%
  mutate(n_surveys=n())%>%
  ungroup()%>%
  mutate(survey_effort=n_surveys*max_duration/60)%>%
  full_join(w10)

# hist(w11$n_surveys)
# summary(w11$n_surveys)
# hist(w11$survey_effort)
# summary(w11$survey_effort)


# Filter point counts to only those that occurred in June

w12 <- w11 %>% filter(month == 6)

# Filter point counts to only include those occuring half an hour before sunrise to 5 hours after sunrise

w13 <- w12 %>% filter(hssr >= -0.5 & hssr <= 5)

### Extract max distance band 

# Function to extract the max distance band
extract_max_dist <- function(method) {
  if (is.na(method) || method == "") return(NA)  # Handle missing values
  
  # Split the string by "-" to extract individual distance bands
  bands <- str_split(method, "-", simplify = TRUE) %>% as.vector()
  
  # Remove any non-distance values (e.g., "ARU")
  bands <- bands[!bands %in% c("ARU")]
  
  # Remove "m" from numeric distances
  cleaned_bands <- str_replace_all(bands, "m", "")  # Remove "m"
  
  # Separate INF from numeric distances
  has_inf <- "INF" %in% cleaned_bands  # Check if INF is present
  numeric_bands <- suppressWarnings(as.numeric(cleaned_bands))  # Convert valid numbers
  
  # Remove any NAs that result from failed conversions (e.g., empty strings)
  numeric_bands <- numeric_bands[!is.na(numeric_bands)]
  
  # Get the max numeric value if present, otherwise return INF
  if (has_inf) {
    return("INF")  # If INF was in the original string, return it
  } else if (length(numeric_bands) > 0) {
    return(as.character(max(numeric_bands)))  # Return max numeric distance
  } else {
    return(NA)  # If nothing valid remains, return NA
  }
}

# Apply function to create new column
w14 <- w13 %>%
  mutate(max_dist_band = sapply(survey_distance_method, extract_max_dist))

# Remove all sites where distance bands is less than 100 m
w15 <- w14 %>%
  filter(as.numeric(max_dist_band) >= 100 | max_dist_band == "INF")

# Remove all entries where individual_count is 0 (i.e., PC where 0 had to be entered for non-detection)
w16 <- w15 %>%
  mutate(.ic = suppressWarnings(as.numeric(as.character(individual_count)))) %>%
  filter(survey_type != "PC" | is.na(.ic) | .ic != 0) %>%
  select(-.ic)

### Add coordinates in UTM

# Load in AB boundary shp file to extract PCS

AB_boundary_10TM <- st_read("Input/Spatial Data/Alberta/AB_boundary.shp")

# Set target CRS

target_crs <- st_crs(AB_boundary_10TM)

##////////////////////////////////////////////////////////////////
# Add XY coordinates that are in meters----
w17<-w16%>%
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE)%>%
  st_transform(crs=target_crs)%>%#transform to projection that uses meters
  dplyr::mutate(x_AEP10TM = sf::st_coordinates(.)[,1], #extract coordinates and add as new columns
                y_AEP10TM = sf::st_coordinates(.)[,2])%>%
  as.data.frame()%>%
  # full_join(w11)%>% #bring back in original lat and lon
  dplyr::select(-geometry)   

##////////////////////////////////////////////////////////////////
# Select and order columns----
## define columns of interest
colnms<- c("organization", "project", "project_id", "survey_type", "location", "location_buffer_m", 
           "lon", "lat", "x_AEP10TM", "y_AEP10TM", "date", "ordinalDay", "year", "month", "start_time", 
           "tz", "sunrise", "hssr", "task_method", "survey_distance_method", "max_dist_band", "survey_duration_method",
           "max_duration", "n_surveys", "survey_effort", "detection_time", "detection_distance", 
           "species_code", "individual_order", "individual_count","vocalization") 

wildtrax_cleaned<-w17%>%
  dplyr::select(all_of(colnms))

# Save
save(wildtrax_cleaned, file=paste0("Output/R Data/wildtrax_cleaned_", Sys.Date(), ".rData"))


# Clear environment
rm(list = ls())  # Removes all objects from the environment





