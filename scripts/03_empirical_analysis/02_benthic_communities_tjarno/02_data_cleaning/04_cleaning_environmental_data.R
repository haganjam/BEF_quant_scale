#'
#' @title: Clean environmental data
#' 
#' @description: Cleaning of abiotic field measurements
#' In the end a combined dataset for abiotic measures, flow speed, 
#' light temperature, depth will be created
#' 
#' @authors: Lara Martins, Benedikt Schrofner-Brunner
#' 

# load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)
library(stringr)
library(lubridate)

# load cleaning functions
source(here("scripts/03_empirical_analysis/02_benthic_communities_tjarno/02_data_cleaning/cleaning_functions.r"))

# check if the correct directories exist
check.dirs() 

# 1. load the abiotic measurement data
abiotic_measurements <- 
  read_csv(here("data/benthic_communities_tjarno_data/ResearchBox 843/Data/abiotic_field_measurements.csv"), 
           col_types = cols(Date = col_date(format = "%d/%m/%Y"), 
                            Depth = col_number(), Temp. = col_number(), 
                            Salinity = col_number(), Oxygen = col_number(), 
                            Secchi = col_number())
           )

# rename the columns
names(abiotic_measurements) = c("date","cluster_id","site_id","panel_depth",
                                "temperature_C_field", "salinity_ppt_field",
                                "oxygen","secchi_depth_m","weather_comment","comment")

# write t1, t2 and t3 in the table
abiotic_measurements$time <- "" 

# add t1, t2 and t3 labels
abiotic_measurements$time[abiotic_measurements$date <= '2022-08-16'] <- 't1'
abiotic_measurements$time[(abiotic_measurements$date >= '2022-08-17') & (abiotic_measurements$date <= '2022-08-31')] <- 't2'
abiotic_measurements$time[abiotic_measurements$date >= '2022-09-01'] <- 't3'

# view the abiotic measurements data
print(abiotic_measurements)

# summarise the abiotic measurements
abiotic_measurements <- 
  abiotic_measurements %>%
  rename(depth = panel_depth) %>%
  group_by(site_id, depth, time) %>%
  summarise(temp_field_C_m = mean(temperature_C_field, na.rm = TRUE),
            sal_field_ppt_m = mean(salinity_ppt_field, na.rm = TRUE),
            oxy_field_m = mean(oxygen, na.rm = TRUE),
            secchi_field_m = mean(secchi_depth_m, na.rm = TRUE),
            .groups = "drop")

# 2. load the cleaned logger data

# load the abiotic measurement data
hobo_dat <- 
  read_csv(here("data/benthic_communities_tjarno_data/data_clean/light_temperature_logger_clean.csv"))

# reorganise the columns
hobo_dat <- 
  hobo_dat %>%
  select(site, id, time, date_time, depth, temp_C, lux.corrected) %>%
  rename(site_id = site)

# summarise the temperature and light data
hobo_dat <- 
  hobo_dat %>%
  group_by(site_id, depth, time) %>%
  summarise(temp_C_m = mean(temp_C, na.rm = TRUE),
            temp_C_cv = (sd(temp_C, na.rm = TRUE)/mean(temp_C, na.rm = TRUE))*100,
            temp_C_max = min(temp_C, na.rm = TRUE),
            temp_C_min = max(temp_C, na.rm = TRUE),
            lux_m = mean(lux.corrected, na.rm = TRUE),
            lux_cv = (sd(lux.corrected, na.rm = TRUE)/mean(temp_C, na.rm = TRUE))*100,
            lux_max = max(lux.corrected, na.rm = TRUE),
            lux_min = min(lux.corrected, na.rm = TRUE),
            .groups = "drop")

# 3. load the gypsum data
gypsum <- 
  read_csv(here("data/benthic_communities_tjarno_data/ResearchBox 843/Data/water_velocity_gypsum.csv"),
           col_types = cols(cluster = col_character(),
                            site = col_character(),
                            date_deploy = col_character(),
                            time_deploy = col_character(),
                            date_retrieve = col_character(),
                            time_retrieve = col_character()
                            ))

# add the year to date deploy and retrieve columns
gypsum$date_deploy <- as.character(dmy(paste0(gypsum$date_deploy, "-2022")))
gypsum$date_retrieve <- as.character(dmy(paste0(gypsum$date_retrieve, "-2022")))

# remove the NAs in the gypsum data
gypsum <- gypsum[complete.cases(gypsum$weight_after), ]

# get date time columns
gypsum <- 
  gypsum  %>%
  mutate(date_time_deploy = ymd_hm(paste(date_deploy, time_deploy, sep = " "), tz = "Europe/Stockholm"),
         date_time_retrieve = ymd_hm(paste(date_retrieve, time_retrieve, sep = " "), tz = "Europe/Stockholm"))

# calculate the duration that the gympsum blocks were out for
gypsum <- 
  gypsum %>%
  mutate(duration_hours = as.numeric(date_time_retrieve - date_time_deploy))

# get relevant columns
names(gypsum)

gypsum <- 
  gypsum %>%
  select(site, cluster, gypsum_id, panel_depth_m, date_deploy, weight_before, weight_difference, duration_hours) %>%
  rename(site_id = site)

# calculate the reduction per hour
gypsum <- 
  gypsum %>%
  mutate(reduction_g_hour = (weight_difference/weight_before)/duration_hours)

# add time variables
gypsum$time <- ""

# convert the date deploy column to simply a date
gypsum$date_deploy <- as.Date(gypsum$date_deploy)

# add labels based on the time points
gypsum$time[gypsum$date_deploy <='2022-07-21'] <- 't1' #chose the column, search for the specific date, write the "names" of that
gypsum$time[(gypsum$date_deploy>='2022-08-02') & (gypsum$date_deploy <= '2022-08-29')] <- 't2'
gypsum$time[gypsum$date_deploy>='2022-08-29'] <- 't3'

# summarise the gypsum data
gypsum <- 
  gypsum %>%
  rename(depth = panel_depth_m) %>%
  group_by(site_id, depth, time) %>%
  summarise(reduction_g_hour_m = mean(reduction_g_hour, na.rm = TRUE),
            .groups = "drop")

# 4. load the site data
site_dat <- read_csv(file = here("data/benthic_communities_tjarno_data/ResearchBox 843/Data/site_data.csv"))

# rename the depth variable
site_dat <- 
  site_dat %>%
  rename(depth = panel_depth_m)

# 5. combine the site data with all the other abiotic measurements

# merge the datasets
env_dat <- 
  left_join( 
  left_join(hobo_dat, abiotic_measurements, by = c("site_id", "depth", "time")),
  gypsum,
  by = c("site_id", "depth", "time"))

# join the env_dat to the site data
site_dat <- right_join(site_dat, env_dat, by = c("site_id", "depth"))

# export this file into the data_clean folder
write_csv(site_dat, here("data/benthic_communities_tjarno_data/data_clean/site_env_data.csv"))

### END
