#'
#' @title: Clean the biomass data
#' 
#' @description: Script to clean the biomass data and put it into the correct
#' format for applying the Isbell et al. (2018) partition.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# next steps: figure out how the monocultures and the mixtures match up

# load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(here)

# load the biomass data
bio_dat <- read_csv(here("data/benthic_communities_tjarno_data/ResearchBox 843/Data/biotic_data.csv"))

# check the variable structure
str(bio_dat)

# replace these excel errors with NA
bio_dat[bio_dat == "#VALUE!"] <- NA

# convert dry_weight_g and wet_weight_g to numeric variables
bio_dat$dry_weight_g <- as.numeric(bio_dat$dry_weight_g)
bio_dat$wet_weight_g <- as.numeric(bio_dat$wet_weight_g)
bio_dat$ash_weight_g <- as.numeric(bio_dat$ash_weight_g)

# if wet_weight_g is zero, then dry_weight_g should also be zero
bio_dat <- 
  bio_dat %>%
  mutate(dry_weight_g_cl = if_else(near(wet_weight_g, 0), 0, dry_weight_g),
         ash_weight_g_cl = if_else(near(wet_weight_g, 0), 0, ash_weight_g))

# are there any dry weights less than 0
bio_dat %>%
  filter(dry_weight_g_cl < 0) %>%
  View()

# make these zero
bio_dat <- 
  bio_dat %>%
  mutate(dry_weight_g_cl = if_else(dry_weight_g_cl < 0, 0, dry_weight_g_cl))

# select relevant rows
bio_dat <- 
  bio_dat %>%
  select(laser_id, panel_id, panel_treatment, OTU, 
         measurement_no, measurement_id, date,
         wet_weight_g, dry_weight_g_cl, ash_weight_g_cl, comment)

# check the summary of the different variables
summary(bio_dat)

# what do these NA's mean? There are no mixtures here so it's probably not that bad
bio_dat %>%
  filter(is.na(dry_weight_g_cl)) %>%
  View()

# remove cases with NAs for dry biomass
bio_dat <- 
  bio_dat %>%
  filter(!is.na(dry_weight_g_cl))

# convert the date column into an actual date
bio_dat$date <- gsub(pattern = "_", replacement = "-", bio_dat$date)
bio_dat$date <- lubridate::as_date(bio_dat$date)

# convert the date into a time point
bio_dat$time <- NA
bio_dat$time[bio_dat$date <= '2022-08-13'] <- 't1' 
bio_dat$time[(bio_dat$date >= '2022-08-13') & (bio_dat$date <= '2022-09-01')] <- 't2'
bio_dat$time[bio_dat$date >= '2022-09-01'] <- 't3'

# reorder the columns and arrange the data
bio_dat <- 
  bio_dat %>%
  select(time, date, laser_id, panel_id,
         measurement_no, measurement_id, 
         panel_treatment, OTU,
         dry_weight_g_cl, comment) %>%
  arrange(time, date, panel_treatment, panel_id)



