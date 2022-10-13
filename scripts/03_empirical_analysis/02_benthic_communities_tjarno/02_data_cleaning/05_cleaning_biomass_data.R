#'
#' @title: Clean the biomass data
#' 
#' @description: Script to clean the biomass data and put it into the correct
#' format for applying the Isbell et al. (2018) partition.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

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

# add a cluster ID column
bio_dat$cluster_id <- substr(bio_dat$panel_id, 1, 1)

# add a buoy ID column
bio_dat$buoy_id <- gsub(pattern = "(-)(.*)", replacement = "", x = bio_dat$panel_id)

# check the I4 buoy
bio_dat %>%
  filter(buoy_id == "I4") %>%
  filter(panel_treatment == "mixture") %>%
  filter(time == "t2") %>%
  View()

# one of the panels was sampled too early
bio_dat <- 
  bio_dat %>%
  filter(laser_id != 519)

# reorder the columns and arrange the data
bio_dat <- 
  bio_dat %>%
  select(cluster_id, buoy_id, time,
         panel_treatment, OTU,
         dry_weight_g_cl) %>%
  arrange(cluster_id, buoy_id, time, panel_treatment, OTU)

# split the data into mixtures and monocultures
mix_dat <- 
  bio_dat %>%
  filter(panel_treatment == "mixture") %>%
  select(-panel_treatment)

mono_dat <- 
  bio_dat %>%
  filter(panel_treatment == "mono") %>%
  select(-panel_treatment)

# process the mixture data
OTU_names <- unique(mix_dat$OTU)
OTU_names <- OTU_names[!is.na(OTU_names)]

# why would there be an NA in the OTU column?
mix_dat %>%
  filter(is.na(OTU))

mix_dat <- 
  
  lapply( split(mix_dat, paste(mix_dat$buoy_id, mix_dat$time, sep = "_")), function(x) {
  
  # get any missing OTUs
  y <- OTU_names[!(OTU_names %in% x$OTU)]
  
  # create a data.frame with missing OTUs with zero abundance
  z <- 
    tibble(cluster_id = x$cluster_id[1],
           buoy_id = x$buoy_id[1],
           time = x$time[1],
           OTU = y,
           dry_weight_g_cl = rep(0, times = length(y)))
  
  # bind the missing OTU data to the full dataset
  df <- bind_rows(x, z)
  
  # return the data.frame
  return(df)
  
})

# bind back into a data.frame
mix_dat <- bind_rows(mix_dat)

# check the summary of the data
summary(mix_dat)

# check how complete are the data?
# probably need to remove I5 and J4 (not the end of the world)
mix_dat %>%
  group_by(cluster_id) %>%
  summarise(n = length(unique(buoy_id)))
  
mix_dat %>%
  group_by(cluster_id, buoy_id) %>%
  summarise(n = length(unique(time))) %>%
  View()

# rename the dry weight column
mix_dat <- 
  mix_dat %>%
  rename(Y = dry_weight_g_cl)

# process the monoculture data
mono_dat <- 
  
  lapply( split(mono_dat, paste(mono_dat$buoy_id, mono_dat$time, sep = "_")), function(x) {
    
    # get any missing OTUs
    y <- OTU_names[!(OTU_names %in% x$OTU)]
    
    # create a data.frame with missing OTUs with zero abundance
    z <- 
      tibble(cluster_id = x$cluster_id[1],
             buoy_id = x$buoy_id[1],
             time = x$time[1],
             OTU = y,
             dry_weight_g_cl = rep(NA, times = length(y)))
    
    # bind the missing OTU data to the full dataset
    df <- bind_rows(x, z)
    
    # return the data.frame
    return(df)
    
  })

# bind back into a data.frame
mono_dat <- bind_rows(mono_dat)

# check the summary of the data
summary(mono_dat)

# rename the dry weight g column to M
mono_dat <- 
  mono_dat %>%
  rename(M = dry_weight_g_cl)

# remove the incomplete data from the mix_dat data.frame
mix_dat <- 
  mix_dat %>%
  filter( !(buoy_id %in% c("I5", "J4")) )

# join the monoculture data
bio_dat <- left_join(mix_dat, mono_dat, by = c("cluster_id", "buoy_id", "time", "OTU"))

# write this into a .csv file
write_csv(x = bio_dat, file = here("data/benthic_communities_tjarno_data/data_clean/biomass_data.csv"))

### END
