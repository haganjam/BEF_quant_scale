#'
#' @title: Link cleaned biomass data and environmental data and calculate BEF effects
#' 
#' @description: In this script, we link the clean biomass data with the clean
#' environmental data and process it in order to calculate the BEF effects.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(here)
library(vegan)

# load plotting theme
source(here("scripts/Function_plotting_theme.R"))

# load the cleaned environmental data
env_dat <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/site_env_data.csv"))
head(env_dat)

# select the relevant variables
names(env_dat)

env_dat <- 
  env_dat %>%
  select(cluster_id, site_id, time, 
         panel_depth_m_measured, distance_between_panel_and_seabeed_m,
         temp_C_m, temp_C_cv, temp_C_max, temp_C_min,
         lux_m, lux_cv, lux_max, lux_min) %>%
  rename(buoy_id = site_id, depth_m = panel_depth_m_measured, seabed_dist_m = distance_between_panel_and_seabeed_m)

# load the cleaned biomass data
bio_dat <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/biomass_data.csv"))
head(bio_dat)

# join the environmental data to the biomass data
bio_env <- full_join(env_dat, bio_dat, by = c("cluster_id", "buoy_id", "time"))

# only keep complete-cases
bio_env <- bio_env[complete.cases(bio_env[,names(bio_env) != "M"]), ]

# check at the distribution of the samples
bio_env %>%
  group_by(cluster_id) %>%
  summarise(n = length(unique(buoy_id)))

incomplete <- 
  bio_env %>%
  group_by(cluster_id, buoy_id) %>%
  summarise(n = length(unique(time))) %>%
  filter(n < 3) %>%
  pull(buoy_id)
print(incomplete)

# remove the incomplete buoys
bio_env <- 
  bio_env %>%
  filter( !(buoy_id %in% incomplete) )

# rename the variables
bio_env

