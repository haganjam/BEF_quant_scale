#create lm to correct light sensor with light sensor data

library(readr)
library(dplyr)

#Load cleaning functions
source(here("scripts/03_empirical_analysis/02_benthic_communities_tjarno/02_data_cleaning/cleaning_functions.r"))

check.dirs()
files=get.data.filenames("light_sensor_correction")




rm(list = ls())
