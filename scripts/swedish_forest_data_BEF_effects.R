
# Plymouth rock-pool data

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# set script to call partition functions from
source(here("scripts/isbell_2018_partition.R"))

# read in the data
swe_dat <- read_tsv( here("data/Swedish_forest_data.txt") )
summary(swe_dat)





