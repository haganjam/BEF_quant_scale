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



