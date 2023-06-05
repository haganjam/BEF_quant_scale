#'
#' @title: Extract the monoculture predictions
#' 
#' @description: This script uses the samples from the posterior distributions of
#' the best monoculture models of each species to fill in the missing monoculture data.
#' An object containing these monoculture predictions is then outputted as a .rds file.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(tidyr)

# load plotting theme
source("scripts/Function_plotting_theme.R")
source("scripts/03_empirical_analysis/helper_functions.R")

# load the analysis data
df_obs <- read_csv("data/case_study_2/data_clean/biomass_env_analysis_data.csv")

# load the model predictions
m0_pred <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/04_m0_predictions.rds")

df_imp <- vector("list", length = nrow(m0_pred))
for(i in 1:nrow(m0_pred)) {
  
  # add predictions from one sample from the posterior distribution
  df_obs[["M_imp"]] <- m0_pred[i,]
  
  # write to a list
  df_imp[[i]] <- df_obs
  
}

# for each sample, add the samples to the NA values
df_imp <- 
  
  lapply(df_imp, function(x) {
  
  # add the imputations to the missing data only
  x$M[which(is.na(x$M))] <- x$M_imp[which(is.na(x$M))]
  
  # remove the M_imp column
  x <- dplyr::select(x, -M_imp)
  
  return(x)
  
})

# remove the unnecessary columns in the data data.frame and put data in correct format
df_imp <- 
  
  lapply(df_imp, function(x) {
  
  x %>%
    select(cluster_id, buoy_id, time, OTU, M, Y) %>%
    rename(place = buoy_id, species = OTU) %>%
    mutate(place = as.integer(as.factor(place)),
           time = as.integer(as.factor(time)),
           species = as.integer(as.factor(species))) %>%
    group_by(cluster_id) %>%
    mutate(sample = as.integer(as.factor(paste0(place, time))) ) %>%
    ungroup() %>%
    select(cluster_id, sample, place, time, species, M, Y)
  
} )

# save the monoculture predictions as a .rds file
saveRDS(object = df_imp, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/05_complete_data.rds")

# generate 100 different Dirichlet distribution
start_RA <- sapply(1:100, function(x) gtools::rdirichlet(n = 1, rep(3, length(unique(df_imp[[1]]$species))) ) )

# save the Dirichlet distribution
saveRDS(object = start_RA, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/05_start_RYE.rds")

### END
