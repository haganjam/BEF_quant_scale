#'
#' @title: Clean the the cover-correction data
#' 
#' @description: Script to clean the cover-correction data. This cover-correction
#' data is used to check whether our weights taken at different time points
#' are representative samples of the same communities.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load required libraries
library(readr)
library(here)

# load the data
cov_dat <- read_csv(here("data/benthic_communities_tjarno_data/ResearchBox 843/Data/mixtures_coverage_comparison.csv"))

# remove the NAs
cov_dat <- cov_dat[complete.cases(cov_dat), ]

# arrange the data
cov_dat <- 
  cov_dat %>%
  arrange(buoy_id, sampling_time)

# replace the t1, t2, t3 etc. with just tfirst and tsecond
length(unique(cov_dat$buoy_id))*2 == nrow(cov_dat)
cov_dat$time <- rep(c("T1", "T2"), each = length(unique(cov_dat$buoy_id)))

# reorder the columns
cov_dat <- 
  cov_dat %>%
  select(buoy_id, sampling_time, date_ymd, time, laser_id, height_cm, width_cm,
         contains("percent"))

# write this into the cleaned data folder
write_csv(cov_dat, here("data/benthic_communities_tjarno_data/data_clean/mixture_coverage_data_clean.csv"))

### END
