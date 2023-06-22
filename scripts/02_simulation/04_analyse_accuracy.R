#'
#' @title: Simulate metacommunities to test the analytical pipeline
#' 
#' @description: Analyse accuracy of assuming different RYEs
#' 
#' @details: This script analyses the simulated data to test how often the
#' estimates of the biodiversity that incorporate RYE capture the true, 
#' observed biodiversity effect in the simulated metacommunities.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(dplyr)
library(ggplot2)

# load the data
output_df <- readRDS("results/mc_sim_test.rds")

# get the subset of values that are affected by the RYEs
output_rye <- 
  output_df %>%
  filter( !(Beff %in% c( "LS", "LC", "TC", "NO" )) )

# get a table of the range of the different BEF effects
output_rye %>%
  group_by(Beff) %>%
  summarise(min_BEF = min(BEF_obs),
            max_BEF = max(BEF_obs))

# does the observed value lie in the 90% percentile interval
output_rye <- 
  output_rye %>%
  mutate(within = ifelse(BEF_obs <= PI95_high & BEF_obs >= PI95_low, 1, 0))

# check proportion within the interval
output_rye %>%
  group_by(Beff) %>%
  summarise(prop_within = sum(within)/n())

# calculate the percentage prediction error
output_rye <- 
  output_rye %>%
  mutate(abs_dev = abs(BEF_obs - mean_BEF) )
  
output_rye %>%
  group_by(Beff) %>%
  summarise(median_abs_dev = median(abs_dev),
            mean_abs_dev = mean(abs_dev),
            PI95_low = quantile(abs_dev, 0.05),
            PI95_high = quantile(abs_dev, 0.95))

### END
