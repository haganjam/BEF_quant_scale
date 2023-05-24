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
output_df <- readRDS("results/MC_sims_post.rds")

# get the subset of values that are affected by the RYEs
output_rye <- 
  output_df %>%
  filter(Beff %in% c("AS", "IT", "LS", "NBE", "TS"))

# get a table of the range of the different BEF effects
output_rye %>%
  group_by(Beff) %>%
  summarise(min_BEF = min(BEF_obs),
            max_BEF = max(BEF_obs))

# does the observed value lie in the 90% percentile interval
output_rye %>%
  mutate(within = ifelse(BEF_obs <= PI90_high & BEF_obs >= PI90_low, 1, 0)) %>%
  group_by(Beff) %>%
  summarise(prop_within = sum(within)/100)


