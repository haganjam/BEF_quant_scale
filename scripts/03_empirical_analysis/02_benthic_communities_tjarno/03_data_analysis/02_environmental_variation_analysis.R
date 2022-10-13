#'
#' @title: Quantify the difference in environmental heterogeneity between clusters
#' 
#' @description: Script that uses the environmental variables that we collected and
#' quantifies the level of environmental heterogeneity in the different clusters,
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(vegan)

# load the cleaned environmental data
env_dat <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/site_env_data.csv"))

# choose the relevant columns
names(env_dat)

env_dat <- 
  env_dat %>%
  select(cluster_id, site_id, time, euclidian_dispersion_layer,
         panel_depth_m_measured, distance_between_panel_and_seabeed_m,
         temp_C_m, temp_C_cv, temp_C_max, temp_C_min, temp_field_C_m,
         lux_m, lux_cv, lux_max,
         sal_field_ppt_m, oxy_field_m, secchi_field_m,
         reduction_g_hour_m)

# summarise across the time points
env_dat_m <- 
  env_dat %>%
  group_by(cluster_id, site_id) %>%
  summarise(across(.cols = names(env_dat)[-c(1:3)], ~mean(., na.rm = TRUE)),
            .groups = "drop")

# generate a Euclidean distance matrix
env_d <- dist(apply(env_dat_m[,-c(1,2)], 2, scale), method = "euclidean")

# multivariate dispersion
env_d <- betadisper(d = env_d, group = env_dat_m$cluster_id)

# add the distance to centroid to the data
env_dat_m$distance_centroid <- env_d$distances

# does this correlate with the originally implemented centroid distances
env_dispersion <- 
  env_dat_m %>%
  group_by(cluster_id) %>%
  summarise(GIS_dispersion = mean(euclidian_dispersion_layer),
            field_dispersion = mean(distance_centroid))

# plot the relationship
plot(env_dispersion$GIS_dispersion, env_dispersion$field_dispersion)
abline(a = 0, b = 1)

# examine the spearman correlation which tests for monotonic relationships
cor(env_dispersion$GIS_dispersion, env_dispersion$field_dispersion, method = "spearman")

