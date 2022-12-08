#'
#' @title: Merge biomass and environmental data
#' 
#' @description: In this script, we link the clean biomass data with the clean
#' environmental data and generate PC axes that we will use to predict the
#' missing monocultures.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(here)

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

# check the abundances of the different species in mixture
bio_env_RA <- 
  bio_env %>%
  select(cluster_id, buoy_id, time, OTU, Y) %>%
  group_by(cluster_id, buoy_id, time) %>%
  mutate(RA = Y/sum(Y)) %>%
  ungroup()

ggplot(data = bio_env_RA, 
       mapping = aes(x = RA)) +
  geom_histogram() +
  facet_wrap(~OTU, scales = "free") +
  theme_bw()

bio_env_RA %>%
  filter(OTU %in% c("Bumpi", "Asci", "Hydro") , RA > 0.1) %>%
  View()

bio_env_RA %>%
  group_by(OTU) %>%
  summarise(RA_mean = mean(RA, na.rm = TRUE),
            RA_sd = sd(RA, na.rm = TRUE),
            RA_min = min(RA, na.rm = TRUE),
            RA_max = max(RA, na.rm = TRUE),
            n_0 = sum(RA > 0, na.rm = TRUE))

# can we remove Asci?
# yes, only 7 cases where it is non-zero and has a mean relative abundance of < 0.05%

bio_env <- 
  bio_env %>%
  filter(OTU != "Asci")

# check the distribution of the monoculture data
bio_env %>%
  ggplot(data = .,
         mapping = aes(x = M)) +
  geom_histogram() +
  facet_wrap(~OTU, scales = "free") +
  theme_bw()

# how are the variables correlated?
names(bio_env)
pairs(bio_env[, c(4, 5, 6, 10, 15, 16)])

ggplot(data = bio_env,
       mapping = aes(x = Y, y = M, colour = OTU)) +
  geom_point() +
  facet_wrap(~OTU, scales = "free") +
  theme_bw()

# run a PCA on the environmental variables
pca.x <- prcomp(reformulate(names(bio_env)[-c(1:3, 13:16)]), 
                data = bio_env, scale = TRUE)

# check the PCA summary
summary(pca.x)

# bind the PCA to identifier variables
pca_env <- bind_cols(bio_env[,c(1:3)], as_tibble(pca.x$x)[,c(1:2)] )

# bind the monoculture and mixture data to the pca_env data
pca_env <- bind_cols(pca_env, bio_env[,c(14:16)])

# calculate the total number of measurements left over after cleaning the data
pca_env %>%
  pull(buoy_id) %>%
  unique() %>%
  length()
  
pca_env %>%
  group_by(cluster_id, buoy_id, time) %>%
  summarise(n = n()) %>%
  nrow()

pca_env %>%
  group_by(cluster_id, buoy_id, time) %>%
  summarise(n = sum(!is.na(M))) %>%
  pull(n) %>%
  sum()

# write out the pca_env_data
write_csv(x = pca_env, here("data/benthic_communities_tjarno_data/data_clean/biomass_env_analysis_data.csv"))

### END
