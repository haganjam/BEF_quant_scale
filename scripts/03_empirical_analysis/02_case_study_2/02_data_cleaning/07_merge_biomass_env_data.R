#'
#' @title: Merge biomass and environmental data
#' 
#' @description: In this script, we link the clean biomass data with the clean
#' environmental data and generate PC axes that we will use to predict the
#' missing monocultures.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load relevant functions
source("scripts/Function_plotting_theme.R")

# load the required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# load the cleaned environmental data
env_dat <- read_csv("data/case_study_2/data_clean/site_env_data.csv")
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
bio_dat <- read_csv("data/case_study_2/data_clean/biomass_data.csv")
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
pca_x <- prcomp(reformulate(names(bio_env)[-c(1:3, 13:16)]), 
                data = bio_env, scale = TRUE)

# check the PCA summary
summary(pca_x)

# bind the PCA to identifier variables
pca_env <- bind_cols(bio_env[,c(1:3)], as_tibble(pca_x$x)[,c(1:2)] )

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

# make a PCA without cluster F for the extended data
unique(pca_env$cluster_id)
pca_x_sum <- summary(pca_x)
print(pca_x_sum)

# get a colour palette
col_pal <- wesanderson::wes_palette(name = "Darjeeling1", n = 9, type = "continuous")

# plot the PCA
p1 <- 
  ggplot(data = pca_env %>% filter(cluster_id != "F")) +
  geom_point(mapping = aes(x = PC1, y = PC2, colour = cluster_id),
             size = 3, alpha = 0.2) +
  scale_colour_manual(values = col_pal) +
  theme_meta() +
  guides(colour = guide_legend(override.aes = list(alpha = 0.8) )) +
  labs(colour = "Cluster") +
  xlab("PC1 (42%)") +
  ylab("PC2 (32%)") +
  theme(legend.position = "right",
        legend.key = element_rect(fill = NA))
plot(p1)

# plot the loading plot
load_plot <- as_tibble(pca_x_sum$rotation)

# get the names of the variables
row.names(pca_x_sum$rotation)
var_names <- c("Depth (m)", "Seabed dist. (m)", "Mean temp. (C)",
               "CV temp. (C)", "Max temp. (C)", "Min temp. (C)",
               "Mean light (lux)", "CV light (lux)", "Max light (lux)")

# add the variable names to the loadings
load_plot$var_name <- var_names

p2 <- 
  ggplot(data = load_plot,
       mapping = aes(x = var_name, y = PC2)) +
  geom_hline(yintercept = 0) +
  geom_col(width = 0.5, fill = "grey", colour = "black") +
  theme_meta() +
  ylab("PC2 loadings") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust = 0.65),
        plot.margin = margin(0.25, 0.25, 0.25, 1, "cm"))
plot(p2)

p3 <-
  ggplot(data = load_plot,
         mapping = aes(x = var_name, y = PC1)) +
  geom_hline(yintercept = 0) +
  geom_col(width = 0.5, fill = "grey", colour = "black") +
  theme_meta() +
  ylab("PC1 loadings") +
  xlab(NULL) +
  theme(axis.text.x = element_blank(),
        plot.margin = margin(0.25, 0.25, 0.25, 1, "cm"))
plot(p3)

# merge the two loading plots
p23 <- 
  ggarrange(p3, p2, ncol = 1, nrow = 2, heights = c(1, 1.7),
            labels = c("b", "c"),
            font.label = list(face = "plain", size = 11),
            hjust = -4.5)
plot(p23)

# merge with the PCA plot
p123 <- 
  ggarrange(p1, p23, 
          ncol = 2, nrow = 1, widths = c(1, 1.5),
          labels = c("a", "", ""),
          font.label = list(face = "plain", size = 11))

ggsave(filename = "manuscript/figures/app_1_fig_s5.png", p123, dpi = 600,
       unit = "cm", width = 20, height = 10)

# write out the pca_env_data
write_csv(x = pca_env, "data/case_study_2/data_clean/biomass_env_analysis_data.csv")

### END
