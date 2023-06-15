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
library(ggpubr)
library(vegan)

# load plotting theme
source(here("scripts/Function_plotting_theme.R"))

# load the cleaned environmental data
env_dat <- read_csv("data/case_study_2/data_clean/site_env_data.csv")

# change the depth treatment to a heterogeneity variable
env_dat$heterogeneity <-  ifelse(env_dat$depth_treatment == "alternating", "Heterogeneous", "Homogeneous")

# choose the relevant columns
names(env_dat)

env_dat <- 
  env_dat %>%
  select(cluster_id, site_id, time, heterogeneity, euclidian_dispersion_layer,
         exposure_layer, depth_vis_layer,
         panel_depth_m_measured, distance_between_panel_and_seabeed_m,
         temp_C_m, temp_C_cv, temp_C_max, temp_C_min,
         lux_m, lux_cv, lux_max,
         sal_field_ppt_m, 
         secchi_field_m,
         reduction_g_hour_m)

# summarise across the time points
env_dat_m <- 
  env_dat %>%
  group_by(cluster_id, site_id, heterogeneity) %>%
  summarise(across(.cols = names(env_dat)[-c(1:4)], ~mean(., na.rm = TRUE)),
            .groups = "drop")

# remove cluster F because we only had two patches
env_dat_m <- 
  env_dat_m %>%
  filter(cluster_id != "F")

# check the distribution of the different variables
env_dat_m %>%
  pivot_longer(cols = names(env_dat_m)[-c(1:4)],
               names_to = "env_var",
               values_to = "value") %>%
  ggplot(data = .,
         mapping = aes(x = value)) +
  geom_histogram() +
  facet_wrap(~env_var, scales = "free") +
  theme_meta()

# run a PCA on these data
pca.x <- 
  prcomp(reformulate(termlabels = names(env_dat_m[,-c(1:4)])),
         data = as_tibble(apply(env_dat_m[,-c(1:4)], 2, scale)), 
         scale = FALSE, center = FALSE)

# check the summary statistics
summary(pca.x)

# check the biplot
biplot(pca.x)

# plot the first two axes
pca_df <- as_tibble(pca.x$x[,c(1, 2)])

# add the cluster_id information
pca_df$cluster_id <- env_dat_m$cluster_id

# add the heterogeneous or homogeneous treatment
pca_df$heterogeneity <- env_dat_m$heterogeneity

# check how many sites there are per cluster
pca_df %>%
  group_by(heterogeneity, cluster_id) %>%
  summarise(n = n())

# generate a centroid dataset to plot cluster IDs
pca_df_sum <- 
  pca_df %>%
  group_by(heterogeneity, cluster_id) %>%
  summarise(PC1 = mean(PC1),
            PC2 = mean(PC2))

# plot the PC plot
p1 <- 
  ggplot() +
  geom_point(data = pca_df, 
             mapping = aes(x = PC1, y = PC2, colour = cluster_id, shape = heterogeneity),
             size = 3, alpha = 0.6) +
  geom_point(data = pca_df_sum, 
             mapping = aes(x = PC1, y = PC2, colour = cluster_id, shape = heterogeneity),
             position = position_dodge2(width = 0.5),
             size = 4.2, alpha = 1) +
  geom_label(data = pca_df_sum, 
             mapping = aes(x = PC1, y = PC2, label = cluster_id),
             position = position_dodge2(width = 0.5),
             label.size = NA, alpha = 0, size = 3, colour = "white",
             fontface = "bold") +  
  scale_colour_manual(values = wesanderson::wes_palette(name = "Cavalcanti1", n = 9, type = "continuous")) +
  guides(colour = "none") +
  theme_meta() +
  xlab("PC1 (38%)") +
  ylab("PC2 (32%)") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA))
plot(p1)

# generate a Euclidean distance matrix
env_d <- dist(apply(env_dat_m[,-c(1:4)], 2, scale), method = "euclidean")
names(env_dat_m[,-c(1:4)])

# multivariate dispersion
env_d <- betadisper(d = env_d, group = env_dat_m$cluster_id)

# add the distance to centroid to the data
env_dat_m$distance_centroid <- env_d$distances

# does this correlate with the originally implemented centroid distances
env_dispersion <- 
  env_dat_m %>%
  group_by(heterogeneity, cluster_id) %>%
  summarise(GIS_dispersion = mean(euclidian_dispersion_layer),
            field_dispersion = mean(distance_centroid),
            .groups = "drop")

# plot the relationship
plot(env_dispersion$GIS_dispersion, env_dispersion$field_dispersion)

# examine the spearman correlation which tests for monotonic relationships
cor(env_dispersion$GIS_dispersion, env_dispersion$field_dispersion, method = "spearman")

# compare the multivariate dispersion between homogeneous and heterogeneous treatments
env_dispersion_s <- 
  env_dispersion %>%
  group_by(heterogeneity) %>%
  summarise(d_m = mean(field_dispersion),
            d_sd = sd(field_dispersion))

p2 <- 
  ggplot() +
  geom_jitter(data = env_dispersion,
              mapping = aes(x = heterogeneity, y = field_dispersion, shape = heterogeneity), 
              width = 0.05, size = 2, alpha = 0.2) +
  geom_point(data = env_dispersion_s,
             mapping = aes(x = heterogeneity, y = d_m, shape = heterogeneity),
             colour = "black", size = 3) +
  geom_errorbar(data = env_dispersion_s,
             mapping = aes(x = heterogeneity, ymin = d_m-d_sd, ymax = d_m+d_sd ),
             width = 0,
             colour = "black", size = 0.5) +
  guides(shape = "none") +
  xlab("") +
  ylab("Multivariate dispersion") +
  theme_meta() +
  theme(legend.position = "none")
plot(p2)

# do a t.test()
t.test(field_dispersion ~ heterogeneity, data = env_dispersion)

# arrange the plots
p12 <- 
  ggarrange(p2, p1, ncol = 2, nrow = 1,
          widths = c(1, 1.25), labels = c("a", "b"),
          font.label = list(face = "plain", size = 11), common.legend = TRUE)
plot(p12)

ggsave(filename = here("figures/SI2_fig_4.svg"), p12,
       unit = "cm", width = 20, height = 10)

# output the multivariate dispersion into an .rds file
saveRDS(object = env_dispersion, file = here("results/benthic_env_dispersion.rds"))

### END
