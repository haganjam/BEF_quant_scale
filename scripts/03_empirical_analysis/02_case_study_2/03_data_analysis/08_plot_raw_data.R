#'
#' @title: Plot the raw biomass data of mixtures and monocultures for benthic data
#' 
#' @description: Plot the raw biomass data for the mixtures and monocultures of
#' the benthic fouling data: Appendix 1: Figure S8.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(dplyr)
library(readr)
library(ggpubr)
library(cowplot)
library(ggplot2)

# load relevant scripts
source("scripts/Function_plotting_theme.R")
source("scripts/03_empirical_analysis/helper_functions.R")

# load the monoculture prediction data
dat <- readRDS(file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/05_complete_data.rds")
length(dat)

# get the mixture data
dat_Y <- 
  dat[[1]] %>%
  group_by(cluster_id, sample, place, time) %>%
  summarise(biomass_mu = sum(Y), .groups = "drop")

# add a species richness variable and convert species to a character
dat_Y <-
  dat_Y %>%
  mutate(SR = 5,
         biomass_low = biomass_mu,
         biomass_high = biomass_mu,
         species = "Mixture (5 sp.)") %>%
  select(cluster_id, sample, place, time, SR, species, biomass_mu, biomass_low, biomass_high)

# bind into a data.frame
dat_df <- bind_rows(dat, .id = "rep")

# summarise into mean and standard deviations
dat_M <- 
  dat_df %>%
  group_by(cluster_id, sample, place, time, species) %>%
  summarise(biomass_mu = mean(M),
            biomass_low = HPDI(M, 0.95)[1],
            biomass_high = HPDI(M, 0.95)[2], .groups = "drop")
print(dat_M)

# add a species richness variable and convert species to a character
dat_M <-
  dat_M %>%
  mutate(SR = 1,
         species = as.character(species)) %>%
  select(cluster_id, sample, place, time, SR, species, biomass_mu, biomass_low, biomass_high)

# bind the monoculture and mixture data
dat_sum <- 
  bind_rows(dat_M , dat_Y) %>%
  arrange(cluster_id, sample, place, time, SR, species, biomass_mu, biomass_low, biomass_high)
print(dat_sum)

# rename the factor levels
dat_sum$species <- factor(dat_sum$species, 
                          levels = c("1", "2", "3", "4", "5", "Mixture (5 sp.)"))
levels(dat_sum$species) <- c("Barn", "Bryo", "Asci", "Hydro", "Ciona", "Mixture (5 sp.)")

# remove cluster F from the analysis
dat_sum <- 
  dat_sum %>%
  filter(cluster_id != "F")

# calculate the average monoculture functioning value
dat_sum %>%
  group_by(SR) %>%
  summarise(M = mean(biomass_mu))

# illustrate the net biodiversity effect
df1 <- 
  dat_sum %>%
  filter(cluster_id == "B", place == 7, time == 1)
df1 %>%
  group_by(SR) %>%
  summarise(biomass_mu = mean(biomass_mu))

# get a colour palette
col_pal <- wesanderson::wes_palette(name = "Darjeeling1", n = 6, type = "continuous")

# make a legend
legend <- dat_sum[1:6, ]

legend <- 
  get_legend( 
    ggplot(data = legend,
           mapping = aes(x = SR, y = biomass_mu, colour = species, shape = species)) +
      geom_point(size = 3) +
      scale_colour_manual(name = "OTU/mixture",
                          values = col_pal) +   
      scale_shape_manual(name = "OTU/mixture",
                         values = c(16, 16, 16, 16, 16, 8)) +
      theme_bw() +
      theme(legend.position = "right") +
      theme(legend.text = element_text(size = 11),
            legend.title = element_text(size = 12))
  )
plot(legend)


# plot the raw data of the different clusters
dat_sum %>%
  group_by(cluster_id) %>%
  summarise(n = length(unique(place)) )

# get a vector of the unique cluster ids
c_id <- unique(dat_sum$cluster_id)

# ylabs
ylabs <- list("", "", "", "", "Dry biomass (g)", "", "", "", "")

plot_list <- vector("list", length = length(c_id))
for(i in 1:length(c_id)) {
  
  df_sub <- dat_sum %>% filter(cluster_id == c_id[i])
  
  # check the number of places
  N_place <- length(unique(df_sub$place)) 
  
  # choose number of rows and columns based on the number of places
  if (N_place == 5) {
    NROW <- 1
    NCOL <- 5
  } else if(N_place == 4) {
    NROW <- 1
    NCOL <- 4
  } else if(N_place == 3) {
    NROW = 1
    NCOL = 3
  }
  
  p1 <- 
    ggplot(data = df_sub,
           mapping = aes(x = time, y = biomass_mu, colour = species, shape = species)) +
    geom_line(position = position_dodge(width = 0.5), alpha = 0.5,
              show.legend = FALSE) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(mapping = aes(x = time, 
                                ymin = (biomass_low),
                                ymax = (biomass_high),
                                colour = species),
                  width = 0, position = position_dodge(width = 0.5),
                  show.legend = FALSE) +
    scale_colour_manual(name = "OTU/mixture",
                        values = col_pal) +   
    scale_shape_manual(name = "OTU/mixture",
                       values = c(16, 16, 16, 16, 16, 8)) +
    ylab(ylabs[[i]]) +
    xlab(if(i == 9){"2-weeks"}else{NULL}) +
    ggtitle(paste0("Cluster - ", c_id[i])) +
    scale_x_continuous(breaks = c(1, 2, 3)) +
    facet_wrap(~place, nrow = NROW, ncol = NCOL) +
    theme_meta() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  plot_list[[i]] <- p1
  
}

p1 <- 
  cowplot::plot_grid(plot_list[[1]], plot_list[[2]],plot_list[[3]],
                     plot_list[[4]],plot_list[[5]],plot_list[[6]],
                     plot_list[[7]],plot_list[[8]],plot_list[[9]],
                     nrow = 9, ncol = 1, 
                     align = "hv",
                     axis = "l")

p1 <- ggarrange(p1, legend, ncol = 2, nrow = 1,
                widths = c(8, 1))

ggsave(filename = "figures/fig_ED1.svg", p1,
       unit = "cm", width = 35, height = 45)


# calculate mixture relative abundance
dat_cov <- 
  dat_df %>%
  filter(cluster_id != "F") %>% 
  group_by(rep, cluster_id, sample, place, time) %>%
  mutate(sum_Y = sum(Y)) %>%
  mutate(RA = Y/sum(Y)) %>%
  ungroup() %>%
  select(-Y, -sum_Y) %>%
  mutate(species = as.character(species))

# replace the species characters with species names
dat_cov$species <- factor(dat_cov $species, 
                          levels = c("1", "2", "3", "4", "5"))
levels(dat_cov$species) <- c("Barn", "Bryo", "Asci", "Hydro", "Ciona")

# calculate the mean for each place across times
cov_place <- 
  dat_cov %>%
  group_by(cluster_id, place, species) %>%
  summarise(M_mu = mean(M),
            M_sd = sd(M),
            RA_mu = mean(RA),
            RA_sd = sd(RA)) %>%
  ungroup()

# change cluster_id levels
cov_place$cluster_id <- factor(cov_place$cluster_id)
levels(cov_place$cluster_id) <- paste0("Cluster ", c("A", "B", "C", "D", "E", "G", "H", "I", "J"))

p1 <- 
  ggplot(data = cov_place) +
  geom_point(mapping = aes(x = M_mu, y = RA_mu, colour = species)) +
  geom_errorbar(mapping = aes(x = M_mu, 
                              ymin = RA_mu - RA_sd, ymax = RA_mu + RA_sd,
                              colour = species), 
                width = 0, linewidth = 0.2, alpha = 0.5,
                show.legend = FALSE) +
  geom_errorbarh(mapping = aes(y = RA_mu, 
                               xmin = M_mu - M_sd, xmax = M_mu + M_sd,
                               colour = species), 
                 height = 0, linewidth = 0.2, alpha = 0.5,
                 show.legend = FALSE) +
  scale_colour_manual(values = col_pal[1:5]) +
  facet_wrap(~cluster_id, scales = "free") +
  guides(colour = guide_legend(override.aes = list(size = 3.5))) +
  labs(colour = "OTU") +
  xlab("Monoculture dry biomass (g)") +
  ylab("Mixture relative abundance") +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA),
        strip.background = element_rect(fill="white"))
plot(p1)

ggsave(filename = "figures/fig_ED2.svg", p1,
       unit = "cm", width = 18, height = 20)

# calculate the mean for each time across places
cov_time <- 
  dat_cov %>%
  group_by(cluster_id, time, species) %>%
  summarise(M_mu = mean(M),
            M_sd = sd(M),
            RA_mu = mean(RA),
            RA_sd = sd(RA)) %>%
  ungroup()

# change cluster_id levels
cov_time $cluster_id <- factor(cov_time$cluster_id)
levels(cov_time $cluster_id) <- paste0("Cluster ", c("A", "B", "C", "D", "E", "G", "H", "I", "J"))

p2 <- 
  ggplot(data = cov_time) +
  geom_point(mapping = aes(x = M_mu, y = RA_mu, colour = species)) +
  geom_errorbar(mapping = aes(x = M_mu, 
                              ymin = RA_mu - RA_sd, ymax = RA_mu + RA_sd,
                              colour = species), 
                width = 0, linewidth = 0.2, alpha = 0.5,
                show.legend = FALSE) +
  geom_errorbarh(mapping = aes(y = RA_mu, 
                               xmin = M_mu - M_sd, xmax = M_mu + M_sd,
                               colour = species), 
                 height = 0, linewidth = 0.2, alpha = 0.5,
                 show.legend = FALSE) +
  scale_colour_manual(values = col_pal[1:5]) +
  facet_wrap(~cluster_id, scales = "free") +
  guides(colour = guide_legend(override.aes = list(size = 3.5))) +
  labs(colour = "OTU") +
  xlab("Monoculture dry biomass (g)") +
  ylab("Mixture relative abundance") +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA),
        strip.background = element_rect(fill="white"))
plot(p2)

ggsave(filename = "figures/fig_ED3.svg", p2,
       unit = "cm", width = 18, height = 20)

### END
