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
      geom_point(size = 2) +
      scale_colour_manual(name = "OTU/mixture",
                          values = col_pal) +   
      scale_shape_manual(name = "OTU/mixture",
                         values = c(16, 16, 16, 16, 16, 8)) +
      theme_bw() +
      theme(legend.position = "bottom")
  )
plot(legend)


# plot the raw data of the different clusters
dat_sum %>%
  group_by(cluster_id) %>%
  summarise(n = length(unique(place)) )

# get a vector of the unique cluster ids
c_id <- unique(dat_sum$cluster_id)

for(i in 1:length(c_id)) {
  
  df_sub <- dat_sum %>% filter(cluster_id == c_id[i])
  
  # check the number of places
  N_place <- length(unique(df_sub$place)) 
  
  # choose number of rows and columns based on the number of places
  if (N_place == 5) {
    NROW <- 2
    NCOL <- 3
  } else if(N_place == 4) {
    NROW <- 2
    NCOL <- 2
  } else if(N_place == 3) {
    NROW = 1
    NCOL = 3
  }
  
  p1 <- 
    ggplot(data = df_sub,
           mapping = aes(x = time, y = biomass_mu, colour = species, shape = species)) +
    geom_line(position = position_dodge(width = 0.5), alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(mapping = aes(x = time, 
                                ymin = (biomass_low),
                                ymax = (biomass_high),
                                colour = species),
                  width = 0, position = position_dodge(width = 0.5)) +
    scale_colour_manual(name = "Species/mixture",
                        values = col_pal) +   
    scale_shape_manual(name = "Species/mixture",
                       values = c(16, 16, 16, 16, 16, 8)) +
    ylab("Biomass (g)") +
    xlab("2-weeks") +
    ggtitle(paste0("Cluster - ", c_id[i])) +
    scale_x_continuous(breaks = c(1, 2, 3)) +
    facet_wrap(~place, nrow = NROW, ncol = NCOL) +
    theme_meta() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, hjust = 0.5))
  
  # add a legend
  p1 <- ggarrange(p1, legend, heights = c(6, 1),
                  nrow = 2, ncol = 1,
                  labels = c(letters[i], ""),
                  font.label = list(size = 11, face = "plain"))
  
  # save this plot
  if (N_place == 5) {
    
    ggsave(filename = paste0("figures/", "figA1_S8", letters[i], ".svg"), p1,
           unit = "cm", width = 21, height = 16)
    
  } else if (N_place == 4) {
    
    ggsave(filename = paste0("figures/", "figA1_S8", letters[i], ".svg"), p1,
           unit = "cm", width = 14, height = 16)
    
  } else if (N_place == 3) {
    
    ggsave(filename = paste0("figures/", "figA1_S8", letters[i], ".svg"), p1,
           unit = "cm", width = 21, height = 8)
    
  }
  
}

# calculate mixture relative abundance
df_sub <- 
  dat[[1]] %>%
  filter(cluster_id == "A")

df_sub <- 
  df_sub %>%
  group_by(cluster_id, sample, place, time) %>%
  mutate(sum_Y = sum(Y)) %>%
  mutate(RA = Y/sum(Y)) %>%
  ungroup() %>%
  select(-Y, -sum_Y) %>%
  mutate(species = as.character(species))

df_sub %>%
  group_by(place, species) %>%
  summarise(M = mean(M),
            RA = mean(RA)) %>%
  ggplot(data = .,
         mapping = aes(x = M, y = RA, colour = species)) +
  geom_line()

# plot the covariance between monoculture functioning and relative abundance
cov_RA_M <- 
  data_M %>%
  filter(cluster_id != "F") %>%
  mutate(species = as.character(species)) %>%
  rename(biomass_mu = M_mu,
         biomass_sd = M_sd) %>%
  select(cluster_id, sample, place, time, species, Y, biomass_mu, biomass_sd)

# calculate relative abundance in mixture
cov_RA_M <- 
  cov_RA_M %>%
  group_by(cluster_id, place, time, sample) %>%
  mutate(Y_total = sum(Y)) %>%
  ungroup() %>%
  mutate(RA = Y/Y_total) %>%
  select(cluster_id, sample, place, time, species, RA, biomass_mu, biomass_sd)
head(cov_RA_M)
summary(cov_RA_M)

### END
