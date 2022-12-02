
# load the required libraries
library(dplyr)
library(here)
library(ggpubr)
library(ggplot2)

# load relevant scripts
source(here("scripts/Function_plotting_theme.R"))

# load the analysis data
data_M <- readRDS(here("results/benthic_BEF_data.rds"))

# change the relevant names
data_M <- 
  data_M %>%
  rename(M_mu = M1) %>%
  mutate(M_sd = ifelse(!is.na(M_mu), 0, NA))

# load the monoculture prediction data
sp_mono <- readRDS(file = here("results/benthic_mono_pred.rds"))

# summarise into a mean and standard deviation
sp_mu <- 
  lapply(sp_mono, function(data) {
  apply(data, 1, mean)
} )

sp_sd <- 
  lapply(sp_mono, function(data) {
    apply(data, 1, sd)
  } )

# add monoculture predictions from one sample into the data.frame
for(j in 1:length(sp_mono)) {
  
  # get missing rows
  x <- which( (is.na(data_M[["M"]])) & (data_M[["species"]] == j) )
  
  data_M[x, ][["M_mu"]] <- sp_mu[[j]]
  data_M[x, ][["M_sd"]] <- sp_sd[[j]]
  
}


# plotting the raw data mixtures and monocultures
data_M$species <- factor(data_M$species)
levels(data_M$species) <- names(sp_mono)

data_M1 <- 
  data_M %>%
  mutate(SR = 1,
         species = as.character(species)) %>%
  rename(biomass_mu = M_mu,
         biomass_sd = M_sd) %>%
  select(cluster_id, sample, place, time, SR, species, biomass_mu, biomass_sd)

data_M2 <- 
  data_M %>%
  group_by(cluster_id, sample, place, time) %>%
  summarise(biomass_mu = sum(Y), .groups = "drop") %>%
  mutate(SR = 4,
         species = "Mixture (5 sp.)", 
         biomass_sd = 0) %>%
  select(cluster_id, sample, place, time, SR, species, biomass_mu, biomass_sd)

data_comb <- 
  bind_rows(data_M1 , data_M2) %>%
  arrange(cluster_id, sample, place, time, SR, species, biomass_mu, biomass_sd)

data_comb $species <- factor(data_comb $species)
levels(data_comb$species) <- list("Barn" = "Barn",
                                  "Bryo" = "Bryo",
                                  "Bumpi" = "Bumpi",
                                  "Hydro" = "Hydro",
                                  "Seasq" = "Seasq",
                                  "Mixture (5 sp.)" = "Mixture (5 sp.)"
                                 )

# make a legend
legend <- data_comb[1:6, ]

legend <- 
  get_legend( 
    ggplot(data = legend,
           mapping = aes(x = SR, y = biomass_mu, colour = species)) +
      geom_point(size = 2) +
      scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "C") +
      theme_bw() +
      theme(legend.position = "right")
  )
plot(legend)

# plot one cluster 1

c_id <- "A"

p1 <- 
  ggplot(data = data_comb %>% filter(cluster_id == c_id),
         mapping = aes(x = time, y = biomass_mu, colour = species)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(mapping = aes(x = time, 
                              ymin = (biomass_mu - biomass_sd),
                              ymax = (biomass_mu + biomass_sd),
                              colour = species),
                width = 0, position = position_dodge(width = 0.5)) +
  scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "C") +
  ylab("Biomass (g)") +
  xlab("2-weeks") +
  ggtitle(paste0("Cluster - ", c_id)) +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  facet_wrap(~place) +
  theme_meta() +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, hjust = 0.5))
plot(p1)



