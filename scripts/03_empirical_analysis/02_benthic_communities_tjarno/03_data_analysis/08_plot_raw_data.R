
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

# check if any monoculture predictions are zero
lapply(sp_mono, function(x) {any((x < 0) == TRUE)} )
lapply(sp_mono, function(x) nrow(x)) %>% unlist() %>% sum() + 287 # good, it should be 615

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

# remove cluster F from the analysis
data_comb <- 
  data_comb %>%
  filter(cluster_id != "F")

# calculate the average monoculture functioning value
data_comb %>%
  group_by(SR) %>%
  summarise(M = mean(biomass_mu))

# illustrate the net biodiversity effect
df1 <- 
  data_comb %>%
  filter(cluster_id == "A", place == 1, time == 1)
df1.sum <- 
  df1 %>%
  group_by(SR) %>%
  summarise(biomass_mu = mean(biomass_mu))


# make a legend
legend <- data_comb[1:6, ]

legend <- 
  get_legend( 
    ggplot(data = legend,
           mapping = aes(x = SR, y = biomass_mu, colour = species, shape = species)) +
      geom_point(size = 2) +
      scale_colour_manual(name = "Species/mixture",
                          values = viridis(n = 6, begin = 0.1, end = 0.9, option = "C")) +   
      scale_shape_manual(name = "Species/mixture",
                         values = c(16, 16, 16, 16, 16, 8)) +
      theme_bw() +
      theme(legend.position = "bottom")
  )
plot(legend)


# plot the raw data of the different clusters

data_comb %>%
  group_by(cluster_id) %>%
  summarise(n = length(unique(place)) )

# get a vector of the unique cluster ids
c_id <- unique(data_comb$cluster_id)

for(i in 1:length(c_id)) {
  
  df_sub <- data_comb %>% filter(cluster_id == c_id[i])
  
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
                                ymin = (biomass_mu - biomass_sd),
                                ymax = (biomass_mu + biomass_sd),
                                colour = species),
                  width = 0, position = position_dodge(width = 0.5)) +
    scale_colour_manual(name = "Species/mixture",
                        values = viridis(n = 6, begin = 0.1, end = 0.9, option = "C")) +   
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
    
    ggsave(filename = here(paste0("figures/", "figS8", letters[i], ".png")), p1,
           unit = "cm", width = 21, height = 16)
    
  } else if (N_place == 4) {
    
    ggsave(filename = here(paste0("figures/", "figS8", letters[i], ".png")), p1,
           unit = "cm", width = 14, height = 16)
    
  } else if (N_place == 3) {
    
    ggsave(filename = here(paste0("figures/", "figS8", letters[i], ".png")), p1,
           unit = "cm", width = 21, height = 8)
    
  }
  
}

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
