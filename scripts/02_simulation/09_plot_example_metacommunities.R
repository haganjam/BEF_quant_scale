#'
#' @title: Simulate metacommunities as illustrations for the supplementary material
#' 
#' @description: Simulate metacommunities and plot them
#' 
#' @details: This script simulates example metacommunities and plots the mixtures,
#' monocultures along with environments with different levels of environmental
#' autocorrelation (Fig. S1, S2 and S3)
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)

# load relevant packages
library(ggplot2)
library(dplyr)
library(here)

# load the metacommunity simulation and data-processing functions
source(here("scripts/02_simulation/01_mcomsimr_simulate_MC2_function.R"))
source(here("scripts/Function_plotting_theme.R"))

# set a seed
set.seed(134985)

# set the number of species, patches and time-points
species <- 5
patches <- 20
timesteps <- 500

t_sel <- c(100, 200, 300, 400, 500)

# using the number of patches, set-up a random evenly spread landscape
t1 <- round(seq(25, 25*4, length.out = 4), 0)
t2 <- round(seq(25, 25*5, length.out = 5), 0)
t12 <- expand.grid(t1, t2)

# pull these evenly spaced patches into a data.frame
landscape.x <- data.frame(x = t12[[1]], y = t12[[2]])
plot(landscape.x)

# generate a random dispersal matrix using the landscape.x landscape
dispersal.x <- mcomsimr::dispersal_matrix(landscape = landscape.x, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

# get the starting abundances
start_abun <- round(runif(n = species, 0, 30), 0)

# pull species attributes into a data.frame
sp_att <- 
  
  data.frame(species = 1:species,
             optima = seq(0.15, 0.85, length.out = species),
             env_niche_breadth = round(runif(species, 0.05, 0.25), 2),
             max_r = 5
  )

# generate the competition matrices
comp_mat <- 
  
  matrix(
    
    runif(n = species*species, 
          min = runif(1, 0, 0.5),
          max = runif(1, 1, 1.5) ),  
    nrow = species, 
    ncol = species
    
  )

# scale the competition matrix
comp_mat <- comp_mat*0.05

# transpose the matrix so that competition is symmetrical
comp_mat[lower.tri(comp_mat)] = t(comp_mat)[lower.tri(comp_mat)]

# set intra-specific competition to be equal
diag(comp_mat) <- 1*0.05

# use the env_generate() function to generate a landscape of environmental variables
env_var <- 
  
  mcomsimr::env_generate(landscape = landscape.x, 
                         env1Scale = round(runif(1, 100, 500), 0), 
                         timesteps = timesteps, 
                         plot = FALSE
  )

# draw a random dispersal value
dispersal <- runif(1, 0.01, 0.1)

# simulate the metacommunity
MC.x <- 
  
  sim_metacomm_BEF(patches = patches, species = species, 
                   dispersal = dispersal,
                   timesteps = timesteps, 
                   start_abun = start_abun,
                   extirp_prob = 0.00001,
                   landscape = landscape.x, disp_mat = dispersal.x, env.df = env_var, 
                   env_traits.df = sp_att, 
                   int_mat = comp_mat
  )

# clean the data.frame to calculate BEF effects
MC.x.BEF <- 
  
  Isbell_2018_cleaner(mix = MC.x$mixture, 
                      mono = MC.x$monoculture, 
                      t_sel  = t_sel
  )

# get the environmental data
MC.x.env <- 
  MC.x$mixture %>%
  filter(time %in% t_sel) %>%
  select(time, place = patch, env) %>%
  distinct()

MC.x.BEF <- left_join(MC.x.BEF, MC.x.env, by = c("time", "place"))

# plot figures for the supplementary material

# 1. mixture functioning
p1 <- 
  ggplot(data = MC.x.BEF %>% mutate(Species = as.character(species)),
       mapping = aes(x = time, y = Y, colour = Species)) +
  geom_line() +
  scale_colour_viridis_d(option = "C", end = 0.95) +
  facet_wrap(~place, scales = "free", ncol = 4, nrow = 5) +
  ylab("Mixture functioning") +
  xlab("Time-point") +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p1)

ggsave(filename = here("figures/figS1.png"), p1, dpi = 350,
       units = "cm", width = 18, height = 22.5)

# monoculture functioning
p2 <- 
  ggplot(data = MC.x.BEF %>% mutate(Species = as.character(species)),
         mapping = aes(x = time, y = M, colour = Species)) +
  geom_line() +
  scale_colour_viridis_d(option = "C", end = 0.95) +
  facet_wrap(~place, scales = "free", ncol = 4, nrow = 5) +
  ylab("Monoculture functioning") +
  xlab("Time-point") +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p2)

ggsave(filename = here("figures/figS2.png"), p2, dpi = 350,
       units = "cm", width = 18, height = 22.5)

# 3. environmental variation
env_var1 <- 
  
  mcomsimr::env_generate(landscape = landscape.x, 
                         env1Scale = 100, 
                         timesteps = timesteps, 
                         plot = FALSE
  )
env_var1$autocor <- 100

env_var2 <- 
  mcomsimr::env_generate(landscape = landscape.x, 
                         env1Scale = 500, 
                         timesteps = timesteps, 
                         plot = FALSE
  )
env_var2$autocor <- 500

# bind the environmental variables
env_plot <- 
  bind_rows(env_var1, env_var2) %>%
  mutate(Patch = as.character(patch))

p3 <- 
  ggplot(data = env_plot,
       mapping = aes(x = time, y = env1, colour = Patch)) +
  geom_line(alpha = 0.6) +
  scale_colour_viridis_d(option = "C") +
  ylab("Environment (0-1)") +
  xlab("Time-points") +
  facet_wrap(~autocor) +
  theme_meta() +
  theme(legend.position = "none")
plot(p3)

ggsave(filename = here("figures/figS3.png"), p3, dpi = 350,
       units = "cm", width = 18, height = 10)

### END
