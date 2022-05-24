
# Project: Quantifying biodiversity effects across scales in natural ecosystems
# Title: Can we reliably quantify biodiversity effects across scales in natural ecosystems?
# Author: James Hagan
# Date: 2022/04/24

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19043)

# attached base packages:
# stats, graphics, grDevice, utils, datasets, methods, base     

# other attached packages:
# mcomsimr_0.1.0, pbapply_1.5-0, here_1.0.1, dplyr_1.0.7, tidyr_1.1.4, ggplot2_3.3.5 

# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
library(pbapply)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))

# set the model inputs

# set a seed
set.seed(54258748)

species <- 5
patches <- 5
timesteps <- 100
dispersal = 0.05
extirp_prob = 0

max_r = 0.5
K_max = 150

int_min = 0.75
int_max = 0.75
intra = 0.75

# landscape parameters
# generate a random landscape
l.1 <- 
  data.frame(x = c(25, 50, 75, 50, 50),
             y = c(50, 50, 50, 75, 25))
l.1 <- landscape_generate(patches = patches, plot = FALSE)

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

# generate species environmental optima
# optima = seq(0.1, 0.9, length.out = species)
optima = seq(0.15, 0.85, length.out = species)
print(optima)
env_niche_breadth = lapply(1:N_REP, function(x) round(runif(species, 0.1, 0.25), 2))
print(env_niche_breadth)

# selected time points
t_sel <- c(10, 20, 30, 40, 50)

# get the starting abundances
start_abun <- round(runif(n = species, 5, 30), 0)

t.1 <- 
  data.frame(species = 1:species,
             optima = optima,
             env_niche_breadth = env_niche_breadth[[1]],
             max_r = max_r,
             K_max = K_max)

# competition matrices
si.1 <- matrix(runif(n = species*species, min = int_min, max = int_max), 
               nrow = species, ncol = species)
si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
diag(si.1) <- intra

# simulate environmental variables
autocorr <- round(runif(1, 100, 500), 0)
e.1 <- mcomsimr::env_generate(landscape = l.1, 
                              env1Scale = autocorr, timesteps = timesteps, plot = FALSE)


# simulate the metacommunity
MC1a <- sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal,
                         timesteps = timesteps, start_abun = start_abun,
                         extirp_prob = extirp_prob,
                         comp = "Lotka_Volterra",
                         landscape = l.1, disp_mat = d.1, env.df = e.1, 
                         env_traits.df = t.1, int_mat = si.1)

# clean the data.frame to calculate BEF effects
MC1b <- Isbell_2018_cleaner(mix = MC1a$mixture, mono = MC1a$monoculture, 
                            t_sel  = t_sel)

# get the environmental data
MC1e <- 
  MC1a$mixture %>%
  filter(time %in% t_sel) %>%
  select(time, place = patch, env) %>%
  distinct()

# assume incomplete monoculture and unknown RYE
MC1_NA <- 
  full_join(MC1b, MC1e, by = c("time", "place")) %>%
  arrange(sample, time, place, species)

p_comb <- combn(unique(MC1_NA$place), m = round(0.30*length(unique(MC1_NA$place)), 0) )
max_env <- 
  apply(p_comb, 2, function(x) {
    
    y <- 
      MC1_NA %>%
      filter(place %in% x) %>%
      pull(env) %>%
      range(.)
    
    return(diff(y))
    
  })

# most likely that we would have say 2 places and all times
MC1_NA$M1 <- ifelse((MC1_NA$place %in% p_comb[, which(max_env == max(max_env))] ), MC1_NA$M, NA)

MC1_NA <- 
  MC1_NA %>%
  group_by(sample, time, place) %>%
  mutate(functioning = sum(Y)) %>%
  mutate(Y = Y/functioning) %>%
  mutate(functioning = rpois(n = 1, functioning)) %>%
  pivot_wider(id_cols = c("sample", "time", "place", "env", "functioning"),
              names_from = "species",
              values_from = "Y")
names(MC1_NA) <- c("sample", "time", "place", "env", "functioning", "sp1", "sp2", "sp3", "sp4", "sp5")

# model the monoculture yields using rethinking()
library(rethinking)

# get a list with the relevant data
MC_dat <- list(
  Func = MC1_NA$functioning,
  sp1 = MC1_NA$sp1,
  sp2 = MC1_NA$sp2,
  sp3 = MC1_NA$sp3,
  sp4 = MC1_NA$sp4,
  sp5 = MC1_NA$sp5,
  E = MC1_NA$env)
str(MC_dat)
length(MC_dat$Func)

# model the full model with E and Y interaction
m1 <- ulam(
  alist(
    Func ~ dpois( lambda ),
    log(lambda) <- a + b1*sp1 + b2*sp2 + b3*sp3 + b4*sp4 + b5*sp5 + b6*E,
    
    # priors
    a ~ dnorm(2, 2),
    b1 ~ normal(0, 1),
    b2 ~ normal(0, 1),
    b3 ~ normal(0, 1),
    b4 ~ normal(0, 1),
    b5 ~ normal(0, 1),
    b6 ~ normal(0, 1)
    
  ), data=MC_dat , chains = 4, iter = 2000, log_lik = FALSE)
traceplot(m1)

MC1_pred <- data.frame(
            sp1 = c(1, 0, 0, 0, 0),
            sp2 = c(0, 1, 0, 0, 0),
            sp3 = c(0, 0, 1, 0, 0),
            sp4 = c(0, 0, 0, 1, 0),
            sp5 = c(0, 0, 0, 0, 1))

MC1_pred <- 
  lapply(MC1_NA$env, function(x) { 
  y <- MC1_pred
  y$E <- x
  return(y)} ) %>%
  bind_rows(.)

# simulate observations from these data
MC1_pred <- sim(m1, data = MC1_pred)
dim(MC1_pred)

MC1_NA






    