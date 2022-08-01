#'
#' @title: Clean the experiment data from the initial measurements
#' 
#' @description: Script to clean and output a cleaned version of the initial experiment
#' data. This uses a combination of the direct measurements and the measurements from the
#' processed images.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com) and Benedikt Schrofner-Brunner (bschrobru(at)gmail.com)
#'

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

# next steps:

# 1. split the code so that we can push the most computationally intensive part of the 
# code i.e. using the posterior distribution and calculating biodiversity effects for each
# sample 100 times with different samples from the Dirichlet distribution

# - currently, this takes five minutes

# if we run 100 models, without that step, it should only take 2 hours and then we don't
# have the hassle of trying to load these obscure packages to the cluster

# once that it is done, we can save the outputs and then, in a new script run 
# the calculates on the Albiorix cluster

# this will save a lot of time


# load the relevant packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
library(pbapply)
library(rethinking)

# load the metacommunity simulation and data-processing functions
source(here("scripts/02_simulation/01_mcomsimr_simulate_MC2_function.R"))

# set the model inputs

# set a seed
set.seed(54258748)

# number of replicate simulations
N_REP <- 2

# set-up basic model inputs
species <- 5
patches <- 20
timesteps <- 100
dispersal <- 0.05
extirp_prob <- 0.00001

# maximum growth parameters
max_r <- 5
K_max <- NA

# set-up parameters for the competition matrices
int_min <- 0.05*0
int_max <- 0.05*1.5
intra <- 0.05*1

# set-up the type of Lotka-Volterra competition
comp <- "Beverton_Holt"

# landscape parameters

# generate a random landscape
n1 <- 5
n2 <- 4
assertthat::assert_that(assertthat::see_if( (n1*n2) == patches ),
                        msg = "change the n1 and n2 values so that their product matches the patches")

# using the number of patches, set-up a random evenly spread landscape
t1 <- round(seq(25, 100, length.out = 5), 0)
t2 <- round(seq(25, 100, length.out = 4), 0)
t12 <- expand.grid(t1, t2)

# pull these evenly spaced patches into a data.frame
landscape.x <- data.frame(x = t12[[1]], y = t12[[2]])

# generate a random dispersal matrix using the landscape.x landscape
dispersal.x <- dispersal_matrix(landscape = landscape.x, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

# generate species environmental optima
optima <- seq(0.15, 0.85, length.out = species)
print(optima)
env_niche_breadth <- lapply(1:N_REP, function(x) round(runif(species, 0.1, 0.25), 2))
print(env_niche_breadth)

# select time-points to calculate biodiversity effects for
t_sel <- c(10, 20, 30, 40, 50)

# set different levels of starting abundances for the different species
start_abun <- round(runif(n = species, 0, 30), 0)


# run the pipeline which does the following:

# 1. simulates a BEF experiment with full mixtures in monocultures for each patch

# 2. uses the known starting abundances and full monoculture data to calculate observed biodiversity effects (Isbell et al. 2018)

# 3. assumes that we only have monoculture data in 30 % of the patches

# 4. uses a simple linear model to model the missing monoculture data using the environment and mixture abundance of each species

MC_sims <- 
  
  pblapply(1:N_REP, function(x) {
  
  # get the starting abundances
  start_abun <- round(runif(n = species, 5, 30), 0)
  
  # pull species attributes into a data.frame
  sp_att <- 
    
    data.frame(species = 1:species,
               optima = optima,
               env_niche_breadth = env_niche_breadth[[x]],
               max_r = max_r,
               K_max = K_max
               )
  
  # generate the competition matrices
  comp_mat <- 
    
    matrix(
      
      runif(n = species*species, min = int_min, max = int_max),  
      nrow = species, 
      ncol = species
      
      )
  
  # transpose the matrix so that competition is symmetrical
  comp_mat[lower.tri(comp_mat)] = t(comp_mat)[lower.tri(comp_mat)]
  
  # set intra-specific competition to be equal
  diag(comp_mat) <- intra
  
  # simulate environmental variables
  
  # randomly draw an autocorrelation value between 250 and 500
  autocorr <- round(runif(1, 250, 500), 0)
  
  # use the env_generate() function to generate a landscape of environmental variables
  env_var <- 
    
    mcomsimr::env_generate(landscape = landscape.x, 
                           env1Scale = autocorr, 
                           timesteps = timesteps, 
                           plot = FALSE
                           )
  
  # simulate the metacommunity
  MC.x <- 
    
    sim_metacomm_BEF(patches = patches, species = species, 
                     dispersal = dispersal,
                     timesteps = timesteps, 
                     start_abun = start_abun,
                     extirp_prob = extirp_prob,
                     comp = comp,
                     landscape = landscape.x, disp_mat = dispersal.x, env.df = env_var, 
                     env_traits.df = sp_att, 
                     int_mat = comp_mat
                     )
  
  # output the model identifiers
  MC.x.ids <-
    
    data.frame("start_abun" = paste(start_abun, collapse = "_"),
               "inter_comp" = mean(comp_mat[lower.tri(comp_mat)]),
               "optima" = paste(sp_att$optima, collapse = "_"),
               "niche_breadth" = paste(sp_att$env_niche_breadth, collapse = "_"),
               "t_steps" = paste(t_sel, collapse = "_"),
               "dispersal" = dispersal
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
  
  # calculate the observed biodiversity effects
  
  # get the known initial relative abundance values
  RYe_in <- 
    
    if(length(start_abun) == 1) { 
    
    rep(1/start_abun, species) 
      
    }  else { 
        
      start_abun/sum(start_abun) 
      
      }
  
  # use the Isbell_2018_sampler function to calculate the actual biodiversity effects
  BEF_obs <- Isbell_2018_sampler(data = MC.x.BEF, RYe_post = FALSE, RYe = RYe_in )
  
  BEF_obs <- 
    
    bind_rows(
      
      BEF_obs$Beff, 
      rename(BEF_obs$L.Beff, Beff = L.Beff) 
      
      ) %>%
    
    mutate(Value = round(Value, 2))
  
  # assume incomplete monoculture and unknown RYE
  MC.x.NA <- 
    
    full_join(MC.x.BEF, 
              MC.x.env, 
              by = c("time", "place")) %>%
    
    arrange(sample, time, place, species)
  
  # select 30% of patches to have monoculture and mixture data
  
  # get the number of monoculture and mixtures patches which we set at 30% 
  n_MoMi <- round( 0.30 * length( unique(MC.x.NA$place) ), 0 )
  p_comb <- combn(unique(MC.x.NA$place), m = n_MoMi)
  
  # calculate the environmental range among each of the combinations of two places
  max_env <- 
    apply(p_comb, 2, function(combination) {
      
      range_comb <- 
        
        MC.x.NA %>%
        filter(place %in% combination) %>%
        pull(env) %>%
        range(.)
      
      # take the difference in range
      output <- diff(range_comb)
      
      return(output)
      
    })
  
  # most likely that we would have say 2 places and all times: Pick the 2 places with most environmental variation
  MC.x.NA$M1 <- ifelse( MC.x.NA$place %in% p_comb[, sample(which(max_env == max(max_env)), 1)], MC.x.NA$M, NA)
  
  
  # model the monoculture yields using rethinking()
  
  # get the complete cases to fit the model (i.e. training data), (MC1_NAcc)
  v <- MC.x.NA[complete.cases(MC.x.NA), ]
  
  # set-up a list to run through ulam()
  MC.x.train <- list(
    M = v$M ,
    Y = standardize(v$Y ) ,
    E = v$env,
    S = v$species)
  
  # model the monocultures with with the E and Y interaction
  m1 <- ulam(
    alist(
      M ~ dpois( lambda ),
      log(lambda) <- aS[S] + b_yS[S]*Y + b_eS[S]*E + b_yeS[S]*E*Y,
      
      # priors
      aS[S] ~ dnorm(2, 2),
      b_yS[S] ~ normal(0, 1),
      b_eS[S] ~ normal(0, 1),
      b_yeS[S] ~ normal(0, 1)
      
    ), data = MC.x.train , chains = 4, iter = 2000, log_lik = FALSE)
  
  # predict the missing data
  MC.x.pred <- MC.x.NA[is.na(MC.x.NA$M1), ]
  MC.x.pred <- MC.x.pred[, c("species", "Y", "env")]
  names(MC.x.pred) <- c("S", "Y", "E") 
  MC.x.pred$Y <- standardize(MC.x.pred$Y)
  
  # simulate observations from these data
  MC.x.pred <- sim(m1, data = MC.x.pred)
  
  # pull the relevant outputs into a list
  output.list <-
    
    list("MC.x.ids" = MC.x.ids,
         "BEF_obs" = BEF_obs,
         "MC.x.NA" = MC.x.NA,
         "MC.x.pred" = MC.x.pred
         )
  
  return(output.list)
  
  })

# save the MC_sims object
saveRDS(object = MC_sims, here("results/MC_sims.rds"))

# generate and save a set of 100 potential starting relative abundances from the Dirichlet distribution
start_RA <- sapply(1:100, function(x) gtools::rdirichlet(n = 1, rep(3,  species) ) )
saveRDS(object = start_RA, here("results/MC_sims_start_RA.rds"))
  
### END
