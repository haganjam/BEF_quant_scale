#'
#' @title: Simulate metacommunities to test the analytical pipeline
#' 
#' @description: Simulate BEF experiments and calculate biodiversity effects
#' 
#' @details: This script simulates a BEF experiment with full mixtures and monocultures
#' in each patch. It then uses the known starting abundances and full monoculture data
#' to calculate observed biodiversity effects using Isbell et al.'s (2018) partition. The 
#' model outputs and the observed biodiversity effects are outputted as .rds files. In addition, it generates a set of 100 starting relative 
#' abundances from the Dirichlet distribution. These outputs are all saved as .rds 
#' files. This script can be run locally on a regular desktop computer.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)

# load the relevant packages
library(dplyr)
library(here)

# load the metacommunity simulation and data-processing functions
source(here("scripts/02_simulation/01_mcomsimr_simulate_MC2_function.R"))

# set the model inputs

# set a seed
set.seed(54258748)

# number of replicate simulations
N_REP <- 1000

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

# run the model for each of the replicates
MC_sims <- 
  
  lapply(1:N_REP, function(x) {
    
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
    
    
    # pull the relevant outputs into a list
    output.list <-
      
      list("MC.x.ids" = MC.x.ids,
           "BEF_obs" = BEF_obs,
           "MC.x.NA" = MC.x.NA
           )
    
    return(output.list)
    
  })

# save the MC_sims object
saveRDS(object = MC_sims, here("results/MC_sims.rds"))

# generate and save a set of 100 potential starting relative abundances from the Dirichlet distribution
start_RA <- sapply(1:100, function(x) gtools::rdirichlet(n = 1, rep(3,  species) ) )
saveRDS(object = start_RA, here("results/MC_sims_start_RA.rds"))

### END
