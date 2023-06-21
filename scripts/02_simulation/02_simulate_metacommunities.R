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

# load the metacommunity simulation and data-processing functions
source("scripts/02_simulation/01_mcomsimr_simulate_MC2_function.R")

# set the model inputs

# set a seed
set.seed(54258748)

# number of replicate simulations
N_REP <- 100

# set-up basic model inputs
species <- 5
patches <- 12
timesteps <- 500
dispersal_min <- 0.01
dispersal_max <- 0.1
extirp_prob <- 0.00001

# maximum growth parameters
max_r <- 5

# set-up parameters for the competition matrices
int_min1 <- 0.05*0
int_min2 <- 0.05*0.5
int_max1 <- 0.05*1
int_max2 <- 0.05*1.5
intra <- 0.05*1

# landscape parameters

# generate a random landscape
n1 <- 4
n2 <- 3
assertthat::assert_that(assertthat::see_if( (n1*n2) == patches ),
                        msg = "change the n1 and n2 values so that their product matches the patches")

# using the number of patches, set-up a random evenly spread landscape
t1 <- round(seq(25, 25*4, length.out = n1), 0)
t2 <- round(seq(25, 25*5, length.out = n2), 0)
t12 <- expand.grid(t1, t2)

# pull these evenly spaced patches into a data.frame
landscape_x <- data.frame(x = t12[[1]], y = t12[[2]])
plot(landscape_x)

# generate a random dispersal matrix using the landscape.x landscape
dispersal_x <- dispersal_matrix(landscape = landscape_x, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

# generate species environmental optima
optima <- seq(0.15, 0.85, length.out = species)
print(optima)
env_niche_breadth <- lapply(1:N_REP, function(x) round(runif(species, 0.05, 0.25), 2))
print(env_niche_breadth)

# select time-points to calculate biodiversity effects for
t_sel <- c(100, 200, 300, 400, 500)

# set the min and maximum starting abundances
min_SA <- 5
max_SA <- 150

# run the model for each of the replicates
MC_sims <- 
  
  lapply(1:N_REP, function(x) {
    
    # get the starting abundances
    start_abun <- sim_start_abun(species = species, patches = patches,
                                 min_SA = min_SA, max_SA = max_SA)
    
    # pull species attributes into a data.frame
    sp_att <- 
      
      data.frame(species = 1:species,
                 optima = optima,
                 env_niche_breadth = env_niche_breadth[[x]],
                 max_r = max_r
      )
    
    # generate the competition matrices
    comp_mat <- 
      
      matrix(
        
        runif(n = species*species, 
              min = runif(1, int_min1, int_min2),
              max = runif(1, int_max1, int_max2) ),  
        nrow = species, 
        ncol = species
        
      )
    
    # transpose the matrix so that competition is symmetrical
    comp_mat[lower.tri(comp_mat)] = t(comp_mat)[lower.tri(comp_mat)]
    
    # set intra-specific competition to be equal
    diag(comp_mat) <- intra
    
    # simulate environmental variables
    
    # randomly draw an autocorrelation value between 300 and 500
    autocorr <- round(runif(1, 100, 500), 0)
    
    # use the env_generate() function to generate a landscape of environmental variables
    env_var <- 
      
      env_generate(landscape = landscape_x, 
                   env1Scale = autocorr, 
                   timesteps = timesteps, 
                   plot = FALSE
                  )
    
    # draw a random dispersal value
    dispersal <- runif(1, dispersal_min, dispersal_max)
    
    # simulate the metacommunity
    MC_x <- 
      
      sim_metacomm_BEF(patches = patches, species = species, 
                       dispersal = dispersal,
                       timesteps = timesteps, 
                       start_abun = start_abun,
                       extirp_prob = extirp_prob,
                       landscape = landscape_x, disp_mat = dispersal_x, env.df = env_var, 
                       env_traits.df = sp_att, 
                       int_mat = comp_mat
      )
    
    # output the model identifiers
    MC_x_ids <-
      
      data.frame("start_abun" = paste(start_abun, collapse = "_"),
                 "inter_comp" = mean(comp_mat[lower.tri(comp_mat)]),
                 "optima" = paste(sp_att$optima, collapse = "_"),
                 "niche_breadth" = paste(sp_att$env_niche_breadth, collapse = "_"),
                 "t_steps" = paste(t_sel, collapse = "_"),
                 "dispersal" = dispersal
      )
    
    # clean the data.frame to calculate BEF effects
    MC_x_BEF <- 
      
      Isbell_2018_cleaner(mix = MC_x$mixture, 
                          mono = MC_x$monoculture, 
                          t_sel  = t_sel
      )
    
    # convert starting abundances to a matrix reorder
    SA_mat <- matrix(start_abun, nrow = species, ncol = patches)
    SA_mat <- SA_mat[, as.integer(unique(MC_x_BEF$place))]
    
    # convert to proportions
    SA_mat <- apply(SA_mat, 2, function(x) round(x/sum(x), 2) )
    
    # convert the rows to a list
    RYE <- vector("list", length = ncol(SA_mat))
    for(i in 1:ncol(SA_mat)) {
      RYE[[i]] <- SA_mat[,i]
    }
    
    # replicate the rows across the time-points
    RYE_list <- RYE
    for(i in 2:length(t_sel)) {
      RYE_list <- c(RYE_list, RYE)
    }
    
    # get the environmental data
    MC_x_env <- 
      MC_x$mixture %>%
      filter(time %in% t_sel) %>%
      select(time, place = patch, env) %>%
      distinct()
    
    # calculate the observed biodiversity effects
    BEF_obs <- Isbell_2018_part(data = MC_x_BEF, RYe_equal = FALSE, RYe = RYE_list )
    
    BEF_obs <- 
      
      bind_rows(
        
        BEF_obs$Beff, 
        rename(BEF_obs$L.Beff, Beff = L.Beff) 
        
      ) %>%
      
      mutate(Value = round(Value, 2))
    
    # assume incomplete monoculture and unknown RYE
    MC_dat <- 
      
      full_join(MC_x_BEF, 
                MC_x_env, 
                by = c("time", "place")) %>%
      
      arrange(sample, time, place, species)
    
    
    # pull the relevant outputs into a list
    output_list <-
      
      list("MC_x_ids" = MC_x_ids,
           "BEF_obs" = BEF_obs,
           "MC_dat" = MC_dat
           )
    
    return(output_list)
    
  })

# save the MC_sims object
saveRDS(object = MC_sims, "results/MC_sims.rds")

# get an example dataset to work with for the dimensions
df_ex <- MC_sims[[1]]$MC_dat

# get the a vector of cluster names
samples <- length(unique(df_ex$sample))

# set the number of unique RYe values
Ns <- 100

start_RA <- vector("list", length = Ns)
for(i in 1:Ns) {

  # for each sample, get a simplex from the Dirichlet distribution
  rye_mat <- round(gtools::rdirichlet(n = samples, rep(3, species)), 2)
  
  # write into a list
  start_RA[[i]] <- rye_mat
  
}

# save the Dirichlet distribution
saveRDS(object = start_RA, "results/MC_sims_start_RA.rds")

### END
