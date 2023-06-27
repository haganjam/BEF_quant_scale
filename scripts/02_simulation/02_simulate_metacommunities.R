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

# load the metacommunity simulation and data-processing functions
source("scripts/02_simulation/01_mcomsimr_simulate_MC2_function.R")
source("scripts/01_partition_functions/01_isbell_2018_partition.R")

# set the model inputs

# load the simulated environments
env_SI <- readRDS("scripts/02_simulation/02_generated_environment_SI.rds")
env_TI <- readRDS("scripts/02_simulation/02_generated_environment_SI.rds")
env_com <- readRDS("scripts/02_simulation/02_generated_environment_com.rds")

# pull the simulated environments into a list
env_list <- list(env_SI, env_TI, env_com)

# set a seed
set.seed(54258748)

# number of replicate simulations
N_REP <- 500

assertthat::assert_that(
  N_REP == length(env_SI) && N_REP == length(env_SI) && N_REP == length(env_com)
)

# set-up basic model inputs

# number of species
species <- 3

# number of patches
patches <- 5

assertthat::assert_that(
  patches == length(unique(env_SI[[1]]$patch))
)

# number of time-steps
timesteps <- 300

assertthat::assert_that(
  timesteps == length(unique(env_SI[[1]]$time))
)

# set the minimum and maximum dispersal rate
dispersal_min <- 0.01
dispersal_max <- 0.1

# set the stochastic extinction probability
extirp_prob <- 0.00001

# select time-points to calculate biodiversity effects for
t_sel <- c(100, 200, 300)

# species niche attributes
sp_att <- data.frame(species = 1:species,
                     optima = seq(0.3, 0.7, length.out = species),
                     env_niche_breadth = c(0.3, 0.3, 0.3),
                     max_r = 5
)

# set-up parameters for the competition matrices
int_min1 <- 0.05*0
int_min2 <- 0.05*0.5
int_max1 <- 0.05*1
int_max2 <- 0.05*1.5
intra <- 0.05*1

# landscape parameters

# pull these evenly spaced patches into a data.frame
landscape_x <- data.frame(x = round(seq(25, 25*5, length.out = patches), 0), 
                          y = 25)
plot(landscape_x)

# generate a random dispersal matrix using the landscape.x landscape
dispersal_x <- dispersal_matrix(landscape = landscape_x, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

# apply over the three different simulated environmental datasets
env_sim_list <- 
  lapply(env_list, function(env_dat) {
  
  # run N_REP simulations
  sim_list <- vector("list", length = N_REP)
  for(i in 1:N_REP) {
    
    # get the starting abundances
    start_abun <- sim_start_abun(species = species, patches = patches,
                                 min_SA = 20, max_SA = 130, tot_N = 150)
    
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
    
    # get the environmental variables
    env_var <- env_dat[[i]]
    
    # draw a random dispersal value
    dispersal <- runif(1, dispersal_min, dispersal_max)
    
    # output the model identifiers
    mc_sim_pars <- data.frame("start_abun" = paste(start_abun, collapse = "_"),
                              "inter_comp" = mean(comp_mat[lower.tri(comp_mat)]),
                              "optima" = paste(sp_att$optima, collapse = "_"),
                              "niche_breadth" = paste(sp_att$env_niche_breadth, collapse = "_"),
                              "t_steps" = paste(t_sel, collapse = "_"),
                              "dispersal" = dispersal
    )
    
    # simulate the metacommunity
    mc_sim <- 
      sim_metacomm_BEF(patches = patches, species = species, 
                       dispersal = dispersal,
                       timesteps = timesteps, 
                       start_abun = start_abun,
                       extirp_prob = extirp_prob,
                       landscape = landscape_x, disp_mat = dispersal_x, env.df = env_var, 
                       env_traits.df = sp_att, 
                       int_mat = comp_mat
      )
    
    # get the environmental data
    mc_env <- 
      mc_sim$mixture |>
      dplyr::filter(time %in% t_sel) |>
      dplyr::select(time, patch, env) |>
      dplyr::rename(place = patch)
    
    # remove duplicates using distinct
    mc_env <- dplyr::distinct(mc_env)
    
    # clean the data.frame to calculate BEF effects
    mc_sim_BEF <- isbell_2018_cleaner(mix = mc_sim$mixture, 
                                      mono = mc_sim$monoculture, 
                                      t_sel  = t_sel)
    
    # join the monoculture mixture data and environmental data
    mc_dat <- dplyr::full_join(mc_sim_BEF, mc_env, by = c("time", "place"))
    
    # arrange the columns appropriately  
    mc_dat <- dplyr::arrange(mc_dat, sample, time, place, species)
    
    # convert starting abundances to a matrix and reorder
    SA_mat <- matrix(start_abun, nrow = species, ncol = patches)
    SA_mat <- SA_mat[, as.integer(unique(mc_sim_BEF$place))]
    
    # convert to proportions by dividing by sum across species
    SA_mat <- apply(SA_mat, 2, function(x) round(x/sum(x), 2) )
    
    # convert the rows to a list
    RYE <- vector("list", length = ncol(SA_mat))
    for(j in 1:ncol(SA_mat)) {
      RYE[[j]] <- SA_mat[,j]
    }
    
    # replicate the rows across the time-points
    RYE_list <- RYE
    for(k in 2:length(t_sel)) {
      RYE_list <- c(RYE_list, RYE)
    }
    
    # calculate the observed biodiversity effects
    BEF_obs <- isbell_2018_part(data = mc_sim_BEF, RYe_equal = FALSE, RYe = RYE_list )
    
    # bind local and regional biodiversity effects
    BEF_obs <- dplyr::bind_rows(BEF_obs$Beff, 
                                dplyr::rename(BEF_obs$L.Beff, Beff = L.Beff))
    
    # round off the BEF effects
    BEF_obs <- dplyr::mutate(BEF_obs, Value = round(Value, 2))
    
    # pull the relevant outputs into a list
    output_list <- list("pars" = mc_sim_pars,
                        "BEF_obs" = BEF_obs,
                        "mc_dat" = mc_dat)
    
    # add to the list collecting the simulations
    sim_list[[i]] <- output_list
    
  }
  
  return(sim_list)
  
} )

# link the lists together
env_sim_list <- do.call("c", env_sim_list)

# check that the list is the correct length
length(env_sim_list) == 3*N_REP

# save the MC_sims object
saveRDS(object = env_sim_list, "results/mc_sim_list.rds")

# get an example dataset to work with for the dimensions
df_ex <- env_sim_list[[1]]$mc_dat

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
saveRDS(object = start_RA, "results/mc_sim_list_start_RA.rds")

### END
