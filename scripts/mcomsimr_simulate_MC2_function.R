
# Function to simulate BEF experiments using metacommunity models

# load the mcomsimr package
library(mcomsimr)

# load key functions
library(here)
source(here("scripts/isbell_2018_partition.R"))


#' Modified version of the simulate_MC function from the mcomsimr package (Thompson et al. 2020, Ecology Letters)
# https://github.com/plthompson/mcomsimr

# args

#' @param patches number of patches to simulate
#' @param species number of species to simulate
#' @param dispersal dispersal probability between 0 and 1
#' @param timesteps number of time-steps to run the model for
#' @param start_abun starting abundance of the community (each species starts with start_abun/species in each patch)
#' @param extirp_prob probability of local extirpation for each population in each time step (should be small e.g. 0.001)
#' @param landscape dataframe with x and y coordinates for each patch
#' @param disp_mat matrix with each column specifying the probability that an individual disperses to each other patch (row)
#' @param env.df optional dataframe with environmental conditions with columns: env1, patch, time
#' @param env_optima optional values of environmental optima, should be a vector of length species
#' @param int_mat optional externally generated competition matrix

simulate_MC2 <- function(patches, species, dispersal = 0.01, timesteps = 1200,
                         start_abun = 150,
                         extirp_prob = 0,
                         landscape, disp_mat, env.df, env_traits.df, int_mat, meas_error = 5){
  
  # load the dplyr library
  library(dplyr)
  
  dynamics.df <- data.frame()
  N <- matrix(rep(round(start_abun/species, 0), species*patches), nrow = patches, ncol = species)
  pb <- txtProgressBar(min = 0, max = timesteps, style = 3)
  for(i in 1:(timesteps)){
    
    # get the first environmental condition
    env <- env.df$env1[env.df$time == (i)]
    
    # we use the equation 3 to determine the realised growth rate
    # of each species (col) in each patch (row)
    r <- env_traits.df$max_r*exp(-(t((env_traits.df$optima - matrix(rep(env, each = species), nrow = species, ncol = patches))/(2*env_traits.df$env_niche_breadth)))^2)
    r <- r + rnorm(n = length(r), mean = 0, sd = 0.01)
    r <- ifelse(r < 0, 0, r)
    
    # we use equation 3 to determine the realised carrying capacity
    # of each species (col) in each patch (row)
    K <- env_traits.df$K_max*exp(-(t((env_traits.df$optima - matrix(rep(env, each = species), nrow = species, ncol = patches))/(2*env_traits.df$env_niche_breadth)))^2)
    K <- K + rnorm(n = length(K), mean = 0, sd = 2)
    K <- ifelse(K <= 0, 1, K)
    
    # here we implement the difference equation
    N_hat <- N + ((N*r) * (1 - ((N %*% int_mat)/K) ))
    
    # add noise from a normal distribution to these data
    N_hat <- N_hat + rnorm(n = length(N_hat), mean = 0, sd = meas_error)
    N_hat[N_hat < 0] <- 0
    N_hat <- round(N_hat, 0)
    
    E <- matrix(rbinom(n = patches * species, size = N_hat, prob = rep(dispersal, each = patches)), nrow = patches, ncol = species)
    dispSP <- colSums(E)
    I_hat_raw <- disp_mat%*%E
    I_hat <- t(t(I_hat_raw)/colSums(I_hat_raw))
    I_hat[is.nan(I_hat)] <- 1
    I <- sapply(1:species, function(x) {
      if(dispSP[x]>0){
        table(factor(sample(x = patches, size = dispSP[x], replace = TRUE, prob = I_hat[,x]), levels = 1:patches))
      } else {rep(0, patches)}
    })
    
    N <- N_hat - E + I
    N[rbinom(n = species * patches, size = 1, prob = extirp_prob)>0] <- 0
    
    dynamics.df <- rbind(dynamics.df, data.frame(N = c(N), patch = 1:patches, species = rep(1:species, each = patches), env = env, time = i))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # join the environmental traits to the dynamics
  dynamics.df <- dplyr::left_join(dynamics.df, env_traits.df, by = "species")
  
  # reorganise the dynamics.df
  dynamics.df <- 
    dynamics.df %>%
    filter(time >= 0) %>%
    select(time, patch, env, species, N) %>%
    arrange(time, patch, species)
  
  return(dynamics.df)
  
}


#' Function to use the simulate_MC2 function to simulate BEF experiments

# args

#' @param patches number of patches to simulate
#' @param species number of species to simulate
#' @param dispersal dispersal probability between 0 and 1
#' @param timesteps number of time-steps to run the model for
#' @param start_abun starting abundance of the community (each species starts with start_abun/species in each patch)
#' @param extirp_prob probability of local extirpation for each population in each time step (should be small e.g. 0.001)
#' @param landscape dataframe with x and y coordinates for each patch
#' @param disp_mat matrix with each column specifying the probability that an individual disperses to each other patch (row)
#' @param env.df optional dataframe with environmental conditions with columns: env1, patch, time
#' @param env_optima optional values of environmental optima, should be a vector of length species
#' @param int_mat optional externally generated competition matrix

# simulate monocultures and mixtures in identical environmental conditions
sim_metacomm_BEF <- function(species = 5, patches = 10,
                             dispersal = 0.2,
                             timesteps = 10, 
                             start_abun = 150,
                             extirp_prob = 0.0001,
                             landscape, 
                             disp_mat, 
                             env.df, 
                             env_traits.df, 
                             int_mat,
                             meas_error = 5
                             ) {
  
  # simulate the mixture of species
  mix <- simulate_MC2(species = species, 
                      patches = patches, 
                      dispersal = dispersal, 
                      start_abun = start_abun,
                      timesteps = timesteps,
                      extirp_prob = extirp_prob,
                      landscape = landscape, 
                      disp_mat = disp_mat, 
                      env.df = env.df, 
                      env_traits.df = env_traits.df, 
                      int_mat = int_mat,
                      meas_error = meas_error
  )
  
  # add a mixture column
  mix$mono_mix <- "mixture"
  
  # add a column for each unique sample
  sample <- unique(with(mix, paste(time, patch)))
  mix$sample <- rep(sort(as.integer(as.factor(sample)) ), each = species)
  
  # reorder the columns
  mix <- mix[, c("mono_mix", "sample", "time", "patch", "env", "species", "N")]
  
  # simulate each monoculture over all times and places
  mono <- vector("list", length = species)
  for (i in 1:species) {
    
    # simulate each species
    x <- simulate_MC2(species = 1, patches = patches, 
                      dispersal = dispersal, start_abun = start_abun,
                      timesteps = timesteps,
                      extirp_prob = extirp_prob,
                      landscape = landscape, 
                      disp_mat = disp_mat, 
                      env.df = env.df, 
                      env_traits.df = env_traits.df[i,], 
                      int_mat = int_mat[i,i],
                      meas_error = meas_error)
    
    # rename the species column
    x$species <- i
    
    # add a mono_mix variable
    x$mono_mix <- "monoculture"
    
    # add a column for each unique sample
    sample <- unique(with(x, paste(time, patch)))
    x$sample <- sort(as.integer(as.factor(sample)) )
    
    # write the dynamics data.frame into a list
    mono[[i]] <- x
    
  }
  
  # bind the monocultures into a data.frame
  mono <- bind_rows(mono)
  
  # reorder the columns
  mono <- mono[, c("mono_mix", "sample", "time", "patch", "env", "species", "N")]
  mono <- 
    mono %>%
    arrange(mono_mix, time, patch, species)
  
  # return the list with the mixture and monoculture data
  return(list("mixture" = mix, "monoculture" = mono))
  
}


# Function to process mixture and monoculture data from the simulations and calculate Isbell et al.'s (2018, Ecology Letters) partition

# args

#' @param mix mixture data from the sim_metacomm_BEF function
#' @param mono monoculture data from the sim_metacomm_BEF function
#' @param from_last how many time points from the last one to include

Isbell_2018_cleaner <- function(mix, mono, t_sel) {
  
  library(dplyr)
  
  # if no time-points are provided, then use the five in the simulation
  if(any(is.na(t_sel)) | missing(t_sel)) {
    nt <- length(unique(mix$time))
    t_sel <- (nt-5):nt
  }
  
  # process the mixture data
  mix <- 
    mix %>%
    filter( time %in% t_sel ) %>%
    select(-mono_mix, -env) %>%
    rename(place = patch, Y = N) %>%
    mutate(sample = as.integer(as.factor(sample)))
  
  # process the monoculture data
  mono <- 
    mono %>%
    filter( time %in% t_sel ) %>%
    select(-mono_mix, -env, -sample) %>%
    rename(place = patch, M = N)
  
  # join the mixture and monoculture data to match the partition format
  mix.mono <- 
    full_join(mono, mix, by = c("time", "place", "species")) %>%
    select(sample, time, place, species, M, Y) %>%
    arrange(sample, time, place, species)
  head(mix.mono)
  
  return(mix.mono)
  
}

### END
