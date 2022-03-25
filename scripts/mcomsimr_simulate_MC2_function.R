
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
                         landscape, disp_mat, env.df, env_traits.df, int_mat){
  
  # load the dplyr library
  library(dplyr)
  
  # check if starting abundance is reasonable
  if (any(env_traits.df$K_max > start_abun)) {
    stop("Starting abundance cannot exceed the maximum carrying capacity, K")
  }
  
  dynamics.df <- data.frame()
  N <- matrix(rep(round(start_abun/species, 0), species*patches), nrow = patches, ncol = species)
  pb <- txtProgressBar(min = 0, max = timesteps, style = 3)
  for(i in 1:(timesteps)){
    
    # get the first environmental condition
    env <- env.df$env1[env.df$time == (i)]
    
    # we use the equation 3 to determine the realised growth rate
    # of each species (col) in each patch (row)
    r <- env_traits.df$max_r*exp(-(t((env_traits.df$optima - matrix(rep(env, each = species), nrow = species, ncol = patches))/(2*env_traits.df$env_niche_breadth)))^2)
    
    # we use equation 3 to determine the realised carrying capacity
    # of each species (col) in each patch (row)
    K <- env_traits.df$K_max*exp(-(t((env_traits.df$optima - matrix(rep(env, each = species), nrow = species, ncol = patches))/(2*env_traits.df$env_niche_breadth)))^2)
    
    # here we implement the difference equation
    N_hat <- N + ((N*r) * (1 - ((N %*% int_mat)/K) ))
    
    N_hat[N_hat < 0] <- 0
    N_hat <- matrix(rpois(n = species*patches, lambda = N_hat), ncol = species, nrow = patches)
    
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
    select(time, patch, env, species, N)
  
  return(dynamics.df)
  
}

### END
