
# Modified version of the simulate_MC function

simulate_MC2 <- function(patches, species, dispersal = 0.01, timesteps = 1200,
                         start_abun = 500,
                         extirp_prob = 0,
                         landscape, disp_mat, env.df, env_traits.df, int_mat){
  
  library(dplyr)
  library(ggplot2)
  
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
  dynamics.df <- dplyr::left_join(dynamics.df, env_traits.df, by = "species")
  dynamics.df <- 
    dynamics.df %>%
    filter(time >= 0) %>%
    select(time, patch, env, species, N)
  
  return(dynamics.df)
  
}

### END
