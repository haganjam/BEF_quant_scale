
# Modified version of the simulate_MC function

simulate_MC2 <- function(patches, species, dispersal = 0.01,
                         plot = TRUE, timesteps = 1200, burn_in = 800, initialization = 200,
                         extirp_prob = 0,
                         landscape, disp_mat, env.df, env_traits.df, int_mat){
  
  library(dplyr)
  library(ggplot2)
  
  dynamics.df <- data.frame()
  N <- matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
  pb <- txtProgressBar(min = 0, max = initialization + burn_in + timesteps, style = 3)
  for(i in 1:(initialization + burn_in + timesteps)){
    if(i <= initialization){
      if(i %in% seq(10,100, by = 10)){
        N <- N + matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
      }
      env <- env.df$env1[env.df$time == 1]
    } else {
      env <- env.df$env1[env.df$time == (i-initialization)]
    }
    r <- env_traits.df$max_r[1]*exp(-(t((env_traits.df$optima - matrix(rep(env, each = species), nrow = species, ncol = patches))/(2*env_traits.df$env_niche_breadth)))^2)
    N_hat <- N*r/(1+N%*%int_mat)
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
    
    dynamics.df <- rbind(dynamics.df, data.frame(N = c(N), patch = 1:patches, species = rep(1:species, each = patches), env = env, time = i-initialization-burn_in))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dynamics.df <- dplyr::left_join(dynamics.df, env_traits.df, by = "species")
  env.df$time_run <- env.df$time - burn_in
  
  env.df_init <- data.frame(env1 = env.df$env1[env.df$time == 1], patch = 1:patches, time = NA, time_run = rep(seq(-(burn_in + initialization), -burn_in, by = 1), each = patches))
  env.df <- rbind(env.df_init,env.df)
  
  if(plot == TRUE){
    sample_patches <- sample(1:patches, size = min(c(patches,6)), replace = FALSE)
    g <- dynamics.df %>%
      filter(time %in% seq(min(dynamics.df$time),max(dynamics.df$time), by =10)) %>%
      filter(patch %in% sample_patches) %>%
      ggplot(aes(x = time, y = N, group = species, color = optima))+
      geom_line()+
      facet_wrap(~patch)+
      scale_color_viridis_c()+
      geom_path(data = filter(env.df, patch %in% sample_patches), aes(y = -5, x = time_run, color = env1, group = NULL), size = 3)
    
    print(g)
  }
  
  return(list(dynamics.df = dynamics.df, landscape = landscape, env.df = env.df, env_traits.df = env_traits.df, disp_mat = disp_mat, int_mat = int_mat))
}
