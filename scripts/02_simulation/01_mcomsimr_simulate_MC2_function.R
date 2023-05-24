
#' @title Functions to simulate BEF experiments using metacommunity models

#' @description Modified version of the simulate_MC function from the mcomsimr package (Thompson et al. 2020, Ecology Letters)
# https://github.com/plthompson/mcomsimr

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
  
  # check the starting abundances
  if(length(start_abun) == 1) {
    SA <- rep(round(start_abun/species, 0), species*patches)
  } else if (length(start_abun) == species) {
    SA <- rep(start_abun, each = patches)
  } else {
    stop("add appropriate starting abundances")
  }
  
  dynamics.df <- data.frame()
  N <- matrix(round(SA, 0), nrow = patches, ncol = species)
  pb <- txtProgressBar(min = 0, max = timesteps, style = 3)
  for(i in 1:(timesteps)){
    
    # get the first environmental condition
    env <- env.df$env1[env.df$time == (i)]
    
    # we use the equation 3 to determine the realised growth rate
    # of each species (col) in each patch (row)
    r <- env_traits.df$max_r*exp(-(t((env_traits.df$optima - matrix(rep(env, each = species), nrow = species, ncol = patches))/(2*env_traits.df$env_niche_breadth)))^2)
    r <- ifelse(r < 0, 0, r)
    
    # here we implement the difference equation
    N_hat <- N*r/(1+N%*%int_mat)
    
    # add noise from a poisson distribution to the data
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
                             comp = "Beverton_Holt",
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
  mix$sample <- rep((as.integer(as.factor(sample)) ), each = species)
  
  # reorder the columns
  mix <- mix[, c("mono_mix", "sample", "time", "patch", "env", "species", "N")]
  
  # simulate each monoculture over all times and places
  mono <- vector("list", length = species)
  for (j in 1:species) {
    
    # simulate each species
    x <- simulate_MC2(species = 1, patches = patches, 
                      dispersal = dispersal, start_abun = sum(start_abun),
                      timesteps = timesteps,
                      extirp_prob = extirp_prob,
                      landscape = landscape, 
                      disp_mat = disp_mat, 
                      env.df = env.df, 
                      env_traits.df = env_traits.df[j,], 
                      int_mat = int_mat[j,j],
                      meas_error = meas_error)
    
    # rename the species column
    x$species <- j
    
    # add a mono_mix variable
    x$mono_mix <- "monoculture"
    
    # add a column for each unique sample
    sample <- unique(with(x, paste(time, patch)))
    x$sample <- (as.integer(as.factor(sample)) )
    
    # write the dynamics data.frame into a list
    mono[[j]] <- x
    
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


#' Generate Dispersal Matrix
#'
#' Generates dispersal matrix for metacommunity simulations (code taken directly)
#' from https://github.com/plthompson/mcomsimr
#'
#' @param landscape landscape generated by landscape_generate()
#' @param torus whether to model the landscape as a torus
#' @param disp_mat optional matrix with each column specifying the probability that an individual disperses to each other patch (row)
#' @param kernel_exp the exponential rate at which dispersal decreases as a function of the distance between patches
#' @param plot option to show plot of environmental variation
#'
#' @return matrix with dispersal probabilities
#'
#' @author Patrick L. Thompson, \email{patrick.thompson@@zoology.ubc.ca}
#'
#' @examples
#' dispersal_matrix(landscape_generate())
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom som.nn dist.torus
#'
#' @export
#'
dispersal_matrix <- function(landscape, torus = TRUE, disp_mat, kernel_exp = 0.1, plot = TRUE){
  if (missing(disp_mat)){
    if(torus == TRUE){
      dist_mat <- as.matrix(som.nn::dist.torus(coors = landscape))
    } else {
      dist_mat <- as.matrix(dist(landscape))
    }
    
    disp_mat <- exp(-kernel_exp * dist_mat)
    diag(disp_mat) <- 0
    disp_mat <- apply(disp_mat, 1, function(x) x / sum(x))
  } else {
    disp_mat <- disp_mat
    rownames(disp_mat) <- 1:nrow(disp_mat)
    colnames(disp_mat) <- 1:ncol(disp_mat)
    if (is.matrix(disp_mat) == FALSE) stop ("disp_mat is not a matrix")
    if (nrow(disp_mat) != nrow(landscape) | ncol(disp_mat) != nrow(landscape)) stop ("disp_mat does not have a row and column for each patch in landscape")
  }
  
  if (sum(colSums(disp_mat) > 1.001) > 0) warning ("dispersal from a patch to all others exceeds 100%. Make sure the rowSums(disp_mat) <= 1")
  if (sum(colSums(disp_mat) < 0.999) > 0) warning ("dispersal from a patch to all others is less than 100%. Some dispersing individuals will be lost from the metacommunity")
  
  if (plot == TRUE){
    g <- as.data.frame(disp_mat) %>%
      dplyr::mutate(to.patch = rownames(disp_mat)) %>%
      tidyr::gather(key = from.patch, value = dispersal, -to.patch) %>%
      dplyr::mutate(from.patch = as.numeric(as.character(from.patch)),
                    to.patch = as.numeric(as.character(to.patch))) %>%
      ggplot2::ggplot(ggplot2::aes(x = from.patch, y = to.patch, fill = dispersal))+
      ggplot2::geom_tile()+
      scale_fill_viridis_c()
    
    print(g)
  }
  
  return (disp_mat)
}

#' Generate Environment
#'
#' Generates density independent environmental conditions for metacommunity simulations
#' code taken directly from https://github.com/plthompson/mcomsimr
#'
#' @param landscape landscape generated by landscape_generate()
#' @param env.df optional dataframe with environmental conditions with columns: env1, patch, time
#' @param env1Scale scale of temporal environmental autocorrelation between -2 (anticorrelated) and 2 (correlated), default is 2
#' @param timesteps number of timesteps to simulate
#' @param plot option to show plot of environmental variation
#'
#' @return dataframe with x and y coordinates, time, and environmental conditions
#'
#' @author Patrick L. Thompson, \email{patrick.thompson@@zoology.ubc.ca}
#'
#' @examples
#' env_generate(landscape_generate())
#'
#' @importFrom synchrony phase.partnered
#' @import ggplot2
#'
#' @export
#'
env_generate <- function(landscape, env.df, env1Scale = 2, timesteps = 1000, plot = TRUE){
  if (missing(env.df)){
    repeat {
      env.df <- data.frame()
      for(i in 1:nrow(landscape)){
        env1 = synchrony::phase.partnered(n = timesteps, gamma = env1Scale, mu = 0.5, sigma = 0.25)$timeseries[,1]
        env.df <- rbind(env.df, data.frame(env1 = vegan::decostand(env1,method = "range"), patch = i, time = 1:timesteps))
      }
      env.initial <- env.df[env.df$time == 1,]
      if((max(env.initial$env1)-min(env.initial$env1)) > 0.6) {break}
    }
  } else {
    if(all.equal(names(env.df), c("env1", "patch", "time")) != TRUE) stop("env.df must be a dataframe with columns: env1, patch, time")
  }
  
  if(plot == TRUE){
    g<-ggplot2::ggplot(env.df, aes(x = time, y = env1, group = patch, color = factor(patch)))+
      ggplot2::geom_line()+
      scale_color_viridis_d(guide = "none")
    print(g)
  }
  
  return(env.df)
  print("This version differs from Thompson et al. 2020 in that it does not produce spatially autocorrelated environmental variables.")
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
