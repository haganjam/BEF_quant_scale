
# Simulate a test of our hypothesis

# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))

# Generate high and low heterogeneity clusters

# simulate a set of metacommunities with high environmental heterogeneity among patches

Simulate_heterogeneity <- function(patches, timesteps, 
                                   het = "low", env_min, env_max, env,
                                   temp_fluc) {
  
  # set the starting enviromental heterogeneity among patches
  if (het == "high") {
    start_het <- seq(env_min, env_max, length.out = patches)
  } else if (het == "low") {
    start_het <- rep(env, patches)
  } else {
    stop("Choose appropriate het value: high or low")
  }
  
  # simulate an environment with high heterogeneity
  env1.sim <- sapply(start_het, function(x) {
    
    # choose a random number between 0 and 0.025 and then decide if it is negative or positive
    y <- runif(n = timesteps-1, 0, temp_fluc)*sample(x = c(-1,1), size = timesteps-1, replace = TRUE)
    
    # cumulatively sum the environmental fluctuations
    z <- cumsum(c(x, y))
    
    # make sure they remain between zero and one
    return(ifelse( (z < 0), 0, ifelse( (z > 1), 1, z)))
    
  })
  
  # bind into a data.frame
  df <- data.frame(env1 = unlist(apply(env1.sim, 1, function(x) list(x) ) ),
                   patch = rep(1:patches, timesteps),
                   time = rep(1:timesteps, each = patches)
                  )
  
  return(df)
  
}

# test the Simulate_heterogeneity function
Simulate_heterogeneity(patches = 3, timesteps = 50, 
                       het = "low", env_min = 0.2, env_max = 0.8, env = 0.5,
                       temp_fluc = 0.025)


# Environmental heterogeneity experiment function
Simulate_hetero_exp <- function(n_rep = 5, temp_fluc = 0.025,
                        env_min = 0.2, env_max = 0.8, t_sel,
                        patches, species, dispersal = 0.05, 
                        timesteps, start_abun = 60,
                        extirp_prob = 0,
                        landscape, disp_mat, 
                        env_traits.df, int_mat) {
  
  MC_high <- lapply(1:n_rep, function(x) {
    
    # test the Simulate_heterogeneity function
    e.1 <- 
      Simulate_heterogeneity(patches = patches, timesteps = timesteps, 
                             het = "high", env_min = env_min, env_max = env_max,
                             temp_fluc = temp_fluc)
    
    # simulate a metacommunity in this cluster
    MC1 <- sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal, 
                            timesteps = ts, start_abun = start_abun,
                            extirp_prob = extirp_prob,
                            landscape = landscape, disp_mat = disp_mat, env.df = e.1, 
                            env_traits.df = env_traits.df, int_mat = int_mat)
    
    MC1 <- Isbell_2018_cleaner(mix = MC1$mixture, mono = MC1$monoculture, 
                               t_sel  = t_sel )
    
    # add a low heterogeneity variable
    MC1$env_het <- "high"
    
    # convert the time variable to an integer
    MC1$time <- as.integer(as.factor(MC1$time))
    
    return(MC1)
    
  } )
  

  MC_low <- lapply(seq(env_min, env_max, length.out = n_rep), function(x) {
    
    # test the Simulate_heterogeneity function
    e.2 <- 
      Simulate_heterogeneity(patches = patches, timesteps = timesteps, 
                             het = "low", env = x,
                             temp_fluc = temp_fluc)
    
    # simulate a metacommunity in this cluster
    MC2 <- sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal, 
                            timesteps = ts, start_abun = start_abun,
                            extirp_prob = extirp_prob,
                            landscape = landscape, disp_mat = disp_mat, env.df = e.2, 
                            env_traits.df = env_traits.df, int_mat = int_mat
                            )
    
    
    MC2 <- Isbell_2018_cleaner(mix = MC2$mixture, mono = MC2$monoculture, 
                               t_sel  = t_sel )
    
    # add a low heterogeneity variable
    MC2$env_het <- "low"
    
    # convert the time variable to an integer
    MC2$time <- as.integer(as.factor(MC2$time))
    
    return(MC2)
    
  } )

  return(c(MC_high, MC_low))
  
}


## Explore model behaviour

# set-up fixed inputs
species <- 3
patches <- 3
ts <- 20
n_rep <- 5

# landscape parameters
# generate a random landscape
l.1 <- 
  data.frame(x = c(25, 50, 75),
             y = c(50, 50, 50))

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

# Exp 1: Narrow niches, neutral competition

# generate species environmental optima
t.1 <- 
  data.frame(species = 1:species,
             optima = seq(0.33, 0.66, length.out = species),
             env_niche_breadth = c(0.1, 0.2, 0.4),
             max_r = 0.5,
             K_max = 150)

# competition matrices
si.1 <- matrix(runif(n = species*species, min = 0.8, max = 0.8), 
               nrow = species, ncol = species)
si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
diag(si.1) <- 1
head(si.1)

# simulate the experiment
Exp1 <- Simulate_hetero_exp(n_rep = n_rep, 
                            temp_fluc = 0.025,
                            env_min = 0.33, env_max = 0.66, t_sel = round(seq(5, ts, length.out = 3),0),
                            patches, species, dispersal = 0.05, 
                            timesteps = ts, start_abun = 60,
                            extirp_prob = 0,
                            landscape = l.1, disp_mat = d.1, 
                            env_traits.df = t.1, int_mat = si.1
                            )

# calculate biodiversity effects for each of these replicate clusters
df_obs <- lapply(Exp1, function(x) {
  
  y <- Isbell_2018_sampler(data = x[, c("sample", "time", "place", "species", "M", "Y")],
                           RYe_post = FALSE, RYe = rep(1/species, species)
                          )
  return(y$Beff)
  
})

# bind rows and convert to wide format
df_obs <- 
  df_obs %>%
  bind_rows(., .id = "ID") %>%
  mutate(exp_het = rep(c("high", "low"), each = n_rep*9),
         exp_rep = rep(rep(c(1:n_rep), each = 9) , 2)) %>%
  select(ID, exp_het, exp_rep, Beff, Value) %>%
  pivot_wider(id_cols = c("ID", "exp_het", "exp_rep"),
              names_from = Beff,
              values_from = Value)
(df_obs)

# plot the data: Net biodiversity effect
ggplot() +
  geom_jitter(data = df_obs,
             mapping = aes(x = exp_het, y = NBE), width = 0.1) +
  geom_point(data = df_obs %>% group_by(exp_het) %>% summarise(NBE = mean(NBE)),
             mapping = aes(x = exp_het, y = NBE), colour = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

# plot the data: Proportion spatial insurance
ggplot() +
  geom_jitter(data = df_obs %>% mutate(SI_prop = SI/NBE),
              mapping = aes(x = exp_het, y = SI_prop), width = 0.1) +
  geom_point(data = df_obs %>% group_by(exp_het) %>% summarise(SI_prop = mean(SI/NBE)),
             mapping = aes(x = exp_het, y = SI_prop), colour = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

# why is spatial insurance greater than net biodiversity?
Exp1[[1]] %>%
  group_by(sample) %>%
  summarise(M_m = mean(M),
            max_M = max(M),
            Y_m = sum(Y))

x <- 
  Exp1[[1]] %>%
  group_by(sample, time, place) %>%
  mutate(RA = Y/sum(Y)) %>%
  group_by(place, species) %>%
  summarise(RA = mean(RA),
            M = mean(M))

plot(x[x$species == 1,]$RA, x[x$species == 1,]$M) 
plot(x[x$species == 2,]$RA, x[x$species == 2,]$M) 
plot(x[x$species == 3,]$RA, x[x$species == 3,]$M)


## Assume incomplete data

# can we recover the observed values sufficiently to make inference






