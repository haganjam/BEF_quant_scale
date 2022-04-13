
# Simulate a test of our hypothesis

# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
library(pbapply)
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
    names(e.1) <- c("env1", "place", "time")
    
    # simulate a metacommunity in this cluster
    MC1 <- sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal, 
                            timesteps = ts, start_abun = start_abun,
                            extirp_prob = extirp_prob,
                            landscape = landscape, disp_mat = disp_mat, env.df = e.1, 
                            env_traits.df = env_traits.df, int_mat = int_mat)
    
    MC1 <- Isbell_2018_cleaner(mix = MC1$mixture, mono = MC1$monoculture, 
                               t_sel  = t_sel )
    
    # add the environmental variables to these data
    MC1 <- left_join(MC1, e.1, by = c("place", "time"))
    
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
    names(e.2) <- c("env1", "place", "time")
    
    # simulate a metacommunity in this cluster
    MC2 <- sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal, 
                            timesteps = ts, start_abun = start_abun,
                            extirp_prob = extirp_prob,
                            landscape = landscape, disp_mat = disp_mat, env.df = e.2, 
                            env_traits.df = env_traits.df, int_mat = int_mat
                            )
    
    
    MC2 <- Isbell_2018_cleaner(mix = MC2$mixture, mono = MC2$monoculture, 
                               t_sel  = t_sel )
    
    # add the environmental variables to these data
    MC2 <- left_join(MC2, e.2, by = c("place", "time"))
    
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
BEFF_obs <- lapply(Exp1, function(x) {
  
  y <- Isbell_2018_sampler(data = x[, c("sample", "time", "place", "species", "M", "Y")],
                           RYe_post = FALSE, RYe = rep(1/species, species)
                          )
  return(y$Beff)
  
})

# bind rows and convert to wide format
BEFF_obs <- 
  BEFF_obs %>%
  bind_rows(., .id = "ID") %>%
  mutate(exp_het = rep(c("high", "low"), each = n_rep*9)) %>%
  select(ID, exp_het, Beff, Value) %>%
  pivot_wider(id_cols = c("ID", "exp_het"),
              names_from = Beff,
              values_from = Value)
(BEFF_obs)

# plot the data: Net biodiversity effect
ggplot() +
  geom_jitter(data = BEFF_obs,
             mapping = aes(x = exp_het, y = NBE), width = 0.1) +
  geom_point(data = df_obs %>% group_by(exp_het) %>% summarise(NBE = mean(NBE)),
             mapping = aes(x = exp_het, y = NBE), colour = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

# plot the data: Proportion spatial insurance
ggplot() +
  geom_jitter(data = BEFF_obs %>% mutate(SI_prop = SI/NBE),
              mapping = aes(x = exp_het, y = SI_prop), width = 0.1) +
  geom_point(data = BEFF_obs %>% group_by(exp_het) %>% summarise(SI_prop = mean(SI/NBE)),
             mapping = aes(x = exp_het, y = SI_prop), colour = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()


## Assume incomplete data

# merge these data into a single data.frame and 
Exp1_mobs <- 
  Exp1 %>%
  bind_rows(., .id = "ID") %>%
  select(ID, sample, time, place, env1, species, M, Y)

# for each cluster, we pick 3 samples (i.e. time, place combination)
# where we have mixture and monoculture data, only mixture data for all others
Exp1_msim <- 
  lapply(Exp1, function(x) {
  df <- x
  x <- sample(unique(df$sample), 3)
  df$M[!(df$sample %in% x)] <- NA
  return(df)
})

# merge these data into a single data.frame and 
Exp1_msim <- 
  Exp1_msim %>%
  bind_rows(., .id = "ID") %>%
  select(ID, sample, time, place, env1, species, M, Y)

# model the monoculture yields using rethinking()
library(rethinking)

Exp1_msim_cc <- Exp1_msim[complete.cases(Exp1_msim), ]
dat <- list(
  M = Exp1_msim_cc$M ,
  Y = standardize(Exp1_msim_cc$Y ) ,
  E = Exp1_msim_cc$env1,
  S = Exp1_msim_cc$species)
str(dat)

# intercept only
m1 <- ulam(
  alist(
    M ~ dpois( lambda ),
    log(lambda) <- aS[S] + b_yS[S]*Y + b_eS[S]*E,
    
    # priors
    aS[S] ~ normal(3, 1),
    b_yS[S] ~ normal(0, 1),
    b_eS[S] ~ normal(0, 1)
    
  ), data=dat , chains = 4 )

# check the traceplots and the precis output
traceplot(m1)
precis(m1, depth = 3)

# predict new data
msim_pred <- list(
  Y = standardize(Exp1_msim[is.na(Exp1_msim$M), ]$Y ) ,
  E = Exp1_msim[is.na(Exp1_msim$M), ]$env1,
  S = Exp1_msim[is.na(Exp1_msim$M), ]$species)

# simulate observations from these data
m1_pred <- sim(m1, data = msim_pred)

# summarise this posterior distribution
mu_m1 <- apply(m1_pred, 2, function(x) mean(x) )
PI_m1 <- apply(m1_pred, 2, function(x) PI(x, 0.95) )

# pull into a data.frame and plot
data.frame(mono_obs = Exp1_mobs[is.na(Exp1_msim$M), ]$M,
           species = as.character(Exp1_msim[is.na(Exp1_msim$M), ]$species),
           mu_m1 = mu_m1, 
           PI_low = apply(PI_m1, 2, function(x) x[1]),
           PI_high = apply(PI_m1, 2, function(x) x[2]) ) %>%
  ggplot(data = .,
         mapping = aes(x = mono_obs, y = mu_m1, colour = species)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = PI_low, ymax = PI_high), width = 0.1) +
  ylab("mean +- 89% prediction interval") +
  xlab("actual monoculture yield") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_classic()
  
# simulate different initial proportions from the Dirichlet distribution
dr <- sapply(1:100, function(x) gtools::rdirichlet(n = 1, rep(3, length(unique(Exp1_msim$species))) ) )
#dr <- sapply(1:100, function(x) rep(1/3, 3) )

Mono_reps <- 
  pbapply(m1_pred[sample(x = 1:nrow(m1_pred), 50),], 1, function(x) {
  
  # fill in the missing monoculture data with one sample from the posterior
  df <- Exp1_msim
  df[is.na(df$M), ]$M <- x
  
  # split the df data by ID i.e. experimental replicate
  df <- split(df, df$ID)
  
  Exp_reps <- lapply(df, function(y) {
    
    RYe_reps <- 
      apply(dr, 2, function(z) {
      
      a <- Isbell_2018_sampler(data = y[,c("sample", "time", "place", "species", "M", "Y")], RYe = z, RYe_post = FALSE)
      return(bind_rows(a$Beff, rename(a$L.Beff, Beff = L.Beff)))
      
    } )
    
    return(bind_rows(RYe_reps, .id = "RYe_rep"))
    
  } )
  
  return(bind_rows(Exp_reps, .id = "ID"))
  
} )

# posterior Beff effects
post_Beff <- 
  bind_rows(Mono_reps, .id = "mono_rep") %>%
  arrange(ID, mono_rep, RYe_rep, Beff, Value) %>%
  pivot_wider(id_cols = c("ID", "mono_rep", "RYe_rep"),
              names_from = "Beff", 
              values_from = "Value")

# add the heterogeneity treatment to this using the observed data
post_Beff <- full_join(post_Beff, BEFF_obs[, c("ID", "exp_het")], by = "ID")

# view the data
View(post_Beff)

ggplot(data = post_Beff,
       mapping = aes(x = exp_het, y = SI/NBE, colour = ID)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.5), outlier.shape = NA) +
  geom_point(data = BEFF_obs, 
             mapping = aes(x = exp_het, y = SI/NBE, colour = ID),
             size = 3.5, position = position_dodge(width = 0.5)) +
  ylab("Spatial insurance prop. (SI/NBE)") +
  xlab("Heterogeneity") +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_viridis_d(option = "C") +
  theme_classic() +
  theme(legend.position = "none")

post_Beff %>%
  filter((SI/NBE) > 1) %>%
  View()



# Calculate the biodiversity effects using different samples from the posterior
# For each sample from the posterior, also assume a variety of different starting relative yields

# Overall, we should technically only use a set of Dirichlet distributions

# Post sample 1: Dirichlet 1, 2, 3, 4, 5, 6 etc.

# Then, we can directly compare different posterior samples and calculate contrasts
# because posterior samples will be matching

# What if we can't get decent monocultures? Can we do something else?
# Can we use DI models?





  
  
  
  
  
  
  
  


