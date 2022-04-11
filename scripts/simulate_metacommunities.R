
# Simulate a test of our hypothesis

# r versus K, should they be positively correlated? Probably not...

# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))

# Environmental heterogeneity experiment function
env_het_exp <- function(n_rep = 20, temp_fluc = 0.025,
                        env_min = 0.2, env_max = 0.8, t_sel,
                        patches, species, dispersal = 0.05, 
                        timesteps, start_abun = 60,
                        extirp_prob = 0,
                        landscape, disp_mat, 
                        env_traits.df, int_mat) {
  
  # simulate a set of metacommunities with high environmental heterogeneity among patches
  high_out <- vector("list", n_rep)
  high_env_var <- vector("list", length = n_rep)
  for (i in 1:n_rep) {
    
    # simulate an environment with high heterogeneity
    env1.high <- 
      sapply(seq(env_min, env_max, length.out = patches), function(x) {
        
        # choose a random number between 0 and 0.025 and then decide if it is negative or positive
        z <- 
          sapply(runif(n = ts-1, min = 0, temp_fluc), function(y) {
            y*sample(x = c(-1,1), size = 1)
          })
        
        # cumulatively sum the environmental fluctuations and make sure they remain between zero and one
        u <- cumsum(c(x, z))
        ifelse( (u < 0), 0, ifelse( (u > 1), 1, u))
      })
    
    # bind into a data.frame
    e.1 <- data.frame(env1 = unlist(apply(env1.high, 1, function(x) list(x) ) ),
                      patch = rep(1:patches, ts),
                      time = rep(1:ts, each = patches))
    
    # simulate the metacommunities
    x <- 
      sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal, 
                       timesteps = ts, start_abun = start_abun,
                       extirp_prob = extirp_prob,
                       landscape = landscape, disp_mat = disp_mat, env.df = e.1, 
                       env_traits.df = env_traits.df, int_mat = int_mat)
    
    high_out[[i]] <- x
    
    x <- Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, 
                             t_sel  = t_sel)
    
    x$env_het <- "high"
    
    high_env_var[[i]] <- x
    
  }
  
  # simulate a set of metacommunities with low environmental heterogeneity among patches
  low_out <- vector("list", n_rep)
  low_env_var <- vector("list", length = n_rep)
  low_var <- seq(env_min, env_max, length.out = n_rep)
  for (i in 1:n_rep) {
    
    # simulate low heterogeneity
    e.2 <- 
      sapply(rep(low_var[i], patches), function(x) {
        z <- 
          sapply(runif(n = ts-1, min = 0, temp_fluc), function(y) {
            y*sample(x = c(-1,1), size = 1)
          })
        u <- cumsum(c(x, z))
        ifelse( (u < 0), 0, ifelse( (u > 1), 1, u))
      })
    
    e.2 <- 
      data.frame(env1 = unlist(apply(e.2, 1, function(x) list(x) ) ),
                 patch = rep(1:patches, ts),
                 time = rep(1:ts, each = patches))
    
    x <- 
      sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal, 
                       timesteps = ts, start_abun = start_abun,
                       extirp_prob = extirp_prob,
                       landscape = landscape, disp_mat = disp_mat, env.df = e.2, 
                       env_traits.df = env_traits.df, int_mat = int_mat)
    low_out[[i]] <- x
    
    x <- Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, 
                             t_sel  = t_sel )
    
    x$env_het <- "low"
    
    low_env_var[[i]] <- x
    
  }
  
  # bind these data.frames together
  bind_rows( bind_rows(high_env_var, .id = "rep"), bind_rows(low_env_var, .id = "rep")) %>%
    arrange(rep, env_het)
  
}

# function to simulate different monocultures, calculate biodiversity effects and plot
calculate_bio_effects <- function(obs_data, mono_n = 20, mono_sd = 10, 
                                  alpha_par = 2, N = 100) {
  
  # simulate several of these datasets with uncertainty in the monoculture data
  mono_n <- 20
  df_unc <- 
    lapply(1:mono_n, function(x) {
      
      df <- obs_data
      y <- rnorm(n = length(df$M), mean = df$M, sd = mono_sd )
      df$M <- ifelse(y < 0, 0, y)
      
      return(df)
      
    } )
  
  # apply the isbell sampler to each rep x env_het combination
  part_exp <- 
    lapply(df_unc, function(z) {
      
      # first, we calculate the BEF effects for the high heterogeneity treatments
      x <- filter(z, env_het == "high")
      x <- split(x, x$rep)
      
      part_high <- 
        lapply(x, function(y) {
          
          Isbell_2018_sampler(data = y, 
                              RYe_post = TRUE, alpha_par = alpha_par, N = N)
          
        })
      part_high <- bind_rows(part_high, .id = "exp_rep")
      part_high$env_het <- "high"
      
      # second, we do the same for the low heterogeneity treatments
      x <- filter(z, env_het == "low")
      x <- split(x, x$rep)
      
      part_low <- 
        lapply(x, function(y) {
          
          Isbell_2018_sampler(data = y, 
                              RYe_post = TRUE, alpha_par = alpha_par, N = N)
          
        })
      part_low <- bind_rows(part_low, .id = "exp_rep")
      part_low$env_het <- "low"
      
      # bind these data together
      bind_rows(part_high, part_low)
      
    } )
  
  # view the dataset
  part_exp <- 
    bind_rows(part_exp, .id = "mono_rep") %>%
    arrange(env_het, exp_rep, mono_rep, sample) %>%
    select(env_het, exp_rep, mono_rep, sample, Beff, Value)
  
  # convert these data to the wide format
  part_exp_wide <- 
    part_exp %>%
    pivot_wider(id_cols = c("env_het", "exp_rep", "mono_rep", "sample"),
                names_from = "Beff",
                values_from = "Value") 
  
  # plot these data
  p1 <- 
    part_exp_wide %>%
    rename(`cluster` = exp_rep) %>%
    filter(NBE > quantile(NBE, 0.05), NBE < quantile(NBE, 0.95)) %>%
    ggplot(data = .,
           mapping = aes(x = env_het, y = NBE, colour = `cluster`)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.05, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.5), outlier.shape = NA) +
    scale_colour_viridis_d(option = "C") +
    scale_fill_viridis_d(option = "C", alpha = 0.5) +
    theme_classic() + 
    ylab("Net biodiversity effect") +
    xlab("Environmental heterogeneity") +
    theme_classic() +
    theme(legend.position = "bottom")
  
  # proportion of the spatial insurance effect
  SI_prop <- 
    part_exp_wide %>%
    mutate(SI_prop = (SI/NBE) ) %>%
    rename(`cluster` = exp_rep) %>%
    filter(SI_prop > quantile(SI_prop, 0.05), SI_prop < quantile(SI_prop, 0.95))
  
  p2 <- 
    ggplot(data = SI_prop,
           mapping = aes(x = env_het, y = SI_prop, colour = `cluster`, fill = `cluster`) ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.05, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.5), outlier.shape = NA) +
    scale_colour_viridis_d(option = "C") +
    scale_fill_viridis_d(option = "C", alpha = 0.5) +
    # scale_y_continuous(limits = c(-0.05, quantile(SI_prop$SI_prop, 0.95) + 0.05 )) +
    theme_classic() + 
    ylab("Spatial insurance/Net biodiversity effect") +
    xlab("Environmental heterogeneity") +
    theme(legend.position = "bottom")
  
  list(df_unc, part_exp_wide, p1, p2)
  
}

## Set-up fixed inputs
species <- 3
patches <- 3
ts <- 100

## Landscape parameters

# generate a random landscape
l.1 <- 
  data.frame(x = c(25, 50, 75),
             y = c(50, 50, 50))
plot(l.1)

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)


## Exp 1: Narrow niches, neutral competition

# generate species environmental optima
t.1 <- 
  data.frame(species = 1:species,
             optima = seq(0.2, 0.8, length.out = species),
             env_niche_breadth = 0.2,
             max_r = 0.5,
             K_max = 150)
head(t.1)


## Competition matrices
si.1 <- matrix(runif(n = species*species, min = 0.1, max = 0.75), 
               nrow = species, ncol = species)
si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
diag(si.1) <- 1
head(si.1)

# simulate the observed dataset for the heterogeneity experiment
df_obs <- env_het_exp(n_rep = 5, temp_fluc = 0.025,
                      env_min = 0.3, env_max = 0.7, t_sel = seq(30, 90, length.out = 3),
                      patches = patches, species = species, dispersal = 0.05, 
                      timesteps = ts, start_abun = 60,
                      extirp_prob = 0,
                      landscape = l.1, disp_mat = d.1,
                      env_traits.df = t.1, int_mat = si.1 )
head(df_obs)

df_obs %>% View()

# add uncertainty to this
exp1 <- calculate_bio_effects(obs_data = df_obs, 
                              mono_n = 20, mono_sd = 10, alpha_par = 2, N = 50
                              )

exp1[[3]]
exp1[[4]]

exp1[[1]]



