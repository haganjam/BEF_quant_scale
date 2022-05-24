
# Project: Quantifying biodiversity effects across scales in natural ecosystems
# Title: Can we reliably quantify biodiversity effects across scales in natural ecosystems?
# Author: James Hagan
# Date: 2022/04/24

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19043)

# attached base packages:
# stats, graphics, grDevice, utils, datasets, methods, base     

# other attached packages:
# mcomsimr_0.1.0, pbapply_1.5-0, here_1.0.1, dplyr_1.0.7, tidyr_1.1.4, ggplot2_3.3.5 

# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
library(pbapply)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))

# set the model inputs

# number of replicate simulations
N_REP <- 10

# set a seed
set.seed(54258748)

species <- 5
patches <- 5
timesteps <- 100
dispersal = 0.05
extirp_prob = 0

max_r = 5
# K_max = 150

int_min = 0.05
int_max = 0.05
intra = 0.05

# landscape parameters
# generate a random landscape
l.1 <- 
  data.frame(x = c(25, 50, 75, 50, 50),
             y = c(50, 50, 50, 75, 25))

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

# generate species environmental optima
# optima = seq(0.1, 0.9, length.out = species)
optima = seq(0.15, 0.85, length.out = species)
print(optima)
env_niche_breadth = lapply(1:N_REP, function(x) round(runif(species, 0.1, 0.25), 2))
print(env_niche_breadth)

# selected time points
t_sel <- c(10, 20, 30, 40, 50)

# simulate different initial proportions from the Dirichlet distribution
dr <- sapply(1:100, function(x) gtools::rdirichlet(n = 1, rep(3,  species) ) )
print(dr)

# start loop
M_test <- 
  lapply(1:N_REP, function(a) {
    
    # get the starting abundances
    start_abun <- round(runif(n = species, 5, 30), 0)
    
    t.1 <- 
      data.frame(species = 1:species,
                 optima = optima,
                 env_niche_breadth = env_niche_breadth[[a]],
                 max_r = max_r,
                 K_max = K_max)
    
    # competition matrices
    si.1 <- matrix(runif(n = species*species, min = int_min, max = int_max), 
                   nrow = species, ncol = species)
    si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
    diag(si.1) <- intra
    
    # simulate environmental variables
    autocorr <- round(runif(1, 100, 500), 0)
    e.1 <- mcomsimr::env_generate(landscape = l.1, 
                                  env1Scale = autocorr, timesteps = timesteps, plot = FALSE)
    
    
    # simulate the metacommunity
    MC1a <- sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal,
                             timesteps = timesteps, start_abun = start_abun,
                             extirp_prob = extirp_prob,
                             comp = "Beverton_Holt",
                             landscape = l.1, disp_mat = d.1, env.df = e.1, 
                             env_traits.df = t.1, int_mat = si.1)
    
    # ggplot(data = MC1a$mixture %>% mutate(species = as.character(species)),
    # mapping = aes(x = time, y = N, colour = species)) +
    # geom_line() +
    # facet_wrap(~patch) +
    # theme_classic()
    
    # clean the data.frame to calculate BEF effects
    MC1b <- Isbell_2018_cleaner(mix = MC1a$mixture, mono = MC1a$monoculture, 
                                t_sel  = t_sel)
    
    # get the environmental data
    MC1e <- 
      MC1a$mixture %>%
      filter(time %in% t_sel) %>%
      select(time, place = patch, env) %>%
      distinct()
    
    # calculate the observed biodiversity effects
    RYe_in <- if(length(start_abun) == 1) rep(1/start_abun, species) else start_abun/sum(start_abun)
    BEF1_obs <- Isbell_2018_sampler(data = MC1b, RYe_post = FALSE, RYe = RYe_in )
    BEF1_obs <- 
      bind_rows(BEF1_obs$Beff, rename(BEF1_obs$L.Beff, Beff = L.Beff)) %>%
      mutate(Value = round(Value, 2))
    
    # assume incomplete monoculture and unknown RYE
    MC1_NA <- 
      full_join(MC1b, MC1e, by = c("time", "place")) %>%
      arrange(sample, time, place, species)
    
    p_comb <- combn(unique(MC1_NA$place), m = round(0.30*length(unique(MC1_NA$place)), 0) )
    max_env <- 
      apply(p_comb, 2, function(x) {
      
      y <- 
        MC1_NA %>%
        filter(place %in% x) %>%
        pull(env) %>%
        range(.)
      
      return(diff(y))
      
    })

    # most likely that we would have say 2 places and all times
    MC1_NA$M1 <- ifelse((MC1_NA$place %in% p_comb[, which(max_env == max(max_env))] ), MC1_NA$M, NA)
    
    # model the monoculture yields using rethinking()
    library(rethinking)
    
    MC1_NAcc <- MC1_NA[complete.cases(MC1_NA), ]
    MC_dat <- list(
      M = MC1_NAcc$M ,
      Y = standardize(MC1_NAcc$Y ) ,
      E = MC1_NAcc$env,
      S = MC1_NAcc$species)
    str(MC_dat)
    
    # model the full model with E and Y interaction
    m1 <- ulam(
      alist(
        M ~ dpois( lambda ),
        log(lambda) <- aS[S] + b_yS[S]*Y + b_eS[S]*E + b_yeS[S]*E*Y,
        
        # priors
        aS[S] ~ dnorm(2, 2),
        b_yS[S] ~ normal(0, 1),
        b_eS[S] ~ normal(0, 1),
        b_yeS[S] ~ normal(0, 1)
        
      ), data=MC_dat , chains = 4, iter = 2000, log_lik = FALSE)
    
    # check the traceplots and the precis output
    # traceplot(m1)
    # precis(m1, depth = 3)
    
    # predict the missing data
    MC1_pred <- MC1_NA[is.na(MC1_NA$M1), ]
    MC1_pred <- MC1_pred[, c("species", "Y", "env")]
    names(MC1_pred) <- c("S", "Y", "E") 
    MC1_pred$Y <- standardize(MC1_pred$Y)
    
    # simulate observations from these data
    MC1_pred <- sim(m1, data = MC1_pred)
    dim(MC1_pred)
    
    # summarise the posterior distribution
    mu_m1 <- apply(MC1_pred, 2, function(x) mean(x) )
    PI_m1 <- apply(MC1_pred, 2, function(x) PI(x, 0.95) )
    
    # pull into a data.frame and plot
    df_plot <- 
      data.frame(mono_obs = MC1_NA$M[is.na(MC1_NA$M1)],
                 species = as.character(MC1_NA$species[is.na(MC1_NA$M1)]),
                 mu_m1 = mu_m1, 
                 PI_low = apply(PI_m1, 2, function(x) x[1]),
                 PI_high = apply(PI_m1, 2, function(x) x[2]) )
    
    ggplot(data = df_plot,
           mapping = aes(x = mono_obs, y = mu_m1, colour = species)) +
      geom_point() +
      geom_errorbar(mapping = aes(ymin = PI_low, ymax = PI_high), width = 0.1) +
      ylab("mean +- 89% prediction interval") +
      xlab("actual monoculture yield") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      theme_classic()
    
    Mono_reps <- 
      pbapply(MC1_pred[sample(x = 1:nrow(MC1_pred), 100),], 1, function(x) {
        
        # fill in the missing monoculture data with one sample from the posterior
        df <- MC1_NA
        df[is.na(df$M1), ]$M1 <- x
        df <- df[,c("sample", "time", "place", "species", "M1", "Y")]
        df <- rename(df, M = M1)
        
        RYe_reps <- 
          apply(dr, 2, function(z) {
            
            a <- Isbell_2018_sampler(data = df, RYe = z, RYe_post = FALSE)
            return(bind_rows(a$Beff, rename(a$L.Beff, Beff = L.Beff)))
            
          } )
        
        return(RYe_reps)
        
      } )
    
    df_unc <- bind_rows(Mono_reps, .id = "ID")
    
    # plot the data
    # ggplot() +
    # geom_density(data = df_unc,
    # mapping = aes(x = Value, colour = Beff, fill = Beff),
    # alpha = 0.3) +
    # geom_vline(xintercept = 0) +
    # geom_vline(data = BEF1_obs,
    # mapping = aes(xintercept = Value), colour = "red") +
    # facet_wrap(~Beff, scales = "free") +
    # theme_classic()
    
    # calculate accuracy metrics
    df_unc_sum <- 
      df_unc %>%
      group_by(Beff) %>%
      summarise(PI_low = PI(Value, prob = 0.90)[1],
                PI_high = PI(Value, prob = 0.90)[2],
                HPDI_low = HPDI(Value, prob = 0.90)[1],
                HPDI_high = HPDI(Value, prob = 0.90)[2],
                mu = mean(Value), .groups="drop") %>%
      mutate(PI_int = PI_high - PI_low,
             HPDI_int = HPDI_high - HPDI_low)
    
    # join the observed values to this data.frame
    df_unc_sum <- full_join(df_unc_sum, rename(BEF1_obs, Value_obs = Value ), by = "Beff")
    
    # calculate deviation from the mean in percent and percent of the width
    df_unc_sum <- 
      df_unc_sum %>%
      mutate(mu_deviation = abs((abs(mu - Value_obs)/Value_obs)*100) ) %>%
      select(Effect = Beff,
             starts_with("PI"), 
             starts_with("HPDI"),
             mu,
             mu_deviation,
             Value_obs)
    
    # add covariates regarding the parameters
    df_unc_sum$start_abun <- paste(start_abun, collapse = "_")
    df_unc_sum$inter_comp <- mean(si.1[lower.tri(si.1)])
    df_unc_sum$optima <- paste(t.1$optima, collapse = "_")
    df_unc_sum$niche_breadth <- paste(t.1$env_niche_breadth, collapse = "_")  
    df_unc_sum$t_steps <- paste(t_sel, collapse = "_")
    df_unc_sum$dispersal <- dispersal
    df_unc_sum$mono_cor <- cor(mu_m1, MC1_NA$M[is.na(MC1_NA$M1)])
    df_unc_sum$bad_mono_n <- sum( apply(PI_m1, 2, diff)/mean(mu_m1) > 2 )/length(mu_m1) 
    
    # reorder the columns
    df_unc_sum <- 
      df_unc_sum %>%
      select(t_steps, dispersal, start_abun, optima, niche_breadth, inter_comp, mono_cor, bad_mono_n,
             Effect, Value_obs, mu, mu_deviation, starts_with("PI"),
             starts_with("HPDI"))
    
    return(list(bind_rows(MC1a), df_unc_sum) )
    
  } )

saveRDS(M_test, here("results/test_sim.rds"))

# load the test data
M_test1 <- readRDS(file = here("results/test_sim.rds"))
head(M_test1)
length(M_test1)

M_test1[[3]][[1]] %>%
  pull(mono_mix) %>%
  unique()

M_test1[[3]][[1]] %>%
  pivot_wider(id_cols = c("sample", "time", "patch", "env", "species"),
              names_from = "mono_mix",
              values_from = "N") %>%
  ggplot(data = .,
         mapping = aes(x = mixture, y = monoculture)) +
  geom_point()

# bind these data
M_test <- bind_rows(lapply(M_test1, function(x) x[[2]]), .id = "ID")
head(M_test)
dim(M_test)
View(M_test)

# check the relationship between the true value and the mean with interval
ggplot(data = M_test %>%
         filter(HPDI_int != max(HPDI_int)),
       mapping = aes(x = Value_obs, y = mu, colour = mono_cor)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = HPDI_low, ymax = HPDI_high)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  scale_colour_viridis_c() +
  facet_wrap(~Effect, scales = "free") +
  theme_classic()

# relationship between interval width and correlation
ggplot(data = M_test %>% filter(mono_cor != min(mono_cor)),
       mapping = aes(x = mono_cor, y = HPDI_int, colour = bad_mono_n)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  scale_colour_viridis_c() +
  facet_wrap(~Effect, scales = "free") +
  theme_classic()

hist(M_test$mono_cor)
hist(M_test$bad_mono_n)

# deviation of the mean from the true value and the correlation 
ggplot(data = M_test,
       mapping = aes(x = mono_cor, y = abs(mu - Value_obs))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  scale_colour_viridis_c() +
  facet_wrap(~Effect, scales = "free") +
  theme_classic()

### END
