
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

# next steps:

# 1. split the code so that we can push the most computationally intensive part of the 
# code i.e. using the posterior distribution and calculating biodiversity effects for each
# sample 100 times with different samples from the Dirichlet distribution

# - currently, this takes five minutes

# if we run 100 models, without that step, it should only take 2 hours and then we don't
# have the hassle of trying to load these obscure packages to the cluster

# once that it is done, we can save the outputs and then, in a new script run 
# the calculates on the Albiorix cluster

# this will save a lot of time


# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
library(pbapply)
source(here("scripts/02_simulation/01_mcomsimr_simulate_MC2_function.R"))

# set the model inputs

# number of replicate simulations
N_REP <- 10

# set a seed
set.seed(54258748)

species <- 5
patches <- 20
timesteps <- 100
dispersal = 0.05
extirp_prob = 0.00001

max_r = 5
K_max = NA

int_min = 0.05*0
int_max = 0.05*1.5
intra = 0.05*1

comp <- "Beverton_Holt"

# landscape parameters
# generate a random landscape
n1 <- 5
n2 <- 4
assertthat::assert_that(assertthat::see_if( (n1*n2) == patches ),
                        msg = "change the n1 and n2 values so that their product matches the patches")

t1 <- round(seq(25, 100, length.out = 5), 0)
t2 <- round(seq(25, 100, length.out = 4), 0)
t12 <- expand.grid(t1, t2)

l.1 <- data.frame(x = t12[[1]], y = t12[[2]])
plot(l.1)

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

# run the pipeline which does the following:

# 1. simulates a BEF experiment with full mixtures in monocultures for each patch

# 2. uses the known starting abundances and full monoculture data to calculate observed biodiversity effects (Isbell et al. 2018)

# 3. assumes that we only have monoculture data in 30 % of the patches

# 4. uses a simple linear model to model the missing monoculture data using the environment and mixture abundance of each species

# 5. calculates Isbell et al.'s (2018) biodiversity effects using all the uncertainty present in the
# - modelled monocultures by drawing different samples from the posterior distribution
# - using starting relative abundances from a Dirichlet distribution

# 6. summarises the distribution for each biodiversity effect and compares it to the observed value
# - we assume that for an effect to be meaningful, the 95% HPDI interval of an effect must be
# within the interval of obs value +- 0.5*observed value. This means that the 95% HPDI interval
# contains the true value within a reasonable distance of the true value


M_test <- 
  
  pblapply(1:N_REP, function(a) {
  
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
  print(autocorr)
  e.1 <- mcomsimr::env_generate(landscape = l.1, 
                                env1Scale = autocorr, timesteps = timesteps, plot = FALSE)
  
  
  # simulate the metacommunity
  MC1a <- sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal,
                           timesteps = timesteps, start_abun = start_abun,
                           extirp_prob = extirp_prob,
                           comp = comp,
                           landscape = l.1, disp_mat = d.1, env.df = e.1, 
                           env_traits.df = t.1, int_mat = si.1)
  
  ggplot(data = MC1a$mixture %>% mutate(species = as.character(species)),
         mapping = aes(x = time, y = N, colour = species)) +
    geom_line() +
    facet_wrap(~patch) +
    theme_classic()
  
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
  
  # get the known initial relative abundance values
  RYe_in <- if(length(start_abun) == 1) rep(1/start_abun, species) else start_abun/sum(start_abun)
  
  # use the Isbell_2018_sampler function to calculate the actual biodiversity effects
  BEF1_obs <- Isbell_2018_sampler(data = MC1b, RYe_post = FALSE, RYe = RYe_in )
  
  BEF1_obs <- 
    bind_rows(BEF1_obs$Beff, rename(BEF1_obs$L.Beff, Beff = L.Beff)) %>%
    mutate(Value = round(Value, 2))
  
  # assume incomplete monoculture and unknown RYE
  MC1_NA <- 
    full_join(MC1b, MC1e, by = c("time", "place")) %>%
    arrange(sample, time, place, species)
  
  # get all combinations of 30% of the places (e.g. 5 places, then we get two places)
  n_NA <- round(0.30*length(unique(MC1_NA$place)), 0)
  p_comb <- combn(unique(MC1_NA$place), m = n_NA)
  
  # calculate the environmental range among each of the combinations of two places
  max_env <- 
    apply(p_comb, 2, function(x) {
      
      y <- 
        MC1_NA %>%
        filter(place %in% x) %>%
        pull(env) %>%
        range(.)
      
      return(diff(y))
      
    })
  
  # most likely that we would have say 2 places and all times: Pick the 2 places with most environmental variation
  MC1_NA$M1 <- ifelse( MC1_NA$place %in% p_comb[, sample(which(max_env == max(max_env)), 1)], MC1_NA$M, NA)
  
  # how many NAs are there compared to all the rows
  sum(is.na(MC1_NA$M1))/nrow(MC1_NA)
  
  
  # model the monoculture yields using rethinking(): Write a proper Stan model!
  library(rethinking)
  
  # get the complete cases to fit the model (i.e. training data)
  MC1_NAcc <- MC1_NA[complete.cases(MC1_NA), ]
  MC_dat <- list(
    M = MC1_NAcc$M ,
    Y = standardize(MC1_NAcc$Y ) ,
    E = MC1_NAcc$env,
    S = MC1_NAcc$species)
  str(MC_dat)
  
  # model the monocultures with with the E and Y interaction
  m1 <- ulam(
    alist(
      M ~ dpois( lambda ),
      log(lambda) <- aS[S] + b_yS[S]*Y + b_eS[S]*E + b_yeS[S]*E*Y,
      
      # priors
      aS[S] ~ dnorm(2, 2),
      b_yS[S] ~ normal(0, 1),
      b_eS[S] ~ normal(0, 1),
      b_yeS[S] ~ normal(0, 1)
      
    ), data = MC_dat , chains = 4, iter = 2000, log_lik = FALSE)
  
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
  
  # fill in the monocultures with one sample from the posterior
  # for that one sample from the posterior, calculate biodiversity effects with 100 samples from the Dirichlet distribution
  Mono_reps <- 
    pbapply(MC1_pred[sample(x = 1:nrow(MC1_pred), 100), ], 1, function(x) {
      
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
  ggplot() + 
    geom_density(data = df_unc,
                 mapping = aes(x = Value, colour = Beff, fill = Beff),
                 alpha = 0.3) +
    geom_vline(xintercept = 0) +
    geom_vline(data = BEF1_obs,
               mapping = aes(xintercept = Value), colour = "red") +
    facet_wrap(~Beff, scales = "free") +
    theme_classic() +
    theme(legend.position = "none")
  
  
  # calculate accuracy metrics
  
  # what percentage of the posterior falls between the 50% mean +-
  
  # set the percentage threshold
  thresh <- 0.5
  
  df_unc_sum <- 
    
    full_join(
      
      df_unc %>% 
        group_by(Beff) %>%
        summarise(mu = mean(Value, na.rm = TRUE),
                  HPDI_low = HPDI(Value, prob = 0.95)[1],
                  HPDI_high = HPDI(Value, prob = 0.95)[2],),
      
      rename(BEF1_obs, Value_obs = Value ), by = "Beff"
      
    ) %>%
    mutate(Value_obs_min = (1-thresh)*Value_obs,
           Value_obs_max = (1+thresh)*Value_obs ) %>%
    mutate(HPDI_true = if_else( (HPDI_low > Value_obs_min) & (HPDI_high < Value_obs_max), TRUE, FALSE ) ) %>%
    mutate(mu_deviation = abs((abs(mu - Value_obs)/Value_obs)*100) ) %>%
    mutate(mu_threshold = thresh)
  
  # check if the different effects are consistent with the NBE
  t1 <- sum(df_unc_sum[df_unc_sum$Beff %in% c("AS", "TC", "NO", "TI", "SI", "ST"), ][["mu"]])
  t2 <- df_unc_sum[df_unc_sum$Beff %in% c("NBE"), ][["mu"]]
  
  assertthat::assert_that(assertthat::see_if(near(t1, t2)),
                          msg = "biodiversity effects do not sum to NBD")

  # add covariates regarding the parameters
  df_unc_sum$start_abun <- paste(start_abun, collapse = "_")
  df_unc_sum$inter_comp <- mean(si.1[lower.tri(si.1)])
  df_unc_sum$optima <- paste(t.1$optima, collapse = "_")
  df_unc_sum$niche_breadth <- paste(t.1$env_niche_breadth, collapse = "_")  
  df_unc_sum$t_steps <- paste(t_sel, collapse = "_")
  df_unc_sum$dispersal <- dispersal
  df_unc_sum$mono_cor <- cor(mu_m1, MC1_NA$M[is.na(MC1_NA$M1)])
  
  # definition of a bad monoculture
  df_unc_sum$bad_mono_n <- sum( apply(PI_m1, 2, diff)/mean(mu_m1) > 1 )/length(mu_m1) 
  
  # reorder the columns
  
  # important: When HPDI_true is TRUE, it means that the 95% HPDI interval is within
  # thresh*obs value +- obs value which means that 95% HPDI interval is narrow enough
  # to be meaningful
  df_unc_sum <- 
    df_unc_sum %>%
    select(t_steps, dispersal, start_abun, optima, niche_breadth, inter_comp, mono_cor, bad_mono_n,
           Beff, Value_obs, mu, mu_deviation,
           starts_with("HPDI"), mu_threshold,
           HPDI_true)
  
  return(list(bind_rows(MC1a), df_unc_sum) )
  
} )

saveRDS(M_test, here("results/test_analytical_pipeline.rds"))

### END
