#'
#' @title: Model missing monocultures
#' 
#' @description: Data simulation to test the model formulation along with
#' prior predictive simulation to choose appropriate priors.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# import helper functions
source("scripts/03_empirical_analysis/helper_functions.R")
source("scripts/Function_plotting_theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(readr)

# check the data

# load the analysis data
data <- read_csv("data/case_study_2/data_clean/biomass_env_analysis_data.csv")

# standardise the predictor variables in the overall data
data$Y <- with(data, (Y - min(Y))/(max(Y)-min(Y)) )
data$PC1 <- with(data, (PC1 - min(PC1))/(max(PC1)-min(PC1)) )
data$PC2 <- with(data, (PC2 - min(PC2))/(max(PC2)-min(PC2)) )
data$cluster_id <- as.integer(as.factor(data$cluster_id))
data$time <- as.integer(as.factor(data$time))/3
data$OTU <- as.integer(as.factor(data$OTU))

ggplot(data = data %>% filter(!is.na(M)),
       mapping = aes(x = as.character(OTU), y = M, colour = as.character(OTU) )) +
  geom_point() +
  facet_wrap(~cluster_id) +
  theme(legend.position = "bottom")


# simulate a simple log-normal hurdle model

# example from: https://dominicmagirr.github.io/post/longitudinal-hurdle-models/
# https://dominicmagirr.github.io/post/longitudinal-hurdle-models-2/


# simulate a dataset to fit the model

# get combinations of S, C and T with three replicates
df_preds <- expand.grid(REP = 1:5, 
                        S = 1:5,
                        C = 1:10,
                        T = (1:3)/3)

# check the data
head(df_preds)
dim(df_preds)

# remove the replicate column
df_preds <- df_preds[, -1]

# set the number of species and clusters
Ns <- 5
Nc <- 10

# perform prior predictive simulation
N_rep <- 100
pps_list <- vector("list", length = N_rep)
par_list <- vector("list", length = N_rep)
for(j in 1:N_rep) {
  
  # simulate Y, PC1 and PC2 variables from the uniform distribution
  df_preds$PC1 <- runif(n = nrow(df_preds), 0, 1)
  df_preds$PC2 <- runif(n = nrow(df_preds), 0, 1)
  df_preds$Y <- runif(n = nrow(df_preds), 0, 1) * rbinom(n = nrow(df_preds), 1, 0.80)
  
  # log-normal model
  
  # get the log-normal sigma value
  sigma <- round(rexp(n = 1, rate = 5), 3)
  
  # sample the alpha bar parameter
  abar <- round(runif(n = 1, min = -2.5, max = 2.5), 1)
  
  # sample the b1 bar parameter
  b1bar <- round(rnorm(n = 1, mean = 0, sd = 1), 1)
  
  # sample the b2 bar parameter
  b2bar <- round(rnorm(n = 1, mean = 0, sd = 1), 1)
  
  # sample the sigma parameter of the multivariate normal
  sigma_v <- round(rexp(n = 3, rate = 4), 3)
  
  # sample the correlation matrix
  x <- round(rlkjcorr(n = 1, K = 3, eta = 2), 2)
  Rho_v <- matrix(data = x, nrow = 3, ncol = 3)
  
  # get a matrix of z-scores
  V <- matrix(rnorm(n = 3*Ns, mean = 0, sd = 1), nrow = 3, ncol = Ns)
  
  # obtain the Cholesky factor from this matrix
  L_v <- chol(Rho_v)
  
  # get the offsets from the z-scores
  v <- t(diag(sigma_v) %*% L_v %*% V )
  
  # convert offsets to standard parameter values
  a <- abar + v[,1]
  b1 <- b1bar + v[,2]
  b2 <- b2bar + v[,3]
  
  # logistic regression model
  
  # sample the alpha bar parameter
  abar_hu <- round(runif(n = 1, min = -2.5, max = 2.5), 1)
  
  # sample the b1 bar parameter
  b1bar_hu <- round(rnorm(n = 1, mean = 0, sd = 1), 1)
  
  # sample the sigma parameter of the multivariate normal
  sigma_z <- round(rexp(n = 2, rate = 4), 3)
  
  # sample the correlation matrix
  x <- round(rlkjcorr(n = 1, K = 2, eta = 2), 2)
  Rho_z <- matrix(data = x, nrow = 2, ncol = 2)
  
  # get a matrix of z-scores
  Z <- matrix(rnorm(n = 2*Ns, mean = 0, sd = 1), nrow = 2, ncol = Ns)
  
  # obtain the Cholesky factor from this matrix
  L_z <- chol(Rho_z)
  
  # get the offsets from the z-scores
  z <- t(diag(sigma_z) %*% L_z %*% Z )
  
  # convert offsets to standard parameter values
  a_hu <- abar_hu + z[,1]
  b1_hu <- b1bar_hu + z[,2]
  
  # simulate observations from these parameters
  preds <- vector(length = nrow(df_preds))
  for(i in 1:nrow(df_preds)) {
    
    # get the mean prediction from the lognormal model on the natural scale
    x <- 
      with(df_preds, 
           (a[S[i]]) + (b1[S[i]] * Y[i])+ (b2[S[i]] * PC1[i]) )
    
    # draw from the log-normal distribution
    x <- rlnorm(n = length(x), x, sigma)
    
    # get the probability of 0
    y <- 
      with(df_preds, 
           a_hu[S[i]] + b1_hu[S[i]]*Y[i])
    y <- 1 - plogis(y)
    
    # draw from the binomial distribution
    y <- rbinom(n = length(y), size = 1, prob = y)
    
    preds[i] <- (x*y)
    
  }
  
  # add the simulated values to the df_sim data
  df_preds[["M"]] <- preds
  
  # add the prior predictive simulation dataset to the pps_list
  pps_list[[j]] <- df_preds
  
  # add the data on the parameter values
  par_list[[j]] <- list(sigma = sigma,
                        vbar = c(abar, b1bar, b2bar),
                        sigma_v = sigma_v,
                        Rho_v = Rho_v,
                        a = a,
                        b1 = b1,
                        b2 = b2,
                        zbar = c(abar_hu, b1bar_hu),
                        sigma_z = sigma_z,
                        Rho_z = Rho_z,
                        a_hu = a_hu,
                        b1_hu = b1_hu)
  
}

# bind the pps_list into a data.frame
pps_df <- bind_rows(pps_list, .id = "rep")

# calculate the mean and sd for each rep for each species
pps_df <- 
  pps_df %>%
  group_by(rep, S, C) %>%
  summarise(mean_M = mean(M),
            min_M = min(M),
            max_M = max(M))

# plot the prior predictive simulations
p1 <- 
  ggplot() +
  geom_line(data = pps_df,
            mapping = aes(x = S, y = mean_M, 
                          colour = as.character(S), group = rep),
            position = position_dodge(0.2), alpha = 0.2) +
  geom_errorbar(data = pps_df,
                mapping = aes(x = S, ymin = min_M,  ymax = max_M, 
                              colour = as.character(S), group = rep),
                width = 0, position = position_dodge(0.2),
                alpha = 0.2) +
  facet_wrap(~C) +
  theme(legend.position = "none")
plot(p1)

# what is the maximum predicted value?
max(pps_df$max_M)


# test if we can recover simulated model parameters

# set-up the data list
id <- sample(1:length(pps_list), 1)
df_sim <- pps_list[[id ]]
dat <- 
  list(N = length(df_sim$M),
       S_N = length(unique(df_sim$S)),
       C_N = length(unique(df_sim$C)),
       M = df_sim$M,
       Y = df_sim$Y,
       PC1 = df_sim$PC1,
       S = df_sim$S)

# compile the model
m_sim <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model0.stan",
                        verbose = TRUE)

# sample the stan model
m_sim_fit <- rstan::sampling(m_sim, data = dat, 
                             iter = 1000, chains = 4, algorithm = c("NUTS"),
                             control = list(adapt_delta = 0.99),
                             cores = 4,
                             seed = 485749)

# check the model output
print(m_sim_fit)

# extract the diagnostic parameters
diag1 <- as.data.frame(rstan::summary(m_sim_fit)$summary)

# compare to custom likelihood model

# compile the model
m_cust <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/lognormal_custom_likelihood_test.stan",
                           verbose = TRUE)

# sample the stan model
m_cust_fit <- rstan::sampling(m_cust, data = dat, 
                             iter = 1000, chains = 4, algorithm = c("NUTS"),
                             control = list(adapt_delta = 0.99),
                             cores = 4,
                             seed = 485749)

# check the output
print(m_cust_fit)

# extract the diagnostic parameters
diag2 <- as.data.frame(rstan::summary(m_cust_fit)$summary)

# compare the two
plot(diag1[-which(diag2$mean == min(diag2$mean)), ]$mean, 
     diag2[-which(diag2$mean == min(diag2$mean)), ]$mean)
abline(0, 1)

# compare a few random rows completely
row_id <- sample(1:nrow(diag1), 5)

diag1[row_id,1:5] %>% round(2)
diag2[row_id,1:5] %>% round(2)

# compare the loo scores
rstan::loo(m_sim_fit)
rstan::loo(m_cust_fit)

# check the model parameters
pars <- m_sim_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars) | grepl("V", pars) | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(m_sim_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- as.data.frame(rstan::summary(m_sim_fit)$summary)

# parameter row names
pars_rows <- gsub(pattern = "\\[(.*?)\\]", replace = "", row.names(diag))

# compare estimated and simulated parameter values

# simulated parameter values
par_list[[id]]

# estimated parameter values
diag[pars_rows %in% names(x),]

# check parameter by parameter
par_vec <- names(par_list[[id]])

# which number?
n <- 6

# simulated value
par_list[[id]][par_vec == par_vec[n]]

# estimated value
diag[par_vec[n] == pars_rows, ]$mean

### END
