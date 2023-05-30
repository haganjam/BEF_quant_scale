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
library(brms)
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

# number of samples
N <- 1000

# simulate two predictor variables
df_t <- data.frame(X1 = rnorm(n = N, mean =  0, sd = 1),
                   X2 = rnorm(n = N, mean = 0, sd = 1))

# simulate a log-normally distributed response variable
Y_ln <- rlnorm(n = N, 0.5 + 0.75*df_t$X1, sd = 1)

# simulate a binomial variable
Y_lg <- rbinom(n = N, size = 1, 1 - plogis(0.5 + 0.2*df_t$X2))

# add a Y-variable to these data
df_t$Y <- Y_ln*Y_lg

# fit the log-normal hurdle model using brms()
m_t <- brm(
  bf(Y ~ X1,
     hu ~ X2),
  data = df_t,
  family = hurdle_lognormal(),
  chains = 4, iter = 1000, warmup = 500, cores = 4
)

# check the stancode for this model
stancode(m_t)


# simulate a dataset to fit the model

# get combinations of S, C and T with three replicates
df_preds <- expand.grid(REP = 1:12, 
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
for(j in 1:N_rep) {
  
  # simulate Y, PC1 and PC2 variables from the uniform distribution
  df_preds$PC1 <- runif(n = nrow(df_preds), 0, 1)
  df_preds$PC2 <- runif(n = nrow(df_preds), 0, 1)
  df_preds$Y <- runif(n = nrow(df_preds), 0, 1) * rbinom(n = nrow(df_preds), 1, 0.80)
  
  # log-normal model
  
  # get the log-normal sigma value
  sigma <- round(rexp(n = 1, rate = 5), 3)
  
  # sample the alpha bar parameters as means of the multivariate normal
  abar <- round(runif(n = Ns, min = -2.5, max = 2.5), 3)
  
  # sample the sigma parameter of the multivariate normal
  sigma_a <- round(rexp(n = Ns, rate = 4), 3)
  
  # sample the correlation matrix
  x <- round(rlkjcorr(n = 1, K = Ns, eta = 2), 2)
  Rho_a <- matrix(data = x, nrow = Ns, ncol = Ns)
  
  # get a matrix of z-scores
  z <- matrix(rnorm(n = Ns*Nc, mean = 0, sd = 1), nrow = Ns, ncol = Nc)
  
  # obtain the Cholesky factor from this matrix
  L <- chol(Rho_a)
  
  # get the offsets from the z-scores
  a <- t(diag(sigma_a) %*% L %*% z )
  
  # sample the beta parameters
  b1bar <- 0
  sigma_b1 <- 0.5
  b1 <- round(rnorm(n = Ns, mean = b1bar, sd = sigma_b1), 2)
  
  # logistic regression model
  
  # sample the alpha parameters
  abar_hu <- 0
  sigma_a_hu <- 1
  a_hu <- round(rnorm(n = Ns, mean = abar_hu, sd = sigma_a_hu), 2)
  
  # sample the beta parameters
  b1bar_hu <- 0
  sigma_b1_hu <- 1
  b1_hu <- round(rnorm(n = Ns, mean = b1bar_hu, sd = sigma_b1_hu), 2)
  
  # simulate observations from these parameters
  preds <- vector(length = nrow(df_preds))
  for(i in 1:nrow(df_preds)) {
    
    # get the mean prediction from the lognormal model on the natural scale
    x <- 
      with(df_preds, 
           (abar[S[i]] + a[C[i], S[i]]) + (b1[S[i]] * Y[i]))
    
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
df_sim <- pps_list[[sample(1:length(pps_list)), 1]]
dat <- 
  list(N = length(df_sim$M),
       S_N = length(unique(df_sim$S)),
       C_N = length(unique(df_sim$C)),
       M = df_sim$M,
       Y = df_sim$Y,
       C = df_sim$C,
       S = df_sim$S)

# compile the model
m_sim <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model4.stan",
                        verbose = TRUE)

# sample the stan model
m_sim_fit <- rstan::sampling(m_sim, data = dat, 
                             iter = 1000, chains = 4, algorithm = c("NUTS"),
                             control = list(adapt_delta = 0.99, max_treedepth = 12),
                             cores = 4)

# check the model output
print(m_sim_fit)

# check the model parameters
pars <- m_sim_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars) | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
rstan::traceplot(m_sim_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- as.data.frame(rstan::summary(m_sim_fit)$summary)

# parameter row names
pars_rows <- gsub(pattern = "\\[(.*?)\\]", replace = "", row.names(diag))

# check simulated parameters
n <- 13

pars[n]
eval(parse(text = pars[n]))
pars_rows[pars[n] == pars_rows]
diag[pars[n] == pars_rows, ]$mean

# check the relationship between matrix parameters
m_true <- eval(parse(text = pars[n]))
m_est <- t(matrix(diag[pars[n] == pars_rows, ]$mean, ncol = 10, nrow = 5))

col <- 5
plot(m_true[,col], m_est[,col])
abline(0, 1)

### END
