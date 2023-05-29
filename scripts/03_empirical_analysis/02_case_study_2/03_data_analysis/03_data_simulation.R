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

# simulate a dataset to fit the model

# simulate all combinations of S, C and T with three replicates
df_preds <- expand.grid(REP = 1:8, 
                        S = 1:5,
                        C = 1:10,
                        T = (1:3)/3)

# check the data
head(df_preds)
dim(df_preds)

# remove the replicate column
df_preds <- df_preds[, -1]

# simulate Y, PC1 and PC2 variables from the uniform distribution
df_preds$PC1 <- runif(n = nrow(df_preds), 0, 1)
df_preds$PC2 <- runif(n = nrow(df_preds), 0, 1)
df_preds$Y <- runif(n = nrow(df_preds), 0, 1)

# set the number of species and clusters
Ns <- 5
Nc <- 10

# log-normal model

# get the log-normal sigma value
sigma <- rexp(n = 1, rate = 3)

# sample the alpha bar parameters as means of the multivariate normal
abar <- round(rnorm(n = Ns, mean = 0, sd = 1), 2)

# sample the sigma parameter of the multivariate normal
sigma_a <- round(rexp(n = Ns, rate = 3), 2)

# sample the correlation matrix
x <- round(rlkjcorr(n = 1, K = Ns, eta = 3), 2)
Rho_a <- matrix(data = x, nrow = Ns, ncol = Ns)

# get a matrix of z-scores
z <- matrix(rnorm(n = Ns*Nc, mean = 0, sd = 1), nrow = Ns, ncol = Nc)
  
# obtain the Cholesky factor from this matrix
L <- chol(Rho_a)
  
# get the offsets from the z-scores
a <- t(diag(sigma_a) %*% L %*% z )
  
# sample the beta parameters
b1bar <- 0
sigma_b1 <- 1
b1 <- round(rnorm(n = Ns, mean = b1bar, sd = sigma_b1), 2)

# logistic regression model

# sample the alpha parameters
abar_hu <- 0.5
sigma_a_hu <- 1
a_hu <- round(rnorm(n = Ns, mean = abar_hu, sd = sigma_a_hu), 2)

# sample the beta parameters
b1bar_hu <- 0
sigma_b1_hu <- 1
b1_hu <- round(rnorm(n = Ns, mean = b1bar_hu, sd = sigma_b1_hu), 2)

# simulate observations from these parameters
preds <- vector(length = nrow(df_preds))
for(i in 1:nrow(df_sim)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_preds, 
         (abar[S[i]] + a[C[i], S[i]]) + (b1[S[i]] * Y[i]))
  
  # draw from the log-normal distribution
  x <- rlnorm(n = length(x), x, sigma)
  
  # get the probability of 0
  y <- 
    with(df_preds, 
         (a_hu[S[i]] + b1_hu[S[i]]*Y[i]))
  y <- plogis(y)
  y <- 1-y
  
  # draw from the binomial distribution
  y <- rbinom(n = length(y), size = 1, prob = y)
  
  preds[i] <- (x*y)
  
}

# add the simulated values to the df_sim data
df_preds[["M"]] <- preds

# set-up the data list
dat <- 
  list(N = length(df_preds$M),
       S_N = length(unique(df_preds$S)),
       C_N = length(unique(df_preds$C)),
       M = df_preds$M,
       Y = df_preds$Y,
       C = df_preds$C,
       S = df_preds$S)

# compile the model
m_sim <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model4.stan",
                        verbose = TRUE)

# sample the stan model
m_sim_fit <- rstan::sampling(m_sim, data = dat, 
                             iter = 1500, chains = 4, algorithm = c("NUTS"),
                             control = list(adapt_delta = 0.99, max_treedepth = 12),
                             seed = 54856, cores = 4)

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

# check simulated parameter
n <- 13

pars[n]
eval(parse(text = pars[n]))
pars_rows[pars[n] == pars_rows]
diag[pars[n] == pars_rows, ]$mean

# check the relationship between these parameters
m_true <- eval(parse(text = pars[n]))
m_est <- t(matrix(diag[pars[n] == pars_rows, ]$mean, ncol = 10, nrow = 5))

col <- 3
plot(m_true[,col], m_est[,col])
abline(0, 1)




