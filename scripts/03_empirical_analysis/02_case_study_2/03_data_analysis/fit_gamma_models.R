
# gamma model reference:

# general gamma:
# https://datascienceplus.com/bayesian-regression-with-stan-beyond-normality/
# specific example:
# https://github.com/rmcelreath/cchunts/blob/master/R/models.R
# https://link.springer.com/article/10.1007/s12110-014-9193-4
# https://discourse.mc-stan.org/t/zero-inflated-gamma-model/6788

# stan example models
# https://github.com/stan-dev/example-models

# load the required libraries
library(readr)
library(rstan)
library(loo)

# load the analysis data
df_obs <- read_csv("data/case_study_2/data_clean/biomass_env_analysis_data.csv")

# standardise the predictor variables in the overall data
df_obs$Y <- with(df_obs, (Y - min(Y))/(max(Y)-min(Y)) )
df_obs$PC1 <- with(df_obs, (PC1 - min(PC1))/(max(PC1)-min(PC1)) )
df_obs$PC2 <- with(df_obs, (PC2 - min(PC2))/(max(PC2)-min(PC2)) )
df_obs$cluster_id <- as.integer(as.factor(df_obs$cluster_id))
df_obs$time <- as.integer(as.factor(df_obs$time))/3
df_obs$OTU <- as.integer(as.factor(df_obs$OTU))

# remove any NAs
df_m_obs <- df_obs[complete.cases(df_obs),]

# make a data.list with the training data
df_m_obs <- 
  list(N = length(df_m_obs$M),
       S_N = length(unique(df_m_obs$OTU)),
       C_N = length(unique(df_m_obs$cluster_id)),
       M = df_m_obs$M,
       Y = df_m_obs$Y,
       PC1 = df_m_obs$PC1,
       PC2 = df_m_obs$PC2,
       C = df_m_obs$cluster_id,
       T = df_m_obs$time,
       S = df_m_obs$OTU)

# model 0

# compile the model
m0 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_gamma_model0.stan",
                        verbose = TRUE)

# sample the stan model
m0_fit <- rstan::sampling(m0, data = df_m_obs, 
                          iter = 3000, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99,
                                         stepsize = 0.05,
                                         max_treedepth = 12),
                          cores = 4)

# check the model
m0_fit@model_pars
m0_fit

# check some traceplots
rstan::traceplot(m0_fit, pars = "a")

# plot the model predictions
post <- rstan::extract(m0_fit)

mu0 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_obs, 
         post$a[, , S[i]][, C[i]] + post$b1[, S[i]]* Y[i])
  x <- exp(x)
  x <- rgamma(n = length(x), shape = (x*x)/post$phi, rate = (x)/post$phi)
  
  # get the probability of 0
  y <- 
    with(df_obs, 
         post$a_hu[, S[i]] + post$b1_hu[, S[i]]* Y[i] )
  y <- 1 - plogis(y)
  y <- rbinom(n = length(y), size = 1, prob = y)
  
  mu0[[i]] <- (x*y)
  
}

# bind into a matrix
mu0 <- do.call("cbind", mu0)

apply(mu0, 2, max) %>%
  max()

# plot the predicted values
plot(apply(mu0[, n_obs_mono], 2, mean), df_obs$M[n_obs_mono])
abline(0, 1)

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_0 <- loo::extract_log_lik(m0_fit, merge_chains = F)
r_eff_0 <- loo::relative_eff(log_lik_0)

# calculate the loocv estimating using PSIS
loo_0 <- rstan::loo(log_lik_0, r_eff = r_eff_0)
print(loo_0)

# check individual points
k_high0 <- which(pareto_k_influence_values(loo_0) > 0.7)

# check the data with high k-values
View(df_obs[k_high0, ])

# calculate the r2 value
# calculate the r2 value
r <- apply(mu0[,n_obs_mono], 2, mean) - df_obs[n_obs_mono,]$M
r <- 1 - (var2(r)/var2(df_obs[n_obs_mono,]$M))
print(r)

max(mu0)

# model 2

# model 0

# compile the model
m0 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_gamma_model2.stan",
                        verbose = TRUE)

# sample the stan model
m0_fit <- rstan::sampling(m0, data = df_m_obs, 
                          iter = 2000, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 12),
                          cores = 4)

# check the model
m0_fit@model_pars
m0_f


# check some traceplots
rstan::traceplot(m0_fit, pars = "a")

# plot the model predictions
post <- rstan::extract(m0_fit)

post$lp__

mu <- post$mu[1,]
mu

y <- rbinom(n = length(post$hu[1,]), size = 1, prob = plogis(post$hu[1,]))

phi <- post$phi[1]

x <- round(rgamma(n = length(mu), shape = (exp(mu)*exp(mu))/phi, rate = exp(mu)/phi  ), 2)

z <- x*y

z

mu0 <- vector("list", length = df_m_obs$N)
for(i in 1:df_m_obs$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(post, 
         rgamma(n = length(exp(mu[,i]) ), shape = (exp(mu[,i]) * exp(mu[, i]) )/phi, rate = exp(mu[,i])/phi ))
  
  y <- 1- rbinom(n = length(post$hu[,i]), size = 1, prob = plogis(post$hu[,i]))
  
  mu0[[i]] <- (x*y)
  
}

# bind into a matrix
mu0 <- do.call("cbind", mu0)

# plot the predicted values
plot(apply(mu0, 2, mean), df_m_obs$M)
abline(0, 1)




