#'
#' @title: Model missing monocultures
#' 
#' @description: Evaluate the log-normal hurdle generalised linear models in Stan 
#' that we will use to impute the missing monoculture data .
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstan)
library(loo)

# load plotting theme
source("scripts/Function_plotting_theme.R")
source("scripts/03_empirical_analysis/helper_functions.R")

# load the analysis data
df_obs <- read_csv("data/case_study_2/data_clean/biomass_env_analysis_data.csv")

# standardise the predictor variables in the overall data
df_obs$Y <- with(df_obs, (Y - min(Y))/(max(Y)-min(Y)) )
df_obs$PC1 <- with(df_obs, (PC1 - min(PC1))/(max(PC1)-min(PC1)) )
df_obs$PC2 <- with(df_obs, (PC2 - min(PC2))/(max(PC2)-min(PC2)) )
df_obs$C <- as.integer(as.factor(df_obs$cluster_id))
df_obs$T <- as.integer(as.factor(df_obs$time))/3
df_obs$S <- as.integer(as.factor(df_obs$OTU))

# which rows have observed monoculture
n_obs_mono <- which(!is.na(df_obs$M))

# model 0
m0_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model0_fit.rds")

# check the traceplots
pars <- m0_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars) | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(m0_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- rstan::summary(m0_fit)
diag$summary[grepl(par_sel, row.names(diag$summary)), ]

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
View(v[k_high0, ])

# plot the model predictions
post <- rstan::extract(m0_fit)

# use the posterior predictive distribution
ppd = TRUE

mu0 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_obs, 
         post$a[, S[i]] + post$b1[, S[i]]* Y[i] + post$b2[, S[i]]* PC1[i])
  
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_obs, 
         post$a_hu[, S[i]] + post$b1_hu[, S[i]]* Y[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu0[[i]] <- (x*y)
  
}

# bind into a matrix
mu0 <- do.call("cbind", mu0)

# plot the predicted values
plot(apply(mu0[, n_obs_mono], 2, mean), df_obs[n_obs_mono,]$M)
points(apply(mu0[, n_obs_mono], 2, mean)[k_high0], df_obs[n_obs_mono,]$M[k_high0], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu0[,n_obs_mono], 2, mean) - df_obs[n_obs_mono,]$M
r <- 1 - (var2(r)/var2(df_obs[n_obs_mono,]$M))
print(r)

# check the overall predicted distributions
max(apply(mu0, 2, mean))
max(apply(mu0, 2, max))

# model 1
m1_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model1_fit.rds")

# check the traceplots
pars <- m1_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars) | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(m1_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- rstan::summary(m1_fit)
diag$summary[grepl(par_sel, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_1 <- loo::extract_log_lik(m1_fit, merge_chains = F)
r_eff_1 <- loo::relative_eff(log_lik_1)

# calculate the loocv estimating using PSIS
loo_1 <- rstan::loo(log_lik_1, r_eff = r_eff_1)
print(loo_1)

# check individual points
k_high1 <- which(pareto_k_influence_values(loo_1) > 0.7)

# check the data with high k-values
View(v[k_high1, ])

# plot the model predictions
post <- extract(m1_fit)

# use the posterior predictive distribution
ppd <- TRUE

mu1 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {

  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_obs, 
         post$a[, S[i]] + post$b1[, S[i]]* Y[i] + post$b2[, S[i]]* PC1[i])
  
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_obs, 
         post$a_hu + post$b1_hu * Y[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu1[[i]] <- (x*y)
  
}

# bind into a matrix
mu1 <- do.call("cbind", mu1)

# plot the predicted values
plot(apply(mu1[, n_obs_mono], 2, mean), df_obs[n_obs_mono,]$M)
points(apply(mu1[, n_obs_mono], 2, mean)[k_high1], df_obs[n_obs_mono,]$M[k_high1], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu1[,n_obs_mono], 2, mean) - df_obs[n_obs_mono,]$M
r <- 1 - (var2(r)/var2(df_obs[n_obs_mono,]$M))
print(r)

# check the overall predicted distributions
max(apply(mu1, 2, mean))
max(apply(mu1, 2, max))

# model 2
m2_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model2_fit.rds")

# check the traceplots
pars <- m2_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars) | pars == "v" | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
traceplot(m2_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- rstan::summary(m2_fit)
diag$summary[grepl(par_sel, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_2 <- loo::extract_log_lik(m2_fit, merge_chains = F)
r_eff_2 <- loo::relative_eff(log_lik_2)

# calculate the loocv estimating using PSIS
loo_2 <- rstan::loo(log_lik_2, r_eff = r_eff_2)
print(loo_2)

# check individual points
k_high2 <- which(pareto_k_influence_values(loo_2) > 0.7)

# check the data with high k-values
View(v[k_high2, ])

# plot the model predictions
post <- extract(m2_fit)

# use the posterior predictive distribution
ppd <- TRUE

mu2 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_obs, 
         post$a[, S[i]] + post$b1[, S[i]]* Y[i])
  
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_obs, 
         post$a_hu + post$b1_hu * Y[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu2[[i]] <- (x*y)
  
}

# bind into a matrix
mu2 <- do.call("cbind", mu2)

# plot the predicted values
plot(apply(mu2[, n_obs_mono], 2, mean), df_obs[n_obs_mono,]$M)
points(apply(mu2[, n_obs_mono], 2, mean)[k_high2], df_obs[n_obs_mono,]$M[k_high2], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu2[,n_obs_mono], 2, mean) - df_obs[n_obs_mono,]$M
r <- 1 - (var2(r)/var2(df_obs[n_obs_mono,]$M))
print(r)

# check the overall predicted distributions
max(apply(mu2, 2, mean))
max(apply(mu2, 2, max))

# model 3
m3_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model3_fit.rds")

# check the traceplots
pars <- m3_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars)  | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
traceplot(m3_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- rstan::summary(m3_fit)
diag$summary[grepl(par_sel, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_3 <- loo::extract_log_lik(m3_fit, merge_chains = F)
r_eff_3 <- loo::relative_eff(log_lik_3)

# calculate the loocv estimating using PSIS
loo_3 <- rstan::loo(log_lik_3, r_eff = r_eff_3)
print(loo_3)

# check individual points
k_high3 <- which(pareto_k_influence_values(loo_3) > 0.7)

# check the data with high k-values
View(v[k_high3, ])

# plot the model predictions
post <- extract(m3_fit)

# use the posterior predictive distribution
ppd <- TRUE

mu3 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_obs, 
         post$a[, S[i]] + post$b1[, S[i]]* Y[i])
  
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_obs, 
         post$a_hu + post$b1_hu * Y[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu3[[i]] <- (x*y)
  
}

# bind into a matrix
mu3 <- do.call("cbind", mu3)

# plot the predicted values
plot(apply(mu3[, n_obs_mono], 2, mean), df_obs[n_obs_mono,]$M)
points(apply(mu3[, n_obs_mono], 2, mean)[k_high3], df_obs[n_obs_mono,]$M[k_high3], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu3[,n_obs_mono], 2, mean) - df_obs[n_obs_mono,]$M
r <- 1 - (var2(r)/var2(df_obs[n_obs_mono,]$M))
print(r)

# check the overall predicted distributions
max(apply(mu3, 2, mean))
max(apply(mu3, 2, max))

# model 4
m4_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model4_fit.rds")

# check the traceplots
pars <- m4_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars)  | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
traceplot(m4_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- rstan::summary(m4_fit)
diag$summary[grepl(par_sel, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_4 <- loo::extract_log_lik(m4_fit, merge_chains = F)
r_eff_4 <- loo::relative_eff(log_lik_4)

# calculate the loocv estimating using PSIS
loo_4 <- rstan::loo(log_lik_4, r_eff = r_eff_4)
print(loo_4)

# check individual points
k_high4 <- which(pareto_k_influence_values(loo_4) > 0.7)

# check the data with high k-values
View(v[k_high4, ])

# plot the model predictions
post <- extract(m4_fit)

# use the posterior predictive distribution
ppd <- TRUE

mu4 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_obs, 
         post$a[, S[i]] + post$b1* Y[i] + post$b2* PC1[i])
  
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_obs, 
         post$a_hu + post$b1_hu * Y[i] + post$b2_hu* PC1[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu4[[i]] <- (x*y)
  
}

# bind into a matrix
mu4 <- do.call("cbind", mu4)

# plot the predicted values
plot(apply(mu4[, n_obs_mono], 2, mean), df_obs[n_obs_mono,]$M)
points(apply(mu4[, n_obs_mono], 2, mean)[k_high4], df_obs[n_obs_mono,]$M[k_high4], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu4[,n_obs_mono], 2, mean) - df_obs[n_obs_mono,]$M
r <- 1 - (var2(r)/var2(df_obs[n_obs_mono,]$M))
print(r)

# check the overall predicted distributions
max(apply(mu4, 2, mean))
max(apply(mu4, 2, max))

# model 5
m5_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model5_fit.rds")

# check the traceplots
pars <- m5_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars)  | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
traceplot(m5_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- rstan::summary(m5_fit)
diag$summary[grepl(par_sel, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_5 <- loo::extract_log_lik(m5_fit, merge_chains = F)
r_eff_5 <- loo::relative_eff(log_lik_5)

# calculate the loocv estimating using PSIS
loo_5 <- rstan::loo(log_lik_5, r_eff = r_eff_5)
print(loo_5)

# check individual points
k_high5 <- which(pareto_k_influence_values(loo_5) > 0.7)

# check the data with high k-values
View(v[k_high5, ])

# plot the model predictions
post <- extract(m5_fit)

# use the posterior predictive distribution
ppd <- TRUE

mu5 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_obs, 
         post$a[, S[i]] + post$b1* Y[i] + post$b2* PC1[i])
  
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_obs, 
         post$a_hu + post$b1_hu * Y[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu5[[i]] <- (x*y)
  
}

# bind into a matrix
mu5 <- do.call("cbind", mu5)

# plot the predicted values
plot(apply(mu5[, n_obs_mono], 2, mean), df_obs[n_obs_mono,]$M)
points(apply(mu5[, n_obs_mono], 2, mean)[k_high5], df_obs[n_obs_mono,]$M[k_high5], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu5[,n_obs_mono], 2, mean) - df_obs[n_obs_mono,]$M
r <- 1 - (var2(r)/var2(df_obs[n_obs_mono,]$M))
print(r)

# check the overall predicted distributions
max(apply(mu5, 2, mean))
max(apply(mu5, 2, max))

# model 6
m6_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model6_fit.rds")

# check the traceplots
pars <- m6_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars)  | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
traceplot(m6_fit, pars = par_sel)

# extract the diagnostic parameters
diag <- rstan::summary(m6_fit)
diag$summary[grepl(par_sel, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_6 <- loo::extract_log_lik(m6_fit, merge_chains = F)
r_eff_6 <- loo::relative_eff(log_lik_6)

# calculate the loocv estimating using PSIS
loo_6 <- rstan::loo(log_lik_6, r_eff = r_eff_6)
print(loo_6)

# check individual points
k_high6 <- which(pareto_k_influence_values(loo_6) > 0.7)

# check the data with high k-values
View(v[k_high6, ])

# plot the model predictions
post <- extract(m6_fit)

# use the posterior predictive distribution
ppd <- TRUE

mu6 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_obs, 
         post$a[, S[i]] + post$b1* Y[i])
  
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_obs, 
         post$a_hu + post$b1_hu * Y[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu6[[i]] <- (x*y)
  
}

# bind into a matrix
mu6 <- do.call("cbind", mu6)

# plot the predicted values
plot(apply(mu6[, n_obs_mono], 2, mean), df_obs[n_obs_mono,]$M)
points(apply(mu6[, n_obs_mono], 2, mean)[k_high6], df_obs[n_obs_mono,]$M[k_high6], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu6[,n_obs_mono], 2, mean) - df_obs[n_obs_mono,]$M
r <- 1 - (var2(r)/var2(df_obs[n_obs_mono,]$M))
print(r)

# check the overall predicted distributions
max(apply(mu6, 2, mean))
max(apply(mu6, 2, max))

# compare the different models
loo_compare(list( "m0" = loo_0, 
                  "m1" = loo_1, 
                  "m2" = loo_2, 
                  "m3" = loo_3, 
                  "m4" = loo_4, 
                  "m5" = loo_5, 
                  "m6" = loo_6) )

# are similar points leading to the high values
k_mult <- table(c(k_high0, k_high1, k_high2, k_high3, k_high4, k_high5, k_high6))
df_obs[as.integer(names(k_mult[k_mult>1])),]

# make a fit to sample plot of the best fitting model: m0

# pull into a data.frame for plotting

# rename the cluster variable
C <- factor(df_obs$cluster_id)
levels(C) <- paste0("Cluster ", LETTERS[1:10])

# add an OTU column
OTU <- factor(df_obs$S)
levels(OTU) <- c("Barn", "Bryo", "Asci", "Hydro", "Ciona")

# make a data.frame for plotting
df_plot <- data.frame(M_obs = df_obs$M,
                      Obs_pred = ifelse(is.na(df_obs$M), "Predicted", "Observed"),
                      S = as.character(df_obs$S),
                      OTU = OTU,
                      C = C,
                      T = as.character(df_obs$time),
                      Y = df_obs$Y,
                      PC1 = df_obs$PC1,
                      PC2 = df_obs$PC2,
                      M_pred_mu = apply(mu0, 2, mean),
                      M_pred_PIlow = apply(mu0, 2, HPDI, 0.90)[1,],
                      M_pred_PIhigh = apply(mu0, 2, HPDI, 0.90)[2,])


# check the min and max
summary(df_plot)

# plot fit-to-sample
p1 <- 
  ggplot(data = df_plot %>% filter(!is.na(M_obs))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  geom_point(mapping = aes(x = M_obs, y = M_pred_mu, colour = OTU),
             shape = 16, alpha = 0.75) + 
  geom_errorbar(mapping = aes(x = M_obs, ymin = M_pred_PIlow, ymax = M_pred_PIhigh, 
                              colour = OTU),
                width = 0, alpha = 0.5, size = 0.25, show.legend = FALSE) +
  ylab("Predicted monoculture (g)") +
  xlab("Observed monoculture (g)") +
  facet_wrap(~C, scales = "free", nrow = 5, ncol = 2) +
  scale_y_continuous(limits = c(0, 31)) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_colour_manual(values = viridis(n = 5, begin = 0.1, end = 0.9, option = "C")) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p1)

# check the overlap between the predicted values and the observed values
n <- 25
samples <- mu0[sample(1:nrow(mu0), n), ]

df_samples <- data.frame(V1 = samples[1,])
for(i in 2:nrow(samples)) {
  df_samples[[i]] <- samples[i,]
}

# add identifier columns
df_samples$Obs_pred <- df_plot$Obs_pred
df_samples$OTU <- df_plot$OTU

# filter the predicted samples
df_samples <- 
  df_samples %>%
  filter(Obs_pred == "Predicted")

# pull into the long-format
df_samples <- 
  df_samples %>%
  pivot_longer(cols = paste0("V", 1:n),
               names_to = "sample",
               values_to = "M_pred")

# plot the overlap between the observed and the predicted data
ggplot() +
  geom_density(data = df_plot %>% filter(Obs_pred == "Observed"),
               mapping = aes(x = M_obs), fill = "red", colour = "white", alpha = 1,
               n = 64) +
  geom_density(data = df_samples,
               mapping = aes(x = M_pred, group = sample), alpha = 0.05,
               linewidth = 0.25, n = 64, fill = "grey", colour = "grey") +
  facet_wrap(~ OTU, scales = "free") +
  theme_meta()

# calculate the MESS index: Zurell et al. (2012)

# create a reference dataset
df_ref <- 
  df_plot %>% 
  filter(Obs_pred == "Observed") %>%
  dplyr::select(Y, PC1)

# create a test dataset
df_test <- 
  df_plot %>% 
  filter(Obs_pred == "Predicted") %>%
  dplyr::select(Y, PC1)

# calculate the MESS index: Zurell et al. (2012)
MESS <- eo.mask(traindata = df_ref, 
                newdata = df_test, nbin = 20, type="EO")

# add the MESS index to the test data
df_test$MESS <- MESS

# check how many data point are predicting out of the range
(sum(df_test$MESS)/nrow(df_test))*100


### END
