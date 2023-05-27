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
library(ggplot2)
library(rstan)
library(loo)

# load plotting theme
source("scripts/Function_plotting_theme.R")
source("scripts/03_empirical_analysis/helper_functions.R")

# load the analysis data
data <- read_csv("data/case_study_2/data_clean/biomass_env_analysis_data.csv")

# remove any NAs
v <- data[complete.cases(data),]

# model 1
m1_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model1_fit.rds")

# check the stan output
print(m1_fit)

# check the traceplots
pars <- m1_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | pars == "Z" | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
traceplot(m1_fit, pars = par_sel)

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

mu1 <- vector(length = spp$N)
for(i in 1:spp$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(spp, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]) + 
           (post$b1bar[, S[i]] + post$b1[, , S[i]][, C[i]] * T[i]) + 
           (post$b2bar[, S[i]] + post$b2[, , S[i]][, C[i]] * Y[i]) + 
           (post$b3bar[, S[i]] + post$b3[, , S[i]][, C[i]] * PC1[i]) + 
           (post$b4bar[, S[i]] + post$b4[, , S[i]][, C[i]] * PC2[i]) + 
           (post$b5bar[, S[i]] + post$b5[, , S[i]][, C[i]] * Y[i] * PC1[i]))
  x <- exp(x + (0.5*(post$sigma^2)))
  
  # get the probability of 0
  y <- 
    with(spp, 
         (post$abar_hu[, S[i]] + post$a_hu[, , S[i]][, C[i]]) + 
           (post$b1bar_hu[, S[i]] + post$b1_hu[, , S[i]][, C[i]] * T[i]) + 
           (post$b2bar_hu[, S[i]] + post$b2_hu[, , S[i]][, C[i]] * Y[i]) + 
           (post$b3bar_hu[, S[i]] + post$b3_hu[, , S[i]][, C[i]] * PC1[i]) + 
           (post$b4bar_hu[, S[i]] + post$b4_hu[, , S[i]][, C[i]] * PC2[i]) + 
           (post$b5bar_hu[, S[i]] + post$b5_hu[, , S[i]][, C[i]] * Y[i] * PC1[i]) )
  y <- plogis(y)
  y <- 1-y
  
  mu1[i] <- mean((x*y))
  
}

# plot the predicted values
plot(mu1, spp$M)
points(mu1[k_high1], spp$M[k_high1], col = "red")
abline(0, 1)

# calculate the r2 value
r <- mu1 - spp$M
r <- 1 - (var2(r)/var2(spp$M))
print(r)





