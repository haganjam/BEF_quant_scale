#'
#' @title: Model the missing monocultures
#' 
#' @description: Fit log-normal hurdle generalised linear models in Stan that we 
#' will use to impute the missing monoculture data .
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'
#'
#' do the standardisation before we predict the new data because this
#' will influence the later predictions
#'
#' for reference, this is the simplest centered model coded in ulam()
#' 
#' mx <- 
#' ulam(
#'    alist(M ~ dlnorm(u, sigma),
#'          
#'          u <- a[C, S],
#'          
#'          # priors
#'          matrix[10, 5]:a ~ multi_normal(abar, Rho_a, sigma_a),
#'          vector[5]:abar ~ dlnorm(0, 2),
#'          vector[5]:sigma_a ~ dexp(1),
#'          corr_matrix[5]:Rho_a ~ lkj_corr(2),
#'          
#'          sigma ~ dexp(1)
#'          
#'    ),
#'    data = spp, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))
#'    

# load the required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(rstan)
library(loo)

# load plotting theme
source("scripts/Function_plotting_theme.R")

# load the analysis data
data <- read_csv("data/case_study_2/data_clean/biomass_env_analysis_data.csv")

# standardise the predictor variables in the overall data
data$Y <- with(data, (Y - min(Y))/(max(Y)-min(Y)) )
data$PC1 <- with(data, (PC1 - min(PC1))/(max(PC1)-min(PC1)) )
data$PC2 <- with(data, (PC2 - min(PC2))/(max(PC2)-min(PC2)) )
data$cluster_id <- as.integer(as.factor(data$cluster_id))
data$time <- as.integer(as.factor(data$time))/3
data$OTU <- as.integer(as.factor(data$OTU))

# remove any NAs
v <- data[complete.cases(data),]

# make a data.list with the training data
spp <- 
  list(N = length(v$M),
       S_N = length(unique(v$OTU)),
       C_N = length(unique(v$cluster_id)),
       M = v$M,
       Y = with(v, (Y - min(Y))/(max(Y)-min(Y)) ),
       PC1 = with(v, (PC1 - min(PC1))/(max(PC1)-min(PC1)) ),
       PC2 = with(v, (PC2 - min(PC2))/(max(PC2)-min(PC2)) ),
       C = as.integer(as.factor(v$cluster_id)),
       T = as.integer(as.factor(v$time))/3,
       S = as.integer(as.factor(v$OTU))
  )

# model 1

# compile the model
m1 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model1.stan",
                        verbose = TRUE)

# sample the stan model
m1_fit <- rstan::sampling(m1, data = spp, 
                          iter = 1500, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.95,
                                         max_treedepth = 12),
                          seed = 54856)

# save the stan model fit object
m1_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m1_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model1_fit.rds")

# model 2

# compile the model
m2 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model2.stan",
                        verbose = TRUE)

# sample the stan model
m2_fit <- rstan::sampling(m2, data = spp, 
                          iter = 1500, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 12),
                          seed = 54856)

# save the stan model fit object
m2_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m2_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model2_fit.rds")

# check the stan output
print(m2_fit)

# check the traceplots
pars <- m2_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | pars == "Z" | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
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

mu2 <- vector(length = spp$N)
for(i in 1:spp$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(spp, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]) + 
           (post$b1bar[, S[i]] + post$b1[, , S[i]][, C[i]] * T[i]) + 
           (post$b2bar[, S[i]] + post$b2[, , S[i]][, C[i]] * Y[i]) + 
           (post$b3bar[, S[i]] + post$b3[, , S[i]][, C[i]] * PC1[i]) + 
           (post$b4bar[, S[i]] + post$b4[, , S[i]][, C[i]] * PC2[i]))
  x <- exp(x + (0.5*(post$sigma^2)))
  
  # get the probability of 0
  y <- 
    with(spp, 
         (post$abar_hu[, S[i]] + post$a_hu[, , S[i]][, C[i]]) + 
           (post$b1bar_hu[, S[i]] + post$b1_hu[, , S[i]][, C[i]] * T[i]) + 
           (post$b2bar_hu[, S[i]] + post$b2_hu[, , S[i]][, C[i]] * Y[i]) + 
           (post$b3bar_hu[, S[i]] + post$b3_hu[, , S[i]][, C[i]] * PC1[i]) + 
           (post$b4bar_hu[, S[i]] + post$b4_hu[, , S[i]][, C[i]] * PC2[i]) )
  y <- plogis(y)
  y <- 1-y
  
  mu[i] <- mean((x*y))
  
}

# plot the predicted values
plot(mu2, spp$M)
points(mu2[k_high2], spp$M[k_high2], col = "red")
abline(0, 1)

# model 3

# compile the model
m3 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model3.stan",
                        verbose = TRUE)

# sample the stan model
m3_fit <- rstan::sampling(m3, data = spp, 
                          iter = 1500, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 12),
                          seed = 54856)

# save the stan model fit object
m3_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m3_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model3_fit.rds")

# check the stan output
print(m3_fit)

# check the traceplots
pars <- m3_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | pars == "Z" | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
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

mu3 <- vector(length = spp$N)
for(i in 1:spp$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(spp, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]) + 
           (post$b1bar[, S[i]] + post$b1[, , S[i]][, C[i]] * T[i]) + 
           (post$b2bar[, S[i]] + post$b2[, , S[i]][, C[i]] * Y[i]))
  x <- exp(x + (0.5*(post$sigma^2)))
  
  # get the probability of 0
  y <- 
    with(spp, 
         (post$abar_hu[, S[i]] + post$a_hu[, , S[i]][, C[i]]) + 
           (post$b1bar_hu[, S[i]] + post$b1_hu[, , S[i]][, C[i]] * T[i]) + 
           (post$b2bar_hu[, S[i]] + post$b2_hu[, , S[i]][, C[i]] * Y[i]) )
  y <- plogis(y)
  y <- 1-y
  
  mu3[i] <- mean((x*y))
  
}

# plot the predicted values
plot(mu3, spp$M)
points(mu3[k_high3], spp$M[k_high3], col = "red")
abline(0, 1)

# model 4

# compile the model
m4 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model4.stan",
                        verbose = TRUE)

# sample the stan model
m4_fit <- rstan::sampling(m4, data = spp, 
                          iter = 1500, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 12),
                          seed = 54856)

# save the stan model fit object
m4_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m4_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model4_fit.rds")

# check the traceplots
pars <- m4_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | pars == "Z" | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

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

mu4 <- vector(length = spp$N)
for(i in 1:spp$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(spp, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]) + 
           (post$b1bar[, S[i]] + post$b1[, , S[i]][, C[i]] * T[i]))
  x <- exp(x + (0.5*(post$sigma^2)))
  
  # get the probability of 0
  y <- 
    with(spp, 
         (post$abar_hu[, S[i]] + post$a_hu[, , S[i]][, C[i]]) + 
           (post$b1bar_hu[, S[i]] + post$b1_hu[, , S[i]][, C[i]] * T[i]) )
  y <- plogis(y)
  y <- 1-y
  
  mu4[i] <- mean((x*y))
  
}

# plot the predicted values
plot(mu4, spp$M)
points(mu4[k_high4], spp$M[k_high4], col = "red")
abline(0, 1)

# model 5

# compile the model
m5 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model5.stan",
                        verbose = TRUE)

# sample the stan model
m5_fit <- rstan::sampling(m5, data = spp, 
                          iter = 1500, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 12),
                          seed = 54856)

# save the stan model fit object
m5_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m5_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model5_fit.rds")

# check the stan output
print(m5_fit)

# check the traceplots
pars <- m5_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | pars == "Z" | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# check the traceplots
par_sel <- sample(pars, 1)
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

mu5 <- vector(length = spp$N)
for(i in 1:spp$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(spp, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]))
  x <- exp(x + (0.5*(post$sigma^2)))
  
  # get the probability of 0
  y <- 
    with(spp, 
         (post$abar_hu[, S[i]] + post$a_hu[, , S[i]][, C[i]]))
  y <- plogis(y)
  y <- 1-y
  
  mu5[i] <- mean((x*y))
  
}

# plot the predicted values
plot(mu5, spp$M)
points(mu5[k_high5], spp$M[k_high5], col = "red")
abline(0, 1)

# compare the different models
loo_compare(loo_1, loo_2, loo_3, loo_4, loo_5)


# check the best predictive model
post <- rethinking::extract.samples(m3)

# use the posterior predictive distribution
pred <- sim(m3)
pred <- pred + min_M

# use the average prediction
pred <- link(m3)
pred <- 
  apply(pred, 2, function(x) {
    
    exp(x + (0.5*(post$sigma^2)))
    
  } )

# pull into a data.frame for plotting
df.pred <- data.frame(M_obs = spp$M,
                      M_S = as.character(spp$S),
                      M_C = as.character(spp$C),
                      M_T = as.character(round(spp$T, 1)),
                      M_pred_mu = apply(pred, 2, mean),
                      M_pred_PIlow = apply(pred, 2, HPDI, 0.90)[1,],
                      M_pred_PIhigh = apply(pred, 2, HPDI, 0.90)[2,])

p1 <- 
  ggplot(data = df.pred,
         mapping = aes(x = M_obs, y = M_pred_mu, colour = M_S, fill = M_S)) +
  geom_point(shape = 16, alpha = 0.5) +
  geom_errorbar(mapping = aes(x = M_obs, ymin = M_pred_PIlow, ymax = M_pred_PIhigh),
                width = 0, alpha = 0.5, size = 0.25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  ylab("Predicted monoculture (g)") +
  xlab("Observed monoculture (g)") +
  facet_wrap(~M_C, scales = "free") +
  # scale_y_continuous(limits = c(0, 17)) +
  # scale_x_continuous(limits = c(0, 17)) +
  theme_meta()
plot(p1)


### END
