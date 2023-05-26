#'
#' @title: Model the missing monocultures
#' 
#' @description: This script attempts to impute the missing monoculture
#' data using Bayesian generalised linear models.
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

# remove any NAs
v <- data[complete.cases(data),]

# calculate the minimum observed biomass
# min_M <- min(v[v$M>0,]$M)

# add this minimum to we can model without zero values
# v$M <- v$M + min_M

# make a data.list with the training data
spp <- 
  list(M = v$M,
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
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 12),
                          seed = 54856)

# check the stan output
print(m1_fit)

# check the traceplots
traceplot(m1_fit, pars = c("abar[1]", "abar[2]", "abar[3]", "abar[4]"))

# extract the diagnostic parameters
diag <- rstan::summary(m1_fit)
par <- "bar"
diag$summary[grepl(par, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_1 <- loo::extract_log_lik(m1_fit, merge_chains = F)
r_eff_1 <- loo::relative_eff(log_lik_1)

# calculate the loocv estimating using PSIS
loo_1 <- rstan::loo(log_lik_1, r_eff = r_eff_1)
print(loo_1)

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

# check the stan output
print(m2_fit)

# check the traceplots
traceplot(m2_fit, pars = c("abar[1]", "abar[2]", "abar[3]", "abar[4]"))

# extract the diagnostic parameters
diag <- rstan::summary(m2_fit)
par <- "bar"
diag$summary[grepl(par, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_2 <- loo::extract_log_lik(m2_fit, merge_chains = F)
r_eff_2 <- loo::relative_eff(log_lik_2)

# calculate the loocv estimating using PSIS
loo_2 <- rstan::loo(log_lik_2, r_eff = r_eff_2)

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

# check the stan output
print(m3_fit)

# check the traceplots
traceplot(m3_fit, pars = c("abar[1]", "abar[2]", "abar[3]", "abar[4]"))

# extract the diagnostic parameters
diag <- rstan::summary(m3_fit)
par <- "bar"
diag$summary[grepl(par, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_3 <- loo::extract_log_lik(m3_fit, merge_chains = F)
r_eff_3 <- loo::relative_eff(log_lik_3)

# calculate the loocv estimating using PSIS
loo_3 <- rstan::loo(log_lik_3, r_eff = r_eff_3)

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

# check the stan output
print(m4_fit)

# check the traceplots
traceplot(m4_fit, pars = c("abar[1]", "abar[2]", "abar[3]", "abar[4]"))

# extract the diagnostic parameters
diag <- rstan::summary(m4_fit)
par <- "bar"
diag$summary[grepl(par, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_4 <- loo::extract_log_lik(m4_fit, merge_chains = F)
r_eff_4 <- loo::relative_eff(log_lik_4)

# calculate the loocv estimating using PSIS
loo_4 <- rstan::loo(log_lik_4, r_eff = r_eff_4)

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

# check the stan output
print(m5_fit)

# check the traceplots
traceplot(m5_fit, pars = c("abar[1]", "abar[2]", "abar[3]", "abar[4]"))

# extract the diagnostic parameters
diag <- rstan::summary(m5_fit)
par <- "bar"
diag$summary[grepl(par, row.names(diag$summary)), ]

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_5 <- loo::extract_log_lik(m5_fit, merge_chains = F)
r_eff_5 <- loo::relative_eff(log_lik_5)

# calculate the loocv estimating using PSIS
loo_5 <- rstan::loo(log_lik_5, r_eff = r_eff_5)


# try the gamma distribution

# compile the model
mx <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_modelgamma.stan",
                        verbose = TRUE)

# sample the stan model
mx_fit <- rstan::sampling(mx, data = spp, 
                          iter = 1000, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 12),
                          seed = 54856)

# check the stan output
print(mx_fit)

# check the traceplots
traceplot(mx_fit, pars = c("abar[1]", "abar[2]", "abar[3]", "abar[4]"))

# calculate the PSIS loocv estimate
# ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
log_lik_x <- loo::extract_log_lik(mx_fit, merge_chains = F)
r_eff_x <- loo::relative_eff(log_lik_x)

# calculate the loocv estimating using PSIS
loo_x <- rstan::loo(log_lik_x, r_eff = r_eff_x)
loo_x


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
