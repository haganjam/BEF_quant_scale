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

# set a seed for reproducibility
set.seed(4597)

# load the analysis data
df_obs <- read_csv("data/case_study_2/data_clean/biomass_env_analysis_data.csv")

# standardise the predictor variables in the overall data
df_obs$Y <- with(df_obs, (Y - min(Y))/(max(Y)-min(Y)) )
df_obs$PC1 <- with(df_obs, (PC1 - min(PC1))/(max(PC1)-min(PC1)) )
df_obs$PC2 <- with(df_obs, (PC2 - min(PC2))/(max(PC2)-min(PC2)) )
df_obs$C <- as.integer(as.factor(df_obs$cluster_id))
df_obs$T <- as.integer(as.factor(df_obs$time))/3
df_obs$S <- as.integer(as.factor(df_obs$OTU))

# subset out the data with missing monocultures
df_m_obs <- 
  df_obs %>%
  filter(!is.na(M))

# function to get a list of model parameters
# stan_object: stan model object
stan_pars <- function(stan_object) {
  
  # check the traceplots
  pars <- stan_object@model_pars
  pars <- pars[!(grepl("Rho", pars) | grepl("Z", pars) | 
                   grepl("V", pars) | grepl("v", pars) | 
                   grepl("z", pars)  | pars == "mu" | 
                   pars == "hu" | pars == "log_lik" | pars == "lp__")]
  
  return(pars)
  
}

# function to get the loocv estimate using PSIS
# stan_object: stan model object
loo_est <- function(stan_object) {
  
  # calculate the PSIS loocv estimate
  # ref: http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/
  log_lik <- loo::extract_log_lik(stan_object, merge_chains = F)
  r_eff <- loo::relative_eff(log_lik)
  
  # calculate the loocv estimating using PSIS
  loo_score <- rstan::loo(log_lik, r_eff = r_eff)
  
  return(loo_score)
  
}

# function to obtain model predictions
# samples: samples from the posterior distribution
# data: dataset to predict from
# mu_dist: lognormal or gamma
stan_predict <- function(samples, data, mu_dist = "lognormal") {
  
  pred_list <- vector("list", length = nrow(data))
  for(i in 1:nrow(data)) {
    
    # get the mean prediction from the lognormal model on the natural scale
    x <- samples[["mu"]][,i]
    
    if(mu_dist == "lognormal") {
      x <- rlnorm(n = length(x), x, samples[["sigma"]])
    } else if(mu_dist == "gamma") {
      x <- exp(x)
      x <- rgamma(n = length(x), shape = (x*x)/samples[["phi"]], rate = (x)/samples[["phi"]])
    }
    
    # get the probability of 0
    y <- samples[["hu"]][,i]
    y <- 1 - plogis(y)
    y <- rbinom(n = length(y), size = 1, prob = y)
    
    pred_list[[i]] <- (x*y)
    
  }
  
  # bind into a matrix
  pred_mat <- do.call("cbind", pred_list)
  
  return(pred_mat)
  
}


# Lognormal hurdle models

# model 1
ln1_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model1_fit.rds")

# get the parameter names
pars <- stan_pars(stan_object = ln1_fit)

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(ln1_fit, pars = par_sel)

# get the loo estimate
loo_est(stan_object = ln1_fit)

# plot the model predictions
post <- rstan::extract(ln1_fit)

# get the predicted values from the stan model object
pred_ln1 <- stan_predict(samples = post, data = df_m_obs, mu_dist = "lognormal")

# plot the predicted values
plot(apply(pred_ln1, 2, mean), df_m_obs$M)
abline(0, 1)

# check the maximum prediction
max(apply(pred_ln1, 2, max))

# model 2
ln2_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model2_fit.rds")

# get the parameter names
pars <- stan_pars(stan_object = ln2_fit)

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(ln2_fit, pars = par_sel)

# get the loo estimate
loo_est(stan_object = ln2_fit)

# plot the model predictions
post <- rstan::extract(ln2_fit)

# get the predicted values from the stan model object
pred_ln2 <- stan_predict(samples = post, data = df_m_obs, mu_dist = "lognormal")

# plot the predicted values
plot(apply(pred_ln2, 2, mean), df_m_obs$M)
abline(0, 1)

# check the maximum prediction
max(apply(pred_ln2, 2, max))

# model 3
ln3_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model3_fit.rds")

# get the parameter names
pars <- stan_pars(stan_object = ln3_fit)

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(ln3_fit, pars = par_sel)

# get the loo estimate
loo_est(stan_object = ln3_fit)

# plot the model predictions
post <- rstan::extract(ln3_fit)

# get the predicted values from the stan model object
pred_ln3 <- stan_predict(samples = post, data = df_m_obs, mu_dist = "lognormal")

# plot the predicted values
plot(apply(pred_ln3, 2, mean), df_m_obs$M)
abline(0, 1)

# check the maximum prediction
max(apply(pred_ln3, 2, max))

# model 4
ln4_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model4_fit.rds")

# get the parameter names
pars <- stan_pars(stan_object = ln4_fit)

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(ln4_fit, pars = par_sel)

# get the loo estimate
loo_est(stan_object = ln4_fit)

# plot the model predictions
post <- rstan::extract(ln4_fit)

# get the predicted values from the stan model object
pred_ln4 <- stan_predict(samples = post, data = df_m_obs, mu_dist = "lognormal")

# plot the predicted values
plot(apply(pred_ln4, 2, mean), df_m_obs$M)
abline(0, 1)

# check the maximum prediction
max(apply(pred_ln4, 2, max))

# model 5
ln5_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model5_fit.rds")

# get the parameter names
pars <- stan_pars(stan_object = ln5_fit)

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(ln5_fit, pars = par_sel)

# get the loo estimate
loo_est(stan_object = ln5_fit)

# plot the model predictions
post <- rstan::extract(ln5_fit)

# get the predicted values from the stan model object
pred_ln5 <- stan_predict(samples = post, data = df_m_obs, mu_dist = "lognormal")

# plot the predicted values
plot(apply(pred_ln5, 2, mean), df_m_obs$M)
abline(0, 1)

# check the maximum prediction
max(apply(pred_ln5,2, max))

# model 6
ln6_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model6_fit.rds")

# get the parameter names
pars <- stan_pars(stan_object = ln6_fit)

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(ln6_fit, pars = par_sel)

# get the loo estimate
loo_est(stan_object = ln6_fit)

# plot the model predictions
post <- rstan::extract(ln6_fit)

# get the predicted values from the stan model object
pred_ln6 <- stan_predict(samples = post, data = df_m_obs, mu_dist = "lognormal")

# plot the predicted values
plot(apply(pred_ln6, 2, mean), df_m_obs$M)
abline(0, 1)

# check the maximum prediction
max(apply(pred_ln6,2, max))

# model 7
ln7_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model7_fit.rds")

# get the parameter names
pars <- stan_pars(stan_object = ln7_fit)

# check the traceplots
par_sel <- sample(pars, 1)
print(par_sel)
rstan::traceplot(ln7_fit, pars = par_sel)

# get the loo estimate
loo_est(stan_object = ln7_fit)

# plot the model predictions
post <- rstan::extract(ln7_fit)

# get the predicted values from the stan model object
pred_ln7 <- stan_predict(samples = post, data = df_m_obs, mu_dist = "lognormal")

# plot the predicted values
plot(apply(pred_ln7, 2, mean), df_m_obs$M)
abline(0, 1)

# check the maximum prediction
max(apply(pred_ln7, 2, max))

# null model
ln_null_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model_null_fit.rds")

# compare the different models by approximating leave-one-out CV

# get a vector of model names
mod_names <- c(paste0("ln", 1:7, "_fit"), "ln_null_fit")

loo_fit <- vector("list", length = length(mod_names))
for(i in 1:length(mod_names)) {
  
  # extract the model fit
  mod_fit <- eval(parse(text = mod_names[i]))
  
  # put the loo estimate into a list
  loo_fit[[i]] <- loo_est(mod_fit)
  
}
names(loo_fit) <- mod_names

# compare the different models using the loo score
loo_compare(loo_fit)

loo_fit[[1]]$estimates
loo_fit[[2]]$estimates

# which model is the best fit?
# ln1

# extract the posterior distribution
post <- rstan::extract(ln1_fit)

# compute the total number of samples
N <- (ln1_fit@stan_args[[1]]$iter - ln1_fit@stan_args[[1]]$warmup)*length(ln1_fit@stan_args)

# get 1000 samples from the posterior distribution
id_samp <- sample(1:N, 1000)

# make a fit to sample plot of the best fitting model
pred_ln1 <- vector("list", length = nrow(df_obs))
for(i in 1:nrow(df_obs)) {
  
  # get predictions for all values using the ln0 model
  x <- 
    with(df_obs, 
         post$a[id_samp, S[i]] + 
           post$b1[id_samp, S[i]]* Y[i] + 
           post$b2[id_samp, S[i]]* PC1[i])
  
  x <- rlnorm(n = length(x), x, post$sigma)
  
  # get the probability of zero for each observation
  y <- 
    with(df_obs, 
         post$a_hu[id_samp, S[i]] + 
           post$b1_hu[id_samp, S[i]]* Y[i] )
  
  y <- 1 - plogis(y)
  y <- rbinom(n = length(y), size = 1, prob = y)
  
  # multiply the model predictions together
  pred_ln1[[i]] <- (x*y)
  
}

# bind into a matrix
pred_ln1 <- do.call("cbind", pred_ln1)

# save this output
saveRDS(pred_ln1, "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/04_m1_predictions.rds")

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
                      M_pred_mu = apply(pred_ln1, 2, mean),
                      M_pred_PIlow = apply(pred_ln1, 2, HPDI, 0.90)[1,],
                      M_pred_PIhigh = apply(pred_ln1, 2, HPDI, 0.90)[2,])


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
  ylab("Predicted monoculture dry biomass (g)") +
  xlab("Observed monoculture dry biomass (g)") +
  facet_wrap(~C, scales = "free", nrow = 4, ncol = 3) +
  scale_y_continuous(limits = c(0, 31)) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_colour_manual(values = viridis(n = 5, begin = 0.1, end = 0.9, option = "C")) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p1)

ggsave(filename = "figures/figA2_SX.png", p1, dpi = 300,
       units = "cm", width = 18, height = 24)

# check the overlap between the predicted values and the observed values
n <- 25
samples <- pred_ln1[sample(1:nrow(pred_ln1), n), ]

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
p2 <- 
  ggplot() +
  geom_density(data = df_plot %>% filter(Obs_pred == "Observed"),
               mapping = aes(x = M_obs), fill = "red", colour = "white", alpha = 1) +
  geom_density(data = df_samples,
               mapping = aes(x = M_pred, group = sample), alpha = 0.05,
               linewidth = 0.25, fill = "grey", colour = "grey") +
  ylab("Density") +
  xlab("Monoculture dry biomass (g)") +
  facet_wrap(~ OTU, scales = "free", nrow = 3, ncol = 2) +
  theme_meta()
plot(p2)

ggsave(filename = "figures/figA2_SY.png", p2, dpi = 300,
       units = "cm", width = 12, height = 18)

# calculate how many zeros are predicted on the observed data
sum(df_plot$M_obs[!is.na(df_plot$M_obs)] == 0)/sum(!is.na(df_plot$M_obs))

pred_zero <- 
  apply(pred_ln1[,!is.na(df_plot$M_obs)], 2, function(x) {
  
  sum(x == 0)/length(x)
  
})
mean(pred_zero)

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

# how how many zeros there are, on average, for each sample from the posterior
mean(apply(pred_ln1, 1, function(x) sum(x == 0)))/ncol(pred_ln1)
sum(df_m_obs$M == 0)/nrow(df_m_obs)

### END
