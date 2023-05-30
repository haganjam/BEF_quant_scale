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

# use the posterior predictive distribution
ppd <- TRUE

mu1 <- vector("list", length = spp$N)
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
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
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
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu1[[i]] <- (x*y)
  
}

# bind into a matrix
mu1 <- do.call("cbind", mu1)

# plot the predicted values
plot(apply(mu1, 2, mean), spp$M)
points(apply(mu1, 2, mean)[k_high1], spp$M[k_high1], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu1, 2, mean) - spp$M
r <- 1 - (var2(r)/var2(spp$M))
print(r)

# model 2
m2_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model2_fit.rds")

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

# use the posterior predictive distribution
ppd <- TRUE

mu2 <- vector("list", length = spp$N)
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
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  }
  
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
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu2[[i]] <- (x*y)
  
}

# bind into a matrix
mu2 <- do.call("cbind", mu2)

# plot the predicted values
plot(apply(mu2, 2, mean), spp$M)
points(apply(mu2, 2, mean)[k_high2], spp$M[k_high2], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu2, 2, mean) - spp$M
r <- 1 - (var2(r)/var2(spp$M))
print(r)

# model 3
m3_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model3_fit.rds")

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

# use the posterior predictive distribution
ppd <- TRUE

mu3 <- vector("list", length = spp$N)
for(i in 1:spp$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(spp, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]) + 
           (post$b1bar[, S[i]] + post$b1[, , S[i]][, C[i]] * T[i]) + 
           (post$b2bar[, S[i]] + post$b2[, , S[i]][, C[i]] * Y[i]))
  x <- exp(x + (0.5*(post$sigma^2)))
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  }
  
  # get the probability of 0
  y <- 
    with(spp, 
         (post$abar_hu[, S[i]] + post$a_hu[, , S[i]][, C[i]]) + 
           (post$b1bar_hu[, S[i]] + post$b1_hu[, , S[i]][, C[i]] * T[i]) + 
           (post$b2bar_hu[, S[i]] + post$b2_hu[, , S[i]][, C[i]] * Y[i]) )
  y <- plogis(y)
  y <- 1-y
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu3[[i]] <- (x*y)
  
}

# bind into a matrix
mu3 <- do.call("cbind", mu3)

# plot the predicted values
plot(apply(mu3, 2, mean), spp$M)
points(apply(mu3, 2, mean)[k_high3], spp$M[k_high3], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu3, 2, mean) - spp$M
r <- 1 - (var2(r)/var2(spp$M))
print(r)

# model 4
m4_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model4_fit.rds")

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

# use the posterior predictive distribution
ppd <- TRUE

mu4 <- vector("list", length = spp$N)
for(i in 1:spp$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(spp, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]) + 
           (post$b1[, S[i]] * Y[i]))
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(spp, 
         post$a_hu[, S[i]] * Y[i] + post$b1_hu[, S[i]] * Y[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu4[[i]] <- (x*y)
  
}

# bind into a matrix
mu4 <- do.call("cbind", mu4)

# plot the predicted values
plot(apply(mu4, 2, mean), spp$M)
points(apply(mu4, 2, mean)[k_high4], spp$M[k_high4], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu4, 2, mean) - spp$M
r <- 1 - (var2(r)/var2(spp$M))
print(r)

# model 5
m5_fit <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model5_fit.rds")

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

# use the posterior predictive distribution
ppd <- TRUE

mu5 <- vector("list", length = spp$N)
for(i in 1:spp$N) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(spp, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]))
  x <- exp(x + (0.5*(post$sigma^2)))
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  }
  
  # get the probability of 0
  y <- 
    with(spp, 
         (post$abar_hu[, S[i]] + post$a_hu[, , S[i]][, C[i]]))
  y <- plogis(y)
  y <- 1-y
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  mu5[[i]] <- (x*y)
  
}

# bind into a matrix
mu5 <- do.call("cbind", mu5)

# plot the predicted values
plot(apply(mu5, 2, mean), spp$M)
points(apply(mu5, 2, mean)[k_high5], spp$M[k_high5], col = "red")
abline(0, 1)

# calculate the r2 value
r <- apply(mu5, 2, mean) - spp$M
r <- 1 - (var2(r)/var2(spp$M))
print(r)

# compare the different models
loo_compare(loo_1, loo_2, loo_3, loo_4, loo_5)

# are similar points leading to the high values
k_mult <- table(c(k_high1, k_high2, k_high3, k_high4, k_high5))

v[as.integer(names(k_mult[k_mult>1])),]


# make a fit to sample plot of the best fitting model: m4

# extract samples from the posterior distribution
post <- rstan::extract(m4_fit)
id_high <- which(post$sigma > quantile(post$sigma, 0.95))
id_low <- which(post$sigma < quantile(post$sigma, 0.05))
id <- c(id_high, id_low)

# pull into a data.frame for plotting
df_pred <- data.frame(M = data$M,
                      S = data$OTU,
                      C = data$cluster_id,
                      T = data$time,
                      Y = data$Y,
                      PC1 = data$PC1, 
                      PC2 = data$PC2)

# use the posterior predictive distribution
ppd <- TRUE

pred_mod <- vector("list", length = nrow(df_pred))
for(i in 1:nrow(df_pred)) {

  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_pred, 
         (post$abar[, S[i]] + post$a[, , S[i]][, C[i]]) + 
           (post$b1[, S[i]] * Y[i]))
  if(ppd) {
    x <- rlnorm(n = length(x), x, post$sigma)
  } else {
    x <- exp(x + (0.5*(post$sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_pred, 
         post$a_hu[, S[i]] * Y[i] + post$b1_hu[, S[i]] * Y[i] )
  y <- 1 - plogis(y)
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  pred_mod[[i]] <- (x*y)
  
}

# bind into a matrix
pred_mod <- do.call("cbind", pred_mod)

# pull into a data.frame for plotting
C <- factor(data$cluster_id)
levels(C) <- paste0("Cluster ", LETTERS[1:10])
df_plot <- data.frame(M_obs = data$M,
                      Obs_pred = ifelse(is.na(data$M), "Predicted", "Observed"),
                      S = as.character(df_pred$S),
                      C = C,
                      T = as.character(data$time),
                      Y = data$Y,
                      PC1 = data$PC1,
                      PC2 = data$PC2,
                      M_pred_mu = apply(pred_mod, 2, mean),
                      M_pred_PIlow = apply(pred_mod, 2, HPDI, 0.90)[1,],
                      M_pred_PIhigh = apply(pred_mod, 2, HPDI, 0.90)[2,])


# check the min and max
summary(df_plot)

# factor OTU order: Barn Bryo Bumpi Hydro Seasq
p1 <- 
  ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  geom_point(data = df_plot,
         mapping = aes(x = M_obs, y = M_pred_mu, colour = S))
  geom_errorbar(mapping = aes(x = M_obs, 
                              ymin = M_pred_PIlow, ymax = M_pred_PIhigh,
                              colour = S),
                width = 0, alpha = 0.5, size = 0.25, show.legend = FALSE) +
  geom_point(shape = 16, alpha = 0.75) +
  ylab("Predicted monoculture (g)") +
  xlab("Observed monoculture (g)") +
  facet_wrap(~C, scales = "free") +
  scale_y_continuous(limits = c(0, 31)) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_colour_manual(name = "OTU",
                      labels = c("Barn", "Bryo", "Asci", "Hydro", "Ciona"),
                      values = viridis(n = 5, begin = 0.1, end = 0.9, option = "C")) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p1)

### END
