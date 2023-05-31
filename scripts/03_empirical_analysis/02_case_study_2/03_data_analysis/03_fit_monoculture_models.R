#'
#' @title: Model the missing monocultures
#' 
#' @description: Fit log-normal hurdle generalised linear models in Stan that we 
#' will use to impute the missing monoculture data .
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'
#' reference on divergent transitions: https://www.martinmodrak.cz/2018/02/19/taming-divergences-in-stan-models/
#'
#' for a conceptual reference, this is the simplest centered model coded in ulam()
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

# model 1

# compile the model
m1 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model1.stan",
                        verbose = TRUE)

# sample the stan model
m1_fit <- rstan::sampling(m1, data = df_m_obs, 
                          iter = 3000, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99),
                          cores = 4)

# save the stan model fit object
m1_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m1_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model1_fit.rds")

# model 2

# compile the model
m2 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model2.stan",
                        verbose = TRUE)

# sample the stan model
m2_fit <- rstan::sampling(m2, data = df_m_obs, 
                          iter = 3000, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99),
                          cores = 4)

# save the stan model fit object
m2_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m2_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model2_fit.rds")

# model 3

# compile the model
m3 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model3.stan",
                        verbose = TRUE)

# sample the stan model
m3_fit <- rstan::sampling(m3, data = df_m_obs, 
                          iter = 3000, chains = 4, algorithm = c("NUTS"), 
                          control = list(adapt_delta = 0.99),
                          cores = 4)

# save the stan model fit object
m3_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m3_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model3_fit.rds")

# model 4

# compile the model
m4 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model4.stan",
                        verbose = TRUE)

# sample the stan model
m4_fit <- rstan::sampling(m4, data = df_m_obs, 
                          iter = 3000, chains = 4, algorithm = c("NUTS"),
                          cores = 4)

# check the model
pairs(m4_fit, pars = c("abar"))

# save the stan model fit object
m4_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m4_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model4_fit.rds")

# model 5

# compile the model
m5 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model5.stan",
                        verbose = TRUE)

# sample the stan model
m5_fit <- rstan::sampling(m5, data = df_m_obs, 
                          iter = 2000, chains = 4, algorithm = c("NUTS"),
                          seed = 54856, cores = 4)

# save the stan model fit object
m5_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m5_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model5_fit.rds")

# model 6

# compile the model
m6 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model6.stan",
                        verbose = TRUE)

# sample the stan model
m6_fit <- rstan::sampling(m6, data = df_m_obs, 
                          iter = 3000, chains = 4, algorithm = c("NUTS"),
                          cores = 4)

# check for divergent transitions
pairs(m6_fit, pars = "a")

pairs(m6_fit, pars = "a")
pairs(m6_fit, pars = c("sigma", "sigma_a", "b1", "a_hu", "b1_hu"))

# save the stan model fit object
m6_fit@stanmodel@dso <- new("cxxdso")
saveRDS(m5_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model6_fit.rds")

### END
