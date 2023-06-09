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


# Lognormal hurdle models

# model 1

# compile the model
ln1 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model1.stan",
                        verbose = TRUE)

# sample the stan model
ln1_fit <- rstan::sampling(ln1, data = df_m_obs, 
                          iter = 3000, warmup = 1000, 
                          chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99),
                          cores = 4,
                          seed = 58455)

# get the number of unconstrained parameters: 40
get_num_upars(ln1_fit)

# save the stan model fit object
ln1_fit@stanmodel@dso <- new("cxxdso")
saveRDS(ln1_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model1_fit.rds")

# model 2

# compile the model
ln2 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model2.stan",
                        verbose = TRUE)

# sample the stan model
ln2_fit <- rstan::sampling(ln2, data = df_m_obs, 
                           iter = 3000, warmup = 1000, 
                           chains = 4, algorithm = c("NUTS"),
                           control = list(adapt_delta = 0.99),
                           cores = 4,
                           seed = 53275)

# get the number of unconstrained parameters: 27
get_num_upars(ln2_fit)

# save the stan model fit object
ln2_fit@stanmodel@dso <- new("cxxdso")
saveRDS(ln2_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model2_fit.rds")

# model 3

# compile the model
ln3 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model3.stan",
                        verbose = TRUE)

# sample the stan model
ln3_fit <- rstan::sampling(ln3, data = df_m_obs, 
                           iter = 3000, warmup = 1000, 
                           chains = 4, algorithm = c("NUTS"),
                           control = list(adapt_delta = 0.99),
                           cores = 4,
                           seed = 40975
                           )

# get the number of unconstrained parameters: 18
get_num_upars(ln3_fit)

# save the stan model fit object
ln3_fit@stanmodel@dso <- new("cxxdso")
saveRDS(ln3_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model3_fit.rds")

# model 4

# compile the model
ln4 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model4.stan",
                        verbose = TRUE)

# sample the stan model
ln4_fit <- rstan::sampling(ln4, data = df_m_obs, 
                           iter = 3000, warmup = 1000,
                           chains = 4, algorithm = c("NUTS"),
                           control = list(adapt_delta = 0.99),
                           cores = 4,
                           seed = 50947
                           )

# get the number of unconstrained parameters: 17
get_num_upars(ln4_fit)

# save the stan model fit object
ln4_fit@stanmodel@dso <- new("cxxdso")
saveRDS(ln4_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model4_fit.rds")

# model 5

# compile the model
ln5 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model5.stan",
                        verbose = TRUE)

# sample the stan model
ln5_fit <- rstan::sampling(ln5, data = df_m_obs, 
                           iter = 3000, warmup = 1000, 
                           chains = 4, algorithm = c("NUTS"),
                           control = list(adapt_delta = 0.99),
                           cores = 4,
                           seed = 9840
                           )

# get the number of unconstrained parameters: 14
get_num_upars(ln5_fit)

# save the stan model fit object
ln5_fit@stanmodel@dso <- new("cxxdso")
saveRDS(ln5_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model5_fit.rds")

# model 6

# compile the model
ln6 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model6.stan",
                        verbose = TRUE)

# sample the stan model
ln6_fit <- rstan::sampling(ln6, data = df_m_obs, 
                           iter = 3000, warmup = 1000, 
                           chains = 4, algorithm = c("NUTS"),
                           control = list(adapt_delta = 0.99),
                           cores = 4,
                           seed =  93125
                           )

# get the number of unconstrained parameters: 12
get_num_upars(ln6_fit)

# save the stan model fit object
ln6_fit@stanmodel@dso <- new("cxxdso")
saveRDS(ln6_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model6_fit.rds")

# model 7

# compile the model
ln7 <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model7.stan",
                        verbose = TRUE)

# sample the stan model
ln7_fit <- rstan::sampling(ln7, data = df_m_obs, 
                           iter = 3000, warmup = 1000, 
                           chains = 4, algorithm = c("NUTS"),
                           control = list(adapt_delta = 0.99),
                           cores = 4,
                           seed = 45982
                           )

# get the number of unconstrained parameters: 11
get_num_upars(ln7_fit)

# save the stan model fit object
ln7_fit@stanmodel@dso <- new("cxxdso")
saveRDS(ln7_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model7_fit.rds")

# null model (intercept only)

# compile the model
ln_null <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model_null.stan",
                         verbose = TRUE)

# sample the stan model
ln_null_fit <- rstan::sampling(ln_null, data = df_m_obs, 
                           iter = 3000, warmup = 1000, 
                           chains = 4, algorithm = c("NUTS"),
                           control = list(adapt_delta = 0.99),
                           cores = 4,
                           seed = 45982
)

# get the number of unconstrained parameters: 3
get_num_upars(ln_null_fit)

# save the stan model fit object
ln_null_fit@stanmodel@dso <- new("cxxdso")
saveRDS(ln_null_fit, file = "scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_lognormal_model_null_fit.rds")

### END
