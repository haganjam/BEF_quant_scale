#'
#' @title: Model the missing monocultures
#' 
#' @description: Fit log-normal hurdle generalised linear models in Stan that we 
#' will use to impute the missing monoculture data .
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
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
       Y = v$Y,
       PC1 = v$PC1,
       PC2 = v$PC2,
       C = v$cluster_id,
       T = v$time,
       S = v$OTU)

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

### END
