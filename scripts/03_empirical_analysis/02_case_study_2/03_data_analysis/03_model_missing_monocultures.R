#'
#' @title: Model the missing monocultures
#' 
#' @description: This script attempts to impute the missing monoculture
#' data using Bayesian generalised linera models.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'
#'
#' do the standardisation before we predict the new data because this
#' will influence the later predictions
#'
#' for reference, this is the simplest centered model
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

# load the required libraries
library(readr)
library(dplyr)
library(ggplot2)

# load the rethinking package
library(rethinking)

# load plotting theme
source("scripts/Function_plotting_theme.R")

# load the analysis data
data <- read_csv("data/case_study_2/data_clean/biomass_env_analysis_data.csv")

# remove any NAs
v <- data[complete.cases(data),]

# calculate the minimum observed biomass
min_M <- min(v[v$M>0,]$M)

# add this minimum to we can model without zero values
v$M <- v$M + min_M

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


# model 1: full model 
m1 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- 
            (abar[S] + a[C, S]) + 
            (b1bar[S] + b1[C, S]*T) + 
            (b2bar[S] + b2[C, S]*Y) + 
            (b3bar[S] + b3[C, S]*PC1) + 
            (b4bar[S] + b4[C, S]*PC2) +
            (b5bar[S] + b5[C, S]*Y*PC1),
          
          # define the effects using other parameters
          matrix[5, 10]:Z ~ normal( 0 , 1 ),
          
          transpars> matrix[10,5]:a <-
            compose_noncentered( sigma_a , L_Rho_a , Z ),
          
          transpars> matrix[10,5]:b1 <-
            compose_noncentered( sigma_b1 , L_Rho_b1 , Z ),
          
          transpars> matrix[10,5]:b2 <-
            compose_noncentered( sigma_b2 , L_Rho_b2 , Z  ),
          
          transpars> matrix[10,5]:b3 <-
            compose_noncentered( sigma_b3 , L_Rho_b3 , Z  ),
          
          transpars> matrix[10,5]:b4 <-
            compose_noncentered( sigma_b4 , L_Rho_b4 , Z  ),
          
          transpars> matrix[10,5]:b5 <-
            compose_noncentered( sigma_b5 , L_Rho_b5 , Z  ),
          
          # fixed priors
          cholesky_factor_corr[5]:c(L_Rho_a, L_Rho_b1, L_Rho_b2, L_Rho_b3, L_Rho_b4, L_Rho_b5) ~ lkj_corr_cholesky( 2 ),
          
          vector[5]:c(sigma_a, sigma_b1, sigma_b2, sigma_b3, sigma_b4, sigma_b5) ~ dexp(1),
          vector[5]:c(abar, b1bar, b2bar, b3bar, b4bar, b5bar) ~ dnorm(0, 2),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[5,5]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[5,5]:Rho_b1 <<- Chol_to_Corr(L_Rho_b1),
          gq> matrix[5,5]:Rho_b2 <<- Chol_to_Corr(L_Rho_b2),
          gq> matrix[5,5]:Rho_b3 <<- Chol_to_Corr(L_Rho_b3),
          gq> matrix[5,5]:Rho_b4 <<- Chol_to_Corr(L_Rho_b4),
          gq> matrix[5,5]:Rho_b5 <<- Chol_to_Corr(L_Rho_b5)
          
    ),
    data = spp, chains = 4, log_lik = TRUE)

# check the precis output
precis( m1, depth = 3 )
# traceplot(mx)
# dev.off()
stancode(m1)

# compile the stan growth rate model
mx <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model1.stan",
                        verbose = TRUE)
print(m1)

# sample the stan model: m1
m1_fit <- rstan::sampling(m1, data = dat, 
                          iter = 2500, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99),
                          seed = 54856)

# plot the observed versus predicted data
post <- sim(m1) - min_M
plot(apply(post, 2, mean), spp$M)
abline(0, 1)

r <- apply(post, 2, mean) - spp$M
resid_var <- var2(r)
outcome_var <- var2( spp$M )
1 - resid_var/outcome_var

# model 2
m2 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- 
            (abar[S] + a[C, S]) + 
            (b1bar[S] + b1[C, S]*T) + 
            (b2bar[S] + b2[C, S]*Y) + 
            (b3bar[S] + b3[C, S]*PC1) + 
            (b4bar[S] + b4[C, S]*PC2),
          
          # define the effects using other parameters
          matrix[5, 10]:Z ~ normal( 0 , 1 ),
          
          transpars> matrix[10,5]:a <-
            compose_noncentered( sigma_a , L_Rho_a , Z ),
          
          transpars> matrix[10,5]:b1 <-
            compose_noncentered( sigma_b1 , L_Rho_b1 , Z ),
          
          transpars> matrix[10,5]:b2 <-
            compose_noncentered( sigma_b2 , L_Rho_b2 , Z  ),
          
          transpars> matrix[10,5]:b3 <-
            compose_noncentered( sigma_b3 , L_Rho_b3 , Z  ),
          
          transpars> matrix[10,5]:b4 <-
            compose_noncentered( sigma_b4 , L_Rho_b4 , Z  ),
          
          # fixed priors
          cholesky_factor_corr[5]:c(L_Rho_a, L_Rho_b1, L_Rho_b2, L_Rho_b3, L_Rho_b4) ~ lkj_corr_cholesky( 2 ),
          
          vector[5]:c(sigma_a, sigma_b1, sigma_b2, sigma_b3, sigma_b4) ~ dexp(1),
          vector[5]:c(abar, b1bar, b2bar, b3bar, b4bar) ~ dnorm(0, 2),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[5,5]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[5,5]:Rho_b1 <<- Chol_to_Corr(L_Rho_b1),
          gq> matrix[5,5]:Rho_b2 <<- Chol_to_Corr(L_Rho_b2),
          gq> matrix[5,5]:Rho_b3 <<- Chol_to_Corr(L_Rho_b3),
          gq> matrix[5,5]:Rho_b4 <<- Chol_to_Corr(L_Rho_b4)
          
    ),
    data = spp, chains = 4, log_lik = TRUE)

# , control = list(adapt_delta = 0.95)

# check the precis output
precis(m2, depth = 3)
# traceplot(m2)
# dev.off()

# plot the observed versus predicted data
post <- sim(m2) - min_M
plot(apply(post, 2, mean), spp$M)
abline(0, 1)

r <- apply(post, 2, mean) - spp$M
resid_var <- var2(r)
outcome_var <- var2( spp$M )
1 - resid_var/outcome_var

# model 3
m3 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- 
            (abar[S] + a[C, S]) + 
            (b1bar[S] + b1[C, S]*T) + 
            (b2bar[S] + b2[C, S]*Y),
          
          # define the effects using other parameters
          matrix[5, 10]:Z ~ normal( 0 , 1 ),
          
          transpars> matrix[10,5]:a <-
            compose_noncentered( sigma_a , L_Rho_a , Z ),
          
          transpars> matrix[10,5]:b1 <-
            compose_noncentered( sigma_b1 , L_Rho_b1 , Z ),
          
          transpars> matrix[10,5]:b2 <-
            compose_noncentered( sigma_b2 , L_Rho_b2 , Z  ),
          
          # fixed priors
          cholesky_factor_corr[5]:c(L_Rho_a, L_Rho_b1, L_Rho_b2) ~ lkj_corr_cholesky( 2 ),
          
          vector[5]:c(sigma_a, sigma_b1, sigma_b2) ~ dexp(1),
          vector[5]:c(abar, b1bar, b2bar) ~ dnorm(0, 2),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[5,5]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[5,5]:Rho_b1 <<- Chol_to_Corr(L_Rho_b1),
          gq> matrix[5,5]:Rho_b2 <<- Chol_to_Corr(L_Rho_b2)
          
    ),
    data = spp, chains = 4, log_lik = TRUE)

# check the precis output
precis( m3, depth = 3 )
# traceplot(m3)
# dev.off()

# plot the observed versus predicted data
post <- sim(m3) - min_M
plot(apply(post, 2, mean), spp$M)
abline(0, 1)

r <- apply(post, 2, mean) - spp$M
resid_var <- var2(r)
outcome_var <- var2( spp$M )
1 - resid_var/outcome_var

# model 4
m4 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, S] + b1[C, S]*T,
          
          u <- 
            (abar[S] + a[C, S]) + 
            (b1bar[S] + b1[C, S]*T),
          
          # define the effects using other parameters
          matrix[5, 10]:Z ~ normal( 0 , 1 ),
          
          transpars> matrix[10,5]:a <-
            compose_noncentered( sigma_a , L_Rho_a , Z ),
          
          transpars> matrix[10,5]:b1 <-
            compose_noncentered( sigma_b1 , L_Rho_b1 , Z ),
          
          # fixed priors
          cholesky_factor_corr[5]:c(L_Rho_a, L_Rho_b1) ~ lkj_corr_cholesky( 2 ),
          
          vector[5]:c(sigma_a, sigma_b1) ~ dexp(1),
          vector[5]:c(abar, b1bar) ~ dnorm(0, 2),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[5,5]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[5,5]:Rho_b1 <<- Chol_to_Corr(L_Rho_b1)
    ),
    data = spp, chains = 4, log_lik = TRUE)

# check the precis output
precis( m4, depth = 3 )
# traceplot(m4)
# dev.off()

# plot the observed versus predicted data
post <- sim(m4) - min_M
plot(apply(post, 2, mean), spp$M)
abline(0, 1)

r <- apply(post, 2, mean) - spp$M
resid_var <- var2(r)
outcome_var <- var2( spp$M )
1 - resid_var/outcome_var

# model 5
m5 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- abar[S] + a[C, S],
          
          # parameters
          transpars> matrix[10, 5]:a <-
            compose_noncentered( sigma_par , L_Rho_a , Z),
          matrix[5, 10]:Z ~ normal( 0 , 1 ),
          
          # fixed priors
          vector[5]:abar ~ dnorm(0, 2),
          vector[5]:sigma_par ~ dexp(1),
          cholesky_factor_corr[5]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[5,5]:Rho_a <<- Chol_to_Corr(L_Rho_a)
          
    ),
    data = spp, chains = 4, log_lik = TRUE)

# check the precis output
precis( m5, depth = 3 )
# traceplot(m5)
# dev.off()

# plot the observed versus predicted data
post <- sim(m5) - min_M
plot(apply(post, 2, mean), spp$M)
abline(0, 1)

r <- apply(post, 2, mean) - spp$M
resid_var <- var2(r)
outcome_var <- var2( spp$M )
1 - resid_var/outcome_var

# compare these three models for their predictive abilities: m3 best
compare(m1, m2, m3, m4, m5, func = "PSIS")

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
