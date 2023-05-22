#'
#' @title: Barn: Model the missing monocultures
#' 
#' @description: This script attempts to model the missing monocultures
#' for the Barn species after preliminary data analysis showed that a single
#' monoculture would not be suitable for all species.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(here)

# load the rethinking package
library(rethinking)

# load the analysis data
data <- read_csv(here("data/case_study_2/data_clean/biomass_env_analysis_data.csv"))

# create a v data.frame with the relevant data
v <- data[data$OTU %in% c("Barn"), ]

# remove NAs from the monoculture column
v <- v[complete.cases(v),]

# make a data.list with the training data
barn <- 
  list(M = v$M,
       Y = with(v, (Y - min(Y))/(max(Y)-min(Y)) ),
       PC1 = with(v, (PC1 - min(PC1))/(max(PC1)-min(PC1)) ),
       PC2 = with(v, (PC2 - min(PC2))/(max(PC2)-min(PC2)) ),
       C = as.integer(as.factor(v$cluster_id)),
       T = as.integer(as.factor(v$time))
  )

# full model: non-centred
m1 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, T] + b1[C]*Y + b2[C]*PC1 + b3[C]*PC2 + b4[C]*Y*PC1 + b5[C]*Y*PC2,
          
          # a parameters
          transpars> matrix[C,3]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[3, C]:z_Ca ~ normal( 0 , 1 ),
          
          # b1-3 parameters
          transpars> vector[10]:b1 <<- bbar[1] + z_Cb[,1],
          transpars> vector[10]:b2 <<- bbar[2] + z_Cb[,2],
          transpars> vector[10]:b3 <<- bbar[3] + z_Cb[,3],
          transpars> vector[10]:b4 <<- bbar[4] + z_Cb[,4],
          transpars> vector[10]:b5 <<- bbar[5] + z_Cb[,5],
          transpars> matrix[10, 5]:z_Cb <- 
            compose_noncentered(sigma_b, L_Rho_b, Z),
          matrix[5, 10]:Z ~ normal(0, 1),
          
          # fixed priors
          cholesky_factor_corr[3]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[3]:sigma_a ~ dexp(1),
          
          cholesky_factor_corr[5]:L_Rho_b ~ lkj_corr_cholesky( 2 ),
          vector[5]:bbar ~ normal(0, 1),
          vector[5]:sigma_b ~ exponential(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[3,3]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[5,5]:Rho_b <<- Chol_to_Corr(L_Rho_b)
          
    ),
    data = barn, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m1, depth = 3 )
traceplot(m1)
dev.off()

# plot the observed versus predicted data
post <- sim(m1)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

r <- apply(post, 2, mean) - barn$M
resid_var <- var2(r)
outcome_var <- var2( barn$M )
1 - resid_var/outcome_var

# model 2
m2 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, T] + b1[C]*Y + b2[C]*PC1 + b3[C]*PC2 + b4[C]*Y*PC1,
          
          # a parameters
          transpars> matrix[C,3]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[3, C]:z_Ca ~ normal( 0 , 1 ),
          
          # b1-3 parameters
          transpars> vector[10]:b1 <<- bbar[1] + z_Cb[,1],
          transpars> vector[10]:b2 <<- bbar[2] + z_Cb[,2],
          transpars> vector[10]:b3 <<- bbar[3] + z_Cb[,3],
          transpars> vector[10]:b4 <<- bbar[4] + z_Cb[,4],
          transpars> matrix[10, 4]:z_Cb <- 
            compose_noncentered(sigma_b, L_Rho_b, Z),
          matrix[4, 10]:Z ~ normal(0, 1),
          
          # fixed priors
          cholesky_factor_corr[3]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[3]:sigma_a ~ dexp(1),
          
          cholesky_factor_corr[4]:L_Rho_b ~ lkj_corr_cholesky( 2 ),
          vector[4]:bbar ~ normal(0, 1),
          vector[4]:sigma_b ~ exponential(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[3,3]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[4,4]:Rho_b <<- Chol_to_Corr(L_Rho_b)
          
    ),
    data = barn, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m2, depth = 3 )
traceplot(m2)
dev.off()

# plot the observed versus predicted data
post <- sim(m2)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

r <- apply(post, 2, mean) - barn$M
resid_var <- var2(r)
outcome_var <- var2( barn$M )
1 - resid_var/outcome_var

# model 3
m3 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, T] + b1[C]*Y + b2[C]*PC1 + b3[C]*PC2,
          
          # a parameters
          transpars> matrix[C,3]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[3, C]:z_Ca ~ normal( 0 , 1 ),
          
          # b1-3 parameters
          transpars> vector[10]:b1 <<- bbar[1] + z_Cb[,1],
          transpars> vector[10]:b2 <<- bbar[2] + z_Cb[,2],
          transpars> vector[10]:b3 <<- bbar[3] + z_Cb[,3],
          transpars> matrix[10, 3]:z_Cb <- 
            compose_noncentered(sigma_b, L_Rho_b, Z),
          matrix[3, 10]:Z ~ normal(0, 1),
          
          # fixed priors
          cholesky_factor_corr[3]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[3]:sigma_a ~ dexp(1),
          
          cholesky_factor_corr[3]:L_Rho_b ~ lkj_corr_cholesky( 2 ),
          vector[3]:bbar ~ normal(0, 1),
          vector[3]:sigma_b ~ exponential(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[3,3]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[3,3]:Rho_b <<- Chol_to_Corr(L_Rho_b)
          
    ),
    data = barn, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m3, depth = 3 )
traceplot(m3)
dev.off()

# plot the observed versus predicted data
post <- sim(m3)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

r <- apply(post, 2, mean) - barn$M
resid_var <- var2(r)
outcome_var <- var2( barn$M )
1 - resid_var/outcome_var

# model 4
m4 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, T] + b1[C]*Y + b2[C]*PC1,
          
          # a parameters
          transpars> matrix[C,3]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[3, C]:z_Ca ~ normal( 0 , 1 ),
          
          # b1-3 parameters
          transpars> vector[10]:b1 <<- bbar[1] + z_Cb[,1],
          transpars> vector[10]:b2 <<- bbar[2] + z_Cb[,2],
          transpars> matrix[10, 2]:z_Cb <- 
            compose_noncentered(sigma_b, L_Rho_b, Z),
          matrix[2, 10]:Z ~ normal(0, 1),
          
          # fixed priors
          cholesky_factor_corr[3]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[3]:sigma_a ~ dexp(1),
          
          cholesky_factor_corr[2]:L_Rho_b ~ lkj_corr_cholesky( 2 ),
          vector[2]:bbar ~ normal(0, 1),
          vector[2]:sigma_b ~ exponential(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[3,3]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[2,2]:Rho_b <<- Chol_to_Corr(L_Rho_b)
          
    ),
    data = barn, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m4, depth = 3 )
traceplot(m4)
dev.off()

# plot the observed versus predicted data
post <- sim(m4)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

apply(post, 2, PI)

r <- apply(post, 2, mean) - barn$M
resid_var <- var2(r)
outcome_var <- var2( barn$M )
1 - resid_var/outcome_var

# model 5
m5 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, T] + z[C]*sigma_b*Y,
          
          # a parameters
          transpars> matrix[C,3]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[3, C]:z_Ca ~ normal( 0 , 1 ),
          
          # fixed priors
          cholesky_factor_corr[3]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[3]:sigma_a ~ dexp(1),
          
          z[C] ~ dnorm( 0 , 1 ),
          sigma_b ~ exponential(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[3,3]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> vector[C]:b <<- z*sigma_b
          
    ),
    data = barn, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m5, depth = 3 )
traceplot(m5)
dev.off()

# plot the observed versus predicted data
post <- sim(m5)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

r <- apply(post, 2, mean) - barn$M
resid_var <- var2(r)
outcome_var <- var2( barn$M )
1 - resid_var/outcome_var

# model 6
m6 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, T],
          
          # a parameters
          transpars> matrix[C,3]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[3, C]:z_Ca ~ normal( 0 , 1 ),
          
          # fixed priors
          cholesky_factor_corr[3]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[3]:sigma_a ~ dexp(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[3,3]:Rho_a <<- Chol_to_Corr(L_Rho_a)
          
    ),
    data = barn, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m6, depth = 3 )
traceplot(m6)
dev.off()

# plot the observed versus predicted data
post <- sim(m6)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

r <- apply(post, 2, mean) - barn$M
resid_var <- var2(r)
outcome_var <- var2( barn$M )
1 - resid_var/outcome_var

# compare using PSIS
compare(m1, m2, m3, m4, m5, m6, func = "PSIS")

# sample from the m.A3 model
post <- rethinking::extract.samples(m4)

# write this posterior distribution list as a .rds object
saveRDS(post, here("results/SP_A_monoculture_posterior.rds"))

# save the model object as a .rds object
saveRDS(m4, here("results/SP_A_model_object.rds")) 

# make a plot of the observed and predicted values
source(here("scripts/Function_plotting_theme.R"))
pred <- sim(m4)
pred <- link(m4)
pred <- 
  apply(pred, 2, function(x) {
  
  exp(x + (0.5*(post$sigma^2)))
  
} )
df.pred <- data.frame(M_obs = barn$M,
                      M_C = as.character(barn$C),
                      M_T = as.character(barn$T),
                      M_pred_mu = apply(pred, 2, mean),
                      M_pred_PIlow = apply(pred, 2, HPDI, 0.90)[1,],
                      M_pred_PIhigh = apply(pred, 2, HPDI, 0.90)[2,])

p1 <- 
  ggplot(data = df.pred,
       mapping = aes(x = M_obs, y = M_pred_mu, colour = M_T, fill = M_T)) +
  geom_point(shape = 16, alpha = 0.5) +
  geom_errorbar(mapping = aes(x = M_obs, ymin = M_pred_PIlow, ymax = M_pred_PIhigh),
                width = 0, alpha = 0.5, size = 0.25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  ylab("Predicted monoculture (g)") +
  xlab("Observed monoculture (g)") +
  facet_wrap(~M_C) +
  scale_y_continuous(limits = c(0, 17)) +
  scale_x_continuous(limits = c(0, 17)) +
  theme_meta()
plot(p1)



# save as an object into the results folder
saveRDS(p1, file = here("results/SP_A_monoculture_plot.rds"))

# ggsave(here("figures/fig_barn_mono.png"), p1,
       # width = 8, height = 7, units = "cm")

### END
