
v <- data[complete.cases(data),]
min_M <- min(v[v$M>0,]$M)

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


# full model: non-centred
m1 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, S] + b1[C, S]*T + b2[C, S]*Y + b3[C, S]*PC1 + b4[C, S]*PC2,
          
          # parameters
          transpars> matrix[C,5]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[5, C]:z_Ca ~ normal( 0 , 1 ),

          transpars> matrix[C,5]:b1 <-
            compose_noncentered( sigma_b1 , L_Rho_b1 , z_Cb1 ),
          matrix[5, C]:z_Cb1  ~ normal( 0 , 1 ),
          
          transpars> matrix[C,5]:b2 <-
            compose_noncentered( sigma_b2 , L_Rho_b2 , z_Cb2  ),
          matrix[5, C]:z_Cb2  ~ normal( 0 , 1 ),
          
          transpars> matrix[C,5]:b3 <-
            compose_noncentered( sigma_b3 , L_Rho_b3 , z_Cb3  ),
          matrix[5, C]:z_Cb3  ~ normal( 0 , 1 ),
          
          transpars> matrix[C,5]:b4 <-
            compose_noncentered( sigma_b4 , L_Rho_b4 , z_Cb4  ),
          matrix[5, C]:z_Cb4  ~ normal( 0 , 1 ),
          
          # fixed priors
          cholesky_factor_corr[5]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_a ~ dexp(1),
          
          cholesky_factor_corr[5]:L_Rho_b1 ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_b1 ~ dexp(1),
  
          cholesky_factor_corr[5]:L_Rho_b2 ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_b2 ~ dexp(1),
          
          cholesky_factor_corr[5]:L_Rho_b3 ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_b3 ~ dexp(1),
          
          cholesky_factor_corr[5]:L_Rho_b4 ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_b4 ~ dexp(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[5,5]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[5,5]:Rho_b1 <<- Chol_to_Corr(L_Rho_b1),
          gq> matrix[5,5]:Rho_b2 <<- Chol_to_Corr(L_Rho_b2),
          gq> matrix[5,5]:Rho_b3 <<- Chol_to_Corr(L_Rho_b3),
          gq> matrix[5,5]:Rho_b4 <<- Chol_to_Corr(L_Rho_b4)
          
    ),
    data = spp, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m1, depth = 3 )
# traceplot(mx)
#dev.off()

# plot the observed versus predicted data
post <- sim(m1) - min_M
plot(apply(post, 2, mean), spp$M)
abline(0, 1)

r <- apply(post, 2, mean) - spp$M
resid_var <- var2(r)
outcome_var <- var2( spp$M )
1 - resid_var/outcome_var

# full model: non-centred
m2 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, S] + b1[C, S]*T + b2[C, S]*Y,
          
          # parameters
          transpars> matrix[C,5]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[5, C]:z_Ca ~ normal( 0 , 1 ),
          
          transpars> matrix[C,5]:b1 <-
            compose_noncentered( sigma_b1 , L_Rho_b1 , z_Cb1 ),
          matrix[5, C]:z_Cb1  ~ normal( 0 , 1 ),
          
          transpars> matrix[C,5]:b2 <-
            compose_noncentered( sigma_b2 , L_Rho_b2 , z_Cb2  ),
          matrix[5, C]:z_Cb2  ~ normal( 0 , 1 ),
          
          # fixed priors
          cholesky_factor_corr[5]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_a ~ dexp(1),
          
          cholesky_factor_corr[5]:L_Rho_b1 ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_b1 ~ dexp(1),
          
          cholesky_factor_corr[5]:L_Rho_b2 ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_b2 ~ dexp(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[5,5]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[5,5]:Rho_b1 <<- Chol_to_Corr(L_Rho_b1),
          gq> matrix[5,5]:Rho_b2 <<- Chol_to_Corr(L_Rho_b2)
          
    ),
    data = spp, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m2, depth = 3 )
# traceplot(m2)
#dev.off()

# plot the observed versus predicted data
post <- sim(m2) - min_M
plot(apply(post, 2, mean), spp$M)
abline(0, 1)

r <- apply(post, 2, mean) - spp$M
resid_var <- var2(r)
outcome_var <- var2( spp$M )
1 - resid_var/outcome_var


# full model: non-centred
m3 <- 
  ulam(
    alist(M ~ dlnorm(u, sigma),
          
          u <- a[C, S] + b1[C, S]*T,
          
          # parameters
          transpars> matrix[C,5]:a <-
            compose_noncentered( sigma_a , L_Rho_a , z_Ca ),
          matrix[5, C]:z_Ca ~ normal( 0 , 1 ),
          
          transpars> matrix[C,5]:b1 <-
            compose_noncentered( sigma_b1 , L_Rho_b1 , z_Cb1 ),
          matrix[5, C]:z_Cb1  ~ normal( 0 , 1 ),
          
          # fixed priors
          cholesky_factor_corr[5]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_a ~ dexp(1),
          
          cholesky_factor_corr[5]:L_Rho_b1 ~ lkj_corr_cholesky( 2 ),
          vector[5]:sigma_b1 ~ dexp(1),
          
          sigma ~ exponential(1),
          
          # convert cholesky to corr matrix
          gq> matrix[5,5]:Rho_a <<- Chol_to_Corr(L_Rho_a),
          gq> matrix[5,5]:Rho_b1 <<- Chol_to_Corr(L_Rho_b1)
          
    ),
    data = spp, chains = 4, log_lik = TRUE, control = list(adapt_delta = 0.99))

# check the precis output
precis( m3, depth = 3 )
# traceplot(m2)
#dev.off()

# plot the observed versus predicted data
post <- sim(m3) - min_M
plot(apply(post, 2, mean), spp$M)
abline(0, 1)

r <- apply(post, 2, mean) - spp$M
resid_var <- var2(r)
outcome_var <- var2( spp$M )
1 - resid_var/outcome_var

# compare these three models for their predictive abilities
compare(m1, m2, m3, func = "PSIS")







# check the best predictive model
post <- rethinking::extract.samples(mx)

# make a plot of the observed and predicted values
source(here("scripts/Function_plotting_theme.R"))
pred <- sim(mx)
pred <- pred + min_M
pred <- link(mx)
pred <- 
  apply(pred, 2, function(x) {
    
    exp(x + (0.5*(post$sigma^2)))
    
  } )
df.pred <- data.frame(M_obs = spp$M,
                      M_S = spp$S,
                      M_C = as.character(spp$C),
                      M_T = as.character(spp$T),
                      M_pred_mu = apply(pred, 2, mean),
                      M_pred_PIlow = apply(pred, 2, HPDI, 0.75)[1,],
                      M_pred_PIhigh = apply(pred, 2, HPDI, 0.75)[2,])

p1 <- 
  ggplot(data = df.pred,
         mapping = aes(x = M_obs, y = M_pred_mu, colour = M_C, fill = M_C)) +
  geom_point(shape = 16, alpha = 0.5) +
  geom_errorbar(mapping = aes(x = M_obs, ymin = M_pred_PIlow, ymax = M_pred_PIhigh),
                width = 0, alpha = 0.5, size = 0.25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  ylab("Predicted monoculture (g)") +
  xlab("Observed monoculture (g)") +
  facet_wrap(~M_S, scales = "free") +
  # scale_y_continuous(limits = c(0, 17)) +
  # scale_x_continuous(limits = c(0, 17)) +
  theme_meta()
plot(p1)







