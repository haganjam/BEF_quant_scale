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
data <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/biomass_env_analysis_data.csv"))

# create a v data.frame with the relevant data
v <- data[data$OTU %in% c("Barn"), ]

# remove NAs from the monoculture column
v <- v[complete.cases(v),]

# make a data.list with the training data
barn <- 
  list(M = v$M,
       Y = v$Y,
       Y2 = v$Y^2,
       PC1 = v$PC1,
       PC2 = v$PC2
  )

# model1
m.A1 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bY2*Y2 + bPC1*PC1 + bPC2*PC2 + bPCI*PC1*PC2,
          
          bY ~ dnorm( 0 , 1 ),
          bY2 ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          bPC2 ~ dnorm( 0, 1),
          bPCI ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = barn, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.A1  )

# check the traceplots
traceplot(m.A1 )
dev.off()

# plot the Y vs M comparison
range(barn$Y)
Y_seq <- seq(-1, 3.2, 0.1)
PC1_seq <- mean(barn$PC1)
PC2_seq <- mean(barn$PC2)
pred_dat <- list(Y = Y_seq, Y2 = Y_seq^2, 
                 PC1 = rep(PC1_seq, length(Y_seq)),
                 PC2 = rep(PC2_seq, length(Y_seq)))
post <- link(m.A1, pred_dat)

plot(barn$Y, barn$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.A1)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

# model2
m.A2 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC2*PC2 + bPCI*PC1*PC2,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          bPC2 ~ dnorm( 0, 1),
          bPCI ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = barn, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.A2  )

# check the traceplots
traceplot(m.A2 )
dev.off()

# plot the Y vs M comparison
range(barn$Y)
Y_seq <- seq(-1, 3.2, 0.1)
PC1_seq <- mean(barn$PC1)
PC2_seq <- mean(barn$PC2)
pred_dat <- list(Y = Y_seq, 
                 PC1 = rep(PC1_seq, length(Y_seq)),
                 PC2 = rep(PC2_seq, length(Y_seq)))
post <- link(m.A2, pred_dat)

plot(barn$Y, barn$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.A2)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)


# model3
m.A3 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC2*PC2,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          bPC2 ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = barn, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.A3  )

# check the traceplots
traceplot(m.A3 )
dev.off()

# plot the Y vs M comparison
range(barn$Y)
Y_seq <- seq(-1, 3.2, 0.1)
PC1_seq <- mean(barn$PC1)
PC2_seq <- mean(barn$PC2)
pred_dat <- list(Y = Y_seq, 
                 PC1 = rep(PC1_seq, length(Y_seq)),
                 PC2 = rep(PC2_seq, length(Y_seq)))
post <- link(m.A3, pred_dat)

plot(barn$Y, barn$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.A3)
plot(apply(post, 2, mean), barn$M, ylab = "Predicted", xlab = "Observed")
abline(0, 1)

# model4
m.A4 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = barn, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.A4  )

# check the traceplots
traceplot(m.A4 )
dev.off()

# plot the Y vs M comparison
range(barn$Y)
Y_seq <- seq(-1, 3.2, 0.1)
PC1_seq <- mean(barn$PC1)
pred_dat <- list(Y = Y_seq, 
                 PC1 = rep(PC1_seq, length(Y_seq)))
post <- link(m.A4, pred_dat)

plot(barn$Y, barn$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.A4)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

# model5
m.A5 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bY2*Y2,
          
          bY ~ dnorm( 0 , 1 ),
          bY2 ~ dnorm( 0, 1 ),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = barn, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.A5  )

# check the traceplots
traceplot(m.A5 )
dev.off()

# plot the Y vs M comparison
range(barn$Y)
Y_seq <- seq(-1, 3.2, 0.1)
pred_dat <- list(Y = Y_seq, 
                 Y2 = Y_seq^2)
post <- link(m.A5, pred_dat)

plot(barn$Y, barn$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.A5)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

# model6
m.A6 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPCY*PC1*Y,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          bPCY ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = barn, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.A6  )

# check the traceplots
traceplot(m.A6 )
dev.off()

# plot the Y vs M comparison
range(barn$Y)
Y_seq <- seq(-1, 3.2, 0.1)
PC1_seq <- mean(barn$PC1)
pred_dat <- list(Y = Y_seq,
                 PC1 = rep(PC1_seq, length(Y_seq)))
post <- link(m.A6, pred_dat)

plot(barn$Y, barn$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.A6)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

# model7
m.A7 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y,
          
          bY ~ dnorm( 0 , 1 ),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = barn, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.A7  )

# check the traceplots
traceplot(m.A7 )
dev.off()

# plot the Y vs M comparison
range(barn$Y)
Y_seq <- seq(-1, 3.2, 0.1)
pred_dat <- list(Y = Y_seq)
post <- link(m.A7, pred_dat)

plot(barn$Y, barn$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.A7)
plot(apply(post, 2, mean), barn$M)
abline(0, 1)

# model8
m.A8 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS,
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = barn, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.A8  )

# check the traceplots
traceplot(m.A8 )
dev.off()

# simulate the posterior and check the range
range(sim(m.A8))

# compare the different models
compare(m.A1, m.A2, m.A3, m.A4,
        m.A5, m.A6, m.A7, m.A8,
        func = PSIS)

compare(m.A1, m.A2, m.A3, m.A4,
        m.A5, m.A6, m.A7, m.A8,
        func = WAIC)

# sample from the m.A3 model
post <- rethinking::extract.samples(m.A3)

# write this posterior distribution list as a .rds object
saveRDS(post, here("results/SP_A_monoculture_posterior.rds"))

# save the model object as a .rds object
saveRDS(m.A3, here("results/SP_A_model_object.rds")) 

# make a plot of the observed and predicted values
source(here("scripts/Function_plotting_theme.R"))
pred <- sim(m.A3)
df.pred <- data.frame(M_obs = barn$M,
                      M_pred_mu = apply(pred, 2, mean),
                      M_pred_PIlow = apply(pred, 2, PI)[1,],
                      M_pred_PIhigh = apply(pred, 2, PI)[2,])

p1 <- 
  ggplot(data = df.pred,
       mapping = aes(x = M_obs, y = M_pred_mu)) +
  geom_point(shape = 16, alpha = 0.5) +
  geom_errorbar(mapping = aes(x = M_obs, ymin = M_pred_PIlow, ymax = M_pred_PIhigh),
                width = 0, alpha = 0.5, size = 0.25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  ylab("Predicted monoculture (g)") +
  xlab("Observed monoculture (g)") +
  theme_meta()
plot(p1)

# save as an object into the results folder
save("p1", file = here("results/SP_A_monoculture_plot.RData"))

# ggsave(here("figures/fig_barn_mono.png"), p1,
       # width = 8, height = 7, units = "cm")

### END
