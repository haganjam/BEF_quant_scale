#'
#' @title: Seasq: Model the missing monocultures
#' 
#' @description: This script attemps to model the missing monocultures
#' for the Seasq species after preliminary data analysis showed that a single
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
v <- data[data$OTU %in% c("Seasq"), ]

# remove NAs from the monoculture column
v <- v[complete.cases(v),]

# make a data.list with the training data
seasq <- 
  list(M = v$M,
       Y = (v$Y),
       Y2 = (v$Y)^2,
       PC1 = v$PC1
  )

# model1
m.E1 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bY2*Y2 + bPC1*PC1,
          
          bY ~ dnorm( 0 , 1 ),
          bY2 ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = seasq, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.E1 )

# check the traceplots
traceplot(m.E1)
dev.off()

# plot the Y vs M comparison
range(seasq$Y)
Y_seq <- seq(-0.9, 2.15, 0.05)
PC1_seq <- mean(seasq$PC1)
pred_dat <- list(Y = Y_seq, Y2 = Y_seq^2, PC1 = rep(PC1_seq, length(Y_seq)) )
post <- link(m.E1, pred_dat)

plot(seasq$Y, seasq$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.E1)
plot(apply(post, 2, mean), seasq$M)
abline(0, 1)

# model2
m.E2 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bY2*Y2,
          
          bY ~ normal( 0 , 1 ),
          bY2 ~ normal( 0 , 1 ),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = seasq, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.E2 )

# check the traceplots
traceplot(m.E2)
dev.off()

# plot the Y vs M comparison
Y_seq <- seq(-0.9, 2.15, 0.05)
pred_dat <- list(Y = Y_seq, Y2 = Y_seq^2 )
post <- link(m.E2, pred_dat)

plot(seasq$Y, seasq$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.E2)
plot(apply(post, 2, mean), seasq$M)
abline(0, 1)

# model3: Intercept only
m.E3 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS,
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = seasq, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.E3 )

# check the traceplots
traceplot(m.E3)
dev.off()

# compare the three models
compare(m.E1, m.E2, m.E3, func = PSIS)
compare(m.E1, m.E2, m.E3, func = WAIC)

# sample from the m.E1 model
post <- rethinking::extract.samples(m.E1)

# write this posterior distribution list as a .rds object
saveRDS(post, here("results/SP_E_monoculture_posterior.rds"))

# save the model object as a .rds object
saveRDS(m.E1, here("results/SP_E_model_object.rds")) 

### END
