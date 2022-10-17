#'
#' @title: Bryo: Model the missing monocultures
#' 
#' @description: This script attemps to model the missing monocultures
#' for the Bryo species after preliminary data analysis showed that a single
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
v <- data[data$OTU %in% c("Bryo"), ]

# make a data.list with the training data
bryo <- 
  list(M = v$M,
       Y = scale(v$Y),
       Y2 = scale(v$Y)^2,
       PC1 = v$PC1,
       PC2 = v$PC2
  )

# model1
m.B1 <- 
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
    data = bryo, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.B1  )

# check the traceplots
traceplot(m.B1 )
dev.off()

# plot the Y vs M comparison
range(bryo$Y)
Y_seq <- seq(-0.7, 2.2, 0.1)
PC1_seq <- mean(bryo$PC1)
PC2_seq <- mean(bryo$PC2)
pred_dat <- list(Y = Y_seq, Y2 = Y_seq^2, 
                 PC1 = rep(PC1_seq, length(Y_seq)),
                 PC2 = rep(PC2_seq, length(Y_seq)))
post <- link(m.B1, pred_dat)

plot(bryo$Y, bryo$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.B1)
plot(apply(post, 2, mean), bryo$M)
abline(0, 1)

# model2
m.B2 <- 
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
    data = bryo, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.B2  )

# check the traceplots
traceplot(m.B2 )
dev.off()

# plot the observed versus predicted data
post <- sim(m.B2)
plot(apply(post, 2, mean), bryo$M)
abline(0, 1)

# model3
m.B3 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPCY*PC1*Y,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          bPCY ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = bryo, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.B3  )

# check the traceplots
traceplot(m.B3 )
dev.off()

# plot the observed versus predicted data
post <- sim(m.B3)
plot(apply(post, 2, mean), bryo$M)
abline(0, 1)

# model4
m.B4 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bPC1*PC1 + bPC2*PC2 + bPCI*PC1*PC2,
          
          bPC1 ~ dnorm( 0, 1 ),
          bPC2 ~ dnorm( 0, 1),
          bPCI ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = bryo, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.B4  )

# check the traceplots
traceplot(m.B4 )
dev.off()

# plot the observed versus predicted data
post <- sim(m.B4)
plot(apply(post, 2, mean), bryo$M)
abline(0, 1)

# model5
m.B5 <- 
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
    data = bryo, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.B5  )

# check the traceplots
traceplot(m.B5 )
dev.off()

# plot the observed versus predicted data
post <- sim(m.B5)
plot(apply(post, 2, mean), bryo$M)
abline(0, 1)

# model6
m.B6 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC2*PC2 + bPCI*PC1*PC2 + bPCY*PC1*Y,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          bPC2 ~ dnorm( 0, 1),
          bPCI ~ dnorm( 0, 1),
          bPCY ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = bryo, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.B6  )

# check the traceplots
traceplot(m.B6 )
dev.off()

# plot the observed versus predicted data
post <- sim(m.B6)
plot(apply(post, 2, mean), bryo$M)
abline(0, 1)

# model7
m.B7 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC2*PC2 + bPCI*PC1*PC2 + bPCY*PC2*Y,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          bPC2 ~ dnorm( 0, 1),
          bPCI ~ dnorm( 0, 1),
          bPCY ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = bryo, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.B7)

# check the traceplots
traceplot(m.B7 )
dev.off()

# plot the observed versus predicted data
post <- sim(m.B7)
plot(apply(post, 2, mean), bryo$M)
abline(0, 1)

# model8
m.B8 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC2*PC2 + bPCI*PC1*PC2 + bPC1Y*PC2*Y + bPC2Y*PC2*Y,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          bPC2 ~ dnorm( 0, 1),
          bPCI ~ dnorm( 0, 1),
          bPC1Y ~ dnorm( 0, 1),
          bPC2Y ~ dnorm( 0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = bryo, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.B8)

# check the traceplots
traceplot(m.B8 )
dev.off()

# plot the observed versus predicted data
post <- sim(m.B8)
plot(apply(post, 2, mean), bryo$M)
abline(0, 1)


# compare the models
compare(m.B1, m.B2, m.B3, m.B4, 
        m.B5, m.B6, m.B7, m.B8,
        func = PSIS)

compare(m.B1, m.B2, m.B3,m.B4, 
        m.B5, m.B6, m.B7, m.B8,
        func = WAIC)

### END
