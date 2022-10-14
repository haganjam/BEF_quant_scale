#'
#' @title: Bumpi: Model the missing monocultures
#' 
#' @description: This script attemps to model the missing monocultures
#' for the Bumpi species after preliminary data analysis showed that a single
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
v <- data[data$OTU %in% c("Bumpi"), ]

# make a data.list with the training data
bumpi <- 
  list(M = v$M,
       Y = scale(v$Y),
       Y2 = scale(v$Y)^2,
       PC1 = v$PC1
  )

# model1
m.C1 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bY2*Y2 + bPC1*PC1,
          
          bY ~ dnorm( 0 , 1 ),
          bY2 ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm( 0, 1 ),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = bumpi, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.C1  )

# check the traceplots
traceplot(m.C1 )
dev.off()

# plot the Y vs M comparison
range(bumpi$Y)
Y_seq <- seq(-0.92, 1.4, 0.05)
PC1_seq <- mean(bumpi$PC1)
pred_dat <- list(Y = Y_seq, Y2 = Y_seq^2, PC1 = rep(PC1_seq, length(Y_seq)) )
post <- link(m.C1, pred_dat)

plot(bumpi$Y, bumpi$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.C1)
plot(apply(post, 2, mean), bumpi$M)
abline(0, 1)

# model2
m.C2 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bY2*Y2,
          
          bY ~ dnorm( 0 , 1 ),
          bY2 ~ dnorm( 0 , 1 ),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = bumpi, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.C2  )

# check the traceplots
traceplot(m.C2 )
dev.off()

# plot the Y vs M comparison
range(bumpi$Y)
Y_seq <- seq(-0.92, 1.4, 0.05)
pred_dat <- list(Y = Y_seq, Y2 = Y_seq^2 )
post <- link(m.C2, pred_dat)

plot(bumpi$Y, bumpi$M)
lines(pred_dat$Y, apply(post, 2, mean))

# plot the observed versus predicted data
post <- sim(m.C2)
plot(apply(post, 2, mean), bumpi$M)
abline(0, 1)

# model3: Intercept only
m.C3 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS,
          
          aS ~ dnorm(0,1),
          
          sigma ~ dexp(1)
          
    ),
    data = bumpi, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.C3  )

# check the traceplots
traceplot(m.C3 )
dev.off()

# compare the three models
compare(m.C1, m.C2, m.C3, func = PSIS)
compare(m.C1, m.C2, m.C3, func = WAIC)

### END
