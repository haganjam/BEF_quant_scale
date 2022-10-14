#'
#' @title: Bumpi: Model the missing monocultures
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

# make a data.list with the training data
seasq <- 
  list(M = v$M,
       Y = scale(v$Y),
       Y2 = scale(v$Y)^2,
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

# plot the observed and modeled data with error
post <- sim(m.E2)
mu <- apply(post, 2, mean)
PI_low <- apply(post, 2, function(x) PI(x, 0.90)[1] )
PI_high <- apply(post, 2, function(x) PI(x, 0.90)[2] )

# are any mu values less than zero?
min(mu)

# plot predictions with error
tibble(obs = seasq$M,
       mu = mu,
       PI_low = PI_low,
       PI_high = PI_high) %>%
  ggplot( data = .,
          mapping = aes(x = obs, y = mu)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_errorbar(mapping = aes(x = obs, 
                              ymin = PI_low,
                              ymax = PI_high)) +
  theme_meta()

### END
