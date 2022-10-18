#'
#' @title: Barn: Model the missing monocultures
#' 
#' @description: This script attemps to model the missing monocultures
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
       Y = scale(v$Y),
       Y2 = scale(v$Y)^2,
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
plot(apply(post, 2, mean), barn$M)
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
post

View(data)

# assign input variables to the names
Barn_NA <- 
  data %>%
  filter(OTU == "Barn", is.na(M))

var_names <- names(m.A3@data)
var_names <- var_names[var_names != "M"]

for(i in var_names) {
  
  assign(x = i, Barn_NA[[i]])
  
}

# assign the parameter values to the names
post_names <- names(post)
post_samp <- sapply(post, function(x) sample(x, 1) )
names(post_samp) <- NULL

# write a loop and assign a sample from the posterior distribution to a parameter name
for (j in 1:length(post_names)) {
  
  assign(x = post_names[j], value = post_samp[j])
  
}

# calculate mu: set-up the expression
form <- parse(text = m.A3@formula[[2]][[3]])

# evaluate the expression
u <- eval(form)

# run the 
dist <- parse(text = m.A3@formula[[1]][[3]])
gsub(pattern = "d", replacement = "r", x = m.A3@formula[[1]][[3]] )
eval(dist)

### END
