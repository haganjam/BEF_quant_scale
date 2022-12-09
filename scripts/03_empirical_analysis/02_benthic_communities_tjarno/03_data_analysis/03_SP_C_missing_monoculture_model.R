#'
#' @title: Bumpi: Model the missing monocultures
#' 
#' @description: This script attempts to model the missing monocultures
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

# remove NAs from the monoculture column
v <- v[complete.cases(v),]

# make a data.list with the training data
bumpi <- 
  list(M = v$M,
       Y = (v$Y),
       Y2 = (v$Y)^2,
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
# traceplot(m.C1 )
# dev.off()

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
# precis( m.C2  )

# check the traceplots
# traceplot(m.C2 )
# dev.off()

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
# traceplot(m.C3 )
# dev.off()

# compare the three models
compare(m.C1, m.C2, m.C3, func = PSIS)
compare(m.C1, m.C2, m.C3, func = WAIC)

# sample from the m.C3 model
post <- rethinking::extract.samples(m.C3)

# write this posterior distribution list as a .rds object
saveRDS(post, here("results/SP_C_monoculture_posterior.rds"))

# save the model object as a .rds object
saveRDS(m.C3, here("results/SP_C_model_object.rds")) 

# make a plot of the observed and predicted values
source(here("scripts/Function_plotting_theme.R"))
pred <- sim(m.C3)
df.pred <- data.frame(M_obs = bumpi$M,
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
saveRDS(p1, file = here("results/SP_C_monoculture_plot.rds"))

# ggsave(here("figures/fig_bumpi_mono.png"), p1,
       # width = 8, height = 7, units = "cm")

### END
