#'
#' @title: Hydro: Model the missing monocultures
#' 
#' @description: This script attemps to model the missing monocultures
#' for the Hydro species after preliminary data analysis showed that a single
#' monoculture would not be suitable for all species.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# load the rethinking package
library(rethinking)

# load the analysis data
data <- read_csv(here("data/case_study_2/data_clean/biomass_env_analysis_data.csv"))

# pivot mixture data longer
mix <- 
  data %>%
  select(cluster_id, buoy_id, time, OTU, Y) %>%
  pivot_wider(id_cols = c("cluster_id", "buoy_id", "time"),
               names_from = "OTU",
               values_from = "Y")

# create a v data.frame with the relevant data
v <- data[data$OTU %in% c("Hydro"), ]

# join these data together
v <- full_join(mix, v, by = c("cluster_id", "buoy_id", "time"))

# remove NAs from the monoculture column
v <- v[complete.cases(v),]

# make a data.list with the training data
hydro <- 
  list(M = v$M,
       Y = v[["Hydro"]],
       Barn = v$Barn,
       Bryo = v$Bryo,
       Bumpi = v$Bumpi,
       Seasq = v$Seasq,
       PC1 = v$PC1,
       PC2 = v$PC2
  )

# check the correlations among the variables
pairs(hydro)

data.frame(v$M, v$Bryo, v$Bumpi, v$Seasq) %>%
  arrange(v.M)

# model1
m.D1 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bBr*Bryo + bYbr*Y*Bryo + bPC1*PC1 + bPC2*PC2 + bPCI*PC1*PC2,
          
          bY ~ dnorm( 0 , 1 ),
          bBr ~ dnorm( 0, 1),
          bYbr ~ dnorm( 0, 1),
          bPC1 ~ dnorm(0, 1),
          bPC2 ~ dnorm(0, 1),
          bPCI ~ dnorm(0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = hydro, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.D1  )

# check the traceplots
# traceplot(m.D1, pars = c("aS", "bPC1") )
# dev.off()

# plot the observed versus predicted data
post <- sim(m.D1)
plot(apply(post, 2, mean), hydro$M)
abline(0, 1)

# model2
m.D2 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC1Y*Y*PC1,
          
          bY ~ dnorm( 0 , 1 ),
          bBr ~ dnorm( 0, 1),
          bYbr ~ dnorm( 0, 1),
          bPC1 ~ dnorm(0, 1),
          bPC1Y ~ dnorm(0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = hydro, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.D2  )

# check the traceplots
# traceplot(m.D2, pars = c("aS", "bPC1") )
# dev.off()

# plot the observed versus predicted data
post <- sim(m.D2)
plot(apply(post, 2, mean), hydro$M)
abline(0, 1)

# model3
m.D3 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC1Y*Y*PC1 + bPC2*PC2,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm(0, 1),
          bPC2 ~ dnorm(0, 1),
          bPC1Y ~ dnorm(0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = hydro, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.D3  )

# check the traceplots
# traceplot(m.D3, pars = c("aS", "bPC1") )
# dev.off()

# plot the observed versus predicted data
post <- sim(m.D3)
plot(apply(post, 2, mean), hydro$M)
abline(0, 1)

# model4
m.D4 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC1Y*Y*PC1 + bPC2*PC2 + bPC2Y*Y*PC2,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm(0, 1),
          bPC2 ~ dnorm(0, 1),
          bPC1Y ~ dnorm(0, 1),
          bPC2Y ~ dnorm(0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = hydro, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.D4  )

# check the traceplots
# traceplot(m.D4, pars = c("aS", "bPC1", "bPC2") )
# dev.off()

# plot the observed versus predicted data
post <- sim(m.D4)
plot(apply(post, 2, mean), hydro$M)
abline(0, 1)

# check if any variables can explain some variation within the middle cloud
x <- apply(post, 2, mean)
bind_rows(hydro)[which((x > 0.045) & (x < 0.07) ),] %>%
  pairs()

# model5
m.D5 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC1Y*Y*PC1 + bPC2*PC2 + bPC2Y*Y*PC2 + bSq*Seasq,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm(0, 1),
          bPC2 ~ dnorm(0, 1),
          bPC1Y ~ dnorm(0, 1),
          bPC2Y ~ dnorm(0, 1),
          bSq ~ dnorm(0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = hydro, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.D5  )

# check the traceplots
# traceplot(m.D5, pars = c("aS", "bPC1", "bPC2") )
# dev.off()

# plot the observed versus predicted data
post <- sim(m.D5)
plot(apply(post, 2, mean), hydro$M)
abline(0, 1)

# model6
m.D6 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC1Y*Y*PC1 + bPC2*PC2 + bPC2Y*Y*PC2 + bSq*Seasq + bPC2Sq*PC2*Seasq,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm(0, 1),
          bPC2 ~ dnorm(0, 1),
          bPC1Y ~ dnorm(0, 1),
          bPC2Y ~ dnorm(0, 1),
          bSq ~ dnorm(0, 1),
          bPC2Sq ~ dnorm(0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = hydro, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.D6  )

# check the traceplots
# traceplot(m.D6, pars = c("aS", "bPC1", "bPC2") )
# dev.off()

# plot the observed versus predicted data
post <- sim(m.D6)
plot(apply(post, 2, mean), hydro$M)
abline(0, 1)

# model7
m.D7 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC2Y*Y*PC1 + bPC2*PC2 + bSq*Seasq + bPC2Sq*PC2*Seasq,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm(0, 1),
          bPC2 ~ dnorm(0, 1),
          bPC2Y ~ dnorm(0, 1),
          bSq ~ dnorm(0, 1),
          bPC2Sq ~ dnorm(0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = hydro, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.D7  )

# check the traceplots
# traceplot(m.D7, pars = c("aS", "bPC1", "bPC2") )
# dev.off()

# plot the observed versus predicted data
post <- sim(m.D7)
plot(apply(post, 2, mean), hydro$M)
abline(0, 1)

# model8
m.D8 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS + bY*Y + bPC1*PC1 + bPC1Y*Y*PC1 + bPC2*PC2 + bPC2Y*Y*PC2 + bSq*Seasq + bPC2Sq*PC2*Seasq + bBr*Bryo,
          
          bY ~ dnorm( 0 , 1 ),
          bPC1 ~ dnorm(0, 1),
          bPC2 ~ dnorm(0, 1),
          bPC1Y ~ dnorm(0, 1),
          bPC2Y ~ dnorm(0, 1),
          bSq ~ dnorm(0, 1),
          bPC2Sq ~ dnorm(0, 1),
          bBr ~ dnorm(0, 1),
          
          aS ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = hydro, chains = 4, log_lik = TRUE)

# check the precis output
precis( m.D8  )

# check the traceplots
# traceplot(m.D8, pars = c("aS", "bPC1", "bPC2") )
# dev.off()

# plot the observed versus predicted data
post <- sim(m.D8)
plot(apply(post, 2, mean), hydro$M)
abline(0, 1)

# based on these models, model 6 is the best
compare(m.D1, m.D2, m.D3, m.D4, 
        m.D5, m.D6, m.D7, m.D8, func = PSIS)

# sample from the m.D6 model
post <- rethinking::extract.samples(m.D6)

# write this posterior distribution list as a .rds object
saveRDS(post, here("results/SP_D_monoculture_posterior.rds"))

# save the model object as a .rds object
saveRDS(m.D6, here("results/SP_D_model_object.rds")) 

# make a plot of the observed and predicted values
source(here("scripts/Function_plotting_theme.R"))
pred <- sim(m.D6)
df.pred <- data.frame(M_obs = hydro$M,
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
saveRDS(p1, file = here("results/SP_D_monoculture_plot.rds"))

# ggsave(here("figures/fig_hydro_mono.png"), p1,
       # width = 8, height = 7, units = "cm")

### END
