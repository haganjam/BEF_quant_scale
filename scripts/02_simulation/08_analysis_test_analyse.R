#'
#' @title: Analyse simulated data to test the analytical pipeline
#' 
#' @description: Examines accuracy of the pipeline on simulated metacommunities
#' 
#' @details: This script analyses the summary data generated from our different
#' simulated metacommunities. We use three key metrics of accuracy:
#' 
#' 1. mu_deviation: This measures the absolute deviation of the mean of the posterior
#' from the observed value as a percentage of the observed value.
#' 
#' 2. PI_obs_true: This measures whether the observed biodiversity effects fall
#' within the 90% percentile interview of the posterior distribution
#' 
#' 3. PI_true: The problem with PI_obs_true is that when there is a lot of uncertainty,
#' the PI 90% can be extremely wide. Therefore, whilst the 90% PI may contain the
#' observed value, it is essentially meaningless as the 90% PI is so wide. To solve
#' this problem, we determine which 90% PI ranges are less than the range of the observed
#' value +- thresh*observed value. If the 90% PI range is less than this threshold range
#' and the observed value falls within the range, then PI_true is given as TRUE.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(here)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(rethinking)

# load relevant functions
source(here("scripts/Function_plotting_theme.R"))

# load the data
mono_init_dat <- readRDS(here("results/BEF_output.rds"))
head(mono_init_dat)

init_dat <- readRDS(here("results/BEF_output2.rds"))
head(init_dat)

# remove the models where the max_mono_width is less than 200
mono_init_dat  <- 
  mono_init_dat  %>%
  filter(max_mono_width < 300) %>%
  filter(mono_cor > 0.7)

# how many models are left over?
length(unique(mono_init_dat$model_ID))
length(unique(init_dat$model_ID))

# set the mean range accuracy: 
thresh <- 0.75

# calculate the accuracy metrics
accuracy <- 
  
  lapply(list(mono_init_dat, init_dat), function(BEF_output) {
  
  # set-up variables for testing the accuracy
  x <- 
    
    BEF_output %>%
    mutate(mu_deviation = round( abs((mu - Value_obs)), 5 )  ) %>%
    mutate(mu_deviation_perc = round( abs(mu_deviation/Value_obs)*100, 5) ) %>%
    mutate(mu_threshold = thresh) %>%
    mutate(Value_obs_range = abs(2*thresh*Value_obs) ) %>%
    mutate(PI_range = PI_high-PI_low) %>%
    mutate(PI_range_true = ifelse( PI_range > Value_obs_range, FALSE, TRUE)  ) %>%
    mutate(PI_obs_true = ifelse( (PI_low < Value_obs) & (PI_high > Value_obs), TRUE, FALSE ) ) %>%
    mutate(PI_true = ifelse( (PI_range_true & PI_obs_true), TRUE, FALSE ))
  
  return(x) 
  
  } )
  

## incomplete monocultures unknown initial relative abundance

# accuracy test 1: PI_obs_true

# fit a binomial regression to model the accuracy based on the effect and the monoculture correlation
m.dat <- 
  list(BE = as.integer(as.factor(BEF_output$Beff)),
       INT = ifelse(BEF_output$PI_obs_true == TRUE, 1, 0) )

m1 <- ulam(
  alist(
    INT ~ dbinom( 1 , p ),
    logit(p) <- a[BE],
    
    a[BE] ~ dnorm( 0 , 2)
    
  ) , data = m.dat , chains = 4, cores = 4 )

# check the model outputs: Rhat values are good and traceplots look decent
# precis( m1 , depth = 2 )
# traceplot( m1 )

# set-up a data.frame of data to simulate
m1.pred <- expand.grid(BE = unique(as.integer(as.factor(BEF_output$Beff))) )

# use sim to simulate observations for this data.frame
m1.sim <- link(fit = m1, data = m1.pred)

# check how many samples to use
n <- 500

# set-up the initial values
m1.post <- data.frame(Value = m1.sim[, 1][sample(x = 1:nrow(m1.sim), n)]  )
m1.post <- bind_cols(data.frame(BE = m1.pred[1,]), m1.post)

# loop over the different rows
for(i in 2:nrow(m1.pred)) {
  
  x <- data.frame(Value = m1.sim[, i][sample(x = 1:nrow(m1.sim), n)])
  y <- bind_cols(data.frame(BE = m1.pred[i,]), x)
  m1.post <- bind_rows(m1.post, y)
  
}

# replace the numbers with the effect codes
m1.post$BE <- rep(levels(as.factor(BEF_output$Beff))[m.dat$BE[1:11]], each = n)

# change the order of the effects
m1.post$BE <- factor(m1.post$BE,
                     levels = c("LC", "LS", "TC", "TS", "NBE", "NO", "IT", "AS", "TI", "SI", "ST"))

# change the order of the observed
BEF_output_sum$BE <- factor(BEF_output_sum$BE,
                            levels = c("LC", "LS", "TC", "TS", "NBE", "NO", "IT", "AS", "TI", "SI", "ST"))

# model the mu deviation percentage

# fit a log-normal regression model to the data
m.dat2 <- 
  list(BE = as.integer(as.factor(BEF_output$Beff)),
       MUD =  BEF_output$mu_deviation_perc)

m2 <- ulam(
  alist(
    MUD ~ dlnorm(mu, sigma),
    mu <- a[BE],
    
    a[BE] ~ dnorm( 0 , 3),
    sigma ~ dexp(1)
    
  ) , data = m.dat2 , chains = 4, cores = 4 )

# check the model outputs: Rhat values are good and traceplots look decent
# precis( m2 , depth = 2 )
# traceplot( m2 )

# set-up a data.frame of data to simulate
m2.pred <- expand.grid(BE = unique(as.integer(as.factor(BEF_output$Beff))) )

# use sim to simulate observations for this data.frame
m2.sim <- sim(fit = m2, data = m2.pred)

# check how many samples to use
n <- 500

# set-up the initial values
m2.post <- data.frame(Value = m2.sim[, 1][sample(x = 1:nrow(m2.sim), n)]  )
m2.post <- bind_cols(data.frame(BE = m2.pred[1,]), m2.post)

# loop over the different rows
for(i in 2:nrow(m2.pred)) {
  
  x <- data.frame(Value = m2.sim[, i][sample(x = 1:nrow(m2.sim), n)])
  y <- bind_cols(data.frame(BE = m2.pred[i,]), x)
  m2.post <- bind_rows(m2.post, y)
  
}

# replace the numbers with the effect codes
m2.post$BE <- rep(levels(as.factor(BEF_output$Beff))[m.dat2$BE[1:11]], each = n)

# change the order of the effects
m2.post$BE <- factor(m2.post$BE,
                     levels = c("LC", "LS", "TC", "TS", "NBE", "NO", "IT", "AS", "TI", "SI", "ST"))


  
  
  
  
  
# plot the results
ggplot() +
  geom_quasirandom(data = accuracy_output[[1]]$m1.post,
                   mapping = aes(x = BE, y = Value, colour = BE),
                   alpha = 0.01, width = 0.2) +
  geom_point(data = accuracy_output[[1]]$BEF_output_sum,
             mapping = aes(x = BE, y = PI_obs_true, colour = BE)) +
  scale_y_continuous(limits = c(0.3, 0.9)) +
  scale_colour_manual(values = v_col_BEF()) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_meta() +
  theme(legend.position = "none")

# plot the results
ggplot() +
  geom_quasirandom(data = accuracy_output[[2]]$m1.post,
                   mapping = aes(x = BE, y = Value, colour = BE),
                   alpha = 0.01, width = 0.2) +
  geom_point(data = accuracy_output[[2]]$BEF_output_sum,
             mapping = aes(x = BE, y = PI_obs_true, colour = BE)) +
  scale_y_continuous(limits = c(0.3, 0.9)) +
  scale_colour_manual(values = v_col_BEF()) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_meta() +
  theme(legend.position = "none")

# plot the results
ggplot() +
  geom_quasirandom(data = m2.post,
                   mapping = aes(x = BE, y = log10(Value), colour = BE),
                   alpha = 0.025, width = 0.4) +
  geom_point(data = BEF_output_sum,
             mapping = aes(x = BE, y = (mu_deviation_m), colour = BE ), 
             size = 2) +
  geom_errorbar(data = BEF_output_sum,
                mapping = aes(x = BE, 
                              ymin = mu_deviation_m - mu_deviation_sd,
                              ymax = mu_deviation_m + mu_deviation_sd,
                              colour = BE),
                width = 0) +
  scale_colour_manual(values = v_col_BEF()) +
  geom_hline(yintercept = 1.69, linetype = "dashed") + # 50% absolute deviation from observed
  theme_meta() +
  theme(legend.position = "none")

### END
