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

# load the data
BEF_output <- readRDS(here("results/BEF_output.rds"))
head(BEF_output)
names(BEF_output)

# check the summary statistics and variable structures
summary(BEF_output)
str(BEF_output)

# view the dataset
View(BEF_output)

# examine the distribution of modelled effects
hist(BEF_output$mu)
range(BEF_output$mu)
mean(BEF_output$mu)

# remove the major, major outliers
BEF_output$cond <- (BEF_output$mu < quantile(BEF_output$mu, 0.99)) & (BEF_output$mu > quantile(BEF_output$mu, 0.01))

BEF_output <- 
  BEF_output %>%
  group_by(model_ID) %>%
  mutate(outliers = ifelse(any(cond == FALSE), 1, 0)) %>%
  filter(outliers != 1) %>%
  select(-cond)

# set-up variables for testing the accuracy

# set the mean range accuracy: 
thresh <- 0.75

BEF_output <- 
  
  BEF_output %>%
  mutate(mu_deviation = round( abs((mu - Value_obs)), 5 )  ) %>%
  mutate(mu_deviation_perc = round( abs(mu_deviation/Value_obs)*100, 5) ) %>%
  mutate(mu_threshold = thresh) %>%
  mutate(Value_obs_range = abs(2*thresh*Value_obs) ) %>%
  mutate(PI_range = PI_high-PI_low) %>%
  mutate(PI_range_true = ifelse( PI_range > Value_obs_range, FALSE, TRUE)  ) %>%
  mutate(PI_obs_true = ifelse( (PI_low < Value_obs) & (PI_high > Value_obs), TRUE, FALSE ) ) %>%
  mutate(PI_true = ifelse( (PI_range_true & PI_obs_true), TRUE, FALSE ))
View(BEF_output)  
summary(BEF_output)

# is there a relationship between monoculture correlation and mu deviation
ggplot(data = BEF_output,
       mapping = aes(x = log10(mono_error), y = (mu_deviation) )) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Beff, scales = "free") +
  theme_bw()

ggplot(data = BEF_output,
       mapping = aes(x = (mono_cor), y = log10(mu_deviation) )) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Beff, scales = "free") +
  theme_bw()

ggplot(data = BEF_output,
       mapping = aes(x = (mono_PI_range), y = log10(mu_deviation) )) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Beff, scales = "free") +
  theme_bw()

# calculate accuracy metrics
BEF_output_sum <- 
  BEF_output %>%
  group_by(Beff) %>%
  filter(mu_deviation_perc < quantile(mu_deviation_perc, 0.95)) %>%
  summarise(PI_obs_true = sum(PI_obs_true)/n(),
            PI_true = sum(PI_true)/n(),
            mu_deviation_m = mean(mu_deviation_perc),
            mu_deviation_sd = sd(mu_deviation_perc)) %>%
  rename(BE1 = Beff)
print(BEF_output_sum)

# what is the average monoculture correlation?
mean(BEF_output$mono_cor)
hist(BEF_output$mono_cor)

# accuracy test 1: PI_true

# fit a binomial regression to model the accuracy based on the effect and the monoculture correlation
m.dat <- 
  list(BE = as.integer(as.factor(BEF_output$Beff)),
       C = BEF_output$mono_cor,
       MR = BEF_output$mono_PI_range,
       ME = log10(BEF_output$mono_error),
       INT = ifelse(BEF_output$PI_true == TRUE, 1, 0) )

m1 <- ulam(
  alist(
    INT ~ dbinom( 1 , p ),
    logit(p) <- a[BE] + b[BE]*C + b1[BE]*MR,
    
    a[BE] ~ dnorm( 0 , 2),
    b[BE] ~ dnorm(0, 2),
    b1[BE] ~ dnorm(0, 2)
    
  ) , data = m.dat , chains = 4, cores = 4 )

# check the model outputs: Rhat values are good and traceplots look decent
precis( m1 , depth = 2 )
traceplot( m1 )

saveRDS(m1, file = here("results/stan_model_m1.rds"))

# set-up a data.frame of data to simulate

# choose the correlations
cor.in <- c(0.1, 0.9)
mono_range <- c(5000)

m1.pred <- expand.grid(C = cor.in,
                       ME = log10(mono_range),
                       BE = unique(as.integer(as.factor(BEF_output$Beff))))

# use sim to simulate observations for this data.frame
m1.sim <- link(fit = m1, data = m1.pred)

# summarise these data
m1.pred$mu <- apply(m1.sim, 2, mean)
m1.pred$PI_low <- apply(m1.sim, 2, function(x) PI(samples = x, prob = 0.90)[1] )  
m1.pred$PI_high <- apply(m1.sim, 2, function(x) PI(samples = x, prob = 0.90)[2] )

# add the labels
m1.pred$BE1 <- rep(levels(as.factor(BEF_output$Beff))[m.dat$BE[1:11]], each = length(cor.in)*length(mono_range))

# plot the results
ggplot() +
  geom_errorbar(data = m1.pred %>% mutate(C = as.character(C)), 
                mapping = aes(x = BE1, ymin = PI_low, ymax = PI_high, 
                              colour = C),
                width = 0,
                position = position_dodge(width = 0.75)) +
  geom_point(data = BEF_output_sum,
             mapping = aes(x = BE1, y = PI_true)) +
  theme_bw()
  