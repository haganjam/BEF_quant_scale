#'
#' @title: Analyse simulated data to test the analytical pipeline
#' 
#' @description: Examines accuracy of the pipeline on simulated metacommunities
#' 
#' @details: This script analyses the summary data generated from our different
#' simulated metacommunities.
#' 
#' PI_true: This measures whether the 90% percentile interval of the posterior distribution
#' is contained within 0.50 x (observed value) and 1.5 x (observed value). This quantity
#' shows whether interval generated is reasonably close to the observed value.
#' 
#' PI_mu_true: This simply measures whether the observed value is within the 
#' 90% percentile interval of the distribution.
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

# check the summary statistics and variable structures
summary(BEF_output)
str(BEF_output)

# view the dataset
View(BEF_output)

# how do we want to test the accuracy
BEF_output <- 
  BEF_output %>%
  mutate(INT_true = ifelse( (PI_true == TRUE) & (PI_mu_true) == TRUE, TRUE, FALSE ))

# is there a relationship between monoculture correlation and mu deviation
ggplot(data = BEF_output,
       mapping = aes(x = log10(mono_error), y = log10(mu_deviation) )) +
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
       mapping = aes(x = log10(mono_error), y = log10(PI_high-PI_low) )) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Beff, scales = "free") +
  theme_bw()

ggplot(data = BEF_output,
       mapping = aes(x = (mono_cor), y = log10(PI_high-PI_low) )) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Beff, scales = "free") +
  theme_bw()

# calculate accuracy metrics
BEF_output_sum <- 
  BEF_output %>%
  group_by(Beff) %>%
  summarise(PI_true_accuracy = sum(PI_true)/n(),
            PI_mu_true_accuracy = sum(PI_mu_true)/n(),
            accuracy = (sum(INT_true/n())))
print(BEF_output_sum)

# fit a binomial regression to model the accuracy based on the effect and the monoculture correlation

# check the distribution of correlation
hist(BEF_output$mono_cor)

m.dat <- 
  list(BE = as.integer(as.factor(BEF_output$Beff)),
       C = BEF_output$mono_cor,
       INT = ifelse(BEF_output$PI_mu_true == TRUE, 1, 0) )

m1 <- ulam(
  alist(
    INT ~ dbinom( 1 , p ),
    logit(p) <- a[BE] + b[BE]*C,
    
    a[BE] ~ dnorm( 0 , 2),
    b[BE] ~ dnorm(0, 2)
    
  ) , data = m.dat , chains = 4, cores = 4 )

# check the model outputs: Rhat values are good and traceplots look decent
precis( m1 , depth = 2 )
traceplot(m1)

# set-up a data.frame of data to simulate
m1.pred <- expand.grid(C = c(0.2, 0.8),
                       BE = unique(as.integer(as.factor(BEF_output$Beff))))

# use sim to simulate observations for this data.frame
m1.sim <- link(fit = m1, data = m1.pred)

# summarise these data
m1.pred$mu <- apply(m1.sim, 2, mean)
m1.pred$PI_low <- apply(m1.sim, 2, function(x) PI(samples = x, prob = 0.90)[1] )  
m1.pred$PI_high <- apply(m1.sim, 2, function(x) PI(samples = x, prob = 0.90)[2] )

# add the labels
m1.pred$BE1 <- rep(levels(as.factor(BEF_output$Beff))[m.dat$BE[1:11]], each = 2)

# convert the correlation to character
m1.pred$C <- as.character(m1.pred$C)

# add the observed values
BEF_output_sum <- 
  BEF_output_sum %>%
  rename(BE1 = Beff, mu = PI_mu_true_accuracy)

# plot the results
ggplot() +
  geom_errorbar(data = m1.pred, 
                mapping = aes(x = BE1, ymin = PI_low, ymax = PI_high, 
                              colour = C),
                width = 0,
                position = position_dodge(width = 0.75)) +
  geom_point(data = BEF_output_sum,
             mapping = aes(x = BE1, y = mu))
  
  
  


