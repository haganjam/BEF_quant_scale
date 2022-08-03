#'
#' @title: Analyse simulated data to test the analytical pipeline
#' 
#' @description: Examines accuracy of the pipeline on simulated metacommunities
#' 
#' @details: This script analyses the summary data generated from our different
#' simulated metacommunities.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(here)
library(dplyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# load the data
BEF_output <- readRDS(here("results/BEF_output.rds"))
head(BEF_output)

# check the summary statistics and variable structures
summary(BEF_output)
str(BEF_output)

# view the dataset
View(BEF_output)

# is there a relationship between monoculture correlation and mu deviation
ggplot(data = BEF_output,
       mapping = aes(x = log(mono_error), y = (mu_deviation) )) +
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

BEF_output_sum <- 
  BEF_output %>%
  group_by(Beff) %>%
  summarise(accuracy_value = sum(PI_true)/n(),
            accuracy_interval = sum(PI_mu_true)/n(),
            Value_obs = mean(Value_obs))

# plot the full set of simulations
ggplot(data = BEF_output_sum,
       mapping = aes(x = Beff, y = Value_obs)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_quasirandom(data = BEF_output %>% filter(mu < 10000 & mu > -10000), 
                   mapping = aes(x = Beff, y = mu), alpha = 0.1 ) +
  geom_point(colour = "red", size = 2) +
  theme_bw()
