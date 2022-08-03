
# examine the accuracy

# notes: PI_mu_true
# whether the 90% PI interval is within 50% of the observed value on either side
# PI_mu is whether the observed value is within the 90% interval

# load the relevant libraries
library(here)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# load the data
BEF_output <- readRDS(here("results/BEF_output.rds"))
View(BEF_output)

# analyse the accuracy
View(BEF_output)

# check for outliers
summary(BEF_output)

# check the large mono-error outlier
BEF_output %>%
  filter(mono_error == max(mono_error)) %>%
  View()

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
