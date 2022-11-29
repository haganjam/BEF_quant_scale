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
library(ggpubr)
library(rethinking)

# load relevant functions
source(here("scripts/Function_plotting_theme.R"))

# load the data
mono_init_dat <- readRDS(here("results/BEF_output.rds"))
head(mono_init_dat)

init_dat <- readRDS(here("results/BEF_output2.rds"))
head(init_dat)

# Table S2: calculate the range of observed biodiversity effects observed
mono_init_dat %>%
  group_by(Beff) %>%
  summarise(min_value = min(Value_obs),
            max_value = max(Value_obs),
            mean_value = mean(Value_obs),
            sd_value = sd(Value_obs))

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

# overwrite the objects the new accuracy information  
mono_init_dat <- accuracy[[1]]
init_dat <- accuracy[[2]]


## incomplete monocultures unknown initial relative abundance

# calculate observed summaries of accuracy
mono_init_sum <- 
  mono_init_dat %>%
  group_by(Beff) %>%
  summarise(PI_obs_true = sum(PI_obs_true)/n(),
            PI_true = sum(PI_true)/n(),
            mu_deviation_m = log10( mean(mu_deviation_perc) ),
            mu_deviation_sd = log10( sd(mu_deviation_perc) )) %>%
  rename(BE = Beff)
print(mono_init_sum)

# reorder the levels
mono_init_sum$BE <- factor(mono_init_sum$BE,
                           levels = c("LC", "TC", "LS", "TS", "NBE", "NO", "IT", "AS", "TI", "SI", "ST"))


# accuracy test 1: PI_obs_true

# fit a binomial regression to model the accuracy based on the effect and the monoculture correlation
m.dat <- 
  list(BE = as.integer(as.factor(mono_init_dat$Beff)),
       INT = ifelse(mono_init_dat$PI_obs_true == TRUE, 1, 0) )

m1 <- ulam(
  alist(
    INT ~ dbinom( 1 , p ),
    logit(p) <- a[BE],
    
    a[BE] ~ dnorm( 0 , 3)
    
  ) , data = m.dat , chains = 4, cores = 4 )

# check the model outputs: Rhat values are good and traceplots look decent
precis( m1 , depth = 2 )
traceplot( m1 )

# set-up a data.frame of data to simulate
m1.pred <- expand.grid(BE = unique(as.integer(as.factor(mono_init_dat$Beff))) )

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
m1.post$BE <- rep(levels(as.factor(mono_init_dat$Beff))[m.dat$BE[1:11]], each = n)

# change the order of the effects
m1.post$BE <- factor(m1.post$BE,
                     levels = c("LC", "TC", "LS", "TS", "NBE", "NO", "IT", "AS", "TI", "SI", "ST"))

# plot the results
p1 <- 
  ggplot() +
  geom_quasirandom(data = m1.post,
                   mapping = aes(x = BE, y = Value, colour = BE),
                   alpha = 0.05, width = 0.2) +
  geom_point(data = mono_init_sum,
             mapping = aes(x = BE, y = PI_obs_true),
             size = 2, shape = 21, colour = "black", fill = "white") +
  scale_y_continuous(limits = c(0.3, 0.9)) +
  scale_colour_manual(values = v_col_BEF()) +
  ylab("Prop: PI-5% < obs. < PI-95%") +
  xlab(NULL) +
  ggtitle("Uncertainty: monocultures + initial RA") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_meta() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

plot(p1)


# model the mu deviation percentage

# fit a log-normal regression model to the data
m.dat2 <- 
  list(BE = as.integer(as.factor(mono_init_dat$Beff)),
       MUD =  mono_init_dat$mu_deviation_perc)

m2 <- ulam(
  alist(
    MUD ~ dlnorm(mu, sigma),
    mu <- a[BE],
    
    a[BE] ~ dnorm( 0 , 3),
    sigma ~ dexp(1)
    
  ) , data = m.dat2 , chains = 4, cores = 4 )

# check the model outputs: Rhat values are good and traceplots look decent
precis( m2 , depth = 2 )
traceplot( m2 )

# set-up a data.frame of data to simulate
m2.pred <- expand.grid(BE = unique(as.integer(as.factor(mono_init_dat$Beff))) )

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
m2.post$BE <- rep(levels(as.factor(mono_init_dat$Beff))[m.dat2$BE[1:11]], each = n)

# change the order of the effects
m2.post$BE <- factor(m2.post$BE,
                     levels = c("LC", "TC", "LS", "TS", "NBE", "NO", "IT", "AS", "TI", "SI", "ST"))

# plot the results
p2 <- 
  ggplot() +
  geom_quasirandom(data = m2.post,
                   mapping = aes(x = BE, y = log10(Value), colour = BE),
                   alpha = 0.15, width = 0.4) +
  geom_point(data = mono_init_sum,
             mapping = aes(x = BE, y = (mu_deviation_m) ), 
             size = 2, shape = 21, colour = "black", fill = "white") +
  scale_colour_manual(values = v_col_BEF()) +
  geom_hline(yintercept = 1.69, linetype = "dashed") + # 50% absolute deviation from observed
  ylab("log10 absolute deviation (%)") +
  ggtitle("") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")

plot(p2)


## unknown initial relative abundance

# which biodiversity effects are sensitive to initial RA?
init_BEF <- c("LS", "TS", "NBE", "IT", "AS")

# subset these effects
init_dat <- 
  init_dat %>%
  filter(Beff %in% init_BEF)

# calculate observed summaries of accuracy
init_sum <- 
  init_dat %>%
  group_by(Beff) %>%
  summarise(PI_obs_true = sum(PI_obs_true)/n(),
            PI_true = sum(PI_true)/n(),
            mu_deviation_m = mean( (mu_deviation_perc) ),
            mu_deviation_sd = sd( (mu_deviation_perc) )) %>%
  rename(BE = Beff)
print(init_sum)

# reorder the levels
init_sum$BE <- factor(init_sum$BE,
                           levels = c("LS", "TS", "NBE", "IT", "AS"))


# accuracy test 1: PI_obs_true

# fit a binomial regression to model the accuracy based on the effect and the monoculture correlation
m.dat3 <- 
  list(BE = as.integer(as.factor(init_dat$Beff)),
       INT = ifelse(init_dat$PI_obs_true == TRUE, 1, 0) )

m3 <- ulam(
  alist(
    INT ~ dbinom( 1 , p ),
    logit(p) <- a[BE],
    
    a[BE] ~ dnorm( 0 , 2)
    
  ) , data = m.dat3 , chains = 4, cores = 4 )

# check the model outputs: Rhat values are good and traceplots look decent
precis( m3 , depth = 2 )
traceplot( m3 )

# set-up a data.frame of data to simulate
m3.pred <- expand.grid(BE = unique(as.integer(as.factor(init_dat$Beff))) )

# use sim to simulate observations for this data.frame
m3.sim <- link(fit = m3, data = m3.pred)

# check how many samples to use
n <- 500

# set-up the initial values
m3.post <- data.frame(Value = m3.sim[, 1][sample(x = 1:nrow(m3.sim), n)]  )
m3.post <- bind_cols(data.frame(BE = m3.pred[1,]), m3.post)

# loop over the different rows
for(i in 2:nrow(m3.pred)) {
  
  x <- data.frame(Value = m3.sim[, i][sample(x = 1:nrow(m3.sim), n)])
  y <- bind_cols(data.frame(BE = m3.pred[i,]), x)
  m3.post <- bind_rows(m3.post, y)
  
}

# replace the numbers with the effect codes
m3.post$BE <- rep(levels(as.factor(init_dat$Beff))[m.dat3$BE[1:5]], each = n)

# change the order of the effects
m3.post$BE <- factor(m3.post$BE,
                     levels = c("LS", "TS", "NBE", "IT", "AS"))

# plot the results
p3 <- 
  ggplot() +
  geom_quasirandom(data = m3.post,
                   mapping = aes(x = BE, y = Value, colour = BE),
                   alpha = 0.05, width = 0.2) +
  geom_point(data = init_sum,
             mapping = aes(x = BE, y = PI_obs_true),
             size = 2, shape = 21, colour = "black", fill = "white") +
  scale_y_continuous(limits = c(0.3, 0.9)) +
  scale_colour_manual(values = v_col_BEF(c("LS", "TS", "NBE", "IT", "AS"))) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ggtitle("Uncertainty: initial RA") +
  ylab("Prop: PI-5% < obs. < PI-95%") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

plot(p3)


# fit a log-normal regression model to the data
m.dat4 <- 
  list(BE = as.integer(as.factor(init_dat$Beff)),
       MUD = (init_dat$mu_deviation_perc))

m4 <- ulam(
  alist(
    MUD ~ dlnorm(mu, sigma),
    mu <- a[BE],
    
    a[BE] ~ dnorm( 0 , 3),
    sigma ~ dexp(1)
    
  ) , data = m.dat4 , chains = 4, cores = 4 )

# check the model outputs: Rhat values are good and traceplots look decent
precis( m4 , depth = 2 )
traceplot( m4 )

# set-up a data.frame of data to simulate
m4.pred <- expand.grid(BE = unique(as.integer(as.factor(init_dat$Beff))) )

# use sim to simulate observations for this data.frame
m4.sim <- sim(fit = m4, data = m4.pred)

# check how many samples to use
n <- 500

# set-up the initial values
m4.post <- data.frame(Value = m4.sim[, 1][sample(x = 1:nrow(m4.sim), n)]  )
m4.post <- bind_cols(data.frame(BE = m4.pred[1,]), m4.post)

# loop over the different rows
for(i in 2:nrow(m4.pred)) {
  
  x <- data.frame(Value = m4.sim[, i][sample(x = 1:nrow(m4.sim), n)])
  y <- bind_cols(data.frame(BE = m4.pred[i,]), x)
  m4.post <- bind_rows(m4.post, y)
  
}

# replace the numbers with the effect codes
m4.post$BE <- rep(levels(as.factor(init_dat$Beff))[m.dat4$BE[1:5]], each = n)

# change the order of the effects
m4.post$BE <- factor(m4.post$BE,
                     levels = c("LS", "TS", "NBE", "IT", "AS"))

# plot the results
p4 <- 
  ggplot() +
  geom_quasirandom(data = m4.post,
                   mapping = aes(x = BE, y = (Value), colour = BE),
                   alpha = 0.1, width = 0.4) +
  scale_colour_manual(values = v_col_BEF(c("LS", "TS", "NBE", "IT", "AS"))) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  geom_point(data = init_sum,
             mapping = aes(x = BE, y = (mu_deviation_m)), 
             size = 2, shape = 21, colour = "black", fill = "white") +
  ylab("Absolute deviation (%)") +
  xlab(NULL) +
  ggtitle("") +
  theme_meta() +
  theme(legend.position = "none")

plot(p4)

# arrange the plot
p1234 <- 
  ggarrange(p3, p1, p4, p2, 
            ncol = 2, nrow = 2,
            widths = c(1, 1.7),
            labels = c("a", "b", "c", "d"),
            font.label = list(size = 11, face = "plain"))

plot(p1234)

ggsave(filename = here("figures/fig1.png"), p1234, dpi = 400,
       unit = "cm", width = 20, height = 15)

### END
