#'
#' @title: Link cleaned biomass data and environmental data and calculate BEF effects
#' 
#' @description: In this script, we link the clean biomass data with the clean
#' environmental data and process it in order to calculate the BEF effects.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(here)
library(vegan)

# load plotting theme
source(here("scripts/Function_plotting_theme.R"))

# load the cleaned environmental data
env_dat <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/site_env_data.csv"))
head(env_dat)

# select the relevant variables
names(env_dat)

env_dat <- 
  env_dat %>%
  select(cluster_id, site_id, time, 
         panel_depth_m_measured, distance_between_panel_and_seabeed_m,
         temp_C_m, temp_C_cv, temp_C_max, temp_C_min,
         lux_m, lux_cv, lux_max, lux_min) %>%
  rename(buoy_id = site_id, depth_m = panel_depth_m_measured, seabed_dist_m = distance_between_panel_and_seabeed_m)

# load the cleaned biomass data
bio_dat <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/biomass_data.csv"))
head(bio_dat)

# join the environmental data to the biomass data
bio_env <- full_join(env_dat, bio_dat, by = c("cluster_id", "buoy_id", "time"))

# only keep complete-cases
bio_env <- bio_env[complete.cases(bio_env[,names(bio_env) != "M"]), ]

# check at the distribution of the samples
bio_env %>%
  group_by(cluster_id) %>%
  summarise(n = length(unique(buoy_id)))

incomplete <- 
  bio_env %>%
  group_by(cluster_id, buoy_id) %>%
  summarise(n = length(unique(time))) %>%
  filter(n < 3) %>%
  pull(buoy_id)
print(incomplete)

# remove the incomplete buoys
bio_env <- 
  bio_env %>%
  filter( !(buoy_id %in% incomplete) )

# rename the variables
View(bio_env)

# check the correlation among these variables
pairs(bio_env[,-c(1:3, 14:16)])

# let's try build a model with three variables:
# 1. depth
# 2. temp_C_m
# 3. lux_C_m

# load the rethinking package
library(rethinking)

# make a data.list with the training data
v <- bio_env[complete.cases(bio_env), ]
head(v)
dim(v)

# check the distribution of the response variables
hist(sqrt(v$M))
min(v$M)
v$Y[v$Y >0] %>% min()
hist(log(0.005 + v$M))

hist(log(0.0005 + v$Y))

train.d <- 
  list(M = sqrt(v$M),
       Y = log(1 + v$Y),
       TMP = (v$temp_C_m - mean(v$temp_C_m))/sd(v$temp_C_m),
       LUX = (v$lux_m - mean(v$lux_m))/sd(v$lux_m), 
       D = v$depth_m,
       S = as.integer(as.factor(v$OTU))
       )

# fit the model without interactions among environmental variables
m1 <- 
  ulam(
  alist(M ~ dnorm(u, sigma),
        
        u <- aS[S] + b_YS[S]*Y + b_TS[S]*TMP + b_LS[S]*LUX + b_DS[S]*D,
        
        b_YS[S] ~ normal( 0 , 1 ),
        b_TS[S] ~ normal( 0 , 1 ),
        b_LS[S] ~ normal( 0 , 1 ),
        b_DS[S] ~ normal( 0, 1 ),
        
        aS[S] ~ normal( 0 , 2 ),
        
        sigma ~ dexp(1)
        
        ),
  data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m1 , depth=2 )

# check the traceplots
traceplot(m1)

# plot the observed data versus the modelled data
post <- sim(m1)
dev.off()
plot(sqrt(v$M), apply(post, 2, mean))
abline(a = 0, b = 1)

# fit the model with two additional interaction terms
m2 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS[S] + b_YS[S]*Y + b_TS[S]*TMP + b_LS[S]*LUX + b_DS[S]*D + b_YT*Y*TMP + b_YD*Y*D,
          
          b_YS[S] ~ normal( 0 , 1 ),
          b_TS[S] ~ normal( 0 , 1 ),
          b_LS[S] ~ normal( 0 , 1 ),
          b_DS[S] ~ normal( 0, 1 ),
          b_YT ~ normal(0, 1),
          b_YD ~ normal(0, 1),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m2, depth = 2 )

# check the traceplots
traceplot(m2)

# plot the observed data versus the modelled data
post <- sim(m2)
plot((v$M - mean(v$M))/sd(v$M), apply(post, 2, mean))
abline(a = 0, b = 1)

# fit the model without mixture data
m3 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS[S] + b_TS[S]*TMP + b_LS[S]*LUX + b_DS[S]*D,
          
          b_TS[S] ~ normal( 0 , 1 ),
          b_LS[S] ~ normal( 0 , 1 ),
          b_DS[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m3, depth = 2 )

# check the traceplots
traceplot(m3)

# plot the observed data versus the modelled data
post <- sim(m3)
plot((v$M - mean(v$M))/sd(v$M), apply(post, 2, mean))
abline(a = 0, b = 1)


# temperature and depth only with full interactions
m4 <- 
  ulam(
    alist(M ~ dnorm(u, sigma),
          
          u <- aS[S] + b_YS[S]*Y + b_TS[S]*TMP + b_DS[S]*D + b_YTS[S]*Y*TMP + b_YDS[S]*Y*D + b_TDS[S]*TMP*D,
          
          b_YTS[S] ~ normal( 0 , 1 ),
          b_YDS[S] ~ normal( 0 , 1 ),
          b_TDS[S] ~ normal( 0 , 1 ),
          b_YS[S] ~ normal( 0 , 1 ),
          b_TS[S] ~ normal( 0 , 1 ),
          b_DS[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m4, depth = 2 )

# check the traceplots
traceplot(m4)

# plot the observed data versus the modelled data
post <- sim(m4)
plot((v$M - mean(v$M))/sd(v$M), apply(post, 2, mean))
abline(a = 0, b = 1)


# temperature and depth only with full interactions using a student's t-distribution
m5 <- 
  ulam(
    alist(M ~ dstudent(2 , u, sigma),
          
          u <- aS[S] + b_YS[S]*Y + b_TS[S]*TMP + b_DS[S]*D + b_YTS[S]*Y*TMP + b_YDS[S]*Y*D + b_TDS[S]*TMP*D,
          
          b_YTS[S] ~ normal( 0 , 1 ),
          b_YDS[S] ~ normal( 0 , 1 ),
          b_TDS[S] ~ normal( 0 , 1 ),
          b_YS[S] ~ normal( 0 , 1 ),
          b_TS[S] ~ normal( 0 , 1 ),
          b_DS[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m5, depth = 2 )

# check the traceplots
traceplot(m5)
dev.off()

# plot the observed data versus the modeled data
post <- sim(m5)
plot((v$M - mean(v$M))/sd(v$M), apply(post, 2, mean))
abline(a = 0, b = 1)

# compare the models using PSIS: m5 is the best by far
compare(m1, m2, m3, m4, m5, func = PSIS)

# plot a proper fit to sample plot
mu <- apply(post, 2, mean)
PI <- apply(post, 2, PI, 0.90)
PI_low <- apply(PI, 2, function(x) x[1])
PI_high <- apply(PI, 2, function(x) x[2])

# pull this information into a data.frame
m5_df <- tibble(species = v$OTU,
                obs = (v$M - mean(v$M))/sd(v$M),
                mu = mu,
                PI_low = PI_low,
                PI_high = PI_high)
head(m5_df)

ggplot(data = m5_df,
       mapping = aes(x = obs, y = mu, colour = species)) +
  geom_point(size = 1.5, shape = 1) +
  geom_errorbar(mapping = aes(x = obs, ymin = PI_low, ymax = PI_high)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_meta()


