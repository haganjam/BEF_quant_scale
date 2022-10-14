#'
#' @title: Model the missing monocultures
#' 
#' @description: In this script, we link the clean biomass data with the clean
#' environmental data and try to find models that fit the monoculture data well.
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

# check the correlation among these variables
# pairs(bio_env[,-c(1:3, 14:16)])

# let's try build a model with three variables:
# 1. mixture biomass
# 2. temp_C_m
# 3. lux_C_m
# 4. seabed distance

# load the rethinking package
library(rethinking)

# make a data.list with the training data
train.d <- bio_env[complete.cases(bio_env), ]
head(train.d)
dim(train.d)

# examine the training data
train.d %>%
  group_by(OTU) %>%
  summarise(max_M = max(M))

# subset barnacles and sea-squirts
barn_sq <- FALSE

# if barnacle then subset the barnacle data
if( barn ) {
  
  v <- train.d[train.d$OTU %in% c("Barn", "Seasq"),]
  
} else {
  
  v <- train.d[!(train.d$OTU %in% c("Barn", "Seasq")),]
  
}

train.d <- 
  list(M = v$M,
       M1 = v$M + min(v$M[v$M > 0]),
       Y = sqrt(v$Y),
       T = (v$temp_C_m - mean(v$temp_C_m))/sd(v$temp_C_m),
       L = (v$lux_m - mean(v$lux_m))/sd(v$lux_m),
       B = (v$seabed_dist_m-mean(v$seabed_dist_m))/sd(v$seabed_dist_m),
       S = as.integer(as.factor(v$OTU))
       )

# model1
m1 <- 
  ulam(
  alist(M1 ~ dlnorm(u, sigma),
        
        u <- aS[S] + bY[S]*Y + bT[S]*T + bL[S]*L + bB[S]*B + bTL[S]*T*L + bTB[S]*T*B + bLB[S]*L*B,
        
        bY[S] ~ normal( 0 , 1 ),
        bT[S] ~ normal( 0 , 1 ),
        bL[S] ~ normal( 0 , 1 ),
        bB[S] ~ normal( 0, 1 ),
        bTL[S] ~ normal( 0, 1 ),
        bTB[S] ~ normal( 0, 1 ),
        bLB[S] ~ normal( 0, 1 ),
        
        aS[S] ~ normal( 0 , 2 ),
        
        sigma ~ dexp(1)
        
        ),
  data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m1 , depth=2 )

# check the traceplots
traceplot(m1)
dev.off()

# plot the observed data versus the modelled data
post <- sim(m1)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)

# reduced model2
m2 <- 
  ulam(
    alist(M1 ~ dlnorm(u, sigma),
          
          u <- aS[S] + bY[S]*Y + bT[S]*T + bL[S]*L +bB[S]*B + bTL[S]*T*L + bTB[S]*T*B,
          
          bY[S] ~ normal( 0 , 1 ),
          bT[S] ~ normal( 0 , 1 ),
          bL[S] ~ normal( 0 , 1 ),
          bB[S] ~ normal( 0, 1 ),
          bTL[S] ~ normal( 0, 1 ),
          bTB[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m2 , depth=2 )

# check the traceplots
traceplot(m2)
dev.off()

# plot the observed data versus the modelled data
post <- sim(m2)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)

# reduced model3
m3 <- 
  ulam(
    alist(M1 ~ dlnorm(u, sigma),
          
          u <- aS[S] + bY[S]*Y + bT[S]*T + bL[S]*L + bB[S]*B + bTB[S]*T*B + bLB[S]*L*B,
          
          bY[S] ~ normal( 0 , 1 ),
          bT[S] ~ normal( 0 , 1 ),
          bL[S] ~ normal( 0 , 1 ),
          bB[S] ~ normal( 0, 1 ),
          bTB[S] ~ normal( 0, 1 ),
          bLB[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m3 , depth=2 )

# check the traceplots
traceplot(m3)
dev.off()

# plot the observed data versus the modeled data
post <- sim(m3)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)


# reduced model4
m4 <- 
  ulam(
    alist(M1 ~ dlnorm(u, sigma),
          
          u <- aS[S] + bY[S]*Y + bT[S]*T + bL[S]*L + bB[S]*B + bTL[S]*T*L + bLB[S]*L*B,
          
          bY[S] ~ normal( 0 , 1 ),
          bT[S] ~ normal( 0 , 1 ),
          bL[S] ~ normal( 0 , 1 ),
          bB[S] ~ normal( 0, 1 ),
          bTL[S] ~ normal( 0, 1 ),
          bLB[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m4 , depth=2 )

# check the traceplots
traceplot(m4)
dev.off()

# plot the observed data versus the modelled data
post <- sim(m4)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)


# reduced model5
m5 <- 
  ulam(
    alist(M1 ~ dlnorm(u, sigma),
          
          u <- aS[S] + bY[S]*Y + bT[S]*T + bL[S]*L + bB[S]*B + bLB[S]*L*B,
          
          bY[S] ~ normal( 0 , 1 ),
          bT[S] ~ normal( 0 , 1 ),
          bL[S] ~ normal( 0 , 1 ),
          bB[S] ~ normal( 0, 1 ),
          bLB[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m5 , depth=2 )

# check the traceplots
traceplot(m5)
dev.off()

# plot the observed data versus the modeled data
post <- sim(m5)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)


# reduced model6
m6 <- 
  ulam(
    alist(M1 ~ dlnorm(u, sigma),
          
          u <- aS[S] + bY[S]*Y + bT[S]*T + bL[S]*L + bB[S]*B + bTB[S]*T*B,
          
          bY[S] ~ normal( 0 , 1 ),
          bT[S] ~ normal( 0 , 1 ),
          bL[S] ~ normal( 0 , 1 ),
          bB[S] ~ normal( 0, 1 ),
          bTB[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m6 , depth=2 )

# check the traceplots
traceplot(m6)
dev.off()

# plot the observed data versus the modelled data
post <- sim(m6)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)


# reduced model7
m7 <- 
  ulam(
    alist(M1 ~ dlnorm(u, sigma),
          
          u <- aS[S] + bY[S]*Y + bT[S]*T + bL[S]*L + bB[S]*B + bTL[S]*T*L,
          
          bY[S] ~ normal( 0 , 1 ),
          bT[S] ~ normal( 0 , 1 ),
          bL[S] ~ normal( 0 , 1 ),
          bB[S] ~ normal( 0, 1 ),
          bTL[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m7 , depth=2 )

# check the traceplots
traceplot(m7)
dev.off()

# plot the observed data versus the modelled data
post <- sim(m7)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)


# reduced model8
m8 <- 
  ulam(
    alist(M1 ~ dlnorm(u, sigma),
          
          u <- aS[S] + bT[S]*T + bL[S]*L + bB[S]*B + bTL[S]*T*L,
          
          bT[S] ~ normal( 0 , 1 ),
          bL[S] ~ normal( 0 , 1 ),
          bB[S] ~ normal( 0, 1 ),
          bTL[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m8 , depth=2 )

# check the traceplots
traceplot(m8)
dev.off()

# plot the observed data versus the modeled data
post <- sim(m8)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)


# reduced model9
m9 <- 
  ulam(
    alist(M1 ~ dlnorm(u, sigma),
          
          u <- aS[S] + bY[S]*Y + bT[S]*T + bL[S]*L + bB[S]*B,
          
          bY[S] ~ normal( 0 , 1 ),
          bT[S] ~ normal( 0 , 1 ),
          bL[S] ~ normal( 0 , 1 ),
          bB[S] ~ normal( 0, 1 ),
          bTL[S] ~ normal( 0, 1 ),
          bTB[S] ~ normal( 0, 1 ),
          bLB[S] ~ normal( 0, 1 ),
          
          aS[S] ~ normal( 0 , 2 ),
          
          sigma ~ dexp(1)
          
    ),
    data = train.d, chains = 4, log_lik = TRUE)

# check the precis output
precis( m9 , depth=2 )

# check the traceplots
traceplot(m9)
dev.off()

# plot the observed data versus the modelled data
post <- sim(m9)
plot(train.d$M1, apply(post, 2, mean))
abline(a = 0, b = 1)


# compare the models using PSIS: m5 is the best by far
compare(m1, m2, m3, m4, m5, m6, m7, m8, func = PSIS)


# the best model is reduced model 7

# sample that model
post <- sim(fit = m8)

# back transform the data by subtracting the minimum value
post <- apply(post, 2, function(x) x + min(train.d$M[train.d$M > 0]))

# plot a proper fit to sample plot
mu <- apply(post, 2, mean)
PI <- apply(post, 2, PI, 0.90)
PI_low <- apply(PI, 2, function(x) x[1])
PI_high <- apply(PI, 2, function(x) x[2])

# pull this information into a data.frame
m7_df <- tibble(species = v$OTU,
                obs = train.d$M,
                mu = mu,
                PI_low = PI_low,
                PI_high = PI_high)
head(m7_df)

ggplot(data = m7_df,
       mapping = aes(x = obs, y = mu, colour = species)) +
  geom_point(size = 1.5, shape = 1) +
  geom_errorbar(mapping = aes(x = obs, ymin = PI_low, ymax = PI_high)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1.5)) +
  theme_meta()


