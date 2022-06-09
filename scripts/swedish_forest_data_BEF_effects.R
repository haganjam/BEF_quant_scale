
# Swedish forest data

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# set script to call partition functions from
source(here("scripts/isbell_2018_partition.R"))

# read in the data
swe_dat <- read_tsv( here("data/Swedish_forest_data.txt") )
summary(swe_dat)
names(swe_dat)

# how many are permanent plots?
length(unique(swe_dat$yttyp))

swe_dat <- 
  swe_dat %>%
  select(region2, 
         year = ar, trakt, palslag, 
         SoilMoist.cont, soil_type = JORDMAN2, parent_mat = JORDART2,
         klimat, humid, temp = tsumma,
         stand_age = bestald, 
         # sdtvtall, sdtvgran, sdtvbjork, sdtvek, sdtvasp, sdtvbok,
         # dbiomass02,
         stvtall, stvgran, stvbjork, stvek, stvasp, stvbok,
         tot.bmass) %>%
  group_by(region2, year, trakt) %>%
  summarise(SoilMoist.cont = mean(SoilMoist.cont, na.rm = TRUE),
            soil_type = median(soil_type, na.rm = TRUE),
            parent_mat = median(parent_mat, na.rm = TRUE),
            klimat = median(klimat, na.rm = TRUE),
            humid = mean(humid, na.rm = TRUE),
            temp = mean(temp, na.rm = TRUE),
            stand_age = mean(stand_age, na.rm = TRUE),
            n = n(),
            #sdtvtall = mean(sdtvtall, na.rm = TRUE), 
            #sdtvgran = mean(sdtvgran, na.rm = TRUE), 
            #sdtvbjork = mean(sdtvbjork, na.rm = TRUE), 
            #sdtvek = mean(sdtvek, na.rm = TRUE),
            #sdtvasp = mean(sdtvasp, na.rm = TRUE), 
            #sdtvbok = mean(sdtvbok, na.rm = TRUE),
            # dbiomass02 = mean(dbiomass02, na.rm = TRUE)
            stvtall = mean(stvtall, na.rm = TRUE), 
            stvgran = mean(stvgran, na.rm = TRUE), 
            stvbjork = mean(stvbjork, na.rm = TRUE), 
            stvek = mean(stvek, na.rm = TRUE),
            stvasp = mean(stvasp, na.rm = TRUE), 
            stvbok = mean(stvbok, na.rm = TRUE),
            tot.bmass = mean(tot.bmass, na.rm = TRUE), .groups = "drop") %>%
  rename(place = trakt,
         scots_pine = stvtall, 
         norway_spruce = stvgran,
         birch = stvbjork,
         oak = stvek,
         aspen = stvasp,
         beech = stvbok)

# get complete cases
sum(complete.cases(swe_dat))
nrow(swe_dat)
swe_dat <- swe_dat[complete.cases(swe_dat), ]

# check how much these six species make up in general
swe_dat %>%
  mutate(prop_six_species = (scots_pine + norway_spruce + birch + oak + aspen + beech)/tot.bmass ) %>%
  pull(prop_six_species) %>%
  summary()

# subset the rows where these six species make up between 0.9 and 1
swe_dat <- 
  swe_dat %>%
  mutate(prop_six_species = (scots_pine + norway_spruce + birch + oak + aspen + beech)/tot.bmass ) %>%
  filter(prop_six_species > 0.7, prop_six_species < 1)
nrow(swe_dat)
summary(swe_dat)

# get complete cases only
nrow(swe_dat)
swe_dat <- swe_dat[complete.cases(swe_dat), ]
nrow(swe_dat)
summary(swe_dat)
hist(swe_dat$prop_six_species)

# load the rethinking package
library(rethinking)

# fit the ulam() model
M_mono <- list(
  SP_foc = swe_dat$aspen,
  SP1 = swe_dat$scots_pine/max(swe_dat$scots_pine),
  SP2 = swe_dat$norway_spruce/max(swe_dat$norway_spruce),
  SP3 = swe_dat$birch/max(swe_dat$birch),
  SP4 = swe_dat$oak/max(swe_dat$oak),
  SP5 = swe_dat$beech/max(swe_dat$beech),
  SM = standardize(swe_dat$SoilMoist.cont),
  H = standardize(swe_dat$humid),
  TE = standardize(swe_dat$temp),
  AGE = standardize(swe_dat$stand_age))
str(M_mono)
hist(M_mono$SP_foc)
hist(M_mono$SP5)
M_mono$SP_foc[M_mono$SP_foc> 0] %>% length()

# model the full model with E and Y interaction
m1 <- ulam(
  alist(
    SP_foc ~ dnorm( u, sd ),
    u <- a + b_sm*SM + b_H*H + b_TE*TE + b_AGE*AGE + b_s1*SP1 + b_s2*SP2 + b_s3*SP3 + b_s4*SP4 + b_s5*SP5,
    
    # priors
    a ~ dnorm(0, 2),
    b_sm ~ normal(0, 2),
    b_H ~ normal(0, 2),
    b_TE ~ normal(0, 2),
    b_AGE ~ normal(0, 2),
    b_s1 ~ normal(0, 2),
    b_s2 ~ normal(0, 2),
    b_s3 ~ normal(0, 2),
    b_s4 ~ normal(0, 2),
    b_s5 ~ normal(0, 2),
    sd ~ dexp(0.5)
    
  ), data=M_mono, chains = 4, iter = 2000, log_lik = FALSE)
traceplot(m1)
dev.off()
precis(m1, depth = 2)

# check the prediction value
pred.m1 <- sim(fit = m1, M_mono[-1])
# apply(pred.m1, 2, mean)

plot(M_mono$SP_foc, apply(pred.m1, 2, mean))
abline(a = 0, b = 1)





