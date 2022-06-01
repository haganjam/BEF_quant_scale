
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
         stvtall, stvgran, stvbjork, stvek, stvasp, stvbok,
         dbiomass02) %>%
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
            stvtall = mean(stvtall, na.rm = TRUE), 
            stvgran = mean(stvgran, na.rm = TRUE), 
            stvbjork = mean(stvbjork, na.rm = TRUE), 
            stvek = mean(stvek, na.rm = TRUE),
            stvasp = mean(stvasp, na.rm = TRUE), 
            stvbok = mean(stvbok, na.rm = TRUE),
            dbiomass02 = mean(dbiomass02, na.rm = TRUE), .groups = "drop") %>%
  rename(place = trakt,
         scots_pine = stvtall, 
         norway_spruce = stvgran,
         birch = stvbjork,
         oak = stvek,
         aspen = stvasp,
         beech = stvbok)

# get complete cases
swe_dat <- swe_dat[complete.cases(swe_dat), ]

# check how much these six species make up in general
swe_dat %>%
  mutate(prop_six_species = (scots_pine + norway_spruce + birch + oak + aspen + beech)/dbiomass02 ) %>%
  pull(prop_six_species) %>%
  summary()

# subset the rows where these six species make up between 0.9 and 1
swe_dat <- 
  swe_dat %>%
  mutate(prop_six_species = (scots_pine + norway_spruce + birch + oak + aspen + beech)/dbiomass02 )

# check how many monocultures there are
swe_dat$SR <- 
  apply(swe_dat[, c("scots_pine", "norway_spruce", "birch", "oak", "aspen", "beech")], 1, function(x) {
  sum(x > 0) } )
sum(swe_dat$SR == 1)

swe_dat %>%
  filter(SR == 1) %>%
  View()

# subset the data to only include plots where the six species make up most of the productivity
swe_dat <- 
  swe_dat %>%
  filter(prop_six_species > 0.9, prop_six_species <= 1) %>%
  select(-prop_six_species)
sum(swe_dat$SR == 1)

# exclude negative values
swe_dat %>%
  filter()

swe_dat <- 
  swe_dat[apply(swe_dat[, c("scots_pine", "norway_spruce", "birch", "oak", "aspen", "beech")], 1, function(x) {
  all(x >= 0)
} ), ]

# get the monoculture data
swe_mono <- 
  swe_dat %>%
  filter(SR == 1) %>%
  select(-dbiomass02, - SR)

swe_mono <- 
  swe_mono %>%
  pivot_longer(cols = c("scots_pine", "norway_spruce", "birch", "oak", "aspen", "beech"),
               names_to = "species",
               values_to = "M")

# load the rethinking package
library(rethinking)

# fit the ulam() model
M_mono <- list(
  M = swe_mono$M ,
  SM = standardize(swe_mono$SoilMoist.cont),
  H = standardize(swe_mono$humid),
  TE = standardize(swe_mono$temp),
  AGE = standardize(swe_mono$stand_age),
  S = as.integer(factor(swe_mono$species)))
str(M_mono)

# model the full model with E and Y interaction
m1 <- ulam(
  alist(
    M ~ dnorm( u, sd ),
    u <- aS[S] + b_SM[S]*SM + b_H[S]*H + b_TE[S]*TE + b_AGE[S]*AGE,
    
    # priors
    aS[S] ~ dnorm(0, 2),
    b_SM[S] ~ normal(0, 2),
    b_H[S] ~ normal(0, 2),
    b_TE[S] ~ normal(0, 2),
    b_AGE[S] ~ normal(0, 2),
    sd ~ dexp(0.5)
    
  ), data=M_mono, chains = 4, iter = 2000, log_lik = FALSE)
traceplot(m1)
dev.off()
precis(m1, depth = 2)

# check the prediction value
pred.m1 <- sim(fit = m1, M_mono[-1])
apply(pred.m1, 2, mean)

plot(M_mono$M, apply(pred.m1, 2, mean))
abline(a = 0, b = 1)





