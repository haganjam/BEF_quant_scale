
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

# environmental data
swe_env <- 
  swe_dat %>%
  select(ar, trakt, palslag, SoilMoisture, ParentMaterial, SoilComp, klimat) %>%
  group_by(ar, trakt) %>%
  summarise(SoilMoisture = median(SoilMoisture, na.rm = TRUE),
            ParentMaterial = median(ParentMaterial),
            SoilComp = median(SoilComp),
            klimat = median(klimat), .groups = "drop") %>%
  mutate(ar = as.integer(factor(ar)),
         trakt = as.integer(factor(trakt))) %>%
  mutate(sample = paste(ar, trakt, sep = "")) %>%
  rename(time = ar, place = trakt) %>%
  select(sample, place, time,
         SoilMoisture,
         ParentMaterial,
         SoilComp,
         klimat)
head(swe_env)  

swe_env <- swe_env[complete.cases(swe_env),]
head(swe_env)

# subset the data
swe_dat <- 
  swe_dat %>%
  select(ar, trakt, palslag, sdtvtall, sdtvgran, sdtvbjork, sdtvek,
         sdtvasp, sdtvbok) %>%
  group_by(ar, trakt) %>%
  summarise( n = n(),
             sdtvtall = mean(sdtvtall), 
             sdtvgran = mean(sdtvgran), 
             sdtvbjork = mean(sdtvbjork), 
             sdtvek = mean(sdtvek),
             sdtvasp = mean(sdtvasp), 
             sdtvbok = mean(sdtvbok), .groups = "drop" ) %>%
  mutate(ar = as.integer(factor(ar)),
         trakt = as.integer(factor(trakt))) %>%
  mutate(sample = paste(ar, trakt, sep = "")) %>%
  rename(time = ar, place = trakt,
         scots_pine = sdtvtall, 
         norway_spruce = sdtvgran,
         birch = sdtvbjork,
         oak = sdtvek,
         aspen = sdtvasp,
         beech = sdtvbok) %>%
  select(sample, place, time,
         scots_pine,
         norway_spruce,
         birch,
         oak,
         aspen,
         beech)
head(swe_dat)

# filter out incomplete cases
swe_dat <- swe_dat[complete.cases(swe_dat), ]
head(swe_dat)

# is there a unique sample for each row?
length(unique(swe_dat$sample)) == nrow(swe_dat)

# calculate species richness
swe_dat$SR <- 
  apply(swe_dat[, c("scots_pine", "norway_spruce", "birch", "oak", "aspen", "beech")], 1, function(x) {
  sum(x > 0)
})

swe_mono <- 
  swe_dat %>%
  filter(SR == 1) %>%
  pivot_longer(cols = c("scots_pine", "norway_spruce", "birch", "oak", "aspen", "beech"),
               names_to = "species",
               values_to = "productivity") %>%
  select(-SR)

swe_mono <- left_join(swe_mono, swe_env , by = c("sample", "place", "time"))
swe_mono

# fit a model to try and predict the monocultures and then try and use them to get the unknown monocultures



