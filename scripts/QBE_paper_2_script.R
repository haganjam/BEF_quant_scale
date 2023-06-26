
# QBE paper 2 script:

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load plotting theme
source("scripts/Function_plotting_theme.R")
source("scripts/03_empirical_analysis/helper_functions.R")

# estimates of the different BEF effects
BEF_dat <- readRDS("results/BEF_effects_case_study_2.rds")

# for all 1000 monoculture estimates and all 100 random RYe values,
# we get a BEF effect for each cluster (i.e. 100 000 BEF effects) for each cluster
head(BEF_dat)
dim(BEF_dat)
summary(BEF_dat)

# load in the environmental dispersion data i.e. multivariate dispersion
# for the different clusters
env_disp <- readRDS(file = "results/benthic_env_dispersion.rds")
head(env_disp)

# add the environmental dispersion data to the BEF effects data
BEF_dat <- full_join(env_disp, BEF_dat, by = "cluster_id")
head(BEF_dat)

# remove cluster F because it only has two usable sites
BEF_dat <- 
  BEF_dat %>%
  filter(cluster_id != "F")

# remove the LS and LC effects because these don't really add much information
BEF_dat <- 
  BEF_dat %>%
  filter( !(Beff %in% c("LC", "LS")) )

# refactor the Beff variable to get rid of LS and LC
BEF_dat$Beff <- factor(BEF_dat$Beff)

# what is the magnitude of different biodiversity effects?

# calculate the summary dataset per cluster
BEF_sum <- 
  BEF_dat %>%
  group_by(cluster_id, Beff) %>%
  summarise(Value_m = mean(Value),
            Value_sd = sd(Value),
            n = n(),
            HPDI_low = HPDI(Value, 0.95)[1],
            HPDI_high = HPDI(Value, 0.95)[2], .groups = "drop")

# check the mean and highest posterior density interval for each effect
print(BEF_sum)


# example analysis using all the different estimated biodiversity effects
# we regress environmental heterogeneity on the spatial insurance effect
# across clusters for all 100 000 estimates of the spatial insurance effect

# this means that we have 100 000 regression lines

# we can then look at, for example, the distribution of regression slopes
# and try to make some inference with it


# do insurance effects depend on environmental heterogeneity?

# get a list of the mean and PI
BEF_SI <- 
  full_join(env_disp, 
            filter(BEF_sum, Beff == "SI"), by = "cluster_id")

# convert the different estimates into a list
SI_list <- 
  BEF_dat %>%
  filter(Beff == "SI") %>%
  mutate(Value_id = paste(mono_rep, RYE, sep = "_")) %>%
  split(., .$Value_id)

# check the split
SI_list[[1]]

# run a regression looped over all the different estimates
lm_list <- vector("list", length = length(SI_list))
for(i in 1:length(SI_list)) {
  
  lm_x <- lm(Value ~ field_dispersion, data =  SI_list[[i]])
  lm_est <- coef(lm_x)
  names(lm_est) <- NULL
  df <- tibble(intercept = lm_est[1],
               slope = lm_est[2])
  
  lm_list[[i]] <- df
  
}

# bind into a data.frame
lm_df <- bind_rows(lm_list, .id = "SI_rep")

# check the distribution of slopes
hist(lm_df$slope)

# check the mean regression slope
mean(lm_df$slope)

### END
