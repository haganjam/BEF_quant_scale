#'
#' @title: Analyse the cover correction data
#' 
#' @description: Script to clean the cover-correction data. This cover-correction
#' data is used to check whether our weights taken at different time points
#' are representative samples of the same communities.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# load the cleaned dataset
cov_dat <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/mixture_coverage_data_clean.csv"))

# pivot this data longer
cov_dat <- 
  cov_dat %>%
  pivot_longer(cols = ends_with("percent"),
               names_to = "species",
               values_to = "cover")

# select the relevant columns
names(cov_dat)
cov_dat <- 
  cov_dat %>%
  select(buoy_id, time, species, cover)

# get the T1 measurements
covT1 <- 
  cov_dat %>%
  filter(time == "T1") %>%
  rename(T1_cover = cover) %>%
  select(-time)

# get the T2 measurements
covT2 <- 
  cov_dat %>%
  filter(time == "T2") %>%
  rename(T2_cover = cover) %>%
  select(-time)

# join these two together
covT12 <- full_join(covT1, covT2, by = c("buoy_id", "species"))

# examine the relationship between covT1 and covT2
plot(covT12$T1_cover, covT12$T2_cover)
cor.test(covT12$T1_cover, covT12$T2_cover, method = "spearman")

# correlation within each buoy
cor_dist <- 
  covT12 %>%
  group_by(buoy_id) %>%
  summarise(Spearman_r = cor(T1_cover, T2_cover))

# calculate the mean correlation
mean(cor_dist$Spearman_r, na.rm = TRUE)
hist(cor_dist$Spearman_r)

# calculate the average deviation across all species-buoy combinations
covT12 %>%
  mutate(T12_diff = abs((T1_cover - T2_cover)) ) %>%
  group_by(buoy_id) %>%
  summarise(T12_diff_m = mean(T12_diff, na.rm = TRUE),
            T12_diff_sd = sd(T12_diff, na.rm = TRUE))

# overall mean difference between the two panels is only 3.31
x <- abs(covT12$T1_cover - covT12$T2_cover)
mean(x)
sd(x)

# what about if we only consider cases where both actually have cover values
y <- (covT12$T1_cover > 0) & (covT12$T2_cover > 0)
z <- abs(covT12$T1_cover[y] - covT12$T2_cover[y])
mean(z)
sd(z)

# what does the distribution look like?
hist(z)

# visualise the relationship between the T1 and T2 cover values at each buoy
ggplot(data = covT12,
       mapping = aes(x = T1_cover, y = T2_cover)) +
  geom_jitter(width = 0.05, shape = 1, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  facet_wrap(~buoy_id) +
  theme_bw()

# what about total cover?
tot_cov <- 
  covT12 %>%
  group_by(buoy_id) %>%
  summarise(T1_tot_cover = sum(T1_cover, na.rm = TRUE),
            T2_tot_cover = sum(T2_cover, na.rm = TRUE))

# what is the correlation?
cor(tot_cov$T1_tot_cover, tot_cov$T2_tot_cover)

# plot the relationship
plot(tot_cov$T1_tot_cover, tot_cov$T2_tot_cover, ylab = "T2_cover", xlab = "T1_cover",
     xlim = c(0, 175), ylim = c(0, 175))
abline(a = 0, b = 1)

### END
