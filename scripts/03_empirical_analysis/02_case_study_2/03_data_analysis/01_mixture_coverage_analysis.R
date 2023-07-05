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
library(ggpubr)

# load the plotting theme
source("scripts/Function_plotting_theme.R")

# load the cleaned dataset
cov_dat <- read_csv("data/case_study_2/data_clean/mixture_coverage_data_clean.csv")

# check how many unique buoys there are
length(unique(cov_dat$buoy_id))

# remove the species that we are not using anymore
cov_dat <- 
  cov_dat %>%
  select(-all_of(c("Asci_percent", "Seasq_bumpi_percent", "Hydro_pink_percent")))

# pivot this data longer
cov_dat <- 
  cov_dat %>%
  pivot_longer(cols = ends_with("percent"),
               names_to = "Species",
               values_to = "cover")

# remove the percent in the words
cov_dat$Species <- gsub(pattern = "_percent", replacement = "", x = cov_dat$Species)

# select the relevant columns
cov_dat <- 
  cov_dat %>%
  select(buoy_id, time, Species, cover)

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
covT12 <- full_join(covT1, covT2, by = c("buoy_id", "Species"))

# examine the relationship between covT1 and covT2
plot(covT12$T1_cover, covT12$T2_cover)
abline(a = 0, b = 1)
SP.x <- cor.test(covT12$T1_cover, covT12$T2_cover, method = "spearman")
print(SP.x)

# plot a proper correlation plot

# get a colour palette
col_pal <- wesanderson::wes_palette(name = "Darjeeling1", n = 6, type = "continuous")

# rename the OTUs
covT12$OTU <- factor(covT12$Species)
levels(covT12$OTU) <- c("Barn", "Bryo", "Asci", "Hydro", "Ciona")

p1 <- 
  ggplot(data = covT12,
       mapping = aes(x = T1_cover, y = T2_cover, colour = OTU)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 2, alpha =0.7) +
  scale_colour_manual(values = col_pal[1:5]) +
  xlab("Cover (%)") +
  ylab("Cover (%)") +
  annotate(geom = "text", x = 10, y = 90, label = paste0("rho = ", round(SP.x$estimate, 2))) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 10),
        legend.spacing.x = unit(0.05, 'cm'))
plot(p1)

# what about total cover?
tot_cov <- 
  covT12 %>%
  group_by(buoy_id) %>%
  summarise(T1_tot_cover = sum(T1_cover, na.rm = TRUE),
            T2_tot_cover = sum(T2_cover, na.rm = TRUE))

# what is the correlation?
SP.y <- cor(tot_cov$T1_tot_cover, tot_cov$T2_tot_cover, method = "spearman")
cor.test(tot_cov$T1_tot_cover, tot_cov$T2_tot_cover, method = "spearman")

# plot a proper correlation plot
p2 <- 
  ggplot(data = tot_cov,
       mapping = aes(x = T1_tot_cover, y = T2_tot_cover)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 2, alpha =0.7) +
  scale_colour_viridis_d(option = "C") +
  xlab("Total cover (%)") +
  ylab("Total cover (%)") +
  ggtitle("") +
  annotate(geom = "text", x = 15, y = 160, label = paste0("rho = ", round(SP.y, 2))) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA),
        plot.title = element_text(size = 35))
plot(p2)

p12 <- 
  ggarrange(p1, p2, ncol = 2, nrow = 1,
            labels =c("a", "b"), font.label = list(face = "plain", size = 11))
plot(p12)

ggsave(filename = "figures/SI3_fig_S11.svg", p12,
       unit = "cm", width = 20, height = 10)

# correlation within each buoy
cor_dist <- 
  covT12 %>%
  group_by(buoy_id) %>%
  summarise(Spearman_r = cor(T1_cover, T2_cover, method = "spearman"))

# calculate the mean correlation
mean(cor_dist$Spearman_r, na.rm = TRUE)
sd(cor_dist$Spearman_r, na.rm = TRUE)
hist(cor_dist$Spearman_r)

# make a proper histogram
ggplot(data = cor_dist,
       mapping = aes(x = Spearman_r)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, 1.05)) +
  geom_vline(xintercept = mean(cor_dist$Spearman_r, na.rm = TRUE), 
             linetype = "dashed", colour = "red") +
  theme_meta()

# calculate the average deviation across all species-buoy combinations
cov_diff <- 
  covT12 %>%
  mutate(T12_diff = abs((T1_cover - T2_cover)) )

ggplot(data = cov_diff,
         mapping = aes(x = T12_diff)) +
  geom_histogram() +
  geom_vline(xintercept = mean( cov_diff$T12_diff), 
             linetype = "dashed", colour = "red" ) +
  theme_meta()
  
# overall mean difference between the two panels is only 3.31
x <- cov_diff$T12_diff
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
  theme_meta()

### END
