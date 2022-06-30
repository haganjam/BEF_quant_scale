
# Project: Quantifying biodiversity effects across scales in natural ecosystems
# Title: Randomise the position of the mixtures and the monocultures
# Author: James Hagan
# Date: 2022/04/24

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19043)

# attached base packages:
# stats, graphics, grDevice, utils, datasets, methods, base

# load relevant packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(googlesheets4)

# read a Google sheet
pan_dat <- read_sheet("https://docs.google.com/spreadsheets/d/1fUn4k0SuZMr7KVWW1PPDzviemPAKGBNtelKgJl11lXU/edit#gid=0", 
                       sheet = "Sheet1")
head(pan_dat)

# randomise the position of the mixtures in each cluster
clus_lab <- pan_dat[, c(1, 2, 3, 4, 5)]
head(clus_lab)

mix_mono <- pan_dat[, c(6, 7)]
head(mix_mono)

# split the mix_mono by cluster id
mix_mono <- split(mix_mono, clus_lab$site_id)

# randomise the row positions of the mixtures
set.seed(547976)
mix_mono <- lapply(mix_mono, function(x) x[sample(1:nrow(x), nrow(x), replace = FALSE), ] )

# bind this list back together
mix_mono <- dplyr::bind_rows(mix_mono)

# bind this back to cluster lab
pan_dat <- dplyr::bind_cols(clus_lab, mix_mono)
View(pan_dat)


x <- pan_dat[c(1:12), ]
x <- x[, c(3, 5, 6, 7)]
head(x)

x$row <- rep(c(1, 2, 3), each = 4)
x$col <- rep(c(1, 2, 3, 4), 3)

# replace NAs with blank data
x$time_point <- ifelse(x$time_point == "NA", "", x$time_point)

# add a column for species
x$species <- "species =                         "
x$species <- ifelse(x$time_point != "", "", x$species)

ggplot(data = x,
       mapping = aes(x = col, y = row, fill = treatment)) +
  geom_tile(colour = "black") +
  geom_label(mapping = aes(label = panel_id), 
             label.padding = unit(0.25, "lines"), label.size = NA,
             vjust = -1.75) +
  geom_label(mapping = aes(label = time_point), 
             label.padding = unit(0.25, "lines"), label.size = NA,
             vjust = -0.25) +
  geom_label(mapping = aes(label = treatment), 
             label.padding = unit(0.25, "lines"), label.size = NA,
             vjust = 1.25) +
  geom_label(mapping = aes(label = species), 
             label.padding = unit(0.25, "lines"), label.size = NA,
             vjust = 2.75) +
  scale_fill_manual(values = c("grey", "white")) +
  theme(legend.position = "none") +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle(label = paste("site_id",  x$site_id[1], sep = " - ") ) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))



