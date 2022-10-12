
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
library(here)

# load plotting function
source(here("scripts/Function_plotting_theme.R"))

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

# add laser id to the panel_id
pan_dat$panel_id <- paste(pan_dat$panel_id, " (", pan_dat$laser_id, ")", sep = "")

# plot a grid for each site_id
lapply( split(pan_dat, pan_dat$site_id), function(z) {
  
  x <- z[, c(3, 5, 6, 7)]
  
  x$row <- rep(c(3, 2, 1), each = 4)
  x$col <- rep(c(1, 2, 3, 4), 3)
  
  # add a column for species
  x$species <- "species =                         "
  x$species <- ifelse(x$time_point != "NA", 
                      "", x$species)
  
  # replace NAs with blank data
  x$time_point <- ifelse(x$time_point == "NA", "time =        ", 
                         paste("time =", x$time_point) ) 
  
  # plot the grid
  p1 <- 
    ggplot(data = x,
           mapping = aes(x = col, y = row, fill = treatment)) +
    geom_tile(colour = "black") +
    geom_label(mapping = aes(label = panel_id), 
               label.padding = unit(0.25, "lines"), label.size = NA,
               vjust = -2, size = 5.5) +
    geom_label(mapping = aes(label = time_point), 
               label.padding = unit(0.25, "lines"), label.size = NA,
               vjust = -0.25, size = 5.5) +
    geom_label(mapping = aes(label = treatment), 
               label.padding = unit(0.25, "lines"), label.size = NA,
               vjust = 1.5, size = 5.5) +
    geom_label(mapping = aes(label = species), 
               label.padding = unit(0.25, "lines"), label.size = NA,
               vjust = 3.25, size = 5.5) +
    scale_fill_manual(values = c("grey", "white")) +
    theme(legend.position = "none") +
    ylab(NULL) +
    xlab(NULL) +
    ggtitle(label = paste("site_id",  x$site_id[1], sep = " - ") ) +
    theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    #theme_meta() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank() )
  
  # export the grid to a .csv file
  ggsave(filename = here(paste("scripts/03_empirical_analysis/02_benthic_communities_tjarno/01_preparation/template_", 
                               x$site_id[1], ".png", sep = "")),
         p1, width = 29, height = 20, units = "cm", dpi = 300)
  
} )

### END
