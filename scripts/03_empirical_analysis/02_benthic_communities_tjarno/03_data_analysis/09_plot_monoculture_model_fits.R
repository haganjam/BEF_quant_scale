#'
#' @title: Plot fit to sample of monoculture models
#' 
#' @description: This script plots the fit of the models used to impute the
#' missing monoculture data of the benthic communities.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(dplyr)
library(here)
library(ggpubr)
library(ggplot2)

# load relevant scripts
source(here("scripts/Function_plotting_theme.R"))

# load the plots for the different species
figS5a <- readRDS(here("results/SP_A_monoculture_plot.rds"))
figS5b <- readRDS(here("results/SP_B_monoculture_plot.rds"))
figS5c <- readRDS(here("results/SP_C_monoculture_plot.rds"))
figS5d <- readRDS(here("results/SP_D_monoculture_plot.rds"))
figS5e <- readRDS(here("results/SP_E_monoculture_plot.rds"))

# arrange these plots
figS5 <- ggarrange(figS5a, figS5b, figS5c, figS5d, figS5e,
                   nrow = 2, ncol = 3,
                   labels = letters[1:5],
                   font.label = list(size = 11, face = "plain")
                   )

# save this plot
ggsave(filename = here("figures/figA1_S5.png"), figS5,
       unit = "cm", width = 21, height = 16)

### END
