
# plot the fit-to-sample monoculture models

# load the required libraries
library(dplyr)
library(here)
library(ggpubr)
library(ggplot2)

# load relevant scripts
source(here("scripts/Function_plotting_theme.R"))

# load the plots for the different species
figS9a <- readRDS(here("results/SP_A_monoculture_plot.rds"))
figS9b <- readRDS(here("results/SP_B_monoculture_plot.rds"))
figS9c <- readRDS(here("results/SP_C_monoculture_plot.rds"))
figS9d <- readRDS(here("results/SP_D_monoculture_plot.rds"))
figS9e <- readRDS(here("results/SP_E_monoculture_plot.rds"))

# arrange these plots
figS9 <- ggarrange(figS9a, figS9b, figS9c, figS9d, figS9e,
                   nrow = 2, ncol = 3,
                   labels = letters[1:5],
                   font.label = list(size = 11, face = "plain")
                   )

# save this plot
ggsave(filename = here("figures/figS9.png"), figS9,
       unit = "cm", width = 21, height = 16)

### END