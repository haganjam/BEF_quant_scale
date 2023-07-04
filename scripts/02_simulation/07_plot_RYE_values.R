#'
#' @title Plot RYes from different parameterisations of the Dirichlet distribution
#' 
#' @description This script makes a plot showing how variation in the alpha
#' parameters of the Dirichlet distribution changes the spread of RYe values
#' 

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load plotting theme
source("scripts/Function_plotting_theme.R")

# simulate N RYes with different alpha values
N <- 200

# set the different alpha values
a <- c(1:6)

# set the number of species to simulate
species <- 3

# loop over the different alpha values
dd_list <- vector("list", length = length(a))

for(i in a) {
  
  # get samples from the dirichlet distribution
  dd <- gtools::rdirichlet(n = N, rep(a[i], species))
  dd <- as_tibble(dd)
  names(dd) <- c("1", "2", "3")
  dd <- bind_cols(tibble(rep = 1:N, alpha = paste0("Alpha = ", a[i]) ), dd)
  
  # pull into the long format
  dd <- 
    dd %>%
    pivot_longer(cols = c("1", "2", "3"),
                 names_to = "Species",
                 values_to = "RYe")
  
  dd_list[[i]] <- dd
  
}

# bind into a data.frame
dd_df <- bind_rows(dd_list)

# make a colour variable to highlight alpha = 3
dd_df <- 
  dd_df %>%
  mutate(alpha_3 = ifelse(alpha == "Alpha = 3", "yes", "no"))

# make a grouping variable
dd_df <- 
  dd_df %>%
  mutate(group = paste(rep, alpha, sep = "_"))

p1 <- 
  ggplot(data = dd_df,
       mapping = aes(x = Species, y = RYe, group = group, colour = alpha_3)) +
  geom_line(alpha = 0.1) +
  facet_wrap(~alpha) +
  scale_colour_manual(values = c("grey", "red")) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Species") +
  ylab(bquote("RY"[E])) +
  theme_meta() +
  theme(legend.position = "none", 
        strip.background = element_blank(),
        strip.text = element_text(size = 11))
plot(p1)

# export as a .svg file
ggsave(filename = "figures/SI2_fig_S2.svg", p1,
       unit = "cm", width = 20, height = 14)

# check the range of RYEs for alpha = 3
dd_df %>%
  group_by(alpha) %>%
  summarise(min = min(RYe),
            max = max(RYe))

### END
