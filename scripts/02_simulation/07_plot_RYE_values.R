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
ggsave(filename = "manuscript/figures/app_3_fig_s7.png", p1, dpi = 600,
       unit = "cm", width = 20, height = 14)

# check the range of RYEs for alpha = 3
dd_df %>%
  group_by(alpha) %>%
  summarise(min = min(RYe),
            max = max(RYe))

# plot a single one of these with alpha = 3
p2 <- 
  ggplot(data = dd_df |> dplyr::filter(alpha_3 == "yes"),
       mapping = aes(x = Species, y = RYe, group = group)) +
  geom_line(alpha = 0.1) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Species") +
  ylab(bquote("RY"[E])) +
  theme_meta() +
  theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
)

# export the figure for further modification
ggsave(filename = "manuscript/figures/app_3_fig_s6_part.png", p2, dpi = 600,
       unit = "cm", width = 10, height = 8, bg = "transparent")

### END
