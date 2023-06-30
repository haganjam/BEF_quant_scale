#'
#' @title: Simulate metacommunities to test the analytical pipeline
#' 
#' @description: Analyse accuracy of assuming different RYEs
#' 
#' @details: This script analyses the simulated data to test how often the
#' estimates of the biodiversity that incorporate RYE capture the true, 
#' observed biodiversity effect in the simulated metacommunities.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(dplyr)
library(readr)
library(ggplot2)

# load custom plotting theme
source("scripts/Function_plotting_theme.R")

# load the data
output_df <- readRDS("results/mc_sim_test.rds")

# get the subset of values that are affected by the RYEs
output_rye <- 
  output_df %>%
  filter( !(Beff %in% c( "LS", "LC", "TC", "NO" )) )

# add a variable describing the different levels of variation
output_rye$env_var <- rep(c("Spatial env. variation", "Temporal env. variation", "Combination"), each = 10*7)

# get a table of the range of the different BEF effects
output_rye %>%
  group_by(env_var, Beff) %>%
  summarise(mean_BEF = mean(BEF_obs),
            min_BEF = min(BEF_obs),
            max_BEF = max(BEF_obs)) %>%
  filter(Beff %in% c("NBE", "SI", "TI"))

# does the observed value lie in the 90% percentile interval
output_rye <- 
  output_rye %>%
  mutate(within = ifelse(BEF_obs <= PI95_high & BEF_obs >= PI95_low, 1, 0))

# check proportion within the interval
prop_df <- 
  output_rye %>%
  group_by(Beff) %>%
  summarise(prop_within = sum(within)/n())
print(prop_df)

# calculate the absolute prediction error (APE)
# calculate the percentage prediction error (PPE)
output_rye <- 
  output_rye %>%
  mutate(APE = abs(BEF_obs - mean_BEF),
         PPE = (abs(BEF_obs-mean_BEF)/BEF_obs)*100)
  
# get the absolute prediction error
# get the percentage prediction error
error_df <- 
  output_rye %>%
  group_by(Beff) %>%
  summarise(APE_m = median(APE, na.rm = TRUE),
            PPE_med = median(PPE, na.rm = TRUE))
print(error_df)

# bind the prop and absolute deviation dfs
RYE_df <- full_join(prop_df, abs_dev_df, by = "Beff")
names(RYE_df) <- c("BEF effect", "Prop. within", "Mean absolute deviation", "PI95 low", "PI95 high")

# export the table
write_csv(x = RYE_df, file = "figures/SI2_table1.csv")

# check the correlations
cor_obs <- 
  output_rye %>%
  group_by(Beff) %>%
  summarise(cor_x = cor(BEF_obs, mean_BEF))
print(cor_obs)

# plot out the relationship
plot_rye <- 
  output_rye |> 
  dplyr::mutate(within = factor(within)) |>
  dplyr::rename(`Within PI95%` = within)
levels(plot_rye$`Within PI95%`) <- c("Yes", "No")

# visualisation of the accuracy of the method
p1<- 
  ggplot(data = plot_rye) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_errorbar(mapping = aes(x = BEF_obs, 
                              ymin = PI95_low, ymax = PI95_high,
                              colour = `Within PI95%`),
                width = 0, alpha = 0.2, show.legend = FALSE) +
  geom_point(mapping = aes(x = BEF_obs, y = mean_BEF, colour = `Within PI95%`), shape = 1,
             alpha = 0.7) +
  facet_wrap(~Beff, scales = "free") +
  scale_colour_manual(values = wesanderson::wes_palette(name = "Darjeeling1", n = 2)) +
  ylab("Estimated biodiversity effect (with 95% PI)") +
  xlab("True, simulated biodiversity effect") +
  guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1,
                                                   stroke = 1))) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA),
        strip.background = element_rect(fill="white"))
plot(p1)

ggsave(filename = "figures/SI2_fig_S5.svg", p1,
       unit = "cm", width = 21, height = 21)

### END
