#'
#' @title: Calculate the BEF effects using the modeled monoculture data
#' 
#' @description: This script summarises and plots the different BEF effects calculated
#' on the benthic marine fouling data using uncertainty in both monoculture functioning
#' and initial relative abundances.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(here)

# load plotting theme
source(here("scripts/Function_plotting_theme.R"))

# load the relevant data
BEF_dat <- readRDS(here("results/benthic_BEF_effects.rds"))
head(BEF_dat)
dim(BEF_dat)
summary(BEF_dat)

# reorder the columns
BEF_dat <- 
  BEF_dat %>%
  select(cluster_id, mono_sample, RA, Beff, Value)

# check how many unique cluster_id, mono_samp and RA there are
pre_trim <- length(unique(paste0(BEF_dat$cluster_id, BEF_dat$mono_sample, BEF_dat$RA)))

# trim the distributions to the 99% quantiles to prevent the few gigantic outliers
head(BEF_dat)
BEF_trim <- 
  BEF_dat %>%
  group_by(Beff) %>%
  mutate(q99 = quantile(Value, 0.99),
         q01 = quantile(Value, 0.01)) %>%
  mutate(within_PI = if_else( (Value < q99) & (Value > q01), 1, 0 )) %>%
  ungroup() %>%
  group_by(cluster_id, mono_sample, RA) %>%
  mutate(within_PI = if_else( any(within_PI != 1), 1, 0 )) %>%
  ungroup() %>%
  filter(within_PI != 1) %>%
  select(-q99, -q01, -within_PI)

# unique cluster_id, mono_samp, RA after trimming
post_trim <- length(unique(paste0(BEF_trim$cluster_id, BEF_trim$mono_sample, BEF_trim$RA)))

# check the difference
((pre_trim - post_trim)/pre_trim)*100

# plot histograms of the effects: A few huge outliers that inflate the SD
ggplot(data = BEF_trim,
       mapping = aes(x = Value)) +
  geom_density(fill = "grey") +
  facet_wrap(~Beff, scales = "free") +
  theme_meta()

ggplot(data = BEF_trim,
       mapping = aes(x = Value)) +
  geom_histogram(fill = "grey") +
  facet_wrap(~Beff, scales = "free") +
  theme_meta()

# check the internal consistency of the net biodiversity effects
BEF_trim %>%
  filter(cluster_id == "A") %>%
  filter( !(Beff %in% c("NBE", "LS", "LC", "TS", "IT")) ) %>%
  group_by(mono_sample, RA) %>%
  summarise(s = sum(Value))

BEF_dat %>%
  filter(cluster_id == "A", Beff == "NBE")

# summarise the data into a pooled estimate
BEF_pool <- 
  BEF_trim %>%
  group_by(cluster_id, Beff) %>%
  summarise(Value_m = round(mean(Value), 2),
            PI_low = rethinking::HPDI(Value, 0.90)[1],
            PI_high = rethinking::HPDI(Value, 0.90)[2], .groups = "drop")

# is the median within the percentile interval 90
BEF_pool %>%
  mutate(within_PI = if_else((Value_m > PI_low) & (Value_m < PI_high), 1, 0 )) %>%
  filter(within_PI != 1) %>%
  View()

# compare total to local selection and complementarity

# which effects to plot
eff_in <- c("LC", "LS", "TC", "TS")

# plot the effects
p1 <- 
  ggplot(data = BEF_pool %>% filter(Beff %in% eff_in) ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_point(mapping = aes(x = Beff, y = Value_m, 
                           colour = Beff, fill = Beff, group = cluster_id),
             position = position_dodge(0.75), size = 2.5) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75)) + 
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2) +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in)) +
  ylab("Effect (g time-1)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p1)

# compare net biodiversity effects, total complementarity and total selection
eff_in <- c("NBE", "TC", "NO", "IT")
p2 <- 
  ggplot(data = BEF_pool %>%
           filter(Beff %in% eff_in) %>%
           mutate(Beff = factor(Beff, levels = eff_in))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_point(mapping = aes(x = Beff, y = Value_m, 
                           colour = Beff, fill = Beff, group = cluster_id),
             position = position_dodge(0.75), size = 2) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75))  +
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2) +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab("Effect (g time-1)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p2)

# examine the distribution of the insurance effects
eff_in <- c("AS", "TI", "SI", "ST")
p3 <- 
  ggplot(data = BEF_pool %>%
           filter(Beff %in% eff_in) %>%
           mutate(Beff = factor(Beff, levels = eff_in))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_point(mapping = aes(x = Beff, y = Value_m, 
                           colour = Beff, fill = Beff, group = cluster_id),
             position = position_dodge(0.75), size = 2) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75)) +
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2) +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab("Effect (g (%) time-1)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p3)

# arrange this plot
p123 <- 
  ggarrange(p1, p2, p3, ncol = 3, nrow = 1,
            labels = c("a", "b", "c"),
            font.label = list(size = 11, face = "plain"))
plot(p123)

# save this plot
ggsave(filename = here("figures/ben_fig1.png"), p123,
       unit = "cm", width = 20, height = 7)








