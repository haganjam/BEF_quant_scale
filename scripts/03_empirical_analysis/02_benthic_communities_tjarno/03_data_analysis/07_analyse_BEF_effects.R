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
library(tidyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(betapart)
library(vegan)
library(here)
library(corrplot)

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

# check how many of each effect there are
BEF_dat %>%
  group_by(cluster_id, Beff) %>%
  summarise(n = n())

# check how many unique cluster_id, mono_samp and RA there are
pre_trim <- length(unique(paste0(BEF_dat$cluster_id, BEF_dat$mono_sample, BEF_dat$RA)))

# trim the distributions to the 99% quantiles to prevent the few gigantic outliers
head(BEF_dat)
BEF_trim <- 
  BEF_dat %>%
  group_by(cluster_id, Beff) %>%
  mutate(q99 = quantile(Value, 0.99),
         q01 = quantile(Value, 0.01)) %>%
  mutate(within_PI = if_else( (Value < q99) & (Value > q01), 1, 0 )) %>%
  ungroup() %>%
  group_by(cluster_id, mono_sample, RA) %>%
  mutate(within_PI = if_else( any(within_PI != 1), 1, 0 )) %>%
  ungroup() %>%
  filter(within_PI != 1) %>%
  select(-q99, -q01, -within_PI)

BEF_trim %>%
  group_by(cluster_id, Beff) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(n) %>%
  mean()

# unique cluster_id, mono_samp, RA after trimming
post_trim <- length(unique(paste0(BEF_trim$cluster_id, BEF_trim$mono_sample, BEF_trim$RA)))

# check the difference
((pre_trim - post_trim)/pre_trim)*100

# check the internal consistency of the net biodiversity effects
BEF_trim %>%
  filter(cluster_id == "A") %>%
  filter( !(Beff %in% c("NBE", "LS", "LC", "TS", "IT")) ) %>%
  group_by(mono_sample, RA) %>%
  summarise(s = sum(Value))

BEF_dat %>%
  filter(cluster_id == "A", Beff == "NBE")

# what is the magnitude of different biodiversity effects?

# summarise the data into a pooled estimate
BEF_pool <- 
  BEF_trim %>%
  group_by(cluster_id, Beff) %>%
  summarise(Value_m = round(mean(Value), 2),
            PI_low = rethinking::PI(Value, 0.90)[1],
            PI_high = rethinking::PI(Value, 0.90)[2], .groups = "drop")

# pooled average and PI
BEF_trim %>%
  filter(cluster_id != "F") %>%
  group_by(Beff) %>%
  summarise(Value_m = round(mean(Value), 2),
            PI_low = rethinking::PI(Value, 0.90)[1],
            PI_high = rethinking::PI(Value, 0.90)[2], .groups = "drop")

# is the mean within the percentile interval 90
BEF_pool %>%
  mutate(within_PI = if_else((Value_m > PI_low) & (Value_m < PI_high), 1, 0 )) %>%
  filter(within_PI != 1) %>%
  View()

# remove cluster F because it only had two usable sites
BEF_pool <- 
  BEF_pool %>%
  filter(cluster_id != "F")

# check some numbers of the manuscript
BEF_pool %>%
  filter(Beff == "NBE")

# compare total to local selection and complementarity

# which effects to plot
eff_in <- c("LC", "TC", "LS", "TS")

# plot the effects
p1 <- 
  ggplot(data = BEF_pool %>% 
           filter(Beff %in% eff_in) %>%
           mutate(Beff = factor(Beff, levels = eff_in) )) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_point(mapping = aes(x = Beff, y = Value_m, 
                           colour = Beff, fill = Beff, group = cluster_id),
             position = position_dodge(0.75), size = 3) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75)) + 
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2.75) +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in)) +
  ylab("Effect (g)") +
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
             position = position_dodge(0.75), size = 3) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75))  +
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2.75) +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab(NULL) +
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
             position = position_dodge(0.75), size = 3) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75)) +
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2.75) +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab(NULL) +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p3)

# arrange this plot
p123 <- 
  ggarrange(p1, p2, p3, ncol = 3, nrow = 1,
            labels = c("a", "b", "c"),
            widths = c(1.15, 1, 1),
            hjust = -0.05,
            vjust = 0.9,
            font.label = list(size = 11, face = "plain"))
plot(p123)

# save this plot
ggsave(filename = here("figures/fig6.png"), p123,
       unit = "cm", width = 20, height = 7.5)

# check the correlations among these different effects
BEF_pool %>%
  select(cluster_id, Beff, Value_m) %>%
  pivot_wider(id_cols = "cluster_id",
              names_from = "Beff",
              values_from = "Value_m") %>%
  select(-cluster_id) %>%
  cor() %>%
  corrplot(method = c("ellipse"))


# do insurance effects depend on environmental heterogeneity?

# load in the environmental dispersion data
env_disp <- readRDS(file = here("results/benthic_env_dispersion.rds"))
head(env_disp)

# using the original data, subset NBE and IT and calculate proportion of IT
SI_env <- 
  BEF_pool %>%
  filter(Beff == c("SI")) %>%
  select(-Beff) %>%
  rename(SI = Value_m)

# join the env_disp data to the BEF_pool data
SI_env <- right_join(env_disp, SI_env, by = "cluster_id")

# plot the relationship between environmental dispersion and the insurance effect

# run the simple linear regressions to test this relationship
lm.a <- lm(SI ~ field_dispersion, data = SI_env)
lm.a <- summary(lm.a)
r2 <- round(lm.a$r.squared, 2)
Fst <- round(lm.a$fstatistic[1], 2)
df1 <- round(lm.a$fstatistic[2], 0)
df2 <- round(lm.a$fstatistic[3], 0)
pval <- round(lm.a$coefficients[2, 4], 2) 

p1 <- 
  ggplot(data = SI_env) +
  geom_point(mapping = aes(x = field_dispersion, y = SI), size = 2) +
  geom_errorbar(mapping = aes(x = field_dispersion, 
                              ymin = PI_low,
                              ymax = PI_high), width = 0) +
  geom_smooth(mapping = aes(x = field_dispersion, y = SI), 
              method = "lm", se = TRUE, alpha = 0.25,
              size = 0.5, colour = "black") +
  annotate(geom = "text", x = 1.5, y = 3.24, label = bquote(r^2~" = "~.(r2) ), size = 3) +
  annotate(geom = "text", x = 2.2, y = 3.2, label = bquote(F[.(df1)~","~.(df2)]~" = "~.(Fst) ), size = 3) +
  annotate(geom = "text", x = 2.9, y = 3.24, label = bquote("P = "~.(pval) ), size = 3) +
  ylab("Spatial insurance (g)") +
  xlab("Multivariate dispersion") +
  theme_meta()
plot(p1)


# run the simple linear regressions to test this relationship
lm.b <- lm(SI ~ field_dispersion, data = SI_env %>% filter(cluster_id != "H"))
lm.b <- summary(lm.b)
r2 <- round(lm.b$r.squared, 2)
Fst <- round(lm.b$fstatistic[1], 2)
df1 <- round(lm.b$fstatistic[2], 0)
df2 <- round(lm.b$fstatistic[3], 0)
pval <- round(lm.b$coefficients[2, 4], 2) 

# plot the relationship between environmental dispersion and the insurance effect
p2 <- 
  ggplot(data = SI_env %>% filter(cluster_id != "H")) +
  geom_point(mapping = aes(x = field_dispersion, y = SI), size = 2) +
  geom_errorbar(mapping = aes(x = field_dispersion, 
                              ymin = PI_low,
                              ymax = PI_high), width = 0) +
  geom_smooth(mapping = aes(x = field_dispersion, y = SI), 
              method = "lm", se = TRUE, alpha = 0.25,
              size = 0.5, colour = "red", linetype = "dashed") +
  annotate(geom = "text", x = 1.5, y = -0.46, label = bquote(r^2~" = "~.(r2) ), size = 3) +
  annotate(geom = "text", x = 2.2, y = -0.5, label = bquote(F[.(df1)~","~.(df2)]~" = "~.(Fst) ), size = 3) +
  annotate(geom = "text", x = 2.9, y = -0.46, label = bquote("P = "~.(pval) ), size = 3) +
  ylab("Spatial insurance (g)") +
  xlab("Multivariate dispersion") +
  theme_meta()
plot(p2)

p12 <- ggarrange(p1, p2, labels = c("a", "b"),
                 font.label = list(size = 11, face = "plain")
                 )

ggsave(filename = here("figures/fig7.png"), p12,
       unit = "cm", width = 16, height = 7)


# BES talk figures

# summarise the data into a pooled estimate
BEF_pool2 <- 
  BEF_trim %>%
  filter(cluster_id != "F") %>%
  group_by(Beff) %>%
  summarise(Value_m = round(mean(Value), 2),
            PI_low = rethinking::PI(Value, 0.90)[1],
            PI_high = rethinking::PI(Value, 0.90)[2], .groups = "drop")

# calculate the percentage of total complementarity due to local complementarity
LC <- filter(BEF_pool2, Beff == "LC")[["Value_m"]]
TS <- filter(BEF_pool2, Beff == "TS")[["Value_m"]]
LS <- filter(BEF_pool2, Beff == "LS")[["Value_m"]]

(LC/(LS - TS + LC))*100

# compare net biodiversity effects, total complementarity and total selection
eff_in <- c("NBE", "TC", "TS")
p.BES7 <- 
  ggplot(data = BEF_pool2 %>%
           filter(Beff %in% eff_in) %>%
           mutate(Beff = factor(Beff, levels = eff_in))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_point(mapping = aes(x = Beff, y = Value_m, 
                           colour = Beff, fill = Beff),
             position = position_dodge(0.75), size = 2.5) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff),
                width = 0, position = position_dodge(0.75))  +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab("Effect (g)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p.BES7)

ggsave(filename = here("figures/fig_BES7.png"), p.BES7,
       unit = "cm", width = 7, height = 7)


# compare net biodiversity effects, total complementarity and total selection
eff_in <- c("NBE", "TC", "NO", "IT")
p.BES8 <- 
  ggplot(data = BEF_pool2 %>%
           filter(Beff %in% eff_in) %>%
           mutate(Beff = factor(Beff, levels = eff_in))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_point(mapping = aes(x = Beff, y = Value_m, 
                           colour = Beff, fill = Beff),
             position = position_dodge(0.75), size = 2.5) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff),
                width = 0, position = position_dodge(0.75))  +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab("Effect (g)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p.BES8)

ggsave(filename = here("figures/fig_BES8.png"), p.BES8,
       unit = "cm", width = 7, height = 7)

# compare net biodiversity effects, total complementarity and total selection
eff_in <- c("AS", "TI", "SI", "ST")
p.BES9 <- 
  ggplot(data = BEF_pool2 %>%
           filter(Beff %in% eff_in) %>%
           mutate(Beff = factor(Beff, levels = eff_in))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_point(mapping = aes(x = Beff, y = Value_m, 
                           colour = Beff, fill = Beff),
             position = position_dodge(0.75), size = 2.5) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff),
                width = 0, position = position_dodge(0.75))  +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab("Effect (g)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p.BES9)

ggsave(filename = here("figures/fig_BES9.png"), p.BES9,
       unit = "cm", width = 7, height = 7)

### END
