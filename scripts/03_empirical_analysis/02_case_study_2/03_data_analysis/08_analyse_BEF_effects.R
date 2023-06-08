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
library(vegan)
library(corrplot)

# load plotting theme
source("scripts/Function_plotting_theme.R")
source("scripts/03_empirical_analysis/helper_functions.R")

# estimates of the different BEF effects
BEF_dat <- readRDS("results/BEF_effects_case_study_2.rds")
head(BEF_dat)
dim(BEF_dat)
summary(BEF_dat)

# load in the environmental dispersion data
env_disp <- readRDS(file = "results/benthic_env_dispersion.rds")
head(env_disp)

# add the environmental dispersion data to the BEF effects data
BEF_dat <- full_join(env_disp, BEF_dat, by = "cluster_id")
head(BEF_dat)

# check how many of each effect there are: Should be 100 000
BEF_dat %>%
  group_by(cluster_id, Beff) %>%
  summarise(n = n())

# check the internal consistency of the net biodiversity effects
BEF_dat %>%
  filter(cluster_id == "A") %>%
  filter( !(Beff %in% c("NBE", "LS", "LC", "TS", "IT")) ) %>%
  group_by(mono_rep, RYE) %>%
  summarise(s = sum(Value))

BEF_dat %>%
  filter(cluster_id == "A", Beff == "NBE")

# remove cluster F because it only has two usable sites
BEF_dat <- 
  BEF_dat %>%
  filter(cluster_id != "F")

# remove the LS and LC effects
BEF_dat <- 
  BEF_dat %>%
  filter( !(Beff %in% c("LS", "LS")) )

# what is the magnitude of different biodiversity effects?

# calculate the summary dataset
BEF_sum <- 
  BEF_dat %>%
  group_by(cluster_id, Beff) %>%
  summarise(Value_m = round(mean(Value), 2),
            HPDI_low = HPDI(Value, 0.89)[1],
            HPDI_high = HPDI(Value, 0.89)[2], .groups = "drop")
  

# get a nice colour palette
col_pal <- wesanderson::wes_palette("Darjeeling1", n = 9, type = "continuous")

# get 100 samples for each biodiversity effect
id <- dplyr::distinct(BEF_dat[,c("mono_rep", "RYE")])

# set seed for reproducibility
set.seed(487532)
id <- id[sample(1:nrow(id), 100), ]
id$id_est <- paste(id$mono_rep, id$RYE, sep = "_")
head(id)

# get a subset for plotting
BEF_plot <- 
  BEF_dat %>%
  group_by(cluster_id, Beff) %>%
  filter(paste(mono_rep, RYE, sep = "_") %in% id$id_est) %>%
  ungroup()

# plot the different biodiversity in different clusters

# set-up the parameter lists
BEF_pars <- list(p1 =  c("TS", "TC", "NBE"),
                 p2 = c("IT", "NO"),
                 p3 = c("ST", "TI", "SI",  "AS"))

# set-up the colour lists
BEF_col <- list(p1 = col_pal[3:1],
                p2 = col_pal[5:4],
                p3 = col_pal[9:6])

# set-up some segments for plot comparisons
segments <- data.frame(xstart = 0,
                       xend = 0.15,
                       effect = 30)

# plot the ith plot
plot_list <- vector("list", length = length(BEF_pars))
for(i in 1:length(BEF_pars)) {
  
  x <- 
    BEF_plot %>% 
    filter(Beff %in% BEF_pars[[i]])
  
  x$Beff <- factor(x$Beff, levels = BEF_pars[[i]])
  
  x_sum <-
    BEF_sum %>%
    filter(Beff %in% BEF_pars[[i]])
  
  x_sum$Beff <- factor(x_sum$Beff, levels = BEF_pars[[i]])
  
  p <-
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(data = x, 
               mapping = aes(x = Beff, y = Value, colour = Beff, group = cluster_id), 
               position = position_jitterdodge(jitter.width = 0.15,
                                               dodge.width = 0.9),
               alpha = 0.2, shape = 1, size = 1, stroke = 0.25) +
    geom_errorbar(data = x_sum,
                  mapping = aes(x = Beff, ymin = HPDI_low, ymax = HPDI_high, 
                                colour = Beff, group = cluster_id),
                  position = position_dodge(width = 0.9), width = 0) +
    geom_point(data = x_sum,
               mapping = aes(x = Beff, y = Value_m, fill = Beff, group = cluster_id),
               position = position_dodge(width = 0.9), shape = 21, size = 2.7,
               colour = "black", stroke = 0.1) +
    geom_label(data = x_sum, 
               mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
               colour = "white", position = position_dodge(0.9),
               label.size = NA, alpha = 0, size = 2.2,
               fontface = "bold") +
    # geom_segment(data = segments,
                 # mapping = aes(x = xstart, xend = xend,
                               # y = effect, yend = effect),
                 # size = 1.25, colour = "red") +
    scale_colour_manual(values = BEF_col[[i]]) +
    scale_fill_manual(values = BEF_col[[i]]) +
    theme_meta() +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    coord_flip()
  
  plot_list[[i]] <- p
  
}

# check the plots
plot_list[[1]]
plot_list[[2]]
plot_list[[3]]

ggsave(filename = "figures/fig_4i.png", plot_list[[1]],
       dpi = 500, units = "cm", width = 10, height = 6)
ggsave(filename = "figures/fig_4ii.png", plot_list[[2]],
       dpi = 500, units = "cm", width = 10, height = 4)
ggsave(filename = "figures/fig_4iii.png", plot_list[[3]],
       dpi = 500, units = "cm", width = 10, height = 8)


# check some numbers of the manuscript
BEF_pool %>%
  filter(Beff %in% c("NBE", "IT") ) %>%
  select(cluster_id, Beff, Value_m) %>%
  pivot_wider(id_cols = c("cluster_id"),
              names_from = "Beff",
              values_from = "Value_m") %>%
  mutate(IT_perc = (IT/NBE)*100 ) %>%
  summarise(IT_perc_M = mean(IT_perc),
            IT_perc_SD = sd(IT_perc))

# compare total to local selection and complementarity
x <- BEF_trim %>%
  filter(cluster_id != "F") %>%
  group_by(Beff) %>%
  summarise(Value_m = round(mean(Value), 2),
            PI_low = rethinking::PI(Value, 0.90)[1],
            PI_high = rethinking::PI(Value, 0.90)[2], .groups = "drop")

(x[x$Beff == "LC",]$Value_m)/((x[x$Beff == "LS",]$Value_m) - (x[x$Beff == "TS",]$Value_m) + (x[x$Beff == "LC",]$Value_m) )


# set-up some segments for plot comparisons
segments <- data.frame(xstart = 0,
                       xend = 0.3,
                       effect = 30)

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
  geom_segment(data = segments,
               mapping = aes(x = xstart, xend = xend,
                             y = effect, yend = effect),
               size = 1.25, colour = "red") +
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
  geom_segment(data = segments,
               mapping = aes(x = xstart, xend = xend,
                             y = effect, yend = effect),
               size = 1.25, colour = "red") +
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
  geom_segment(data = segments,
               mapping = aes(x = xstart, xend = xend,
                             y = effect, yend = effect),
               size = 1.25, colour = "red") +
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
env_disp <- readRDS(file = "results/benthic_env_dispersion.rds")
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

# check the summary of the SI_env data
summary(SI_env)

p1 <- 
  ggplot(data = SI_env) +
  geom_point(mapping = aes(x = field_dispersion, y = SI), size = 2) +
  geom_errorbar(mapping = aes(x = field_dispersion, 
                              ymin = PI_low,
                              ymax = PI_high), width = 0) +
  annotate(geom = "text", x = 1.35, y = 3.26, label = bquote(r^2~" = "~.(r2) ), size = 2.8) +
  annotate(geom = "text", x = 2.15, y = 3.2, label = bquote(F[.(df1)~","~.(df2)]~" = "~.(Fst) ), size = 2.8) +
  annotate(geom = "text", x = 2.95, y = 3.24, label = bquote("P = "~.(pval) ), size = 2.8) +
  ylab("Spatial insurance (g)") +
  xlab("Multivariate dispersion") +
  scale_y_continuous(limits = c(-1.1, 3.6)) +
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
  geom_smooth(mapping = aes(x = field_dispersion, y = SI), 
              method = "lm", se = TRUE, alpha = 0.25,
              size = 0.5, colour = "red", linetype = "dashed") +
  geom_point(mapping = aes(x = field_dispersion, y = SI), size = 2) +
  geom_errorbar(mapping = aes(x = field_dispersion, 
                              ymin = PI_low,
                              ymax = PI_high), width = 0) +
  annotate(geom = "text", x = 1.3, y = 3.26, label = bquote(r^2~" = "~.(r2) ), size = 2.8) +
  annotate(geom = "text", x = 2.1, y = 3.2, label = bquote(F[.(df1)~","~.(df2)]~" = "~.(Fst) ), size = 2.8) +
  annotate(geom = "text", x = 2.95, y = 3.24, label = bquote("P = "~.(pval) ), size = 2.8) +
  ylab("") +
  xlab("Multivariate dispersion") +
  scale_y_continuous(limits = c(-1.1, 3.6)) +
  theme_meta() +
  theme(axis.title.y = element_text(size = 1))
plot(p2)

p12 <- ggarrange(p1, p2, labels = c("a", "b"),
                 font.label = list(size = 11, face = "plain")
                 )
plot(p12)

ggsave(filename = here("figures/fig7.png"), p12,
       unit = "cm", width = 14, height = 8)

### END
