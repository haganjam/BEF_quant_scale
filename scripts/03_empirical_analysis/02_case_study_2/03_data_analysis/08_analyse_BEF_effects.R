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
library(metafor)

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
  filter( !(Beff %in% c("LC", "LS")) )

# refactor the Beff variable to get rid of LS and LC
BEF_dat$Beff <- factor(BEF_dat$Beff)

# what is the magnitude of different biodiversity effects?

# calculate the summary dataset
BEF_sum <- 
  BEF_dat %>%
  group_by(cluster_id, Beff) %>%
  summarise(Value_m = mean(Value),
            Value_sd = sd(Value),
            n = n(),
            HPDI_low = HPDI(Value, 0.95)[1],
            HPDI_high = HPDI(Value, 0.95)[2], .groups = "drop")
  
# calculate the pooled effect size across clusters
beff_vec <- unique(BEF_sum$Beff)
meta_list <- vector("list", length = length(beff_vec))
for(i in 1:length(meta_list)) {
  
  x <- 
    BEF_sum %>%
    filter(Beff == beff_vec[i])
  
  # calculate the effect size
  y <- escalc(mi = x$Value_m, sdi = x$Value_sd, ni = x$n, measure = "MN")
  
  # calculate the grand mean
  z <- rma(yi, vi, method = "ML", data = y, slab = x$cluster_id,
           test = "t")
  beta <- z$beta[,1]
  names(beta) <- NULL
  row.names(beta) <- NULL
  
  # pull these into a data.frame
  df <- tibble(Beff = beff_vec[i],
               grand_mean = beta,
               CI95_low = z$ci.lb,
               CI95_high = z$ci.ub,
               SE = z$se,
               tval = z$zval,
               df = z$dfs,
               pval = z$pval)
  
  meta_list[[i]] <- df
  
}

# bind into a data.frame
BEF_grand <- bind_rows(meta_list)

# output this table
BEF_grand$Beff <- factor(BEF_grand$Beff,
                         levels = c("NBE", "TC", "TS", "NO", "IT", "AS", "SI", "TI", "ST"))

# rearrange the effects into the correct order
table2 <- arrange(BEF_grand, Beff)

# output the table as a .csv file
write_csv(x = table2, file = "figures/table_2.csv")

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
  
  x_grand <- 
    BEF_grand %>%
    filter(Beff %in% BEF_pars[[i]])
  
  x_grand$Beff <- factor(x_grand$Beff, levels = BEF_pars[[i]])
  
  p <-
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 30, 
               linetype = "dashed", colour = "red", 
               alpha = 0.5, linewidth = 0.25) +
    geom_point(data = x, 
               mapping = aes(x = Beff, y = Value, colour = Beff, group = cluster_id), 
               position = position_jitterdodge(jitter.width = 0.15,
                                               dodge.width = 0.9),
               alpha = 0.15, shape = 1, size = 1.4, stroke = 0.20) +
    geom_errorbar(data = x_sum,
                  mapping = aes(x = Beff, ymin = HPDI_low, ymax = HPDI_high, 
                                colour = Beff, group = cluster_id),
                  position = position_dodge(width = 0.9), width = 0.01, 
                  linewidth = 0.3, alpha = 1) +
    geom_errorbar(data = x_grand,
                  mapping = aes(x = Beff, ymin = CI95_low, ymax = CI95_high,
                                colour = Beff), 
                  width = 0.05, linewidth = 0.35, colour = "black") +
    geom_point(data = x_grand,
               mapping = aes(x = Beff, y = grand_mean, fill = Beff),
               position = position_dodge(width = 0.9), shape = 23, size = 2.2,
               colour = "black", stroke = 0.5) +
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

# export the plots
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


# do insurance effects depend on environmental heterogeneity?

# get a list of the mean and PI
BEF_SI <- 
  full_join(env_disp, 
            filter(BEF_sum, Beff == "SI"), by = "cluster_id")

# convert the different estimates into a list
SI_list <- 
  BEF_dat %>%
  filter(Beff == "SI") %>%
  mutate(Value_id = paste(mono_rep, RYE, sep = "_")) %>%
  split(., .$Value_id)

# check the split
SI_list[[1]]

# run a regression looped over all the different estimates
lm_list <- vector("list", length = length(SI_list))
for(i in 1:length(SI_list)) {
  
  lm_x <- lm(Value ~ field_dispersion, data =  SI_list[[i]])
  lm_est <- coef(lm_x)
  names(lm_est) <- NULL
  df <- tibble(intercept = lm_est[1],
               slope = lm_est[2])
  
  lm_list[[i]] <- df
  
}

# bind into a data.frame
lm_df <- bind_rows(lm_list, .id = "SI_rep")

# check the distribution
hist(lm_df$slope)
mean(lm_df$slope)

# sample a few of these slopes
set.seed(4908325)
lm_df_plot <- lm_df[sample(1:nrow(lm_df), 200), ]

# get the regression line through the mean
lm_mean <- lm(Value_m ~ field_dispersion, data = BEF_SI)

# plot the relationship between SI and environmental dispersion
p1 <- 
  ggplot() +
  geom_point(data = BEF_plot %>% filter(Beff == "SI"),
             mapping = aes(x = field_dispersion, y = Value),
             position = position_jitter(width = 0.075), shape = 1,
             stroke = 0.1, alpha = 0.3) +
  geom_point(data = BEF_SI,
                mapping = aes(x = field_dispersion, y = Value_m)) +
  geom_errorbar(data = BEF_SI,
                mapping = aes(x = field_dispersion, ymin = HPDI_low, ymax = HPDI_high),
                width = 0, linewidth = 0.3) +
  geom_abline(data = lm_df_plot,
              mapping = aes(intercept = intercept, slope = slope, group = SI_rep),
              linewidth = 0.1, alpha = 0.1, colour = col_pal[1]) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  ylab("Spatial insurance (g)") +
  xlab("Multivariate dispersion") +
  theme_meta()

slope_sum <- 
  lm_df %>%
  summarise(slope_m = mean(slope),
            HPDI_low = HPDI(slope, 0.95)[1],
            HPDI_high = HPDI(slope, 0.95)[2])
  
p2 <- 
  ggplot() +
  geom_line(data = lm_df,
            mapping = aes(x = slope), stat="density", 
            colour = col_pal[1], linewidth = 0.8,
            alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = slope_sum,
             mapping = aes(x = slope_m, y = 2.1), size = 2, colour = col_pal[1]) +
  geom_errorbarh(data = slope_sum,
                 mapping = aes(y = 2.1,xmin = HPDI_low, xmax = HPDI_high),
                 height = 0, colour = col_pal[1], linewidth = 0.75) +
  scale_y_continuous(limits = c(0, 2.3), expand = c(0, 0)) +
  ylab("Density") +
  xlab("Slope estimate") +
  theme_meta()

p12 <- ggarrange(p1, p2, labels = c("a", "b"),
                 font.label = list(size = 11, face = "plain"),
                 nrow = 1, ncol = 2, widths = c(1,0.7))
plot(p12)

ggsave(filename = "figures/fig_5.png", p12,
       unit = "cm", width = 15, height = 8)

### END
