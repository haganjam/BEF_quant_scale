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
library(cowplot)
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
  
# check summary statistics for the manuscript
BEF_sum %>%
  filter(Beff == "NBE")

# calculate the pooled effect size across clusters
beff_vec <- unique(BEF_sum$Beff)
meta_list <- vector("list", length = length(beff_vec))
for(i in 1:length(meta_list)) {
  
  x <- 
    BEF_sum %>%
    filter(Beff == beff_vec[i])
  
  # calculate the effect size: raw mean
  y <- escalc(mi = x$Value_m, sdi = x$Value_sd, ni = x$n, measure = "MN")
  
  # calculate the grand mean
  z <- rma(yi, vi, method = "REML", data = y, slab = x$cluster_id,
           test = "knha")
  beta <- z$beta[,1]
  names(beta) <- NULL
  row.names(beta) <- NULL
  
  # pull these into a data.frame
  df <- tibble(Beff = beff_vec[i],
               grand_mean = beta,
               tau = sqrt(z$tau2),
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
ED_table_2 <- arrange(BEF_grand, Beff)
print(ED_table_2)

# output the table as a .csv file
write_csv(x = ED_table_2, file = "figures/ED_table_2.csv")

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
BEF_col <- list(p1 = v_col_BEF(BEF_pars$p1),
                p2 = v_col_BEF(BEF_pars$p2),
                p3 = v_col_BEF(BEF_pars$p3))

BEF_fullnames <- list(p1 = c("Total selection", "Total complementarity", "Net biodiversity effect"),
                      p2 = c("Total insurance", "Non-random overyielding"),
                      p3 = c("Spatio-temporal insurance", "Temporal insurance", "Spatial insurance", "Average selection"))

# set-up significance annotations
BEF_sig <- 
  list(p1 = tibble(Beff = c("TS", "TC", "NBE"), Value = c(80, 80, 80), sig = "*"),
       p2 = tibble(Beff = c("IT", "NO"), Value = c(85, 85), sig = "*"),
       p3 = tibble(Beff = c("ST", "TI", "SI",  "AS"), Value = c(95, 95, 95, 95), sig = c("", "", "*", "*")))

# plot the ith plot
plot_list <- vector("list", length = length(BEF_pars))
for(i in 1:length(BEF_pars)) {
  
  x <- 
    BEF_plot %>% 
    filter(Beff %in% BEF_pars[[i]])
  
  x$Beff <- factor(x$Beff, levels = BEF_pars[[i]])
  levels(x$Beff) <- BEF_fullnames[[i]]
  
  x_sum <-
    BEF_sum %>%
    filter(Beff %in% BEF_pars[[i]])
  
  x_sum$Beff <- factor(x_sum$Beff, levels = BEF_pars[[i]])
  levels(x_sum$Beff) <- BEF_fullnames[[i]]
  
  x_grand <- 
    BEF_grand %>%
    filter(Beff %in% BEF_pars[[i]])
  
  x_grand$Beff <- factor(x_grand$Beff, levels = BEF_pars[[i]])
  levels(x_grand$Beff) <- BEF_fullnames[[i]]
  
  x_sig <- BEF_sig[[i]]
  
  x_sig$Beff <- factor(x_sig$Beff, levels = BEF_pars[[i]])
  levels(x_sig$Beff) <- BEF_fullnames[[i]]
  
  # make sure the colours have the correct names
  names(BEF_col[[i]]) <- BEF_fullnames[[i]]
  
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
    geom_text(data = x_sig,
              mapping = aes(x = Beff, y = Value, label = sig),
              size = 6) +
    scale_colour_manual(values = BEF_col[[i]]) +
    scale_fill_manual(values = BEF_col[[i]]) +
    theme_meta() +
    xlab(NULL) +
    ylab(NULL) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, 
                                                            width = 17)) +
    theme(legend.position = "none", 
          axis.text.y = element_text(size = 9)) +
    coord_flip()
  
  plot_list[[i]] <- p
  
}

# cluster order (top to bottom): J, I, H, G, E, D, C, B, A 

# arrange the plots
p123 <- 
  plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
            nrow = 3, ncol = 1, align = "v",
            labels = c("a", "b", "c"), label_size = 11,
            label_fontface = "plain",
            rel_heights = c(1.5, 1, 2)) 
# plot(p123)

ggsave(filename = "figures/MAIN_fig_4.svg", 
       p123, units = "cm", width = 13, height = 18)


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
  ggbeeswarm::geom_quasirandom(data = BEF_plot %>% filter(Beff == "SI"),
             mapping = aes(x = field_dispersion, y = Value), shape = 1,
             stroke = 0.1, alpha = 0.6, colour = "darkgrey",
             width = 0.05) +
  geom_point(data = BEF_SI,
                mapping = aes(x = field_dispersion, y = Value_m),
             colour = "#333333") +
  geom_errorbar(data = BEF_SI,
                mapping = aes(x = field_dispersion, ymin = HPDI_low, ymax = HPDI_high),
                width = 0, linewidth = 0.3, colour = "#333333") +
  geom_abline(data = lm_df_plot,
              mapping = aes(intercept = intercept, slope = slope, group = SI_rep),
              linewidth = 0.1, alpha = 0.1, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  ylab("Spatial insurance (g)") +
  xlab("Multivariate dispersion") +
  theme_meta()
plot(p1)

slope_sum <- 
  lm_df %>%
  summarise(slope_m = mean(slope),
            HPDI_low = HPDI(slope, 0.95)[1],
            HPDI_high = HPDI(slope, 0.95)[2])
print(slope_sum)
  
p2 <- 
  ggplot() +
  geom_line(data = lm_df,
            mapping = aes(x = slope), stat="density", 
            colour = "red", linewidth = 0.8,
            alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = slope_sum,
             mapping = aes(x = slope_m, y = 2.1), size = 2, colour = "red") +
  geom_errorbarh(data = slope_sum,
                 mapping = aes(y = 2.1,xmin = HPDI_low, xmax = HPDI_high),
                 height = 0, colour = "red", linewidth = 0.75) +
  scale_y_continuous(limits = c(0, 2.3), expand = c(0, 0)) +
  ylab("Density") +
  xlab("Slope estimate") +
  theme_meta()

p12 <- ggarrange(p1, p2, labels = c("a", "b"),
                 font.label = list(size = 11, face = "plain"),
                 nrow = 1, ncol = 2, widths = c(1,0.7))
plot(p12)

ggsave(filename = "figures/MAIN_fig_5.svg", p12,
       unit = "cm", width = 15, height = 8)


# remove the outlier and rerun the analysis

# get a list of the mean and PI
BEF_SI_out <- 
  full_join(env_disp, 
            filter(BEF_sum, Beff == "SI"), by = "cluster_id")

# remove cluster H
BEF_SI_out <- 
  BEF_SI_out %>%
  filter(cluster_id != "H")

# convert the different estimates into a list
SI_list_out <- 
  BEF_dat %>%
  filter(Beff == "SI") %>%
  filter(cluster_id != "H") %>%
  mutate(Value_id = paste(mono_rep, RYE, sep = "_")) %>%
  split(., .$Value_id)

# run a regression looped over all the different estimates
lm_list_out <- vector("list", length = length(SI_list_out))
for(i in 1:length(SI_list_out)) {
  
  lm_x <- lm(Value ~ field_dispersion, data =  SI_list_out[[i]])
  lm_est <- coef(lm_x)
  names(lm_est) <- NULL
  df <- tibble(intercept = lm_est[1],
               slope = lm_est[2])
  
  lm_list_out[[i]] <- df
  
}

# bind into a data.frame
lm_df_out <- bind_rows(lm_list_out, .id = "SI_rep")

# check the distribution
hist(lm_df_out$slope)
mean(lm_df_out$slope)

# sample a few of these slopes
set.seed(4908325)
lm_df_plot_out <- lm_df_out[sample(1:nrow(lm_df_out), 200), ]

# plot the relationship between SI and environmental dispersion
p3 <- 
  ggplot() +
  ggbeeswarm::geom_quasirandom(data = BEF_plot %>% filter(Beff == "SI", cluster_id != "H"),
                               mapping = aes(x = field_dispersion, y = Value), shape = 1,
                               stroke = 0.1, alpha = 0.6, colour = "darkgrey",
                               width = 0.05) +
  geom_point(data = BEF_SI_out,
             mapping = aes(x = field_dispersion, y = Value_m),
             colour = "#333333") +
  geom_errorbar(data = BEF_SI_out,
                mapping = aes(x = field_dispersion, ymin = HPDI_low, ymax = HPDI_high),
                width = 0, linewidth = 0.3, colour = "#333333") +
  geom_abline(data = lm_df_plot_out,
              mapping = aes(intercept = intercept, slope = slope, group = SI_rep),
              linewidth = 0.1, alpha = 0.1, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  ylab("Spatial insurance (g)") +
  xlab("Multivariate dispersion") +
  theme_meta()
plot(p3)

slope_sum_out <- 
  lm_df_out %>%
  summarise(slope_m = mean(slope),
            HPDI_low = HPDI(slope, 0.95)[1],
            HPDI_high = HPDI(slope, 0.95)[2])
print(slope_sum_out)

p4 <- 
  ggplot() +
  geom_line(data = lm_df_out,
            mapping = aes(x = slope), stat="density", 
            colour = "red", linewidth = 0.8,
            alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = slope_sum_out,
             mapping = aes(x = slope_m, y = 2.1), size = 2, colour = "red") +
  geom_errorbarh(data = slope_sum_out,
                 mapping = aes(y = 2.1,xmin = HPDI_low, xmax = HPDI_high),
                 height = 0, colour = "red", linewidth = 0.75) +
  scale_y_continuous(limits = c(0, 2.3), expand = c(0, 0)) +
  ylab("Density") +
  xlab("Slope estimate") +
  theme_meta()

p34 <- ggarrange(p3, p4, labels = c("a", "b"),
                 font.label = list(size = 11, face = "plain"),
                 nrow = 1, ncol = 2, widths = c(1,0.7))
plot(p34)

ggsave(filename = "figures/ED_fig_4.svg", p34,
       unit = "cm", width = 15, height = 8)

### END
