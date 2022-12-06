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


# what is the magnitude of different biodiversity effects?

# summarise the data into a pooled estimate
BEF_pool <- 
  BEF_trim %>%
  group_by(cluster_id, Beff) %>%
  summarise(Value_m = round(mean(Value), 2),
            PI_low = rethinking::HPDI(Value, 0.90)[1],
            PI_high = rethinking::HPDI(Value, 0.90)[2], .groups = "drop")

# is the mean within the percentile interval 90
BEF_pool %>%
  mutate(within_PI = if_else((Value_m > PI_low) & (Value_m < PI_high), 1, 0 )) %>%
  filter(within_PI != 1) %>%
  View()

# remove cluster F because it only had two usable sites
BEF_pool <- 
  BEF_pool %>%
  filter(cluster_id != "F")

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
             position = position_dodge(0.75), size = 2.5) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75)) + 
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2) +
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
             position = position_dodge(0.75), size = 2.5) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75))  +
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2) +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab("Effect (g)") +
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
             position = position_dodge(0.75), size = 2.5) +
  geom_errorbar(mapping = aes(x = Beff, ymin = PI_low, ymax = PI_high,
                              colour = Beff, group = cluster_id),
                width = 0, position = position_dodge(0.75)) +
  geom_label(mapping = aes(x = Beff, y = Value_m, label = cluster_id, group = cluster_id),
             colour = "white", position = position_dodge(0.75),
             label.size = NA, alpha = 0, size = 2) +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in )) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in )) +
  ylab("Effect (g)") +
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
p.BES6a <- 
  ggplot(data = SI_env) +
  geom_point(mapping = aes(x = field_dispersion, y = SI), size = 2) +
  geom_errorbar(mapping = aes(x = field_dispersion, 
                              ymin = PI_low,
                              ymax = PI_high), width = 0) +
  geom_smooth(mapping = aes(x = field_dispersion, y = SI), 
              method = "lm", se = TRUE, alpha = 0.25,
              size = 0.5, colour = "black") +
  ylab("Spatial insurance (g)") +
  xlab("Multivariate dispersion") +
  theme_meta()
plot(p.BES6a)

ggsave(filename = here("figures/fig_BES6a.png"), p.BES6a,
       unit = "cm", width = 9, height = 7)

lm.a <- lm(SI ~ field_dispersion, data = SI_env)
summary(lm.a)

# plot the relationship between environmental dispersion and the insurance effect
p.BES6b <- 
  ggplot(data = SI_env %>% filter(cluster_id != "H")) +
  geom_point(mapping = aes(x = field_dispersion, y = SI), size = 2) +
  geom_errorbar(mapping = aes(x = field_dispersion, 
                              ymin = PI_low,
                              ymax = PI_high), width = 0) +
  geom_smooth(mapping = aes(x = field_dispersion, y = SI), 
              method = "lm", se = TRUE, alpha = 0.25,
              size = 0.5, colour = "red", linetype = "dashed") +
  ylab("Spatial insurance (g)") +
  xlab("Multivariate dispersion") +
  theme_meta()
plot(p.BES6b)

lm.b <- lm(SI ~ field_dispersion, data = SI_env %>% filter(cluster_id != "H"))
summary(lm.b)

ggsave(filename = here("figures/fig_BES6b.png"), p.BES6b,
       unit = "cm", width = 9, height = 7)


# BES talk figures

# summarise the data into a pooled estimate
BEF_pool2 <- 
  BEF_trim %>%
  group_by(Beff) %>%
  summarise(Value_m = round(mean(Value), 2),
            PI_low = rethinking::HPDI(Value, 0.90)[1],
            PI_high = rethinking::HPDI(Value, 0.90)[2], .groups = "drop")

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



# does spatial or temporal beta-diversity affect the magnitude of insurance effects?

# load the analysis data
data <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/biomass_env_analysis_data.csv"))
head(data)

# quantify community metrics

# convert the data to a wide format
data_wide <- 
  data %>%
  select(-M) %>%
  pivot_wider(id_cols = c("cluster_id", "buoy_id", "time", "PC1", "PC2"),
              names_from = "OTU",
              values_from = "Y")

# get species by site matrix
com <- 
  data_wide %>%
  select(Barn, Bryo, Bumpi, Hydro, Seasq)
  
# get a site only matrix
site <- 
  data_wide %>%
  select(cluster_id, buoy_id, time, PC1, PC2)

# make a list of the community matrix with and without Barnacles
com_list <- list(com, com[, names(com) != "Barn"])

# make an output list
IT_env_beta <- vector("list", length = length(com_list))

for(i in 1:length(com_list)) {
  
  # calculate spatial beta-diversity within clusters
  beta_spat <- 
    
    split(com_list[[i]], paste0(site$cluster_id, site$time)) %>%
    
    lapply(., function(sp_mat) {
      
      # calculate abundance-based multi-site beta-diversity
      x <- beta.multi.abund(x = sp_mat, index.family = "bray")
      
      # calculate presence-absence based multi-site beta-diversity
      y <- beta.multi(x = vegan::decostand(sp_mat, "pa") )
      
      # pull into a data.frame
      df <- tibble(S.beta_bray_bal = x$beta.BRAY.BAL,
                   S.beta_bray_gra = x$beta.BRAY.GRA,
                   S.beta_bray = x$beta.BRAY,
                   S.beta_sor_turn = y$beta.SIM,
                   S.beta_sor_nes = y$beta.SNE,
                   S.beta_sor = y$beta.SOR)
      
      return(df)
      
    } )
  
  # summarise the average spatial beta diversity across time for each cluster
  beta_spat_df <- 
    bind_cols(site %>%
                select(cluster_id, time) %>%
                distinct(), 
              bind_rows(beta_spat)
    ) %>%
    group_by(cluster_id) %>%
    summarise(across(.cols = starts_with("S.beta"), mean, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(cluster_id)
  
  
  # calculate temporal beta-diversity within clusters
  beta_temp <- 
    
    split(com_list[[i]], paste0(site$cluster_id, site$buoy_id)) %>%
    
    lapply(., function(sp_mat) {
      
      # calculate abundance-based multi-site beta-diversity
      x <- beta.multi.abund(x = sp_mat, index.family = "bray")
      
      # calculate presence-absence based multi-site beta-diversity
      y <- beta.multi(x = vegan::decostand(sp_mat, "pa") )
      
      # pull into a data.frame
      df <- tibble(T.beta_bray_bal = x$beta.BRAY.BAL,
                   T.beta_bray_gra = x$beta.BRAY.GRA,
                   T.beta_bray = x$beta.BRAY,
                   T.beta_sor_turn = y$beta.SIM,
                   T.beta_sor_nes = y$beta.SNE,
                   T.beta_sor = y$beta.SOR)
      
      return(df)
      
    } )
  
  # summarise the average spatial beta diversity across time for each cluster
  beta_temp_df <- 
    bind_cols(site %>%
                select(cluster_id, buoy_id) %>%
                distinct(), 
              bind_rows(beta_temp)
    ) %>%
    group_by(cluster_id) %>%
    summarise(across(.cols = starts_with("T.beta"), mean, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(cluster_id)
  
  
  # add spatial and temporal beta diversity to the IT_env data
  beta_env_df <- 
    left_join(IT_env, 
              full_join(beta_spat_df, beta_temp_df, by = "cluster_id"),
              by = "cluster_id"
    )
  
  # add the IT_env_beta to a list
  IT_env_beta[[i]] <- beta_env_df
  
}

# correlation among different beta-diversity components and total insurance effects
IT_env_beta[[1]][, -c(1:2)] %>%
  cor(method = "spearman") %>%
  corrplot(corr = ., method = "ellipse")

IT_env_beta[[2]][, -c(1:2)] %>%
  cor(method = "spearman") %>%
  corrplot(corr = ., method = "ellipse")

# how does spatial beta diversity without barnacles correlate with insurance effects?
plot(IT_env_beta[[2]]$S.beta_sor, IT_env_beta[[2]]$IT_m)
cor.test(IT_env_beta[[2]]$S.beta_sor, IT_env_beta[[2]]$IT_m, method = "spearman")

# how does spatial beta diversity without barnacles correlate with environmental dispersion?
plot(IT_env_beta[[2]]$GIS_dispersion, IT_env_beta[[2]]$S.beta_sor)
cor.test(IT_env_beta[[2]]$GIS_dispersion, IT_env_beta[[2]]$S.beta_sor, method = "spearman")
range(IT_env_beta[[2]]$S.beta_sor)

# how does spatial beta diversity with barnacles correlate with environmental dispersion?
plot(IT_env_beta[[1]]$GIS_dispersion, IT_env_beta[[1]]$S.beta_sor)
cor.test(IT_env_beta[[1]]$GIS_dispersion, IT_env_beta[[1]]$S.beta_sor, method = "spearman")
range(IT_env_beta[[1]]$S.beta_sor)


# what is the beta diversity like when we exclude barnacles?

# compare beta diversity with and without barnacles
compare_beta <- 
  bind_rows(IT_env_beta[[1]] %>%
            select(contains("beta")) %>%
            mutate(barn_pa = "present"),
          IT_env_beta[[2]] %>%
            select(contains("beta")) %>%
            mutate(barn_pa = "absent")) %>%
  pivot_longer(cols = contains("beta"),
               names_to = "beta_metric",
               values_to = "value")

ggplot(data = compare_beta,
       mapping = aes(x = beta_metric, y = value, fill = barn_pa)) +
  geom_col(position = position_dodge(width = 0.5), width = 0.5) +
  theme_meta() +
  theme(axis.text.x = element_text(angle = 90))


# were barnacles present everywhere and how abundant were they?

# barnacles present at x percentage of sites
sum(com$Barn>0)/nrow(com)*100

# on average, they are this abundant
Barn_RA <- com$Barn/apply(com, 1, sum)
mean(Barn_RA, na.rm = TRUE)
sd(Barn_RA, na.rm = TRUE)

# in this proportion of sites-time points, they made up more than 90% of the functioning
sum(Barn_RA > 0.90, na.rm = TRUE)/length(Barn_RA)

### END

