
# Plymouth rock-pool data

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(here)
library(viridis)
library(ggpubr)
library(rethinking)

# set script to call partition functions from
source(here("scripts/01_partition_functions/01_isbell_2018_partition.R"))
source(here("scripts/Function_plotting_theme.R"))

# read in the data
ply_dat <- read_delim( here("data/Plymouth_data.csv"), delim = "," )

# check the data
head(ply_dat)
summary(ply_dat)

# analyse only the focal species
foc_spp <- 
  c("sargassum_muticum", "bifurcaria_bifurcata", "fucus_serratus", "laminaria_digitata")

# how many focal species are there?
tot_spp <- 4

# subset the focal species
ply_dat <- 
  ply_dat %>%
  mutate(species_richness = (tot_spp - removed) ) %>%
  select( month, shore, pool, composition, removed,
          grid_squares, species_richness, totcover, cover_nocrusts, all_of(foc_spp) )

# create a new mixture/monoculture column
# reorder the columns again
ply_dat <- 
  ply_dat %>%
  mutate(mono_mix = if_else(removed == "3", "monoculture",
                            if_else(removed == "4", "control", 
                                    if_else(removed == "0", "max_mixture", "mixture"))))
head(ply_dat)

# remove the composition = "None" treatment as this is a general control
ply_dat <- 
  ply_dat %>% 
  filter( composition %in% c("All", "F", "B", "S", "L") )

# check how many replicates there were of each
ply_dat %>%
  group_by(shore, month, composition) %>%
  summarise(n = n()) %>%
  View()

# add column for total cover of the four focal species
ply_dat$total_focal_cover <- 
  ply_dat %>%
  select(all_of(foc_spp)) %>%
  rowSums(.)

# reorder columns
names(ply_dat)

ply_sum <- 
  ply_dat %>%
  group_by(shore, month, mono_mix, composition) %>%
  summarise(sargassum_muticum = mean(sargassum_muticum, na.rm = TRUE),
            bifurcaria_bifurcata = mean(bifurcaria_bifurcata, na.rm = TRUE),
            fucus_serratus = mean(fucus_serratus, na.rm = TRUE),
            laminaria_digitata = mean(laminaria_digitata, na.rm = TRUE), .groups = "drop")

# sort out the month variables which don't seem to sort properly
ply_sum$month <- factor(ply_sum$month, levels = c("March", "July", "September"))

ply_mix <- 
  ply_sum %>%
  filter(composition == "All") %>%
  select(-composition, -mono_mix) %>%
  pivot_longer(cols = all_of(foc_spp),
               names_to = "species",
               values_to = "Y") %>%
  rename(place = shore, time = month) %>%
  arrange(place, time, species)

ply_mono <- 
  ply_sum %>% 
  filter(composition != "All") %>%
  pivot_longer(cols = all_of(foc_spp),
               names_to = "species",
               values_to = "M") %>%
  mutate(species_comp = toupper(substr(species, 1, 1)) ) %>%
  filter(composition == species_comp) %>%
  select(place = shore, time = month, species, M) %>%
  arrange(place, time, species)

# join the mixtures and the monocultures
ply_part <- full_join(ply_mono, ply_mix, by = c("time", "place", "species"))

# round off the M and Y columns
ply_part <- 
  ply_part %>%
  mutate(M = round(M, 2),
         Y = round(Y, 2))

# how much cover do these focal species make up?
ply_dat %>%
  mutate(proportion_focal = total_focal_cover/cover_nocrusts) %>%
  pull(proportion_focal) %>%
  summary()

# get summary statistics of species cover
ply_part %>%
  group_by(place, time) %>%
  summarise(Y = sum(Y)) %>% 
  pull(Y) %>%
  mean()

mean(ply_part$M)

# calculate the net biodiversity effect as a test
ply_part %>%
  mutate(RYO = ifelse(M == 0, 0, Y/M)) %>%
  mutate(NBE = (RYO-0.25)*M ) %>%
  pull(NBE) %>%
  sum()

# calculate difference in mean and monoculture across times and places
# this approximates the NBE which gives confidence in the explanation
ply_part %>%
  group_by(place, time) %>%
  summarise(M = mean(M),
            Y = sum(Y)) %>%
  mutate(D = Y-M) %>%
  pull(D) %>%
  sum()

# convert the place and times to integers
ply_part <- 
  ply_part %>%
  mutate(place = as.integer(factor(place)),
         time = as.integer(factor(time)))

# make a sample column
ply_part <- 
  ply_part %>%
  mutate(sample = paste(place, time, sep = "")) %>%
  select(sample, place, time, species, M, Y)

# get random relative expected yields
dr <- sapply(1:100, function(x) gtools::rdirichlet(n = 1, rep(3,  4) ) )
print(dr)

RYe_reps <- 
  apply(dr, 2, function(z) {
    
    a <- Isbell_2018_sampler(data = ply_part, RYe = z, RYe_post = FALSE)
    return(bind_rows(a$Beff, rename(a$L.Beff, Beff = L.Beff)))
    
  } )

df_unc <- bind_rows(RYe_reps, .id = "ID")

# check the data
summary(df_unc)

# summarise the df_unc data
df_unc_sum <- 
  df_unc %>%
  group_by(Beff) %>%
  summarise(Value_m = mean(Value),
            PI_low = PI(Value, prob = 0.90)[1],
            PI_high = PI(Value, prob = 0.90)[2])

# calculate the percentage of total complementarity due to local complementarity
LC <- filter(df_unc_sum, Beff == "LC")[["Value_m"]]
TS <- filter(df_unc_sum, Beff == "TS")[["Value_m"]]
LS <- filter(df_unc_sum, Beff == "LS")[["Value_m"]]

(LC/(LS - TS + LC))*100

# set-up some segments for plot comparisons
segments <- data.frame(xstart = 0,
                       xend = 0.3,
                       effect = 30)

# compare total to local selection and complementarity
eff_in <- c("LC", "TC", "LS", "TS")
p1 <- 
  ggplot(data = df_unc_sum %>%
         filter(Beff %in% eff_in)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_col(mapping = aes(x = Beff, y = Value_m, colour = Beff, fill = Beff),
             width = 0.3) +
  geom_errorbar(mapping = aes(x = Beff, 
                              ymin = PI_low,
                              ymax = PI_high),
                width = 0, colour = "black") +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in) ) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in)) +
  geom_segment(data = segments,
               mapping = aes(x = xstart, xend = xend,
                             y = effect, yend = effect),
               size = 1.25, colour = "red") +
  ylab("Effect (%, cover)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p1)
  
# compare net biodiversity effects, total complementarity and total selection
eff_in <- c("NBE", "TC", "NO", "IT")
p2 <- 
  ggplot(data = df_unc_sum %>%
         filter(Beff %in% eff_in) %>%
         mutate(Beff = factor(Beff, levels = eff_in))
       ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_col(mapping = aes(x = Beff, y = Value_m, colour = Beff, fill = Beff),
           width = 0.3) +
  geom_errorbar(mapping = aes(x = Beff, 
                              ymin = PI_low,
                              ymax = PI_high),
                width = 0, colour = "black") +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in)) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in)) +
  geom_segment(data = segments,
               mapping = aes(x = xstart, xend = xend,
                             y = effect, yend = effect),
               size = 1.25, colour = "red") +
  ylab("") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 1))
plot(p2)

# examine the distribution of the insurance effects
eff_in <- c("AS", "TI", "SI", "ST")
p3 <- 
  ggplot(data = df_unc_sum %>%
         filter(Beff %in% eff_in) %>%
         mutate(Beff = factor(Beff, levels = eff_in))
       ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_col(mapping = aes(x = Beff, y = Value_m, colour = Beff, fill = Beff),
           width = 0.3) +
  geom_errorbar(mapping = aes(x = Beff, 
                              ymin = PI_low,
                              ymax = PI_high),
                width = 0, colour = "black") +
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in)) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in)) +
  geom_segment(data = segments,
               mapping = aes(x = xstart, xend = xend,
                             y = effect, yend = effect),
               size = 1.25, colour = "red") +
  ylab("") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 1))
plot(p3)

# arrange this plot
p123 <- 
  ggarrange(p1, p2, p3, ncol = 3, nrow = 1,
            widths = c(1.1, 1.05, 1),
            labels = c("a", "b", "c"),
            hjust = -0.2,
            vjust = 1,
            font.label = list(size = 11, face = "plain"))
plot(p123)

ggsave(filename = here("figures/fig4.png"), p123,
       unit = "cm", width = 20, height = 7.5)

# why does local complementarity not vary with the RYe
# formula uses average change in relative yield and average relative yield is the same
apply(dr, 2, mean)

# can we explain these effects?

# LC and TC are high whereas LS and TS are low 
# this implies that local complementarity is important

# calculate proportion of species overyielding
OY <- 
  apply(dr, 2, function(x) {
   ply_part$Y - ply_part$M*x
} )

OY <- sapply(OY, function(x) x)

# make a histogram of these numbers
OY_df <- tibble(species = rep(ply_part$species, 100), OY = OY)
OY_df <- arrange(OY_df, species)
summary(OY_df)

# only 1% of species across times and places and integrated over random
# expected relative yields underyielded
sum(OY_df$OY < 0)/nrow(OY_df)


# plotting the raw data mixtures and monocultures
ply_part_M <- 
  ply_part %>%
  mutate(SR = 1) %>%
  rename(Cover = M) %>%
  select(sample, place, time, SR, species, Cover)
  
ply_part_Y <- 
  ply_part %>%
  group_by(sample, place, time) %>%
  summarise(Cover = sum(Y), .groups = "drop") %>%
  mutate(SR = 4,
         species = "mixture") %>%
  select(sample, place, time, SR, species, Cover)

ply_part_sum <- 
  bind_rows(ply_part_M, ply_part_Y) %>%
  arrange(sample, place, time, SR, species, Cover)

# convert to a factor
ply_part_sum$species <- factor(ply_part_sum$species)
levels(ply_part_sum$species) <- list("B.bifurcata" = "bifurcaria_bifurcata",
                                     "F. serratus" = "fucus_serratus",
                                     "L. digitata" = "laminaria_digitata",
                                     "S. muticum" = "sargassum_muticum",
                                     "Mixture (4 sp.)" = "mixture"
                                     )

# make a legend
legend2 <- ply_part_sum[1:5, ]
legend2 <- 
  get_legend( 
    ggplot(data = legend2,
           mapping = aes(x = SR, y = Cover, colour = species, shape = species)) +
      geom_point(size = 2) +
      scale_colour_manual(name = "Species/mixture",
                          labels = c("B.bifurcata", "F. serratus", "L. digitata", "S. muticum", "Mixture (4 sp.)"),
                          values = viridis(n = 5, begin = 0.1, end = 0.9, option = "C")) +   
      scale_shape_manual(name = "Species/mixture",
                         labels = c("B.bifurcata", "F. serratus", "L. digitata", "S. muticum", "Mixture (4 sp.)"),
                         values = c(16, 16, 16, 16, 8)) +
      theme_bw() +
      theme(legend.position = "right")
  )
plot(legend2)

# convert place 1 and 2 into actual places
ply_part_sum$place <- ifelse(ply_part_sum$place == 1, "Challaborough", "Kingsand")
ply_part_sum$time2 <- factor(ply_part_sum$time)
levels(ply_part_sum$time2) <- c("Mar", "Jul", "Sep")

p1 <- 
  ggplot() +
  geom_point(data = ply_part_sum %>% filter(place == "Challaborough"),
             mapping = aes(x = time2, y = Cover, colour = species, shape = species),
             size = 2) +
  geom_line(data = ply_part_sum %>% filter(place == "Challaborough"),
            mapping = aes(x = time, y = Cover, colour = species)) +
  scale_colour_manual(name = "Species/mixture",
                      labels = c("B.bifurcata", "F. serratus", "L. digitata", "S. muticum", "Mixture (4 sp.)"),
                      values = viridis(n = 5, begin = 0.1, end = 0.9, option = "C")) +   
  scale_shape_manual(name = "Species/mixture",
                     labels = c("B.bifurcata", "F. serratus", "L. digitata", "S. muticum", "Mixture (4 sp.)"),
                     values = c(16, 16, 16, 16, 8)) +
  scale_y_continuous(limits = c(0, 45)) +
  ggtitle("Place 1: Challaborough") +
  ylab("Cover (%)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))
plot(p1)

p2 <- 
  ggplot() +
  geom_point(data = ply_part_sum %>% filter(place == "Kingsand"),
             mapping = aes(x = time2, y = Cover, colour = species, shape = species),
             size = 2) +
  geom_line(data = ply_part_sum %>% filter(place == "Kingsand"),
            mapping = aes(x = time, y = Cover, colour = species)) +
  scale_colour_manual(name = "Species/mixture",
                      labels = c("B.bifurcata", "F. serratus", "L. digitata", "S. muticum", "Mixture (4 sp.)"),
                      values = viridis(n = 5, begin = 0.1, end = 0.9, option = "C")) +   
  scale_shape_manual(name = "Species/mixture",
                     labels = c("B.bifurcata", "F. serratus", "L. digitata", "S. muticum", "Mixture (4 sp.)"),
                     values = c(16, 16, 16, 16, 8)) +
  scale_y_continuous(limits = c(0, 45)) +
  ggtitle("Place 2: Kingsand") +
  ylab("") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))
plot(p2)


# check the correlation between relative yield and monoculture yields
p3 <- 
  
  ply_part %>%
  group_by(sample) %>%
  mutate(rel_abun = Y/sum(Y)) %>%
  group_by(place, species) %>%
  summarise(M = mean(M),
            rel_abun = mean(rel_abun), .groups = "drop") %>%
  ggplot(data = .,
         mapping = aes(x = M, y = rel_abun, colour = species)) +
  geom_point(size = 2) +
  ggtitle("Across places") +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  scale_colour_viridis_d(begin = 0.1, end = 0.74, option = "C") +
  ylab("Relative abundance") +
  xlab("Monoculture cover (%)") +
  theme_meta() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))
plot(p3)

p4 <- 
  
  ply_part %>%
  group_by(sample) %>%
  mutate(rel_abun = Y/sum(Y)) %>%
  group_by(time, species) %>%
  summarise(M = mean(M),
            rel_abun = mean(rel_abun)) %>%
  ggplot(data = .,
         mapping = aes(x = M, y = rel_abun, colour = species)) +
  geom_point(size = 2) +
  ggtitle("Across times") +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  scale_colour_viridis_d(begin = 0.1, end = 0.74, option = "C") +
  ylab("") +
  xlab("Monoculture cover (%)") +
  theme_meta()  +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))
plot(p4)

# arrange the plots in a 2 x 2 matrix
p1234 <- 
  ggarrange(p1, p2, p3, p4, 
            heights = c(1, 1.1),
            labels = c("a", "b", "c", "d"),
            font.label = list(size = 11, face = "plain"))

# add the legend
p1234 <- ggarrange(p1234, legend2, widths = c(6, 1))

ggsave(filename = here("figures/fig5.png"), p1234,
       unit = "cm", width = 21, height = 14)



# BES talk figures

# plot of the NBE effect and total complementarity and selection
# compare total to local selection and complementarity
eff_in <- c("NBE", "TC", "TS")
p.BES1 <- 
  ggplot(data = df_unc_sum %>%
           filter(Beff %in% eff_in)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_col(mapping = aes(x = Beff, y = Value_m, colour = Beff, fill = Beff), 
           width = 0.3) +
  geom_errorbar(mapping = aes(x = Beff, 
                              ymin = Value_m - Value_sd,
                              ymax = Value_m + Value_sd),
                width = 0) + 
  scale_colour_manual(values = v_col_BEF(eff_in = eff_in) ) +
  scale_fill_manual(values = v_col_BEF(eff_in = eff_in)) +
  ylab("Effect (%, cover)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p.BES1)

ggsave(filename = here("figures/fig_BES1.png"), p.BES1,
       unit = "cm", width = 7, height = 7)

# make a plot of difference in mixture and monocultures across all times and places
psum <- 
  ply_part %>%
  filter(place == 1, time == 1) %>%
  group_by(place, time) %>%
  summarise(M = mean(M),
            Y = sum(Y), .groups = "drop")

p.BES2 <- 
  ggplot(data = ply_part %>% filter(place == 1, time == 1)) +
  geom_point(mapping = aes(x = 1, y = M),
             position = position_dodge2(width = c(0.1, 0, 0.1, 0))) +
  geom_point(mapping = aes(x = 1, psum$M), colour = "red") +
  geom_point(mapping = aes(x = 4, psum$Y), colour = "red") +
  scale_y_continuous(limits = c(-2, 41)) +
  scale_x_continuous(limits = c(0.8, 4.5)) +
  ylab("Species cover (%)") +
  xlab("Species richness") +
  theme_meta()
plot(p.BES2)

ggsave(filename = here("figures/fig_BES2.png"), p.BES2,
       unit = "cm", width = 7, height = 8)

p.BES3 <- 
  ply_part %>%
  group_by(place, time) %>%
  summarise(M = mean(M),
            Y = sum(Y), .groups = "drop") %>%
  mutate(D = (Y - M)) %>%
  ggplot(data = .) +
  geom_segment(mapping = aes(x = time, xend = time,
                             y = 0, yend = D), colour = "red") +
  scale_x_continuous(breaks = c(1, 2, 3), limits = c(0.75, 3.25)) +
  scale_y_continuous(limits = c(-2.5, 40)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~place) +
  ylab("Cover difference (%)") +
  xlab("") +
  theme_meta()

ggsave(filename = here("figures/fig_BES3.png"), p.BES3,
       unit = "cm", width = 11, height = 8)

# plot a histogram of the overyielding of species
p.BES4 <- 
  ggplot(data = OY_df, 
         mapping = aes(x = OY)) +
  geom_histogram(alpha = 0.3, colour = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.75, colour = "red") +
  # scale_fill_viridis_d(begin = 0.1, end = 0.9, option = "C") +
  ylab("Count") +
  xlab("Overyielding") +
  theme_meta() +
  theme(legend.position = "none")
plot(p.BES4)

ggsave(filename = here("figures/fig_BES4.png"), p.BES4,
       unit = "cm", width = 12, height = 7)

# plot covariance between monoculture and relative abundance
p.BES5 <- 
  ply_part %>%
  group_by(sample) %>%
  mutate(rel_abun = Y/sum(Y)) %>%
  group_by(sample, species) %>%
  summarise(M = mean(M),
            rel_abun = mean(rel_abun)) %>%
  ggplot(data = .,
         mapping = aes(x = M, y = rel_abun, colour = species)) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(begin = 0.1, end = 0.74, option = "C") +
  scale_y_continuous(limits = c(0, 0.85)) +
  ylab("Relative abundance") +
  xlab("Monoculture (%, cover)") +
  theme_meta()  +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))
plot(p.BES5)

ggsave(filename = here("figures/fig_BES5.png"), p.BES5,
       unit = "cm", width = 10, height = 8)

### END
