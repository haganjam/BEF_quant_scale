
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

# set script to call partition functions from
source(here("scripts/01_partition_functions/02_isbell_2018_partition.R"))
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
            Value_sd = sd(Value))

# compare total to local selection and complementarity
p1 <- 
  ggplot(data = df_unc_sum %>%
         filter(Beff %in% eff_in)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_col(mapping = aes(x = Beff, y = Value_m, colour = Beff, fill = Beff), width = 0.5) +
  geom_errorbar(mapping = aes(x = Beff, 
                              ymin = Value_m - Value_sd,
                              ymax = Value_m + Value_sd),
                width = 0) + 
  scale_colour_manual(values = v_col_BEF(eff_in = c("LC", "LS", "TC", "TS")) ) +
  scale_fill_manual(values = v_col_BEF(eff_in = c("LC", "LS", "TC", "TS"))) +
  ylab("Effect (cover (%) time-1)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p1)
  
# compare net biodiversity effects, total complementarity and total selection
p2 <- 
  ggplot(data = df_unc_sum %>%
         filter(Beff %in% eff_in) %>%
         mutate(Beff = factor(Beff, levels = eff_in))
       ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_col(mapping = aes(x = Beff, y = Value_m, colour = Beff, fill = Beff), 
           width = 0.5) +
  geom_errorbar(mapping = aes(x = Beff, 
                              ymin = Value_m - Value_sd,
                              ymax = Value_m + Value_sd),
                width = 0) +
  scale_colour_manual(values = v_col_BEF(eff_in = c("NBE", "TC", "NO", "IT"))) +
  scale_fill_manual(values = v_col_BEF(eff_in = c("NBE", "TC", "NO", "IT"))) +
  ylab("Effect (cover (%) time-1)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")
plot(p2)

# examine the distribution of the insurance effects
p3 <- 
  ggplot(data = df_unc_sum %>%
         filter(Beff %in% eff_in) %>%
         mutate(Beff = factor(Beff, levels = eff_in))
       ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_col(mapping = aes(x = Beff, y = Value_m, colour = Beff, fill = Beff), 
           width = 0.5) +
  geom_errorbar(mapping = aes(x = Beff, 
                              ymin = Value_m - Value_sd,
                              ymax = Value_m + Value_sd),
                width = 0) +
  scale_colour_manual(values = v_col_BEF(eff_in = c("AS", "TI", "SI", "ST"))) +
  scale_fill_manual(values = v_col_BEF(eff_in = c("AS", "TI", "SI", "ST"))) +
  ylab("Effect (cover (%) time-1)") +
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

ggsave(filename = here("figures/ply_fig1.png"), p123,
       unit = "cm", width = 20, height = 7)

# why does local complementarity not vary with the RYe
# formula uses average change in relative yield and average relative yield is the same
apply(dr, 2, mean)

# check the correlation between relative yield and monoculture yields

# make a legend
legend <- ply_part[1:4, ]
legend$Species <- c("B. bifurcata", "F. serratus", "L. digitata", "S. muticum")
legend <- 
  get_legend( 
    ggplot(data = legend,
           mapping = aes(x = M, y = Y, colour = Species)) +
      geom_point(size = 2) +
      scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "C") +
      theme_bw() +
      theme(legend.position = "top")
    )
  
p1 <- 
  
  ply_part %>%
  group_by(sample) %>%
  mutate(rel_abun = Y/sum(Y)) %>%
  group_by(place, species) %>%
  summarise(M = mean(M),
            rel_abun = mean(rel_abun), .groups = "drop") %>%
  ggplot(data = .,
         mapping = aes(x = M, y = rel_abun, colour = species)) +
  geom_point() +
  ggtitle("Across places") +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "C") +
  ylab("Relative abundance") +
  xlab("Monoculture cover (%)") +
  theme_meta() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))
plot(p1)

p2 <- 
  
  ply_part %>%
  group_by(sample) %>%
  mutate(rel_abun = Y/sum(Y)) %>%
  group_by(time, species) %>%
  summarise(M = mean(M),
            rel_abun = mean(rel_abun)) %>%
  ggplot(data = .,
         mapping = aes(x = M, y = rel_abun, colour = species)) +
  geom_point() +
  ggtitle("Across times") +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "C") +
  ylab("") +
  xlab("Monoculture cover (%)") +
  theme_meta()  +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))
plot(p2)

p12 <- 
  ggarrange(p1, p2, legend.grob = legend, labels = c("a", "b"),
            font.label = list(size = 11, face = "plain"))

ggsave(filename = here("figures/ply_fig2.png"), p12,
       unit = "cm", width = 13, height = 8)

### END
