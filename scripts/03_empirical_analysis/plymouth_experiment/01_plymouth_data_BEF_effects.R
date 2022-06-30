
# Plymouth rock-pool data

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# set script to call partition functions from
source(here("scripts/isbell_2018_partition.R"))

# read in the data
ply_dat <- read_delim( here("data/Plymouth_data.csv"), delim = "," )

# check the data
View(ply_dat)
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
  select( all_of(id_vars), grid_squares, species_richness, totcover, cover_nocrusts, all_of(foc_spp) )

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
  pivot_longer(cols = foc_spp,
               names_to = "species",
               values_to = "Y") %>%
  rename(place = shore, time = month) %>%
  arrange(place, time, species)


ply_mono <- 
  ply_sum %>% 
  filter(composition != "All") %>%
  pivot_longer(cols = foc_spp,
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

# plot the data
ggplot() +
  geom_jitter(data = df_unc,
               mapping = aes(x = Beff, y = Value, colour = Beff, fill = Beff),
               alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(legend.position = "none")

# why does local complementarity vary with the RYe
# formula uses average change in relative yield and average relative yield is the same
apply(dr, 2, mean)

# check the correlation between relative yield and monoculture yields
ply_part %>%
  group_by(sample) %>%
  mutate(rel_abun = Y/sum(Y)) %>%
  group_by(place, species) %>%
  summarise(M = mean(M),
            rel_abun = mean(rel_abun)) %>%
  ggplot(data = .,
         mapping = aes(x = M, y = rel_abun, colour = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

ply_part %>%
  group_by(sample) %>%
  mutate(rel_abun = Y/sum(Y)) %>%
  group_by(time, species) %>%
  summarise(M = mean(M),
            rel_abun = mean(rel_abun)) %>%
  ggplot(data = .,
         mapping = aes(x = M, y = rel_abun, colour = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()


