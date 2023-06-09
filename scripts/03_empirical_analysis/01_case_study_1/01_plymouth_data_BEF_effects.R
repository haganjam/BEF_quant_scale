#'
#' @title: Analyse the Plymouth experiment (case study 2)
#' 
#' @description: This scripts cleans the data, performs the analysis and 
#' exports the figures for case study 1: Figures 5 and 6.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

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
source("scripts/01_partition_functions/01_isbell_2018_partition.R")
source("scripts/Function_plotting_theme.R")
source("scripts/03_empirical_analysis/helper_functions.R")

# read in the data
ply_dat <- read_delim( "data/case_study_1/Plymouth_data.csv", delim = "," )

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

# set the number of RYe replicates
Ns <- 100

start_RA <- vector("list", length = Ns)
for(i in 1:Ns) {
  
  # get the number of samples required
  N_RYE <- length(unique(ply_part$sample))
  
  # for each sample, get a simplex from the Dirichlet distribution
  rye_mat <- gtools::rdirichlet(n = N_RYE, rep(3, length(unique(ply_part$species))))
  
  start_RA[[i]] <- rye_mat
  
}

# iterate over all possible initial RYE values
RYE_rep <- 
  
  lapply(start_RA, function(RA) {
    
    RYE <- vector("list", length = nrow(RA))
    for(j in 1:nrow(RA)) { RYE[[j]] <- RA[j,] }
    
    # calculate the BEF effects
    BEF_post <- Isbell_2018_part(data = ply_part, RYe = RYE)
    names(BEF_post[["L.Beff"]])[names(BEF_post[["L.Beff"]]) == "L.Beff"] <- "Beff"
    
    # combine the general biodiversity effects and local effects into one data.frame
    BEF_post <- rbind(BEF_post[["Beff"]], BEF_post[["L.Beff"]])
    
    # convert to a data.frame
    BEF_post <- as.data.frame(BEF_post, row.names = NULL)
    
    return(BEF_post)
    
  })

# bind into a data.frame
BEF_dat <- as_tibble(bind_rows(RYE_rep, .id = "RYE"))

# remove the LS and LC effects
BEF_plot <- 
  BEF_dat %>%
  filter( !(Beff %in% c("LC", "LS")) )

# make a Beff column a factor
BEF_plot$Beff <- factor(BEF_plot$Beff)
 
# calculate the summary dataset
BEF_sum <- 
  BEF_plot %>%
  group_by(Beff) %>%
  summarise(Value_m = mean(Value),
            Value_sd = sd(Value),
            n = n(),
            HPDI_low = HPDI(Value, 0.95)[1],
            HPDI_high = HPDI(Value, 0.95)[2], .groups = "drop")

# replace the TC and NO effects with zeros
BEF_plot <- 
  BEF_plot %>%
  mutate(Value = ifelse(Beff %in% c("TC", "NO"), NA, Value))

# plot the different biodiversity in different clusters

# get a nice colour palette
col_pal <- wesanderson::wes_palette("Darjeeling1", n = 9, type = "continuous")

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
    BEF_plot%>% 
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
               mapping = aes(x = Beff, y = Value, colour = Beff), 
               position = position_jitterdodge(jitter.width = 0.4,
                                               dodge.width = 0),
               alpha = 0.25, shape = 1, size = 1.2, stroke = 0.30) +
    geom_errorbar(data = x_sum,
                  mapping = aes(x = Beff, ymin = HPDI_low, ymax = HPDI_high), width = 0.05, 
                  linewidth = 0.3, alpha = 1, colour = "black") +
    geom_point(data = x_sum,
               mapping = aes(x = Beff, y = Value_m, fill = Beff), shape = 23, size = 2.2,
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
ggsave(filename = "figures/fig_3i.png", plot_list[[1]],
       dpi = 500, units = "cm", width = 10, height = 6)
ggsave(filename = "figures/fig_3ii.png", plot_list[[2]],
       dpi = 500, units = "cm", width = 10, height = 4)
ggsave(filename = "figures/fig_3iii.png", plot_list[[3]],
       dpi = 500, units = "cm", width = 10, height = 8)


# calculate the percentage of total complementarity due to local complementarity
# LC <- filter(BEF_sum, Beff == "LC")[["Value_m"]]
# TS <- filter(BEF_sum, Beff == "TS")[["Value_m"]]
# LS <- filter(BEF_sum, Beff == "LS")[["Value_m"]]

# (LC/(LS - TS + LC))*100

# why does local complementarity not vary with the RYe?
# formula uses average change in relative yield and average relative yield is the same
# apply(dr, 2, mean)

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

### END
