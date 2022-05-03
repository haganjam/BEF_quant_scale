
# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
library(pbapply)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))

# set the model inputs

# set a seed
set.seed(54258748)

# fixed parameters
max_r = 0.5
K_max = 150

species <- 5
patches <- 10
timesteps <- 200
extirp_prob = 0.00001

# landscape parameters
# generate a random landscape
# l.1 <- landscape_generate(patches = patches)
l.1 <- data.frame(x = seq(0, 100, length.out = patches),
           y = 50)
plot(l.1)

# generate a random dispersal matrix
# d.1 <- dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)
d.1 <- matrix(runif(n = patches*patches, min = 0.25, max = 0.25), 
               nrow = patches, ncol = patches)
d.1[lower.tri(d.1)] = t(d.1)[lower.tri(d.1)]
diag(d.1) <- 0

# generate species environmental optima
optima = seq(0.1, 0.9, length.out = species)
print(optima)
env_niche_breadth = rep(0.3, species)

# get the starting abundances
start_abun <- rep(30, species)

t.1 <- 
  data.frame(species = 1:species,
             optima = optima,
             env_niche_breadth = env_niche_breadth[[1]],
             max_r = max_r,
             K_max = K_max)

# simulate a fluctuating environment
temp_fluc <- 0.001
spat_het <- seq(0.1, 0.9, length.out = patches)
env1.sim <- 
  sapply(spat_het, function(x) {
  # choose a random number between 0 and 0.025 and then decide if it is negative or positive
  y <- runif(n = timesteps-1, 0, temp_fluc)*sample(x = c(-1,1), size = timesteps-1, replace = TRUE)
  # cumulatively sum the environmental fluctuations
  z <- cumsum(c(x, y))
  # make sure they remain between zero and one
  return(ifelse( (z < 0), 0, ifelse( (z > 1), 1, z)))
})
# bind into a data.frame
e.1 <- data.frame(env1 = unlist(apply(env1.sim, 1, function(x) list(x) ) ),
                  patch = rep(1:patches, timesteps),
                  time = rep(1:timesteps, each = patches))

ggplot(data = e.1,
       mapping = aes(x = time, y = env1, colour = as.character(patch) )) +
  geom_line() +
  theme_classic() +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  theme(legend.position = "bottom")

int_min = 1
int_max = 1
intra = 1

# competition matrices
si.1 <- matrix(runif(n = species*species, min = int_min, max = int_max), 
               nrow = species, ncol = species)
si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
diag(si.1) <- intra
si.1

dispersal <- 0.1

# simulate the metacommunity
MC1a <- sim_metacomm_BEF(patches = patches, species = species, dispersal = dispersal,
                         timesteps = timesteps, start_abun = start_abun,
                         extirp_prob = extirp_prob,
                         landscape = l.1, disp_mat = d.1, env.df = e.1, 
                         env_traits.df = t.1, int_mat = si.1,
                         meas_error = 1)

ggplot(data = MC1a$mixture %>% mutate(species = as.character(species)),
       mapping = aes(x = time, y = N, colour = species)) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  #geom_hline(yintercept = 150, linetype = "dashed") +
  theme_classic()

reg.scale <- 
  bind_rows(
    MC1a$mixture %>%
      filter(time == max(time)) %>%
      summarise(S = length( unique(species[N>0]) ),
                N = sum(N)) %>%
      ungroup(),
    MC1a$monoculture %>%
      filter(time == max(time)) %>%
      group_by(species) %>%
      summarise(S = length( unique(species[N>0]) ),
                N = sum(N)) %>%
      ungroup() %>%
      select(-species)
  )

loc.scale <- 
  bind_rows( 
  MC1a$mixture %>%
    filter(time == max(time)) %>%
    group_by(sample, time, patch) %>%
    summarise(S = sum(N > 0),
              N = sum(N)) %>%
    ungroup() %>%
    mutate(species = "mix") %>%
    select(sample, time, patch, species, S, N), 
  MC1a$monoculture %>%
    filter(time == max(time)) %>%
    group_by(sample, time, patch, species) %>%
    summarise(S = sum(N > 0),
              N = sum(N)) %>%
    ungroup() %>%
    filter(S > 0) %>%
    mutate(species = as.character(species))
  ) %>%
  arrange(sample, time, patch, species)

# local scale relationship
ggplot(data = loc.scale,
       mapping = aes(x = S, y = N)) +
  geom_jitter(width = 0.1) + 
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  theme_classic()

lm(loc.scale$N/max(loc.scale$N) ~ loc.scale$S)

# patches aggregated together
ggplot(data = reg.scale,
         mapping = aes(x = S, y = N)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  theme_classic()

lm(reg.scale$N/max(reg.scale$N) ~ reg.scale$S)


