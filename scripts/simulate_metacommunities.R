
# Simulate a test of our hypothesis

# r versus K, should they be positively correlated? Probably not...

# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))

## initialise outputs

## Set-up fixed inputs
species <- 3
patches <- 3
ts <- 50


## Landscape parameters

# generate a random landscape
l.1 <- 
  data.frame(x = c(25, 50, 75),
             y = c(50, 50, 50))
plot(l.1)

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)


## Species traits

# generate species environmental optima
t.1 <- 
  data.frame(species = 1:species,
             optima = seq(0.25, 0.75, length.out = species),
             env_niche_breadth = 0.3,
             max_r = 0.5,
             K_max = 150)
head(t.1)


## Competition matrices

# generate species interaction matrix randomly

# matrix 1: Intraspecific > interspecific competition
si.1 <- matrix(runif(n = species*species, min = 0.1, max = 0.75), 
               nrow = species, ncol = species)
si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
diag(si.1) <- 1
head(si.1)

# matrix 2: Neutral competition
si.2 <- species_int_mat(species = species, intra = 1, min_inter = 1, 
                        max_inter = 1, comp_scaler = 1, plot = FALSE)
head(si.2)

# simulate a practice environmental matrix
e.test <- env_generate(landscape = l.1, timesteps = ts, plot = FALSE)
e.test <- e.1
e.test <- e.2

### Explore the behaviour of the model a bit
x <- 
  sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.05, 
                   timesteps = ts, start_abun = 200,
                   extirp_prob = 0,
                   landscape = l.1, disp_mat = d.1, env.df = e.1, 
                   env_traits.df = t.1, int_mat = si.2)

# what does the environmental variation look like?
ggplot(data = e.test,
       mapping = aes(x = time, y = env1)) +
  geom_line() +
  scale_y_continuous(limits = c(0.2, 1)) +
  facet_wrap(~patch) +
  theme_classic()

t.1

# mixtures?
x$mixture %>%
  ggplot(data = .,
         mapping = aes(x = time, y = N, colour = as.character(species))) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()

# monocultures?
x$monoculture %>%
  ggplot(data = .,
         mapping = aes(x = time, y = N, colour = as.character(species))) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()

# do different species dominate in different patches?
x$monoculture %>%
  group_by(time, patch) %>%
  summarise(mean_N = mean(N))

x$mixture %>%
  group_by(time, patch) %>%
  summarise(sum_N = sum(N))


# looks like spatial insurance and local complementarity should be present
mix <- x$mixture
mono <- x$monoculture

t_sel <- round(seq(5, ts, length.out = 5), 0)

# process the mixture data
mix <- 
  mix %>%
  filter( time %in% t_sel ) %>%
  select(-mono_mix, -env) %>%
  rename(place = patch, Y = N) %>%
  mutate(sample = as.integer(as.factor(sample)))
head(mix)

# process the monoculture data
mono <- 
  mono %>%
  filter( time %in% t_sel ) %>%
  select(-mono_mix, -env, -sample) %>%
  rename(place = patch, M = N)
head(mono)

# join the mixture and monoculture data to match the partition format
mix.mono <- 
  full_join(mono, mix, by = c("time", "place", "species")) %>%
  select(sample, time, place, species, M, Y) %>%
  arrange(sample, time, place, species)
head(mix.mono)

# apply the Isbell partition to these data
n_species <- length(unique(mix.mono$species))
part.df <- Isbell_2018_sampler(data = mix.mono, RYe = rep(1/n_species, n_species), RYe_post = FALSE)

ggplot(data = part.df$Beff,
       mapping = aes(x = Beff, y = Value)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

ggplot(data = part.df$L.Beff,
       mapping = aes(x = L.Beff, y = Value)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()


### Environmental heterogeneity experiment

# set number of replicates
n_rep <- 1000

# simulate a set of metacommunities with high environmental heterogeneity among patches
high_out <- vector("list", n_rep)
high_env_var <- vector("list", length = n_rep)
for (i in 1:n_rep) {
  
  # simulate an environment with high heterogeneity
  env1.high <- 
    sapply(seq(0.2, 0.8, length.out = patches), function(x) {
      z <- 
        sapply(runif(n = ts-1, min = 0, 0.025), function(y) {
          y*sample(x = c(-1,1), size = 1)
        })
      u <- cumsum(c(x, z))
      ifelse( (u < 0), 0, ifelse( (u > 1), 1, u))
    })
  
  e.1 <- data.frame(env1 = unlist(apply(env1.high, 1, function(x) list(x) ) ),
                    patch = rep(1:patches, ts),
                    time = rep(1:ts, each = patches))
  
  x <- 
    sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.05, 
                     timesteps = ts, start_abun = 200,
                     extirp_prob = 0,
                     landscape = l.1, disp_mat = d.1, env.df = e.1, 
                     env_traits.df = t.1, int_mat = si.2)
  
  high_out[[i]] <- x
  
  x <- Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, 
                           t_sel  = round(seq(5, ts, length.out = 5), 0))
  
  x$env_het <- "high"
  
  high_env_var[[i]] <- x
  
}

# simulate a set of metacommunities with low environmental heterogeneity among patches
low_out <- vector("list", n_rep)
low_env_var <- vector("list", length = n_rep)
low_var <- seq(0.2, 0.8, length.out = n_rep)
for (i in 1:n_rep) {

  # simulate low heterogeneity
  e.2 <- 
    sapply(rep(low_var[i], patches), function(x) {
      z <- 
        sapply(runif(n = ts-1, min = 0, 0.025), function(y) {
          y*sample(x = c(-1,1), size = 1)
        })
      u <- cumsum(c(x, z))
      ifelse( (u < 0), 0, ifelse( (u > 1), 1, u))
    })
  
  e.2 <- 
    data.frame(env1 = unlist(apply(e.2, 1, function(x) list(x) ) ),
               patch = rep(1:patches, ts),
               time = rep(1:ts, each = patches))
  
  x <- 
    sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.05, 
                     timesteps = ts, start_abun = 200,
                     extirp_prob = 0,
                     landscape = l.1, disp_mat = d.1, env.df = e.2, 
                     env_traits.df = t.1, int_mat = si.2)
  low_out[[i]] <- x
  
  x <- Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, 
                           t_sel  = round(seq(5, ts, length.out = 5), 0) )
  
  x$env_het <- "low"
  
  low_env_var[[i]] <- x

}

# bind these data.frames together
df <- bind_rows( bind_rows(high_env_var, .id = "run"), bind_rows(low_env_var, .id = "run"))
# readr::write_csv(x = df, here("data/het_experiment_data.csv"))

# is the environmental heterogeneity greater in the high treatment?
sapply(high_out, function(x) {
  sd(x$mixture$env)/mean(x$mixture$env)
} ) %>%
  mean(.)

sapply(low_out, function(x) {
  sd(x$mixture$env)/mean(x$mixture$env)
} ) %>%
  mean(.)

# is there a difference in the mean environment?
sapply(high_out, function(x) {
  mean(x$mixture$env)
} ) %>%
  mean(.)

sapply(low_out, function(x) {
  mean(x$mixture$env)
} ) %>%
  mean(.)

# plot the high versus low environmental heterogeneity comparison
ggplot(data = df,
       mapping = aes(x = env_het, y = NBE)) +
  geom_jitter(width = 0.2, alpha = 0.05) +
  geom_boxplot(width = 0.1) +
  theme_classic()

df %>%
  group_by(env_het) %>%
  summarise(mean_NBE = mean(NBE))

df_plot <- 
  df %>%
  pivot_longer(cols = c("NBE", "TC", "TS", "NO","AS", "IT", "TI", "SI", "ST"),
               names_to = "BEFF",
               values_to = "Value")

df_plot$BEFF <- factor(df_plot$BEFF,
                       levels =c("NBE", "TC", "TS", "NO","AS", "IT", "TI", "SI", "ST") )

ggplot(data = df_plot,
       mapping = aes(x = BEFF, y = Value, colour = env_het)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, position=position_dodge(1), outlier.shape = NA) +
  theme_classic() +
  scale_y_continuous(limits = c(-250, 550)) +
  theme(legend.position = "bottom")

df %>%
  filter(NO < -100)

high_out[[49]]$mixture %>%
  ggplot(data = .,
         mapping = aes(x = time, y = N, colour = as.character(species) )) +
  geom_line() +
  facet_wrap(~patch) +
  theme_classic()

high_out[[49]]$monoculture %>%
  ggplot(data = .,
         mapping = aes(x = time, y = N, colour = as.character(species) )) +
  geom_line() +
  facet_wrap(~patch) +
  theme_classic()

