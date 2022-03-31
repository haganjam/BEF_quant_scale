
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
species <- 5
patches <- 5
ts <- 50


## Landscape parameters

# generate a random landscape
l.1 <- 
  data.frame(x = c(50, 25, 25, 75, 75),
             y = c(50, 25, 75, 25, 75))
plot(l.1)

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)


## Species traits

# generate species environmental optima
t.1 <- 
  data.frame(species = 1:species,
             optima = seq(0.1, 0.9, length.out = species),
             env_niche_breadth = 0.2,
             max_r = 0.5,
             K_max = 150)
head(t.1)


## Competition matrices

# generate species interaction matrix randomly

# matrix 1: Intraspecific > interspecific competition
si.1 <- matrix(runif(n = species*species, min = 0.1, max = 1), 
               nrow = species, ncol = species)
si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
diag(si.1) <- 1
head(si.1)

# matrix 2: Neutral competition
si.2 <- species_int_mat(species = species, intra = 1, min_inter = 1, 
                        max_inter = 1, comp_scaler = 1, plot = FALSE)
head(si.2)


## Environmental matrices

# generate a random environmental matrix
e.1 <- env_generate(landscape = l.1, env1Scale = 500, timesteps = (ts), plot = FALSE )

# start with some environmental heterogeneity in each patch and fluctuate around this

# simulate high environmental heterogeneity
env1.high <- 
  sapply(seq(0.1, 0.9, length.out = patches), function(x) {
  z <- 
    sapply(runif(n = ts-1, min = 0, 0.05), function(y) {
    y*sample(x = c(-1,1), size = 1)
  })
  u <- cumsum(c(x, z))
  ifelse( (u < 0), 0, ifelse( (u > 1), 1, u))
})

e.1 <- data.frame(env1 = unlist(apply(env1.high, 1, function(x) list(x) ) ),
                  patch = rep(1:patches, ts),
                  time = rep(1:ts, each = patches))

# simulate low heterogeneity
e.2 <- 
  lapply(seq(0.1, 0.9, length.out =patches), function(y) {
    
    a <- sapply(rep(y, patches), function(x) {
      z <- 
        sapply(runif(n = ts-1, min = 0, 0.025), function(y) {
          y*sample(x = c(-1,1), size = 1)
        })
      u <- cumsum(c(x, z))
      ifelse( (u < 0), 0, ifelse( (u > 1), 1, u))
    })
    
    data.frame(env1 = unlist(apply(a, 1, function(u) list(u) ) ),
               patch = rep(1:patches, ts),
               time = rep(1:ts, each = patches))
    
  }  )


### Explore the behaviour of the model a bit
x <- 
  sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.05, 
                   timesteps = ts, start_abun = 150,
                   extirp_prob = 0,
                   landscape = l.1, disp_mat = d.1, env.df = e.1, 
                   env_traits.df = t.1, int_mat = si.1)

ggplot(data = e.1,
       mapping = aes(x = time, y = env1)) +
  geom_line() +
  facet_wrap(~patch) +
  theme_classic()

x$mixture %>%
  ggplot(data = .,
         mapping = aes(x = time, y = N, colour = as.character(species))) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()

x$monoculture %>%
  ggplot(data = .,
         mapping = aes(x = time, y = N, colour = as.character(species))) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()


x$monoculture %>%
  group_by(patch, species) %>%
  summarise(mean_N = mean(N))

x$monoculture %>%
  group_by(patch) %>%
  summarise(env = range(env))







### Environmental heterogeneity experiment

# write a function to scale a variable between any two values
# https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range
normalise_x <- function(x, b, a) {
  y <- (b - a) * ((x - min(x))/(max(x) - min(x))) + a
  return(y)
}

# set number of replicates
n_rep <- 30

# simulate a set of metacommunities with high environmental heterogeneity among patches
high_out <- vector("list", n_rep)
high_env_var <- vector("list", length = n_rep)
for (i in 1:n_rep) {
  
  x <- 
    sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.05, 
                     timesteps = ts, start_abun = 150,
                     extirp_prob = 0,
                     landscape = l.1, disp_mat = d.1, env.df = e.1, 
                     env_traits.df = t.1, int_mat = si.2)
  
  high_out[[i]] <- x
  
  x <- Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, 
                           t_sel  = c(50, 100))
  
  x$env_het <- "high"
  
  high_env_var[[i]] <- x
  
}

# simulate a set of metacommunities with low environmental heterogeneity among patches
low_out <- vector("list", n_rep)
low_env_var <- vector("list", length = n_rep)
for (i in 1:n_rep) {

x <- 
  sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.05, 
                   timesteps = ts, start_abun = 150,
                   extirp_prob = 0,
                   landscape = l.1, disp_mat = d.1, env.df = e.2[[sample(1:5, 1)]], 
                   env_traits.df = t.1, int_mat = si.2)
low_out[[i]] <- x

x <- Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, 
                         t_sel  = c(50, 100))

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
  pivot_longer(cols = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST", "LC", "LS"),
               names_to = "BEFF",
               values_to = "Value") %>%
ggplot(data = .,
       mapping = aes(x = BEFF, y = Value, colour = env_het)) +
  geom_jitter(alpha = 0.1, width = 0.1) +
  geom_boxplot(width = 0.1) +
  theme_classic()

df %>%
  filter(NO == min(NO))

t.1

high_out[[43]]$mixture %>%
  ggplot(data = .,
         mapping = aes(x = time, y = N, colour = as.character(species))) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()

high_out[[43]]$monoculture %>%
  group_by(patch, species) %>%
  summarise(mean_N = mean(N))

high_out[[43]]$monoculture %>%
  group_by(patch) %>%
  summarise(env = range(env))

high_out[[43]]$monoculture %>%
  ggplot(data = .,
         mapping = aes(x = time, y = N, colour = as.character(species))) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()

full_join(bind_rows(lapply(high_out, function(x) x$monoculture %>% select(-mono_mix) )),
          bind_rows(lapply(high_out, function(x) { 
            x$mixture %>%
              group_by(time, patch) %>%
              mutate(RA = N/sum(N)) %>%
              select(sample, time, patch, env, species, RA)} ))) %>% 
  filter(time %in% c(50, 100)) %>%
  ggplot(data = .,
         mapping = aes(x = N, y = RA, colour = as.character(species) )) +
  geom_point() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()
  
full_join(bind_rows(lapply(low_out, function(x) x$monoculture %>% select(-mono_mix) )),
          bind_rows(lapply(low_out, function(x) { 
            x$mixture %>%
              group_by(time, patch) %>%
              mutate(RA = N/sum(N)) %>%
              select(sample, time, patch, env, species, RA)} ))) %>% 
  filter(time %in% c(50, 100)) %>%
  ggplot(data = .,
         mapping = aes(x = N, y = RA, colour = as.character(species) )) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~patch, scales = "free") +
  theme_classic()




ggplot(data = df,
       mapping = aes(x = env_het, y = IT)) +
  geom_point() +
  theme_classic()

# what about the proportion of the total biodiversity effect
df %>%
  mutate(AS = AS/NBE) %>%
  ggplot(data = .,
       mapping = aes(x = env_het, y = AS)) +
  geom_point() +
  theme_classic()

df %>%
  mutate(LS = LS/NBE) %>%
  ggplot(data = .,
         mapping = aes(x = env_het, y = LS)) +
  geom_point() +
  theme_classic()

df %>%
  mutate(SI = SI/NBE) %>%
  ggplot(data = .,
         mapping = aes(x = env_het, y = SI)) +
  geom_point() +
  theme_classic()

df %>%
  mutate(TI = TI/NBE) %>%
  ggplot(data = .,
         mapping = aes(x = env_het, y = TI)) +
  geom_point() +
  theme_classic()

df %>%
  mutate(ST = ST/NBE) %>%
  ggplot(data = .,
         mapping = aes(x = env_het, y = ST)) +
  geom_point() +
  theme_classic()
