
# Simulate a test of our hypothesis

# r versus K, should they be positively correlated? Probably not...

# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))


# initialise outputs
## Set-up fixed inputs
species <- 5
patches <- 5
ts <- 100


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

### Environmental heterogeneity experiment

# write a function to scale a variable between any two values
# https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range
normalise_x <- function(x, b, a) {
  y <- (b - a) * ((x - min(x))/(max(x) - min(x))) + a
  return(y)
}

# set number of replicates: 10
n_rep <- 10

# simulate a set of metacommunities with high environmental heterogeneity among patches
high_env_var <- vector("list", length = n_rep)
for (i in 1:n_rep) {
  
  x <- 
    sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.05, 
                     timesteps = ts, start_abun = 100,
                     extirp_prob = 0.00001,
                     landscape = l.1, disp_mat = d.1, env.df = e.1, 
                     env_traits.df = t.1, int_mat = si.2)
  
  x <- Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, 
                           t_sel  = seq(10, ts, 20))
  
  x$env_het <- "high"
  
  high_env_var[[i]] <- x
  
}

# simulate a set of metacommunities with low environmental heterogeneity among patches
low_env_var <- vector("list", length = n_rep)
for (i in 1:n_rep) {

# modify the environmental parameters
e.1.mod <- e.1

# get the minimum and maximum values
y <- runif(1, 0, 1)
if ( (1-y) > 0.5 ) { 
  a <- y
  b <- y + 0.5
}
if ( (1-y) < 0.5 ) { 
  a <- y - 0.5
  b <- y
}

# squish the environmental heterogeneity
e.1.mod$env1 <- normalise_x(x = e.1.mod$env1, a = a, b = b)

x <- 
  sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.05, 
                   timesteps = ts, start_abun = 100,
                   extirp_prob = 0.00001,
                   landscape = l.1, disp_mat = d.1, env.df = e.1.mod, 
                   env_traits.df = t.1, int_mat = si.2)

x <- Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, 
                         t_sel  = seq(10, ts, 20))

x$env_het <- "low"

low_env_var[[i]] <- x

}

# bind these data.frames together
df <- bind_rows( bind_rows(high_env_var), bind_rows(low_env_var), .id = "run" )

# plot the high versus low environmental heterogeneity comparison
ggplot(data = df,
       mapping = aes(x = env_het, y = NBE)) +
  geom_point() +
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
