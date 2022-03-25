
# Simulate a test of our hypothesis

# Load the functions and packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))



### Environmental heterogeneity experiment

# write a function to scale a variable between any two values
# https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range
normalise_x <- function(x, b, a) {
  y <- (b - a) * ((x - min(x))/(max(x) - min(x))) + a
  return(y)
}

# set number of replicates: 10
n_rep <- 100

# simulate a set of metacommunities with high environmental heterogeneity among patches
high_env_var <- vector("list", length = n_rep)
for (i in 1:n_rep) {
  
  x <- 
    sim_metacomm_BEF(species = species, patches = patches, 
                     dispersal = 0.2, 
                     timesteps = ts, burn_in = bi, initialization = init,
                     extirp_prob = 0.0001,
                     landscape = l.1, 
                     disp_mat = d.1, 
                     env.df = e.1, 
                     env_traits.df = t.1, 
                     int_mat = si.1
                     )
  
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
  sim_metacomm_BEF(species = species, patches = patches, 
                   dispersal = 0.2, 
                   timesteps = ts, burn_in = bi, initialization = init,
                   extirp_prob = 0.0001,
                   landscape = l.1, 
                   disp_mat = d.1, 
                   env.df = e.1.mod, 
                   env_traits.df = t.1, 
                   int_mat = si.1
  )

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
