
## Calculate Isbell et al.'s (2018) effect on simulated data

# Load the functions
library(here)
source(here("scripts/mcomsimr_simulate_MC2_function.R"))

# set the seed
set.seed(5489485)

## Set-up simulation inputs

## set-up fixed inputs
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
             optima = c(0.1, 0.3, 0.5, 0.7, 0.9),
             env_niche_breadth = 0.2,
             max_r = 0.5,
             K_max = 150)
head(t.1)


## Competition matrices

# generate species interaction matrix randomly

# matrix 1: Intraspecific > interspecific competition
si.1 <- species_int_mat(species = species, intra = 1, min_inter = 0, 
                        max_inter = 0.75, comp_scaler = 1, plot = FALSE)
head(si.1)

# matrix 2: Neutral competition
si.2 <- species_int_mat(species = species, intra = 1, min_inter = 1, 
                        max_inter = 1, comp_scaler = 1, plot = FALSE)
head(si.2)


## Environmental matrices

# generate a random environmental matrix
e.1 <- env_generate(landscape = l.1, env1Scale = 900, timesteps = (ts), plot = FALSE )

# only spatial variation
e.2 <- expand.grid(1:patches, 1:(ts))
e.2$env1 <- rep(c(0.1, 0.3, 0.5, 0.7, 0.9), ts ) 
names(e.2) <- c("patch", "time", "env1")
e.2 <- e.2[,c("env1", "patch", "time")] 
head(e.2)

# only temporal variation
e.3 <- expand.grid(1:patches, 1:ts)
e.3$env1 <- rep(sample(x = seq(0.1, 0.9, 0.1), size = ts, replace = TRUE), 
                each = patches ) 
names(e.3) <- c("patch", "time", "env1")
e.3 <- e.3[,c("env1", "patch", "time")] 


## Test the metacommunity simulation function
x <- 
  sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.1, 
                   timesteps = ts, start_abun = 150,
                   extirp_prob = 0.001,
                   landscape = l.1, disp_mat = d.1, env.df = e.2, 
                   env_traits.df = t.1, int_mat = si.2)

ggplot(data = x$mixture %>% mutate(species = as.character(species)),
       mapping = aes(x = time, y = N, colour = species)) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()

Isbell_2018_cleaner(mix = x$mixture, mono = x$monoculture, from_last = 10)









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
