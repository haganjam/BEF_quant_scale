
# Explore theoretically whether we can detect an effect of environmental heterogeneity
# on the magnitude of different types of biodiversity effects

# next steps:

# quantifying biodiversity effects

# replication? i.e. only have one rep per env het?

# simulate some empirical data and see how to model monocultures
# use simulations of relative yield expected
# and incorporate all that uncertainty into our estimated biodiversity effects

# adaptive tracking

# simulate a dataset like this and see if we can quantify it (Rudman et al. 2021)

# what other effects can we explore in this field season?

# Install the Thompson et al. (2020): mcomsimr
# devtools::install_github("plthompson/mcomsimr")
library(mcomsimr)

# load key functions
library(here)
source(here("scripts/isbell_2018_partition.R"))
source(here("scripts/mcomsimr_simulate_MC2_function.R"))


# simulate monocultures and mixtures in identical environmental conditions
sim_metacomm_BEF <- function(species = 5, patches = 10, 
                             start_abun = 150,
                             dispersal = 0.2, 
                             timesteps = 10, burn_in = 50, initialization = 50,
                             extirp_prob = 0.0001,
                             landscape, 
                             disp_mat, 
                             env.df, 
                             env_traits.df, 
                             int_mat
                             ) {
  
  # simulate the mixture of species
  mix <- simulate_MC2(species = species, 
                      patches = patches, 
                      dispersal = dispersal, 
                      start_abun = start_abun,
                      timesteps = timesteps,
                      extirp_prob = extirp_prob,
                      landscape = landscape, 
                      disp_mat = disp_mat, 
                      env.df = env.df, 
                      env_traits.df = env_traits.df, 
                      int_mat = int_mat
                      )
  
  # add a mixture column
  mix$mono_mix <- "mixture"
  
  # add a column for each unique sample
  sample <- unique(with(mix, paste(time, patch)))
  mix$sample <- rep(sort(as.integer(as.factor(sample)) ), each = species)
  
  # reorder the columns
  mix <- mix[, c("mono_mix", "sample", "time", "patch", "env", "species", "N")]
  
  # simulate each monoculture over all times and places
  mono <- vector("list", length = species)
  for (i in 1:species) {
    
    # simulate each species
    x <- simulate_MC2(species = 1, patches = patches, 
                      dispersal = dispersal, start_abun = start_abun,
                      timesteps = timesteps,
                      extirp_prob = extirp_prob,
                      landscape = landscape, 
                      disp_mat = disp_mat, 
                      env.df = env.df, 
                      env_traits.df = env_traits.df[i,], 
                      int_mat = int_mat[i,i])
    
    # rename the species column
    x$species <- i
    
    # add a mono_mix variable
    x$mono_mix <- "monoculture"
    
    # add a column for each unique sample
    sample <- unique(with(x, paste(time, patch)))
    x$sample <- sort(as.integer(as.factor(sample)) )
    
    # write the dynamics data.frame into a list
    mono[[i]] <- x
    
  }
  
  # bind the monocultures into a data.frame
  mono <- bind_rows(mono)
  
  # reorder the columns
  # reorder the columns
  mono <- mono[, c("mono_mix", "sample", "time", "patch", "env", "species", "N")]
  
  # return the list with the mixture and monoculture data
  return(list("mixture" = mix, "monoculture" = mono))
  
}

# Set-up simulations

# set-up fixed inputs
species <- 5
patches <- 5
ts <- 100

# Generate a random landscape
l.1 <- landscape_generate(patches = patches, plot = FALSE)

# Generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1,
                        torus = TRUE,
                        kernel_exp = 0.1, plot = FALSE)

# Generate species environmental optima
t.1 <- env_traits(species = species, max_r = 0.5, 
                  min_env = 0, max_env = 1, 
                  env_niche_breadth = 0.2, 
                  plot = FALSE, optima_spacing = "random")
t.1$K_max <- 150
head(t.1)

## Competition matrices

# generate species interaction matrix randomly
si.1 <- species_int_mat(species = 5, intra = 1, min_inter = 0, 
                        max_inter = 0.75, comp_scaler = 1, plot = FALSE)
head(si.1)

## Environmental matrices

# generate a random environmental matrix
e.1 <- env_generate(landscape = l.1, env1Scale = 900, timesteps = (ts), plot = FALSE )
head(e.1)
hist(e.1$env1)

# environment space
env.1 <- runif(n = patches, min = 0, max = 1)
# env.1 <- c(0.1, 0.5, 0.9)
e.s <- expand.grid(1:patches, 1:(ts))
e.s$env1 <- rep(env.1,(ts) ) 
names(e.s) <- c("patch", "time", "env1")
e.s <- e.s[,c("env1", "patch", "time")] 
head(e.s)

# environment time
env.1 <- runif(n = (ts + bi + init), min = 0, max = 1)
e.t <- expand.grid(1:patches, 1:(ts + bi + init))
e.t$env1 <- rep(env.1, each = patches ) 
names(e.t) <- c("patch", "time", "env1")
e.t <- e.t[,c("env1", "patch", "time")] 

# test the new metacommunity simulation function
x <- 
  sim_metacomm_BEF(patches = patches, species = species, dispersal = 0.5, 
                   timesteps = ts, start_abun = 150,
                   extirp_prob = 0.001,
                   landscape = l.1, disp_mat = d.1, env.df = e.s, 
                   env_traits.df = t.1, int_mat = si.1)
head(x)

ggplot(data = x$mixture %>% mutate(species = as.character(species)),
       mapping = aes(x = time, y = N, colour = species)) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()







# transform the monoculture and mixture data and run the Isbell et al. (2018) partition

# arrange and rename the relevant columns
mono <- 
  mono %>%
  arrange(time, patch, species) %>%
  rename(place = patch, M = N)

# join the mixture and monoculture data to match the partition format
mix.mono <- 
  full_join(mono, mix, by = c("time", "place", "species")) %>%
  select(sample, time, place, species, M, Y) %>%
  arrange(sample, time, place, species)
head(mix.mono)

# apply the Isbell partition to these data
part.df <- Isbell_2018_sampler(data = mix.mono, RYe = rep(1/species, species), RYe_post = FALSE)

# convert to the wide format
part.df <- 
  tidyr::pivot_wider(bind_rows(part.df$Beff, select(part.df$L.Beff, Beff = L.Beff, Value)),
                     names_from = "Beff",
                     values_from = "Value")






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
