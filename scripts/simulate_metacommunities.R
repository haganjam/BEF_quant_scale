
# Simulate test data

# Install the Thompson et al. (2020): mcomsimr
devtools::install_github("plthompson/mcomsimr")
library(mcomsimr)

# set-up some core parameters
np <- 10
ns <- 10
ts <- 20
bi <- 50
init <- 50
ext_prob <- 0.0001
disp_rate <- 0.2

# generate a random landscape
l.1 <- landscape_generate(patches = np, plot = FALSE)
head(l.1)

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1,
                        torus = TRUE,
                        kernel_exp = 0.1, plot = FALSE)
head(d.1)

# generate a random environmental matrix
e.1 <- env_generate(landscape = l.1, env1Scale = 500, timesteps = (ts + bi + init), plot = FALSE )
head(e.1)

# generate random species traits
t.1 <- env_traits(species = ns, max_r = 5, 
                  min_env = 0, max_env = 1, 
                  env_niche_breadth = 0.5, plot = FALSE, optima_spacing = "random")
head(t.1)

# generate species interaction matrix
si.1 <- species_int_mat(species = ns, intra = 1, min_inter = 0, max_inter = 1.5, 
                        comp_scaler = 0.05, plot = FALSE)
head(si.1)

# put these data into the simulate_MC function and simulate the mixture
sim.list <- simulate_MC2(species = ns, patches = np, dispersal = disp_rate, plot = TRUE, 
                    timesteps = ts, burn_in = bi, initialization = init,
                    extirp_prob = ext_prob,
                    landscape = l.1, 
                    disp_mat = d.1, 
                    env.df = e.1, 
                    env_traits.df = t.1, 
                    int_mat = si.1)

# write the dynamics part of the simulation into a data.frame
mix <- sim.list$dynamics.df

# reorganise this data.frame
mix <- mix[, c("time", "patch", "species", "N")]
mix <- mix[mix$time >= 0, ]
mix <- 
  mix %>%
  arrange(time, patch, species) %>%
  rename(place = patch, Y = N)

# add a column for each unique sample
sample <- unique(with(mix, paste(time, place)))
mix$sample <- rep(as.integer(as.factor(sample) ), each = ns)

# simulate each monoculture over all times and places
mono <- vector("list", length = ns)
for (i in 1:ns) {
  
  # simulate each species
  x <- simulate_MC2(species = 1, patches = np, dispersal = disp_rate, plot = FALSE, 
                    timesteps = ts, burn_in = bi, initialization = init,
                    extirp_prob = ext_prob,
                    landscape = l.1, 
                    disp_mat = d.1, 
                    env.df = e.1, 
                    env_traits.df = t.1[t.1$species == i,], 
                    int_mat = si.1[i, i])
  y <- x$dynamics.df[, c("time", "patch", "species", "N")]
  y$species <- i
  
  # write the dynamics data.frame into a list
  mono[[i]] <- y[y$time >= 0, ]
  
}

# bind the monocultures into a data.frame
mono <- bind_rows(mono)

# arrange and rename the relevant columns
mono <- 
  mono %>%
  arrange(time, patch, species) %>%
  rename(place = patch, M = N)

# join the mixture and monoculture data to match the partition format
mix.mono <- 
  full_join(mono, mix, by = c("time", "place", "species")) %>%
  select(sample, time, place, species, M, Y)
head(mix.mono)

