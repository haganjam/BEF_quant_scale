
# Simulate test data

# Install the Thompson et al. (2020): mcomsimr
devtools::install_github("plthompson/mcomsimr")
library(mcomsimr)

# load key functions
library(here)
source(here("scripts/isbell_2018_partition.R"))
source(here("scripts/mcomsimr_simulate_MC2_function.R"))

# set-up some core parameters
np <- 5
ns <- 3
ts <- 10
bi <- 10
init <- 10
ext_prob <- 0.0001
disp_rate <- 0.3

# set-up fixed inputs

# generate a random landscape
l.1 <- landscape_generate(patches = np, plot = FALSE)

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1,
                        torus = TRUE,
                        kernel_exp = 0.1, plot = FALSE)

# generate species interaction matrix randomly
si.1 <- species_int_mat(species = ns, intra = 1, min_inter = 0, max_inter = 1, 
                        comp_scaler = 1, plot = FALSE)

# or do it manually
si.1 <- matrix(rnorm(n = ns*ns, mean = 0.5, sd = 0.25), nrow = ns, ncol = ns)
diag(si.1) <- 1
head(si.1)

# generate a random environmental matrix
e.1 <- env_generate(landscape = l.1, env1Scale = 500, timesteps = (ts + bi + init), plot = FALSE )
head(e.1)
summary(e.1)

# write a function to scale a variable between any two values
# https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range
normalise_x <- function(x, b, a) {
  y <- (b - a) * ((x - min(x))/(max(x) - min(e.1$env1))) + a
  return(y)
}

# generate random species traits
t.1 <- env_traits(species = ns, max_r = 5, 
                  min_env = 0, max_env = 1, 
                  env_niche_breadth = 0.5, 
                  plot = FALSE, optima_spacing = "even")


# choose number of runs and loop over this number
n_runs <- 2
n_sim_out <- vector("list", length = n_runs)

for (j in 1:n_runs) {
  
  # put these data into the simulate_MC function and simulate the mixture
  sim.list <- simulate_MC2(species = ns, patches = np, dispersal = disp_rate, plot = FALSE, 
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
  mix <- mix[mix$time %in% seq(0, ts, 10),]
  mix <- 
    mix %>%
    arrange(time, patch, species) %>%
    rename(place = patch, Y = N)
  
  # add a column for each unique sample
  sample <- unique(with(mix, paste(time, place)))
  mix$sample <- rep(sort(as.integer(as.factor(sample)) ), each = ns)
  
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
    y <- y[y$time >= 0, ]
    y <- y[y$time %in% seq(0, ts, 10),]
    mono[[i]] <- y
    
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
    select(sample, time, place, species, M, Y) %>%
    arrange(sample, time, place, species)
  head(mix.mono)
  
  # quantify environmental heterogeneity
  env.het <- sim.list$env.df
  env.het <- env.het[env.het$time_run >= 0, ]
  
  # spatial heterogeneity
  sp.h <- mean( aggregate(env.het$env1, list(time = env.het$time_run), function(x) sd(x)/mean(x) )$x ) 
  
  # temporal heterogeneity
  t.h <- aggregate(env.het$env1, list(time = env.het$time), function(x) mean(x) )
  t.h <- sd(t.h$x)/mean(t.h$x)
  
  # total heterogeneity
  tot.h <- sd(env.het$env1)/mean(env.het$env1)
  
  # trait variation in environmental optima
  
  # env.optima
  optima.cv <- sd(sim.list$env_traits.df$optima)/mean(sim.list$env_traits.df$optima)
  
  # max.r
  maxr.cv <- sd(sim.list$env_traits.df$max_r)/mean(sim.list$env_traits.df$max_r)
  
  # niche.breadth.cv
  niche.cv <- sd(sim.list$env_traits.df$env_niche_breadth)/mean(sim.list$env_traits.df$env_niche_breadth)
  
  # apply the Isbell partition to these data
  part.df <- Isbell_2018_sampler(data = mix.mono, RYe = rep(1/ns, ns), RYe_post = FALSE)
  
  # convert to the wide format
  part.df <- 
    tidyr::pivot_wider(bind_rows(part.df$Beff, select(part.df$L.Beff, Beff = L.Beff, Value)),
                       names_from = "Beff",
                       values_from = "Value")
  
  # create a data.frame of predictor variables
  preds <- 
    tibble(species = ns,
           places = np,
           times = ts,
           spatial_het = sp.h,
           temporal_het = t.h,
           total_het = tot.h,
           env_optima_cv = optima.cv,
           max_r_cv = maxr.cv,
           niche.cv = niche.cv)
  
  # join the part.df and preds dataframes
  n_sim_out[[j]] <-  bind_cols(preds, part.df)
  
}

# bind the simulated data into a data.frame
n_sim_out <- bind_rows(n_sim_out, .id = "run_id")
names(n_sim_out)

# plot a few comparisons
ggplot(data = n_sim_out,
       mapping = aes(x = spatial_het, y = SI)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  theme_classic()

ggplot(data = n_sim_out,
       mapping = aes(x = temporal_het, y = TI)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  theme_classic()

ggplot(data = n_sim_out,
       mapping = aes(x = total_het, y = IT, colour = env_optima_LH)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  theme_classic()

View(n_sim_out)

df <- mix.mono
df <- t1a

# define expected relative yields
df$RYe <- rep(rep(1/2, 2), n_unique(df$sample))

# define observed relative yields
df$RYo <- ifelse(df$M == 0, 0, (df$Y/df$M))

df$dRY <- (df$RYo - df$RYe)

sum( df$M * df$dRY )
sum(df$dRY)

