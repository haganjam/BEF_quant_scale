
# Simulate test data

# Install the Thompson et al. (2020): mcomsimr
devtools::install_github("plthompson/mcomsimr")
library(mcomsimr)

# set-up some core parameters
np <- 3
ns <- 3
ts <- 100
bi <- 50
init <- 50

# generate a random landscape
l.1 <- landscape_generate(patches = np)
head(l.1)

# generate a random dispersal matrix
d.1 <- dispersal_matrix(landscape = l.1,
                        torus = TRUE,
                        kernel_exp = 0.1)
head(d.1)

# generate a random environmental matrix
e.1 <- env_generate(landscape = l.1, env1Scale = 500, timesteps = (ts + bi + init) )
head(e.1)

# generate random species traits
t.1 <- env_traits(species = ns, max_r = 5, 
                  min_env = 0, max_env = 1, 
                  env_niche_breadth = 0.5, plot = TRUE, optima_spacing = "random")
head(t.1)

# generate species interaction matrix
si.1 <- species_int_mat(species = ns, intra = 1, min_inter = 0, max_inter = 1.5, 
                        comp_scaler = 0.05, plot = TRUE)
head(si.1)

# put these data into the simulate_MC function
simulate_MC2(species = ns, patches = np, dispersal = 0.1, plot = TRUE, 
             timesteps = ts, burn_in = bi, initialization = init,
             extirp_prob = 0.001,
             landscape = l.1, 
             disp_mat = d.1, 
             env.df = e.1, 
             env_traits.df = t.1, 
             int_mat = si.1)

# try with just one species
simulate_MC2(species = 1, patches = np, dispersal = 0.1, plot = TRUE, 
             timesteps = ts, burn_in = bi, initialization = init,
             extirp_prob = 0.001,
             landscape = l.1, 
             disp_mat = d.1, 
             env.df = e.1, 
             env_traits.df = t.1[t.1$species == 1,], 
             int_mat = si.1[1, 1])


