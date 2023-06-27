
#' @title simulate environments using a random walk
#'
#' @description script generates environmental datasets for use in the
#' simulations to test how well using the random RYe values work
#'
#' @author James G. Hagan, \email{james_hagan@@outlook.com}
#'

# load relevant libraries
library(dplyr)
library(ggplot2)

# function to simulate a random walk
random_walk <- function(N, x0, mu, variance) {
  z<-cumsum(rnorm(n=N, mean=0, 
                  sd=sqrt(variance)))
  t<-1:N
  x<-x0+t*mu+z
  return(x)
}

# function to standardise a variable between two values
range_stand <- function(x, min, max) {
  
  ((max-min)*((x - min(x))/(max(x)-min(x)))) + min
  
}

# set the number of replicates
N_REP <- 500

# set the minimum and maximum of different patches
patches <- 5
timesteps <- 300


# environmental variation for the spatial insurance effect

# data.frame of minimum and maxima for standardisation
df_SI <- data.frame(min = c(0, 0.2, 0.4, 0.6, 0.8),
                     max = c(0.2, 0.4, 0.6, 0.8, 1))
var <- 1
init <- 10
mu <- 0

# replicate N_REP times
env_SI <- vector("list", length = N_REP)
for(j in 1:N_REP) {
  
  # loop over the different patches
  env <- vector("list", length = patches)
  for(i in 1:patches) {
    
    x <- random_walk(timesteps, init, mu, var)
    y <- range_stand(x = x, min = df_SI[["min"]][i], max = df_SI[["max"]][i])
    env[[i]] <- y 
    
  }
  
  # wrap into a data.frame
  df <- 
    data.frame(env1 = round(unlist(env), 3),
               patch = rep(1:patches, each = timesteps),
               time = rep(1:timesteps, patches)
               )
  
  # add the data to the spatial insurance environmental list
  env_SI[[j]] <- df
  
}

# check an example
ggplot(data = env_SI[[1]], 
       mapping = aes(x = time, y = env1, colour = as.character(patch) )) +
  geom_line()


# environmental variation for the temporal insurance effect

# maximimising temporal insurance effects
var <- 0.3
init <- 10
mu <- 0

# replicate N_REP times
env_TI <- vector("list", length = N_REP)
for(j in 1:N_REP) {
  
  env <- vector("list", length = patches)
  for(i in 1:patches) {
    
    x <- random_walk(timesteps, init, mu, var)
    env[[i]] <- x
    
  }
  
  # unlist the environmental variable
  env <- round(unlist(env), 3)
  
  # range standardise between 0 and 1
  env <- range_stand(x = env, min = 0, max = 1)
  
  # wrap into a data.frame
  df <- 
    data.frame(env1 = env,
               patch = rep(1:patches, each = timesteps),
               time = rep(1:timesteps, patches)
    )
  
  # add to output list
  env_TI[[j]] <- df
  
}

# plot an example
ggplot(data = env_TI[[1]], 
       mapping = aes(x = time, y = env1, colour = as.character(patch) )) +
  geom_line()


# environmental variation with spatial and temporal insurance

# both temporal and spatial variation
init <- runif(n = patches, 5, 25)
var <- 0.3
mu <- 0

# replicate N_REP times
env_com <- vector("list", length = N_REP)
for(j in 1:N_REP) {
  
  env <- vector("list", length = patches)
  for(i in 1:patches) {
    
    x <- random_walk(timesteps, init[i], mu, var)
    env[[i]] <- x
    
  }
  
  # unlist the environmental variable
  env <- round(unlist(env), 3)
  
  # range standardise between 0 and 1
  env <- range_stand(x = env, min = 0, max = 1)
  
  # wrap into a data.frame
  df <- 
    data.frame(env1 = env,
               patch = rep(1:patches, each = timesteps),
               time = rep(1:timesteps, patches)
    )
  
  env_com[[j]] <- df
  
}

# plot an example
ggplot(data = env_com[[1]], 
       mapping = aes(x = time, y = env1, colour = as.character(patch) )) +
  geom_line()

# quantify the variation in space and time

# space
bind_rows(env_SI, .id = "rep") %>%
  group_by(rep, patch) %>%
  summarise(mean_env1 = mean(env1)) %>%
  pull(mean_env1) %>%
  var()

bind_rows(env_TI, .id = "rep") %>%
  group_by(rep, patch)  %>%
  summarise(mean_env1 = mean(env1)) %>%
  pull(mean_env1) %>%
  var()

bind_rows(env_com, .id = "rep") %>%
  group_by(rep, patch)  %>%
  summarise(mean_env1 = mean(env1)) %>%
  pull(mean_env1) %>%
  var()

# time
bind_rows(env_SI, .id = "rep") %>%
  group_by(rep, time)  %>%
  summarise(mean_env1 = mean(env1)) %>%
  pull(mean_env1) %>%
  var()

bind_rows(env_TI, .id = "rep") %>%
  group_by(rep, time)%>%
  summarise(mean_env1 = mean(env1)) %>%
  pull(mean_env1) %>%
  var()

bind_rows(env_com, .id = "rep") %>%
  group_by(rep, time) %>%
  summarise(mean_env1 = mean(env1)) %>%
  pull(mean_env1) %>%
  var()

# comparing the variances confirms the simulations

# export these data for use in the simulations
saveRDS(env_SI, "scripts/02_simulation/02_generated_environment_SI.rds")
saveRDS(env_TI, "scripts/02_simulation/02_generated_environment_TI.rds")
saveRDS(env_com, "scripts/02_simulation/02_generated_environment_com.rds")

### END
