
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

# function to standardise a variable between two values
range_stand <- function(x, min, max) {
  
  ((max-min)*((x - min(x))/(max(x)-min(x)))) + min
  
}

# set the number of replicates
N_REP <- 500

# set the minimum and maximum of different patches
patches <- 5
timesteps <- 300

# the scale of environmental variation
env1Scale = 4

# environmental variation for the spatial insurance effect

# data.frame of minimum and maxima for standardisation
df_SI <- data.frame(min = c(0, 0.2, 0.4, 0.6, 0.8),
                     max = c(0.2, 0.4, 0.6, 0.8, 1))

# replicate N_REP times
env_SI <- vector("list", length = N_REP)
for(j in 1:N_REP) {
  
  env_df <- data.frame()
  for(i in 1:patches){
    env1 <- synchrony::phase.partnered(n = timesteps, gamma = env1Scale, mu = 0.5, sigma = 0.25)$timeseries[,1]
    env1 <- range_stand(x = env1, min = df_SI[["min"]][i], max = df_SI[["max"]][i])
    env_df <- rbind(env_df, data.frame(env1 = env1, patch = i, time = 1:timesteps))
    
  }
  
  # add the data to the spatial insurance environmental list
  env_SI[[j]] <- env_df
  
}

# check an example
ggplot(data = env_SI[[1]], 
       mapping = aes(x = time, y = env1, colour = as.character(patch) )) +
  geom_line()


# environmental variation for the temporal insurance effect

# replicate N_REP times
env_TI <- vector("list", length = N_REP)
for(j in 1:N_REP) {
  
  env_df <- data.frame()
  for(i in 1:patches){
    env1 <- synchrony::phase.partnered(n = timesteps, gamma = env1Scale, mu = 0.5, sigma = 0.25)$timeseries[,1]
    env1 <- range_stand(env1, 0.1, 0.9) 
    env_df <- rbind(env_df, data.frame(env1 = env1, patch = i, time = 1:timesteps))
    
  }
  
  # add to output list
  env_TI[[j]] <- env_df
  
}

# plot an example
ggplot(data = env_TI[[1]], 
       mapping = aes(x = time, y = env1, colour = as.character(patch) )) +
  geom_line()


# environmental variation with spatial and temporal insurance

# replicate N_REP times
env_com <- vector("list", length = N_REP)
for(j in 1:N_REP) {
  
  env_df <- data.frame()
  for(i in 1:patches){
    env1 <- synchrony::phase.partnered(n = timesteps, gamma = env1Scale, mu = runif(1, 0, 1), sigma = 0.25)$timeseries[,1]
    env_df <- rbind(env_df, data.frame(env1 = env1, patch = i, time = 1:timesteps))
    
  }
  
  env_df$env1 <- range_stand(env_df$env1, 0, 1)
  
  env_com[[j]] <- env_df
  
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
