
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

# load relevant functions
source("scripts/Function_plotting_theme.R")

# get a colour palette
col_pal <- wesanderson::wes_palette(name = "Darjeeling1", 5)

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
env1Scale <- 5

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

# plot out 5 examples of the spatially fluctuating environment
p1 <- 
  bind_rows(env_SI[1], .id = "rep") %>% 
  as_tibble %>%
  mutate(patch = as.character(patch),
         id = paste(rep, patch, sep = "_") ) %>%
  ggplot(data = ., 
         mapping = aes(x = time, y = env1, colour = patch, group = id )) +
  geom_line(alpha = 1) +
  labs(colour = "Patch") +
  ylab("Environment (0-1)") +
  xlab("Time") +
  scale_colour_manual(values = col_pal) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p1)

# environmental variation for the temporal insurance effect

# replicate N_REP times
env_TI <- vector("list", length = N_REP)
for(j in 1:N_REP) {
  
  # get a fluctuating environment
  env1 <- synchrony::phase.partnered(n = timesteps, gamma = env1Scale, mu = 0.5, sigma = 0.25)$timeseries[,1]
  
  env_df <- data.frame()
  for(i in 1:patches){
    env1 <- env1 + rnorm(n = 1, mean = 0, sd = 0.02)
    env_df <- rbind(env_df, data.frame(env1 = env1, patch = i, time = 1:timesteps))
    
  }
  
  env_df$env1 <- range_stand(env_df$env1, 0.1, 0.9) 
  
  # add to output list
  env_TI[[j]] <- env_df
  
}

# plot out 5 examples of the temporally fluctuating environment
p2 <- 
  bind_rows(env_TI[1], .id = "rep") %>% 
  as_tibble %>%
  mutate(patch = as.character(patch),
         id = paste(rep, patch, sep = "_") ) %>%
  ggplot(data = ., 
         mapping = aes(x = time, y = env1, colour = patch, group = id )) +
  geom_line(alpha = 0.5) +
  labs(colour = "Patch") +
  ylab("") +
  xlab("Time") +
  scale_colour_manual(values = col_pal) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p2)
 

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

# plot out 5 examples of the temporally fluctuating environment
p3 <- 
  bind_rows(env_com[1], .id = "rep") %>% 
  as_tibble %>%
  mutate(patch = as.character(patch),
         id = paste(rep, patch, sep = "_") ) %>%
  ggplot(data = ., 
         mapping = aes(x = time, y = env1, colour = patch, group = id )) +
  geom_line(alpha = 0.5) +
  labs(colour = "Patch") +
  ylab("") +
  xlab("Time") +
  scale_colour_manual(values = col_pal) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p3)

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

# make a plot for the supplementary
ggpubr::ggarrange(p1, p2, p3, 
                  labels = c("a", "b", "c"),
                  font.label = list(size = 11, face = "plain"), 
                  common.legend = TRUE,
                  nrow = 1)

# export these data for use in the simulations
saveRDS(env_SI, "scripts/02_simulation/02_generated_environment_SI.rds")
saveRDS(env_TI, "scripts/02_simulation/02_generated_environment_TI.rds")
saveRDS(env_com, "scripts/02_simulation/02_generated_environment_com.rds")

### END
