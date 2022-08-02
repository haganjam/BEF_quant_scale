#'
#' @title: Summarise the model outputs
#' 
#' @description: Compares biodiversity effects from the analytical pipeline to the observed biodiversity effects
#' 
#' @details: Summarises the distribution for each biodiversity effect and compares it to 
#' the observed value. We assume that for an effect to be meaningful, the 95% HPDI interval 
#' of an effect must be within the interval of the observed value +- 0.5*observed value. 
#' This means that the 95% HPDI interval contains the true value within a reasonable distance of the true value
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(here)
library(dplyr)
library(rethinking)
library(ggplot2)
library(ggbeeswarm)

# set the percentage threshold from the actual biodiversity effect (see details above)
thresh <- 0.5

# load the observed BEF values
MC_sims <- readRDS(file = here("results/MC_sims.rds"))

# load the posterior data
BEF_post <- readRDS(file = here("results/BEF_post.rds"))

# how many models were run?
n <- length(MC_sims)
BEF_sum_list <- vector("list", length = n)

for (i in 1:n) {
  
  BEF_sum <- 
    
    full_join(
      
      BEF_post[[i]] %>% 
        group_by(Beff) %>%
        summarise(mu = mean(Value, na.rm = TRUE),
                  PI_low = round( PI(Value, prob = 0.95)[1], 2),
                  PI_high = round( PI(Value, prob = 0.95)[2], 2) ),
      
      rename(MC_sims[[i]] [["BEF_obs"]], Value_obs = Value ), 
      
      by = "Beff"
      
    )
  
  # generate some additional summary variables
  
  BEF_sum <- 
    
    BEF_sum %>%
    mutate(Value_obs_min = (1-thresh)*Value_obs,
           Value_obs_max = (1+thresh)*Value_obs ) %>%
    mutate(PI_true = if_else( ( (Value_obs > PI_low) & (Value_obs < PI_high) ), TRUE, FALSE ) ) %>%
    mutate(PI_mu_true = if_else( (PI_low > Value_obs_min) & (PI_high < Value_obs_max), TRUE, FALSE ) ) %>%
    mutate(mu_deviation = round( abs((abs(mu - Value_obs)/Value_obs)*100), 2 )  ) %>%
    mutate(mu_threshold = thresh)
  
  # add the identifer variables  
  BEF_sum <- bind_cols( MC_sims[[i]] [["MC.x.ids"]], BEF_sum)
  
  # add bad monoculture correlation information
  
  # calculate the overall monoculture error
  
  # summarise the posterior distribution
  mu_m <- apply(MC_sims[[i]] [["MC.x.pred"]], 2, function(x) mean(x) )
  PI_m <- apply(MC_sims[[i]] [["MC.x.pred"]], 2, function(x) PI(x, 0.95) )
  
  # calculate the correlation coefficient
  BEF_sum[["mono_cor"]] <- cor(mu_m, MC_sims[[i]] [["MC.x.NA"]]$M[is.na(MC_sims[[i]] [["MC.x.NA"]]$M1)])
  
  # calculate absolute error between monocultures and mixtures for all posterior samples
    mono_error <- 
      
    apply(MC_sims[[i]] [["MC.x.pred"]], 1, function(x) {
      
      M_obs <- MC_sims[[i]] [["MC.x.NA"]] [is.na(M_obs[["M1"]]), ] [["M"]]
      sum( abs( x - M_obs ) )
      
    })
    
  BEF_sum[["mono_error"]] <- mean(mono_error, na.rm = TRUE)
  
  # reorder the columns
  
  # important: When PI_true is TRUE, it means that the 95% PI interval is within
  # thresh*obs value +- obs value which means that 95% PI interval is narrow enough
  # to be meaningful
  BEF_sum <- 
    
    BEF_sum %>%
    select(t_steps, dispersal, start_abun, optima, niche_breadth, inter_comp, mono_cor, mono_error,
           Beff, Value_obs, mu, mu_deviation, mu_threshold,
           PI_low, PI_high,
           PI_true, PI_mu_true)
  
  BEF_sum_list[[i]] <- BEF_sum
  
}

# bind this list into a data.frame that can be analysed
BEF_output <- bind_rows(BEF_sum_list, .id = "model_ID")

# analyse the accuracy
View(BEF_output)

# check for outliers
summary(BEF_output)

# check the large mono-error outlier
BEF_output %>%
  filter(mono_error == max(mono_error)) %>%
  View()

# is there a relationship between monoculture correlation and mu deviation
ggplot(data = BEF_output,
         mapping = aes(x = log(mono_error), y = (mu_deviation) )) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Beff, scales = "free") +
  theme_bw()

ggplot(data = BEF_output,
       mapping = aes(x = log10(mono_error), y = log10(PI_high-PI_low) )) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Beff, scales = "free") +
  theme_bw()

BEF_output_sum <- 
  BEF_output %>%
  group_by(Beff) %>%
  summarise(accuracy_value = sum(PI_true)/n(),
            accuracy_interval = sum(PI_mu_true)/n(),
            Value_obs = mean(Value_obs))

# plot the full set of simulations
ggplot(data = BEF_output_sum,
       mapping = aes(x = Beff, y = Value_obs)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_quasirandom(data = BEF_output %>% filter(mu < 10000 & mu > -10000), 
                   mapping = aes(x = Beff, y = mu), alpha = 0.1 ) +
  geom_point(colour = "red", size = 2) +
  theme_bw()

### END
