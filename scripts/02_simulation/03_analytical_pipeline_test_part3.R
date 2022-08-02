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

# set the percentage threshold from the actual biodiversity effect (see details above)
thresh <- 0.5

# load the observed BEF values
MC_sims <- readRDS(file = here("results/MC_sims.rds"))

# load the posterior data
BEF_post <- readRDS(file = here("results/BEF_post.rds"))

# how many models were run?
n <- length(MC_sims)
BEF_sum_list <- vector("list", length = n)

for (i in 1:length(MC_sims)) {
  
  BEF_sum <- 
    
    full_join(
      
      BEF_post[[i]] %>% 
        group_by(Beff) %>%
        summarise(mu = mean(Value, na.rm = TRUE),
                  HPDI_low = HPDI(Value, prob = 0.95)[1],
                  HPDI_high = HPDI(Value, prob = 0.95)[2],),
      
      rename(MC_sims[[i]] [["BEF_obs"]], Value_obs = Value ), 
      
      by = "Beff"
      
    )
  
  # generate some additional summary variables
  
  BEF_sum <- 
    
    BEF_sum %>%
    mutate(Value_obs_min = (1-thresh)*Value_obs,
           Value_obs_max = (1+thresh)*Value_obs ) %>%
    mutate(HPDI_true = if_else( (HPDI_low > Value_obs_min) & (HPDI_high < Value_obs_max), TRUE, FALSE ) ) %>%
    mutate(mu_deviation = abs((abs(mu - Value_obs)/Value_obs)*100) ) %>%
    mutate(mu_threshold = thresh)
  
  # add the identifer variables  
  BEF_sum <- bind_cols( MC_sims[[i]] [["MC.x.ids"]], BEF_sum)
  
  # add bad correlation information
  
  # summarise the posterior distribution
  mu_m <- apply(MC_sims[[i]] [["MC.x.pred"]], 2, function(x) mean(x) )
  PI_m <- apply(MC_sims[[i]] [["MC.x.pred"]], 2, function(x) PI(x, 0.95) )
  
  # calculate the correlation coefficient
  BEF_sum[["mono_cor"]] <- cor(mu_m, MC_sims[[i]] [["MC.x.NA"]]$M[is.na(MC_sims[[i]] [["MC.x.NA"]]$M1)])
  
  # definition of a bad monoculture
  BEF_sum[["bad_mono_n"]] <- sum( apply(PI_m, 2, diff)/mean(mu_m) > 1 )/length(mu_m)
  
  # reorder the columns
  
  # important: When HPDI_true is TRUE, it means that the 95% HPDI interval is within
  # thresh*obs value +- obs value which means that 95% HPDI interval is narrow enough
  # to be meaningful
  BEF_sum <- 
    
    BEF_sum %>%
    select(t_steps, dispersal, start_abun, optima, niche_breadth, inter_comp, mono_cor, bad_mono_n,
           Beff, Value_obs, mu, mu_deviation,
           starts_with("HPDI"), mu_threshold,
           HPDI_true)
  
  BEF_sum_list[[i]] <- BEF_sum
  
}

# bind this list into a data.frame that can be analysed
BEF_output <- bind_rows(BEF_sum_list, .id = "model_ID")

### END
