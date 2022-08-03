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
library(foreach)
library(doParallel)

# set the percentage threshold from the actual biodiversity effect (see details above)
thresh <- 0.5

# load the observed BEF values
MC_sims2 <- readRDS(file = here("BEF_quant_scale/results/MC_sims2.rds"))

# load the posterior data
BEF_post <- readRDS(file = here("BEF_quant_scale/results/BEF_post.rds"))

# set-up a parallel for-loop
n.cores <- 10

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

BEF_sum_list <- foreach(
  
  i = 1:10  # length(MC_sims2)
  
) %dopar% {
  
  # load the dplyr package
  library(dplyr)
  
  # set-up the PI function from the rethinking package
  PI <- function (samples, prob = 0.89) {
    x <- sapply(prob, function(p) {
      a <- (1 - p)/2
      quantile(samples, probs = c(a, 1 - a))
    })
    n <- length(prob)
    result <- rep(0, n * 2)
    for (i in 1:n) {
      low_idx <- n + 1 - i
      up_idx <- n + i
      result[low_idx] <- x[1, i]
      result[up_idx] <- x[2, i]
      a <- (1 - prob[i])/2
      names(result)[low_idx] <- concat(round(a * 100, 0), "%")
      names(result)[up_idx] <- concat(round((1 - a) * 100, 
                                            0), "%")
    }
    return(result)
  }
  
  BEF_sum <- 
    
    full_join(
      
      BEF_post[[i]] %>% 
        group_by(Beff) %>%
        summarise(mu = mean(Value, na.rm = TRUE),
                  PI_low = round( PI(Value, prob = 0.95)[1], 2),
                  PI_high = round( PI(Value, prob = 0.95)[2], 2) ),
      
      rename(MC_sims2[[i]] [["BEF_obs"]], Value_obs = Value ), 
      
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
  BEF_sum <- bind_cols( MC_sims2[[i]] [["MC.x.ids"]], BEF_sum)
  
  # add bad monoculture correlation information
  
  # calculate the overall monoculture error
  
  # summarise the posterior distribution
  mu_m <- apply(MC_sims2[[i]] [["MC.x.pred"]], 2, function(x) mean(x) )
  PI_m <- apply(MC_sims2[[i]] [["MC.x.pred"]], 2, function(x) PI(x, 0.95) )
  
  # calculate the correlation coefficient
  BEF_sum[["mono_cor"]] <- cor(mu_m, MC_sims2[[i]] [["MC.x.NA"]]$M[is.na(MC_sims2[[i]] [["MC.x.NA"]]$M1)])
  
  # calculate absolute error between monocultures and mixtures for all posterior samples
  mono_error <- 
    
    apply(MC_sims2[[i]] [["MC.x.pred"]], 1, function(x) {
      
      M_obs <- MC_sims2[[i]] [["MC.x.NA"]] [is.na(MC_sims2[[i]] [["MC.x.NA"]][["M1"]]), ] [["M"]]
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
  
  return(BEF_sum)
  
  }
  
# bind this list into a data.frame that can be analysed
BEF_output <- dplyr::bind_rows(BEF_sum_list, .id = "model_ID")

# output this as a .rds file
# save this as an RDS file
saveRDS(object = BEF_output, file = here("BEF_quant_scale/results/BEF_output.rds"))

### END
