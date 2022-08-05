#'
#' @title: Summarise the model outputs
#' 
#' @description: Compares biodiversity effects from the analytical pipeline to the observed biodiversity effects
#' 
#' @details: Summarises the distribution for each biodiversity effect generated
#' from the uncertainty from the monocultures and the starting relative
#' abundances. This distrbution is then compared to the observed value calculated from
#' known monoculture yields and starting relative abundances from the model. There
#' are two important metrics that are not obvious from the variable titles:
#' 
#' PI_true: This simply measures whether the observed value is within the 
#' 90% percentile interval of the distribution.
#' 
#' PI_mu_true: This measures whether the 90% percentile interval of the posterior distribution
#' is contained within 0.50 x (observed value) and 1.5 x (observed value). This quantity
#' shows whether interval generated is reasonably close to the observed value.
#' 
#' This script is also very computationally intensive because it runs on very large
#' data files generated in the previous steps. Therefore, it was run in parellel on a 
#' computer cluster using 10 cores (Albiorix: http://mtop.github.io/albiorix/).
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(here)
library(dplyr)
library(foreach)
library(doParallel)

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
  
  i = 1:length(MC_sims2)
  
) %dopar% {
  
  # load the dplyr package
  library(dplyr)
  
  # set-up the PI function from the rethinking package
  
  
  BEF_sum <- 
    
    full_join(
      
      BEF_post[[i]] %>% 
        group_by(Beff) %>%
        summarise(mu = mean(Value, na.rm = TRUE),
                  PI_low = round( quantile(Value, prob = 0.05), 2),
                  PI_high = round( quantile(Value, prob = 0.95), 2) ),
      
      rename(MC_sims2[[i]] [["BEF_obs"]], Value_obs = Value ), 
      
      by = "Beff"
      
    )
  
  # add the identifer variables  
  BEF_sum <- bind_cols( MC_sims2[[i]] [["MC.x.ids"]], BEF_sum)
  
  # calculate the overall monoculture error
  
  # summarise the posterior distribution
  mu_m <- apply(MC_sims2[[i]] [["MC.x.pred"]], 2, function(x) mean(x) )
  
  # calculate the correlation coefficient between the mean of the posterior and the true monoculture values
  BEF_sum[["mono_cor"]] <- cor(mu_m, MC_sims2[[i]] [["MC.x.NA"]]$M[is.na(MC_sims2[[i]] [["MC.x.NA"]]$M1)])
  
  # calculate absolute error between monocultures and mixtures for all posterior samples
  mono_error <- 
    
    apply(MC_sims2[[i]] [["MC.x.pred"]], 1, function(x) {
      
      M_obs <- MC_sims2[[i]] [["MC.x.NA"]] [is.na(MC_sims2[[i]] [["MC.x.NA"]][["M1"]]), ] [["M"]]
      sum( abs( x - M_obs ) )
      
    })
  
  BEF_sum[["mono_error"]] <- mean(mono_error, na.rm = TRUE)
  
  # reorder the columns
  BEF_sum <- 
    
    BEF_sum %>%
    select(t_steps, dispersal, start_abun, optima, niche_breadth, inter_comp, mono_cor, mono_error,
           Beff, Value_obs, mu,
           PI_low, PI_high)
  
  return(BEF_sum)
  
  }
  
# bind this list into a data.frame that can be analysed
BEF_output <- dplyr::bind_rows(BEF_sum_list, .id = "model_ID")

# output this as a .rds file
# save this as an RDS file
saveRDS(object = BEF_output, file = here("BEF_quant_scale/results/BEF_output.rds"))

### END
