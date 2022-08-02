#'
#' @title: Simulate metacommunities to test the analytical pipeline
#' 
#' @description: Simulate BEF experiments and calculate biodiversity effects
#' 
#' @details: This script uses the simulated metacommunities and then then assumes that 
#' we only have monoculture data in 30% of the patches and uses a simple linear model 
#' with mixture yield and environmental data to generate posterior predictions for the 
#' missing monoculture data. In addition, it generates a set of 100 starting relative 
#' abundances from the Dirichlet distribution. These outputs are all saved as .rds 
#' files.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)

# load relevant libraries
library(here)
library(rstan)
library(foreach)
library(doParallel)

# read in the MC_sims object
MC_sims <- readRDS(here::here("BEF_quant_scale/results/MC_sims.rds"))

# use lapply to extract the MC.x.NA data
MC_sims_NA <- lapply(MC_sims, function(x) x[["MC.x.NA"]])

# remove it from the MC_sims list
MC_sims <- lapply(MC_sims, function(x) { 
  
  list("MC.x.ids" = x[["MC.x.ids"]],
       "BEF_obs" = x[["BEF_obs"]]) 
  
  }  )

# set-up a parallel for-loop
n.cores <- 10

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

MC.x.pred <- foreach(
  
  i = 1:length(MC_sims_NA)
  
) %dopar% { 

  MC.x.NA <- MC_sims_NA[[i]]
  
  # select 30% of patches to have monoculture and mixture data
  
  # get the number of monoculture and mixtures patches which we set at 30% 
  n_MoMi <- round( 0.30 * length( unique(MC.x.NA$place) ), 0 )
  p_comb <- combn(unique(MC.x.NA$place), m = n_MoMi)
  
  # calculate the environmental range among each of the combinations of two places
  max_env <- 
    apply(p_comb, 2, function(combination) {
      
      range_comb <- 
        
        range( MC.x.NA[MC.x.NA[["place"]] %in% combination, ][["env"]] )
      
      # take the difference in range
      output <- diff(range_comb)
      
      return(output)
      
    })
  
  # most likely that we would have say 2 places and all times: Pick the 2 places with most environmental variation
  MC.x.NA$M1 <- ifelse( MC.x.NA$place %in% p_comb[, sample(which(max_env == max(max_env)), 1)], MC.x.NA$M, NA)
  
  # model the monoculture yields using rstan()
  
  # get the complete cases to fit the model (i.e. training data)
  v <- MC.x.NA[complete.cases(MC.x.NA), ]
  
  # set-up a list to run through ulam()
  MC.x.train <- list(
    M = as.integer(v$M) ,
    Y = (v$Y - mean(v$Y))/sd(v$Y) ,
    E = v$env,
    S = as.integer(v$species))
  
  # fit the model in stan
  m1 <- rstan::stan_model(here::here("BEF_quant_scale/scripts/02_simulation/missing_monoculture_glm.stan"))
  m1.fit <- rstan::sampling(m1, data = MC.x.train)
  
  # extract the posterior distribution
  m1.post <- rstan::extract(m1.fit)
  
  # predict the missing data
  MC.x.pred <- MC.x.NA[is.na(MC.x.NA$M1), ]
  MC.x.pred <- MC.x.pred[, c("species", "Y", "env")]
  names(MC.x.pred) <- c("S", "Y", "E") 
  MC.x.pred$Y <- (MC.x.pred$Y - mean(MC.x.pred$Y))/sd(MC.x.pred$Y)
  
  # loop over each sample from the posterior distribution
  
  # set an output list
  pred.list <- vector("list", length = 2000)
  for (j in 1:2000) {
    
    sample_pred <- 
      
      apply(MC.x.pred, 1, function(obs) {
        
        # set the relevant values for the predictor variables
        S <- obs[1]; Y <- obs[2]; E <- obs[3]
        
        # run through the linear model to get the log-lambda value
        lambda <- m1.post$aS[j, S] + (m1.post$b_yS[j, S]*Y) + (m1.post$b_eS[j, S]*E) + (m1.post$b_yeS[j, S]*E*Y)
        
        # convert the log-lambda value to lambda
        lambda <- exp(lambda)
        
        # run through a Poisson distribution
        return(rpois(n = 1, lambda = lambda))
        
      } )
    
    pred.list[[j]] <- sample_pred
    
  }
  
  # bind this into a matrix
  MC.x.pred <- do.call("rbind", pred.list)
  
  return( list("MC.x.NA" = MC.x.NA, "MC.x.pred" = MC.x.pred) )
  
  }


# attach the model output to the other simulation list
MC_sims <- 
  
  mapply(function(x, y) {
    
    # attach the predictions
    return( c(x, y) )
    
    }, 
  
  MC_sims, 
  MC.x.pred, 
  
  SIMPLIFY = FALSE )

# save the MC_sims object
saveRDS(object = MC_sims, here::here("BEF_quant_scale/results/MC_sims2.rds"))

### END
