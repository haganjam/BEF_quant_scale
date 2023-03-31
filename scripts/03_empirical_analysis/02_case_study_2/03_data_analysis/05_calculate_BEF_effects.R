#'
#' @title: Calculate the BEF effects using the modeled monoculture data
#' 
#' @description: This script calculates Isbell et al.'s (2018) biodiversity effects 
#' on the benthic marine fouling communitiesusing all the variation in the 
#' modelled monocultures along with the uncertainty from the Dirichlet distribution.
#'  This script is very computationally intensive and was run in parellel on a computer cluster using 10
#' cores (Albiorix: http://mtop.github.io/albiorix/).
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(dplyr)
library(here)
library(foreach)
library(doParallel)

# load the relevant functions
source(here("BEF_quant_scale/scripts/01_partition_functions/02_isbell_2018_partition.R"))

# load the analysis data
data_M <- readRDS(here("BEF_quant_scale/results/benthic_BEF_data.rds"))

# load the starting relative abundance data
start_RA <- readRDS(here("BEF_quant_scale/results/benthic_start_RA.rds"))

# load the monoculture prediction data
sp_mono <- readRDS(file = here("BEF_quant_scale/results/benthic_mono_pred.rds"))

# check how many samples from the posterior distribution there are
n_samp <- (sapply(sp_mono, function(x) ncol(x)))
all(n_samp == max(n_samp))
n_samp <- max(n_samp)

# set-up a parallel for-loop
n.cores <- 10

# create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

BEF_eff <- foreach(
  
  i = 1:n_samp
  
) %dopar% { 
  
  # load the dplyr package
  library(dplyr)
  
  # add monoculture predictions from one sample into the data.frame
  for(j in 1:length(sp_mono)) {
    
    data_M[which( (is.na(data_M[["M"]])) & (data_M[["species"]] == j) ), ][["M1"]] <- sp_mono[[j]][,i]
    
  }
  
  # rename the monoculture columns
  data_M2 <- 
    data_M %>%
    select(cluster_id, sample, place, time, species, M1, Y) %>%
    rename(M = M1)
  
  # split the data by cluster_id
  clus_df <- split(data_M2[, -1], data_M2$cluster_id)
  
  BEF.x <- 
    
    lapply(clus_df, function(data_in) {
      
      RYe_reps <- 
        
        apply(
          
          X = start_RA, 
          
          MARGIN = 2, 
          
          FUN = function(RA) {
            
            # calculate the biodiversity effects for each of the potential starting abundances
            BEF_post.x <- Isbell_2018_sampler(data = data_in, RYe = RA, RYe_post = FALSE)
            names(BEF_post.x[["L.Beff"]])[names(BEF_post.x[["L.Beff"]]) == "L.Beff"] <- "Beff"
            
            # combine the general biodiversity effects and local effects into one data.frame
            BEF_post.x <- rbind(BEF_post.x[["Beff"]], BEF_post.x[["L.Beff"]])
            
            # convert to a data.frame
            BEF_post.x <- as.data.frame(BEF_post.x, row.names = NULL)
            
            return(BEF_post.x)
            
          } )
      
      # bind into a data.frame
      return(bind_rows(RYe_reps, .id = "RA"))
      
    } )
  
  # bind the rows from the different clusters
  BEF_eff <- bind_rows(BEF.x, .id = "cluster_id")
  
  # return the BEF_eff
  return(BEF_eff)
  
  }

# bind into a data.frame
BEF_eff <- 
  bind_rows(BEF_eff, .id = "mono_sample") %>%
  as_tibble() %>%
  arrange(cluster_id, mono_sample, RA, Value)

# save this object
saveRDS(object = BEF_eff, file = here("BEF_quant_scale/results/benthic_BEF_effects.rds"))

### END
