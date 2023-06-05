#'
#' @title: Calculate the BEF effects using the modeled monoculture data
#' 
#' @description: This script calculates Isbell et al.'s (2018) biodiversity effects 
#' on the benthic marine fouling communities using all the variation in the 
#' modelled monocultures along with the uncertainty from the Dirichlet distribution.
#' This script is very computationally intensive and was run in parellel on a computer cluster using 10
#' cores (Albiorix: http://mtop.github.io/albiorix/).
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(dplyr)
library(foreach)
library(doParallel)

# load the relevant functions
source("scripts/01_partition_functions/01_isbell_2018_partition.R")

# load the analysis data
df_imp <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/05_complete_data.rds")

# load the starting relative abundance data
start_RA <- readRDS("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/05_start_RYE.rds")

# test some of the code
# extract a single dataset
data <- df_imp[[1]]

# for each cluster we lapply
clust_list <- split(data, data[["cluster_id"]])

lapply()

x <- df_imp[[1]] %>% filter(cluster_id == "A")

y <- 
  
  apply(start_RA, 2, function(RA) {
    
    # extract the cluster id
    cluster_id <- unique(x[["cluster_id"]])
    
    # calculate the BEF effects
    BEF_post <- Isbell_2018_sampler(data = x[, names(x) != "cluster_id"], RYe = RA, RYe_post = FALSE)
    names(BEF_post[["L.Beff"]])[names(BEF_post[["L.Beff"]]) == "L.Beff"] <- "Beff"
    
    # combine the general biodiversity effects and local effects into one data.frame
    BEF_post <- rbind(BEF_post[["Beff"]], BEF_post[["L.Beff"]])
    
    # convert to a data.frame
    BEF_post <- as.data.frame(BEF_post, row.names = NULL)
    
    # add the cluster id variable
    BEF_post[["cluster_id"]] <- cluster_id
    
    return(BEF_post)
    
  })





    
RYe_reps <- 
  
  apply(
    
    X = start_RA, 
    
    MARGIN = 2, 
    
    FUN = function(RA) {
      
      # calculate the biodiversity effects for each of the potential starting abundances
      BEF_post <- Isbell_2018_sampler(data = df_imp[[1]], RYe = RA, RYe_post = FALSE)
      names(BEF_post[["L.Beff"]])[names(BEF_post[["L.Beff"]]) == "L.Beff"] <- "Beff"
      
      # combine the general biodiversity effects and local effects into one data.frame
      BEF_post <- rbind(BEF_post[["Beff"]], BEF_post[["L.Beff"]])
      
      # convert to a data.frame
      BEF_post <- as.data.frame(BEF_post, row.names = NULL)
      
      return(BEF_post)
      
    } )
      
# bind into a data.frame
bind_rows(RYe_reps, .id = "RA")


# bind the rows from the different clusters
BEF_eff <- bind_rows(BEF.x, .id = "cluster_id")








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
