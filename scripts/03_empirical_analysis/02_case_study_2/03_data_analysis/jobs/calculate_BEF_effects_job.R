#'
#' @title: Calculate the BEF effects using the imputed monoculture data
#' 
#' @description: This script calculates Isbell et al.'s (2018) biodiversity effects 
#' on the benthic marine fouling communities using all the variation in the 
#' modelled monocultures along with the uncertainty from the Dirichlet distribution.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(future)
library(future.apply)
library(dplyr)

# number of workers
n.cores <- 8

# set up parallel back‑end
plan(multisession, workers = n.cores)

# load the relevant functions
source(here::here("scripts/01_partition_functions/01_isbell_2018_partition.R"))

# load the analysis data
df_imp <- readRDS(here::here("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/05_complete_data.rds"))

# load the starting relative abundance data
start_RA <- readRDS(here::here("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/05_start_RYE.rds"))

BEF_list <-
  future_lapply(
    X = seq_along(df_imp),
    FUN = function(i) {
      
      # extract the focal data.frame
      data <- df_imp[[i]]
      
      # split the data into a list by cluster
      clust_list <- split(data, data[["cluster_id"]])
      
      # iterate over each cluster
      clust_rep <- 
        
        lapply(clust_list, function(x) {
          
          # extract the cluster id
          cluster_id <- unique(x[["cluster_id"]])
          
          # iterate over all possible initial RYE values
          RYE_rep <- 
            
            lapply(start_RA[[cluster_id]], function(RA) {
              
              RYE <- vector("list", length = nrow(RA))
              for(j in 1:nrow(RA)) { RYE[[j]] <- RA[j,] }
              
              # calculate the BEF effects
              BEF_post <- isbell_2018_part(data = x[, names(x) != "cluster_id"], RYe = RYE)
              names(BEF_post[["L.Beff"]])[names(BEF_post[["L.Beff"]]) == "L.Beff"] <- "Beff"
              
              # combine the general biodiversity effects and local effects into one data.frame
              BEF_post <- rbind(BEF_post[["Beff"]], BEF_post[["L.Beff"]])
              
              # convert to a data.frame
              BEF_post <- as.data.frame(BEF_post, row.names = NULL)
              
              # add the cluster id variable
              BEF_post[["cluster_id"]] <- cluster_id
              
              return(BEF_post)
              
            })
          
          
          # bind into a data.fram
          RYE_rep <- dplyr::bind_rows(RYE_rep, .id = "RYE")
          
          return(RYE_rep)
          
        })
      
      # bind the list into a data.frame
      clust_df <- dplyr::bind_rows(clust_rep)
      
      # convert to a tibble
      clust_df <- dplyr::as_tibble(clust_df)
      
      # reorder the columns
      clust_df <- dplyr::select(clust_df, cluster_id, RYE, Beff, Value)
      
      return(clust_df)
      
    }
  )

# bind into a data.frame
BEF_list <- dplyr::bind_rows(BEF_list, .id = "mono_rep")

# convert to a tibble
BEF_list <- dplyr::as_tibble(BEF_list)

# rearrange the columns
BEF_list <- dplyr::select(BEF_list, cluster_id, mono_rep, RYE, Beff, Value)

# arrange by column
BEF_list <- dplyr::arrange(BEF_list, cluster_id, mono_rep, RYE, Beff, Value)

# save this object
saveRDS(object = BEF_list, file = here::here("results/BEF_effects_case_study_2.rds"))

### END
