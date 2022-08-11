#'
#' @title: Simulate metacommunities to test the analytical pipeline
#' 
#' @description: Calculates biodiversity effects when initial relative abundance data is missing
#' 
#' @details: This script uses a distribution starting relative abundances from 
#' the Dirichlet distribution to generate posterior distributions of biodiversity effects. 
#' This script is very computationally intensive and was run in parellel on a computer cluster using 10
#' cores (Albiorix: http://mtop.github.io/albiorix/).
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(here)
library(dplyr)
library(foreach)
library(doParallel)

# load the relevant functions
source(here("BEF_quant_scale/scripts/01_partition_functions/02_isbell_2018_partition.R"))

# read in the data
MC_sims2 <- readRDS(file = here("BEF_quant_scale/results/MC_sims2.rds"))
start_RA <- readRDS(file = here("BEF_quant_scale/results/MC_sims_start_RA.rds"))

# set-up a parallel for-loop
n.cores <- 10

# create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# get biodiversity effects considering uncertainty in:

# 1. initial relative abundances

BEF_post2 <- foreach(
  
  i = 1:length(MC_sims2)
  
) %dopar% {
  
  RYe_reps <- 
    
    apply(
      
      X = start_RA, 
      
      MARGIN = 2, 
      
      FUN = function(RA) {
        
        # calculate te biodiversity effects for each of the potential starting abundances
        BEF_post.x <- Isbell_2018_sampler(data = MC_sims2[[i]][["MC.x.NA"]], RYe = RA, RYe_post = FALSE)
        names(BEF_post.x[["L.Beff"]])[names(BEF_post.x[["L.Beff"]]) == "L.Beff"] <- "Beff"
        
        # combine the general biodiversity effects and local effects into one data.frame
        BEF_post.x <- rbind(BEF_post.x[["Beff"]], BEF_post.x[["L.Beff"]])
        
        # convert to a data.frame
        BEF_post.x <- as.data.frame(BEF_post.x, row.names = NULL)
        
        return(BEF_post.x)
        
      } )
  
  output <- dplyr::bind_rows(RYe_reps, .id = "rep")
  
  return(output)
  
}

# save this as an RDS file
saveRDS(object = BEF_post2, file = here("BEF_quant_scale/results/BEF_post2.rds"))

### END
