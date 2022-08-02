#'
#' @title: Simulate metacommunities to test the analytical pipeline
#' 
#' @description: Calculates biodiversity effects when monoculture data is missing
#' 
#' @details: This script uses the posterior distributions generated for each missing monoculture
#' and a distribution of starting relative abundances from the Dirichlet distribution to generate
#' posterior distributions of biodiversity effects. This script is to be run on a computer
#' cluster because it is very computationally intensive.
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
MC_sims <- readRDS(file = here("BEF_quant_scale/results/MC_sims2.rds"))
start_RA <- readRDS(file = here("BEF_quant_scale/results/MC_sims_start_RA.rds"))

# set-up a parallel for-loop
n.cores <- 10

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

BEF_post <- foreach(
  
  i = 1:length(MC_sims)
  
) %dopar% {
  
  MC.rep <- MC_sims[[i]]
  
  BEF_mod <-
    
    apply(
      
      X = MC.rep[["MC.x.pred"]] [sample(x = 1:nrow(MC.rep[["MC.x.pred"]]), 10), ], 
      
      MARGIN = 1, 
      
      FUN = function(x) {
        
        # fill in the missing monoculture data with one sample from the posterior
        
        # generate a data.frame called MC.x.post
        MC.x.post <- MC.rep[["MC.x.NA"]]
        
        # for the missing monoculture values, we fill them in with posterior distribution samples
        MC.x.post[is.na(MC.x.post$M1), ]$M1 <- x
        
        # select the relevant columns
        MC.x.post <- MC.x.post[,c("sample", "time", "place", "species", "M1", "Y")]
        
        # rename the M1 column
        names(MC.x.post)[names(MC.x.post) == "M1"] <- "M"
        
        RYe_reps <- 
          
          apply(
            
            X = start_RA, 
            
            MARGIN = 2, 
            
            FUN = function(RA) {
              
              # calculate te biodiversity effects for each of the potential starting abundances
              BEF_post.x <- Isbell_2018_sampler(data = MC.x.post, RYe = RA, RYe_post = FALSE)
              names(BEF_post.x[["L.Beff"]])[names(BEF_post.x[["L.Beff"]]) == "L.Beff"] <- "Beff"
              
              # combine the general biodiversity effects and local effects into one data.frame
              BEF_post.x <- rbind(BEF_post.x[["Beff"]], BEF_post.x[["L.Beff"]])
              
              # convert to a data.frame
              BEF_post.x <- as.data.frame(BEF_post.x, row.names = NULL)
              
              return(BEF_post.x)
              
            } )
        
        return(RYe_reps)
        
      } )
  
  output <- dplyr::bind_rows(BEF_mod, .id = "rep")
  
  return(output)
  
}

# save this as an RDS file
saveRDS(object = BEF_post, file = here("BEF_quant_scale/results/BEF_post.rds"))

### END
