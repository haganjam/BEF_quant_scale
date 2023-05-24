#'
#' @title: Simulate metacommunities to test the analytical pipeline
#' 
#' @description: Calculates biodiversity effects when initial relative abundance data is missing
#' 
#' @details: This script uses a distribution starting relative abundances from 
#' the Dirichlet distribution to generate posterior distributions of biodiversity effects.
#' The code takes around 20 minutes to run on a regular desktop computer. 
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(dplyr)
library(ggplot2)

# load the relevant functions
source("scripts/01_partition_functions/01_isbell_2018_partition.R")

# read in the data
MC_sims <- readRDS(file = "results/MC_sims.rds")
start_RA <- readRDS(file = "results/MC_sims_start_RA.rds")

# check an example dataset
n <- 6

# plot the mixtures
MC_sims[[n]]$MC_dat %>%
  ggplot(data = .,
         mapping = aes(x = time, y = Y, colour = as.character(species))) +
  geom_line() +
  facet_wrap(~place) +
  theme_test() +
  theme(legend.position = "none")

# plot the monocultures
MC_sims[[n]]$MC_dat %>%
  ggplot(data = .,
         mapping = aes(x = time, y = M, colour = as.character(species))) +
  geom_line() +
  facet_wrap(~place) +
  theme_test() +
  theme(legend.position = "none")


# calculate the biodiversity incorporating uncertainty in RYe

# generate output list
output_list <- vector("list", length = length(MC_sims))

# loop over each each simulated dataset and over each sample from the Dirichlet
for(i in 1:length(MC_sims)) {
  
  RYe_reps <- 
    
    apply(
      
      X = start_RA, 
      
      MARGIN = 2, 
      
      FUN = function(RA) {
        
        # calculate te biodiversity effects for each of the potential starting abundances
        BEF_post <- Isbell_2018_sampler(data = MC_sims[[i]][["MC_dat"]], RYe = RA, RYe_post = FALSE)
        names(BEF_post[["L.Beff"]])[names(BEF_post[["L.Beff"]]) == "L.Beff"] <- "Beff"
        
        # combine the general biodiversity effects and local effects into one data.frame
        BEF_post <- rbind(BEF_post[["Beff"]], BEF_post[["L.Beff"]])
        
        # convert to a data.frame
        BEF_post <- as.data.frame(BEF_post, row.names = NULL)
        
        return(BEF_post)
        
      } )
  
  # pull output into a data.frame
  output <- dplyr::bind_rows(RYe_reps, .id = "rep")
  
  # summarise these data
  output <- 
    output %>%
    group_by(Beff) %>%
    summarise(PI90_low = quantile(Value, 0.05),
              PI90_high = quantile(Value, 0.95),
              mean_BEF = mean(Value, na.rm = TRUE), .groups = "drop")
  
  # add the observed values
  names(MC_sims[[i]][["BEF_obs"]]) <- c("Beff", "BEF_obs")
  
  # join these datasets
  output <- full_join(output, MC_sims[[i]][["BEF_obs"]], by = "Beff")
  
  # write output to the list
  output_list[[i]] <- output
  
}

# bind the output list into a data.frame
output_df <- bind_rows(output_list, .id = "sim_rep")

# save this as an RDS file
saveRDS(object = output_df, file = "results/MC_sims_post.rds")

### END
