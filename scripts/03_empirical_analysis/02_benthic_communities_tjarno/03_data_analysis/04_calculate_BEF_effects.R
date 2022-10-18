#'
#' @title: Calculate the BEF effects using the modelled monoculture data
#' 
#' @description: This script uses the samples from the posterior distributions of
#' the best monoculture models of each species to fill in the missing monoculture data
#' and then calculates Isbell et al.'s (2018) biodiversity effects using all the variation
#' in the modelled monocultures along with the uncertainty from the Dirichlet distribution.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# load the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# load the relevant functions
source(here("scripts/01_partition_functions/02_isbell_2018_partition.R"))
source(here("scripts/Function_plotting_theme.R"))

# load the analysis data
data <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/biomass_env_analysis_data.csv"))

# make a Y2 variable
data$Y2 <- data$Y^2

# pivot mixture data longer
mix <- 
  data %>%
  select(cluster_id, buoy_id, time, OTU, Y) %>%
  pivot_wider(id_cols = c("cluster_id", "buoy_id", "time"),
              names_from = "OTU",
              values_from = "Y")

# join these data together
data <- full_join(data, mix, by = c("cluster_id", "buoy_id", "time"))

# set the OTUs
sp <- sort(unique(data$OTU))

# get a list of file names
files <- list.files(here("results/"))

# get a file name list of posterior distributions
files_post <- files[grepl(pattern = "posterior", x = files )]

# get a file name list of model objects
files_mod <- files[grepl(pattern = "model_object", x = files )]

# choose how many samples to draw from the posterior
n_samp <- 100

sp_mono <- vector("list", length = length(sp))
names(sp_mono) <- sp
for (i in 1:length(sp)) {

  # load the posterior distrbution
  post <- readRDS(paste0(here("results"), "/", files_post[i]))
  
  # load the model object
  model_ob <- readRDS(paste0(here("results"), "/", files_mod[i]))
  
  # assign input variables to the names
  data_NA <- 
    data %>%
    filter(OTU == sp[i], is.na(M))
  
  # extract the variable names
  var_names <- names(model_ob@data)
  var_names <- var_names[var_names != "M"]
  
  # extract the parameter names
  post_names <- names(post)
  
  # extract the correct distribution
  dist <- gsub(pattern = "d", replacement = "r", x = model_ob@formula[[1]][[3]] )
  dist <- paste0(dist[1], "(n = length(Y),", dist[2], ",", dist[3], ")")
  
  # get these predictions for n samples
  post_pred <- 
    
    sapply(1:n_samp, function(x) {
      
      for(j in var_names) {
        
        assign(x = j, data_NA[[j]])
        
      }
      
      # take a sample from the posterior distribution
      sample_id <- sample(x = 1:length(post[[1]]), 1)
      
      # assign the parameter values to the names
      post_samp <- sapply(post, function(x) x[sample_id] )
      
      # write a loop and assign a sample from the posterior distribution to a parameter name
      for (k in 1:length(post_names)) {
        
        assign(x = post_names[k], value = post_samp[k])
        
      }
      
      # calculate mu: set-up the expression
      form <- parse(text = model_ob@formula[[2]][[3]])
      
      # evaluate the expression
      u <- eval(form)
      
      # run the u values through the distribution
      dist <- parse(text = dist )
      M1 <- eval(dist)
      
      # if the value is less than zero then set it to zero
      M1 <- ifelse(M1 < 0, 0, M1)
      
      return(M1)
      
    } )
  
  # return the matrix predictions
  sp_mono[[i]] <- post_pred
  
}

# check if the data have the correct dimensions
all(data %>%
      filter(is.na(M)) %>%
      group_by(OTU) %>%
      summarise(n = n()) %>%
      pull(n) == sapply(sp_mono, nrow)
    ) 

# remove the unnecessary columns in the data data.frame and put data in correct format
names(data)

data_M <- 
  data %>%
  select(cluster_id, buoy_id, time, OTU, M, Y) %>%
  rename(place = buoy_id, species = OTU) %>%
  mutate(place = as.integer(as.factor(place)),
         time = as.integer(as.factor(time)),
         species = as.integer(as.factor(species))) %>%
  group_by(cluster_id) %>%
  mutate(sample = as.integer(as.factor(paste0(place, time))) ) %>%
  select(cluster_id, sample, place, time, species, M, Y)

# generate 100 different Dirichlet distribution
start_RA <- sapply(1:100, function(x) gtools::rdirichlet(n = 1, rep(3, length(unique(data_M$species))) ) )


# loop over all the different monoculture predictions

# sp_mono[[j]][,1] <- 1 needs to be an i and we need to loop over the columns in the monoculture predictions

# add monoculture predictions from one sample into the data.frame
for(j in 1:length(sp)) {

  data_M[which( (is.na(data_M[["M"]])) & (data_M[["species"]] == j) ), ][["M"]] <- sp_mono[[j]][,1]
  
}

# split the data by cluster_id
clus_df <- split(data_M[, -1], data_M$cluster_id)

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
bind_rows(BEF.x, .id = "cluster_id")





