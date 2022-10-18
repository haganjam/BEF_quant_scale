
# calculate the BEF effects

# load the required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(here)

# load the analysis data
data <- read_csv(here("data/benthic_communities_tjarno_data/data_clean/biomass_env_analysis_data.csv"))

# set the OTUs
sp <- sort(unique(data$OTU))

# read in the posterior distribution
post <- readRDS(here("results/SP_A_monoculture_posterior.rds"))

# read in the model object
model_ob <- readRDS(here("results/SP_A_model_object.rds"))

# assign input variables to the names
data_NA <- 
  data %>%
  filter(OTU == "Barn", is.na(M))

# extract the variable names
var_names <- names(model_ob@data)
var_names <- var_names[var_names != "M"]

# extract the parameter names
post_names <- names(post)
names(post_samp) <- NULL

# extract the correct distribution
dist <- gsub(pattern = "d", replacement = "r", x = model_ob@formula[[1]][[3]] )
dist <- paste0(dist[1], "(", dist[2], ",", dist[3], ")")

# choose how many samples to draw from the posterior
n_samp <- 100

# get these predictions for n samples
post_pred <- 
  
  sapply(1:n_samp, function(x) {
    
    for(j in var_names) {
      
      assign(x = j, data_NA[[j]])
      
    }
    
    # take a sample from the posterior distribution
    sample_id <- sample(x = 1:length(post$bY), 1)
    
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



# next steps

# add monoculture predictions from one sample into the data.frame
data[which( (is.na(data[["M"]])) & (data[["OTU"]] == "Barn") ), ][["M"]] <- post_pred[,1]
View(data)



