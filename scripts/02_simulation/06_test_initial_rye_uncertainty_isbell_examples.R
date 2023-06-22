#'
#' @title Simulate metacommunities to test the analytical pipeline
#' 
#' @description Here we test whether assuming different RYE proportions allows
#' the accurate estimation of the example effects reported in Isbell et al.
#' (2018)
#' 

# load the functions
source("scripts/01_partition_functions/01_isbell_2018_partition.R")

# load the test data
test_data <- readRDS(file = "data/Isbell_test_data.rds")
ans_data <- readRDS(file = "data/Isbell_test_data_solutions.rds")

# remove the first two examples that don't have temporal data
test_data <- test_data[3:length(test_data)]
ans_data <- ans_data[3:length(ans_data)]

# get the number of unique samples
samples <- length(unique(test_data[[1]]$sample))

# get the number of species
species <- length(unique(test_data[[1]]$species))

# set the number of unique RYe values
Ns <- 100

# get 100 unique RYe values to use in the calculations
start_RA <- vector("list", length = Ns)
for(i in 1:Ns) {
  
  # for each sample, get a simplex from the Dirichlet distribution
  rye_mat <- round(gtools::rdirichlet(n = samples, rep(3, species)), 2)
  
  # write into a list
  start_RA[[i]] <- rye_mat
  
}

# generate output list
output_list <- vector("list", length = length(test_data))

# loop over each each simulated dataset and over each sample from the Dirichlet
for(i in 1:length(test_data)) {
  
  # iterate over all possible initial RYE values
  RYE_rep <- 
    
    lapply(start_RA, function(RA) {
      
      RYE <- vector("list", length = nrow(RA))
      for(j in 1:nrow(RA)) { RYE[[j]] <- RA[j,] }
      
      # calculate the BEF effects
      BEF_post <- isbell_2018_part(data = test_data[[i]], RYe = RYE)
      names(BEF_post[["L.Beff"]])[names(BEF_post[["L.Beff"]]) == "L.Beff"] <- "Beff"
      
      # combine the general biodiversity effects and local effects into one data.frame
      BEF_post <- rbind(BEF_post[["Beff"]], BEF_post[["L.Beff"]])
      
      # convert to a data.frame
      BEF_post <- as.data.frame(BEF_post, row.names = NULL)
      
      return(BEF_post)
      
    })
  
  
  # bind into a data.fram
  RYE_rep <- dplyr::bind_rows(RYE_rep, .id = "RYE")
  
  # summarise these data
  output <- 
    RYE_rep |>
    dplyr::group_by(Beff) |>
    dplyr::summarise(PI95_low = quantile(Value, 0.025),
                     PI95_high = quantile(Value, 0.975),
                     mean_BEF = mean(Value, na.rm = TRUE), .groups = "drop"
    )
  
  # add the observed values
  names(ans_data[[i]]) <- c("Beff", "BEF_obs")
  
  # join these datasets
  output <- dplyr::full_join(output, ans_data[[i]], by = "Beff")
  
  # write output to the list
  output_list[[i]] <- output
  
}

# bind into a data.frame
output_df <- dplyr::bind_rows(output_list, .id = "case")

# get the subset of values that are affected by the RYEs
output_rye <- 
  output_df %>%
  filter( !(Beff %in% c( "LS", "LC", "TC", "NO" )) )

# round of to three decimal places
output_rye <- 
  output_rye %>%
  mutate(PI95_low = round(PI95_low, 3),
         PI95_high = round(PI95_high, 3),
         mean_BEF = round(mean_BEF, 3),
         BEF_obs = round(BEF_obs, 3))

# does the observed value lie in the 90% percentile interval
output_rye <- 
  output_rye %>%
  mutate(within = ifelse(BEF_obs <= PI95_high & BEF_obs >= PI95_low, 1, 0))

# check proportion within the interval
output_rye %>%
  group_by(Beff) %>%
  summarise(prop_within = sum(within)/n())

# calculate the percentage prediction error
output_rye <- 
  output_rye %>%
  mutate(abs_dev = abs(BEF_obs - mean_BEF) )

output_rye %>%
  group_by(Beff) %>%
  summarise(median_abs_dev = median(abs_dev),
            mean_abs_dev = mean(abs_dev),
            PI95_low = quantile(abs_dev, 0.05),
            PI95_high = quantile(abs_dev, 0.95))

output_rye %>%
  filter(Beff == "NBE")

### END
