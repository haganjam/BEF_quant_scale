#'
#' @title Simulate metacommunities to test the analytical pipeline
#' 
#' @description Here we test whether assuming different RYE proportions allows
#' the accurate estimation of the example effects reported in Isbell et al.
#' (2018)
#' 

# load the test data
test_data <- readRDS(file = "data/Isbell_test_data.rds")
ans_data <- readRDS(file = "data/Isbell_test_data_solutions.rds")

# load initial relative abundance data
start_RA <- sapply(1:100, function(x) gtools::rdirichlet(n = 1, rep(3,  2) ) )

# load the functions
source("scripts/01_partition_functions/01_isbell_2018_partition.R")

# create an output list
output_list <- vector("list", length = length(test_data))
for(i in 1:length(test_data)) {
  
  RYe_reps <- 
    
    apply(
      
      X = start_RA, 
      
      MARGIN = 2, 
      
      FUN = function(RA) {
        
        # calculate te biodiversity effects for each of the potential starting abundances
        BEF_post <- Isbell_2018_sampler(data = test_data[[i]], RYe = RA, RYe_post = FALSE)
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
  names(ans_data[[i]]) <- c("Beff", "BEF_obs")
  
  # join these datasets
  output <- full_join(output, ans_data[[i]], by = "Beff")
  
  # write output to the list
  output_list[[i]] <- output
  
}

# bind into a data.frame
output_df <- bind_rows(output_list, .id = "test_dataset")

# get the subset of values that are affected by the RYEs
output_rye <- 
  output_df %>%
  filter(Beff %in% c("AS", "IT", "LS", "NBE", "TS")) %>%
  filter(!is.na(BEF_obs))

# get a table of the range of the different BEF effects
output_rye %>%
  group_by(Beff) %>%
  summarise(min_BEF = min(BEF_obs),
            max_BEF = max(BEF_obs))

# does the observed value lie in the 90% percentile interval
output_rye %>%
  mutate(within = ifelse(BEF_obs <= PI90_high & BEF_obs >= PI90_low, 1, 0)) %>%
  group_by(Beff) %>%
  summarise(prop_within = sum(within),
            n = n())

### END
