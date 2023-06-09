#'
#' @title Test the Isbell_2018_sampler() function on the data from Isbell et al. (2018, Ecology Letters)
#' 
#' @description Here, we use the functions that we wrote for calculating the biodiversity
#' effects proposed by Isbell et al. (2018, Ecology Letters) on data test data proposed
#' in the same data to make sure that our functions accurately calculate the effects
#' 

# load the test data
test.data <- readRDS(file = "data/Isbell_test_data.rds")
ans.data <- readRDS(file = "data/Isbell_test_data_solutions.rds")

# load the functions
source("scripts/01_partition_functions/01_isbell_2018_partition.R")

# run the test
results <- vector(length = length(test.data))
for (i in 1:length(test.data)) {
  
  u <- Isbell_2018_part(data = test.data[[i]], RYe_equal = TRUE, RYe = c(0.5, 0.5))
  v <- u$Beff
  print(v)
  
  w <- ans.data[[i]]
  x <- which(!is.na(w$Value))
  
  y <- dplyr::near(w$Value[x], v$Value[x], tol = 0.1)
  
  results[i] <- any(y != TRUE)
  
}

# check if any of the biodiversity effects were incorrectly calculated
if ( any(results) ) { 
  
  warning("Functions do not correctly calculate biodiversity effects of the test data") 

} else { 
  
  message("Functions correctly calculate biodiversity effects on the test data")

}

# remove the test data objects
rm(test.data, ans.data)

### END
