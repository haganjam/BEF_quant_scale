
# Test the Isbell_2018_sampler() function on the data from Isbell et al. (2018, Ecology Letters)

# load the test data
library(here)
test.data <- readRDS(file = here("data/Isbell_test_data.rds"))
ans.data <- readRDS(file = here("data/Isbell_test_data_solutions.rds"))

# run the test
results <- vector(length = length(test.data))
for (i in 1:length(test.data)) {

u <- Isbell_2018_sampler(data = test.data[[i]], RYe = c(0.5, 0.5))
v <- u$Beff

w <- ans.data[[i]]
x <- which(!is.na(w$Value))

y <- dplyr::near(w$Value[x], v$Value[x], tol = 0.1)

results[i] <- any(y != TRUE)

}

if ( any(results) ) { 
  
  warning("Functions do not correctly calculate biodiversity effects of the test data") 

} else { 
  
  message("Functions correctly calculate biodiversity effects on the test data")

}

# remove the test data objects
rm(test.data, ans.data)