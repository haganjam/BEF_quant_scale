
## Test the Isbell et al. (2018) partition code

# get functions
library(here)
source(here("scripts/isbell_2018_partition.R"))

# table 1A
t1a <- data.frame(sample = c(1,1,2,2),
                  time = c(1,1,1,1), 
                  place = c(1,1,2,2), 
                  species = c(1,2,1,2), 
                  M = c(200,200,100,100),
                  Y = c(200,200,0,0))

t1a.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,    100, NA, NA, NA, NA, NA, NA))

# table 1B
t1b <- data.frame(sample = c(1,1,2,2),
                  time = c(1,1,1,1), 
                  place = c(1,1,2,2), 
                  species = c(1,2,1,2), 
                  M = c(50,350,0.44,1),
                  Y = c(8.15,291.7,0.88,0))

t1b.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100, 100, 0, NA, NA, NA, NA, NA, NA))

# case 1
cs1 <- data.frame(sample = rep(c(1:4), each = 2),
                   time=c(1,1,1,1,2,2,2,2), 
                   place=c(1,1,2,2,1,1,2,2), 
                   species=c(1,2,1,2,1,2,1,2),
                   M=c(100,50,100,50,100,50,100,50), 
                   Y=c(100,0,100,0,100,0,100,0))

cs1.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,  100,   0,   100,   100,   0,   0,   0))

# case 2
cs2 <- data.frame(sample = rep(c(1:4), each = 2), 
                   time=c(1,1,1,1,2,2,2,2), 
                   place=c(1,1,2,2,1,1,2,2), 
                   species=c(1,2,1,2,1,2,1,2),
                   M=c(100,50,100,50,50,100,50,100), 
                   Y=c(100,0,100,0,0,100,0,100))

cs2.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,    100,   0,   100,   0,   100,   0,   0))

# case 3
cs3 <- 
  data.frame(sample = rep(c(1:4), each = 2),
             time=c(1,1,1,1,2,2,2,2), 
             place=c(1,1,2,2,1,1,2,2), 
             species=c(1,2,1,2,1,2,1,2),
             M=c(100,50,50,100,100,50,50,100), 
             Y=c(100,0,0,100,100,0,0,100))

cs3.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,   100,   0,   100,   0,    0,   100,   0))

# case 4
cs4 <- 
  data.frame(sample = rep(c(1:4), each = 2),
             time=c(1,1,1,1,2,2,2,2), 
             place=c(1,1,2,2,1,1,2,2), 
             species=c(1,2,1,2,1,2,1,2),
             M=c(100,50,50,100,50,100,100,50), 
             Y=c(100,0,0,100,0,100,100,0))

cs4.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   0,   100,   0,   100,   0,    0,     0,   100))

# case 5
cs5 <- data.frame(sample = rep(c(1:4), each = 2),
                   time=c(1,1,1,1,2,2,2,2), 
                   place=c(1,1,2,2,1,1,2,2), 
                   species=c(1,2,1,2,1,2,1,2),
                   M=c(75,75,75,75,75,75,75,75), 
                   Y=c(50,50,50,50,50,50,50,50))

cs5.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   100,   0,   0,    0,   0,    0,     0,    0))

# case 6
cs6 <- data.frame(sample = rep(c(1:4), each = 2),
                   time=c(1,1,1,1,2,2,2,2), 
                   place=c(1,1,2,2,1,1,2,2), 
                   species=c(1,2,1,2,1,2,1,2),
                   M=c(100,50,100,50,100,50,100,50), 
                   Y=c(50,50,50,50,50,50,50,50))

cs6.ans <- data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                      Value = c(100,   150,   -50,   -50,    0,   0,    0,     0,    0))

# write the raw data into a list
test.data <- list(t1a, t1b, cs1, cs2, cs3, cs4, cs5, cs6)

# write the answer data into a list
ans.data <- list(t1a.ans, t1b.ans, cs1.ans, cs2.ans, cs3.ans, cs4.ans, cs5.ans, cs6.ans)

# run the test
results <- vector(length = length(test.data))
for (i in 1:length(test.data)) {

  u <- Isbell_2018_sampler(data = test.data[[i]], RYe = c(0.5, 0.5))
  u
  v <- u$Beff
  v
  
  w <- ans.data[[i]]
  w
  x <- which(!is.na(w$Value))
  x
  
  y <- w$Value[x] == v$Value[x]
  y
  
  results[i] <- any(y != TRUE)
  
}

results


### END
