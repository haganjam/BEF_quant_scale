
## Define example data from Isbell et al. (2018) and Fox (2005)

# Isbell et al. (2018)

# table 1A
t1a <- data.frame(sample = c(1,1,2,2),
                  time = c(1,1,1,1), 
                  place = c(1,1,2,2), 
                  species = c(1,2,1,2), 
                  M = c(200,200,100,100),
                  Y = c(200,200,0,0))

# table 1B
t1b <- data.frame(sample = c(1,1,2,2),
                  time = c(1,1,1,1), 
                  place = c(1,1,2,2), 
                  species = c(1,2,1,2), 
                  M = c(50,350,0.44,1),
                  Y = c(8.15,291.7,0.88,0))

# case 1
cs.1 <- data.frame(sample = rep(c(1:4), each = 2),
                   time=c(1,1,1,1,2,2,2,2), 
                   place=c(1,1,2,2,1,1,2,2), 
                   species=c(1,2,1,2,1,2,1,2),
                   M=c(100,50,100,50,100,50,100,50), 
                   Y=c(100,0,100,0,100,0,100,0))

# case 2
cs.2 <- data.frame(sample = rep(c(1:4), each = 2), 
                   time=c(1,1,1,1,2,2,2,2), 
                   place=c(1,1,2,2,1,1,2,2), 
                   species=c(1,2,1,2,1,2,1,2),
                   M=c(100,50,100,50,50,100,50,100), 
                   Y=c(100,0,100,0,0,100,0,100))

# case 3
cs.3 <- 
  data.frame(sample = rep(c(1:4), each = 2),
             time=c(1,1,1,1,2,2,2,2), 
             place=c(1,1,2,2,1,1,2,2), 
             species=c(1,2,1,2,1,2,1,2),
             M=c(100,50,50,100,100,50,50,100), 
             Y=c(100,0,0,100,100,0,0,100))

# case 4
cs.4 <- 
  data.frame(sample = rep(c(1:4), each = 2),
             time=c(1,1,1,1,2,2,2,2), 
             place=c(1,1,2,2,1,1,2,2), 
             species=c(1,2,1,2,1,2,1,2),
             M=c(100,50,50,100,50,100,100,50), 
             Y=c(100,0,0,100,0,100,100,0))

# case 5
cs.5 <- data.frame(sample = rep(c(1:4), each = 2),
                   time=c(1,1,1,1,2,2,2,2), 
                   place=c(1,1,2,2,1,1,2,2), 
                   species=c(1,2,1,2,1,2,1,2),
                   M=c(75,75,75,75,75,75,75,75), 
                   Y=c(50,50,50,50,50,50,50,50))

# case 6
cs.6 <- data.frame(sample = rep(c(1:4), each = 2),
                   time=c(1,1,1,1,2,2,2,2), 
                   place=c(1,1,2,2,1,1,2,2), 
                   species=c(1,2,1,2,1,2,1,2),
                   M=c(100,50,100,50,100,50,100,50), 
                   Y=c(50,50,50,50,50,50,50,50))


# Fox (2005)
f1 <- data.frame(sample = rep(c(1, 2, 3, 4, 5, 6), each = 2),
                 species = rep(c(1, 2), 6),
                 M = rep(c(500, 250), 6),
                 Y = c(300, 100, 330, 110, 360, 120, 390, 130, 420, 140, 450, 150))

### END
