
## Quantifying the effects of biodiversity across times and places


# Define function for number of unique elements in a vector

# args
# x - vector
n_unique <- function(x) length(unique(x))


# Define function to calculate raw covariance

# args
# x - vector
# y - vector
raw_cov <- function(x, y) {
  if(length(x) != length(y)) {
    stop("Vector x and y differ in length")
  }
  c1 <- x - mean(x)
  c2 <- y - mean(y)
  sum( (c1*c2) )/length(c1)
}


# Define a function to perform Isbell et al.'s (2018) partition

# args
# data - data.frame in the format defined by Isbell et al. (2018)
# RYe - expected frequencies for each species

# symbol descriptions
# Mi = monoculture of each species
# Yoi = observed yield of each species in mixture
# Yo = observed mixture yield - sum(Yoi)
# RYei = expected relative yield of each species (usually 1/n species but can be anything)
# RYoi = observed relative yield (Yoi/Mi) i.e. measures the extent to which species i overyields
# dRY = RYoi - RYei
# N = n spp in mixture
# Poi = observed proportion of each species in mixture (i.e. RYoi/sum(RYoi))

Isbell_2018_part <- function(data, RYe) {
  
  ## Prepare
  
  # libraries to load
  library(dplyr)
  
  # get the number of species in the data
  n_sp <- n_unique(data$species)
  
  # get the number of times
  n_t <- n_unique(data$time)
  
  # get the number places
  n_p <- n_unique(data$place)
  
  # sort the data.frame
  df <- 
    data %>%
    arrange(sample, time, place, species)
  
  # define expected relative yields
  df$RYe <- rep(RYe, n_unique(df$sample))
  
  # define observed relative yields
  df$RYo <- ifelse(df$M == 0, 0, (df$Y/df$M))
  
  # define the change in relative yield
  df$dRY <- (df$RYo - df$RYe)
  
  # calculate expected yield for each species
  df$Ye <- (df$M*df$RYe)
  
  # calculate observed proportion of each species in mixture (po,ijk, Isbell et al. 2018)
  df <- 
    df %>%
    group_by(sample) %>%
    mutate(Poi = if_else(is.na(Y/sum(Y)), 0, Y/sum(Y))  ) %>%
    ungroup()
  
  # calculate change in observed proportion relative to the expectation (d.po,ijk, Isbell et al. 2018)
  df$d.Poi <- (df$Poi - df$RYe)
  
  # calculate change in observed proportion (dRYo,ijk Isbell et al. 2018)
  df$d.RYoi <- (df$RYo - df$Poi)
  
  # add means from different hierarchical levels
  # species means for each time across places (pij and Mij single bar, Isbell et al. 2018)
  sm_t <- aggregate(df[, c("d.Poi", "M") ], list(df$time, df$species), mean)
  names(sm_t) <- c("time", "species", "d.Poi.t", "M.t")
  
  # species means for each place across times (pik and Mik single bar, Isbell et al. 2018)
  sm_p <- aggregate(df[, c("d.Poi", "M") ], list(df$place, df$species), mean)
  names(sm_p) <- c("place", "species", "d.Poi.p", "M.p")
  
  # overall species mean across all times and places
  sm_s <- aggregate(df[, c("d.Poi", "M") ], list(df$species), mean)
  names(sm_s) <- c("species", "d.Poi.s", "M.s")
  
  
  ## Calculate
  
  # calculate total number of species, times and places
  N <- n_sp*n_t*n_p
  
  # 1. Net biodiversity (E10)
  NBE <- sum(df$dRY*df$M)
  
  # 2. Total complementarity (E10)
  TC <- N * mean(df$dRY) * mean(df$M)
  
  # 3. Total selection effect (Fig. 1)
  TS <- NBE - TC
  
  # 4. Non-random overyielding (E10)
  NO <- N * raw_cov(df$d.RYoi, df$M)
  
  # 5. Average selection (E9)
  AS <- n_t * n_p * sum((sm_s$d.Poi.s - mean(df$d.Poi)) * (sm_s$M.s - mean(df$M)))
  
  # 6. Temporal insurance (E9)
  TI_df <- merge(sm_t, sm_s)
  TI <- n_p * sum( (TI_df$d.Poi.t - TI_df$d.Poi.s)*(TI_df$M.t - TI_df$M.s) )
  
  # 7. Spatial insurance (E9)
  SI_df <- merge(sm_p, sm_s)
  SI <- n_t * sum( (SI_df$d.Poi.p - SI_df$d.Poi.s)*(SI_df$M.p - SI_df$M.s) )
  
  # 8. Spatio-temporal insurance (Fig. 1)
  ST <- NBE - TC - NO - AS - TI - SI
  
  # 9. Total insurance effect (Fig. 1)
  IT <- AS + TI + SI + ST
  
  # Local complementarity and local selection
  LC_LS <- 
    df %>%
    group_by(sample) %>%
    summarise(dRY_m = mean(dRY),
              M_m = mean(M),
              cov_m = raw_cov(dRY, M),
              n = n()) %>%
    mutate(LC = n*dRY_m*M_m,
           LS = n*cov_m) %>%
    summarise(LC = sum(LC),
              LS = sum(LS))
  
  # 10. Local complementarity
  LC <- LC_LS$LC
  
  # 11. Local selection
  LS <- LC_LS$LS
  
  
  ## Output
  
  list(Beff = data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                         Value = c(NBE, TC, TS, NO, IT, AS, TI, SI, ST)) ,
       L.Beff = data.frame(L.Beff = c("LC", "LS"),
                           Value = c(LC, LS)),
       RYe = RYe)
  
}

# Use the Isbell_2018_part and incorporate RYe uncertainty

# args
# data - data.frame in the format defined by Isbell et al. (2018)
# RYe - expected frequencies for each species
# RYe_post - TRUE or FALSE deciding to calculate effects using samples from the Dirichlet distrbution
# N - number of samples to use from the Dirichlet distribution
# alpha_par - alpha parameter from the Dirichlet distribution

Isbell_2018_sampler <- function(data, RYe, RYe_post = FALSE, N = 100, alpha_par = 4) {
  
  # libraries to load
  library(dplyr)
  
  if (!RYe_post) {
    
    # test for RYe values
    if (any(!is.na(RYe)) & !dplyr::near(sum(RYe), 1) ) {
      stop("Expected relative yield values do not sum to 1")
    } 
    
    if (n_unique(data$species) != length(RYe) | any(is.na(RYe))) {
      stop("Expected relative yield values are missing")
    }
    
    Isbell_2018_part(data = data, RYe = RYe)
    
    
  } else if (RYe_post) {
    
    Beff_post <- vector("list", length = N)
    for (i in 1:N) {
      
      RYe <- gtools::rdirichlet(n = 1, alpha = rep(alpha_par, n_unique(data$species)))
      RYe <- round(sapply(RYe, function(x) x), 1)
      
      x <- Isbell_2018_part(data = data, RYe = RYe)
      
      Beff_post[[i]] <- x$Beff
      
    }
    
    Beff_post <- bind_rows(Beff_post, .id = "sample")
    
    return(Beff_post)
    
  }

}

### END
