
#'
#' @title Function to implement Isbell et al.'s (2018, Ecology Letters) biodiversity effect partition
#' 
#' @description This script contains a set of functions that are used to implement
#' Isbell et al.'s (2018, Ecology Letters) partition of biodiversity effects. The script includes
#' helper functions to calculate the number of unique elements, n_unique(), and
#' a function to calculate raw covariance, raw_cov(). The function to calculate the
#' different biodiversity effects is Isbell_2018_part(). In addition, we implement an
#' extension that uses samples from the Dirichlet distribution to estimate expected
#' relative yields ()
#'

# function to check if a package (x) is installed
install_if <- function(x) {
  
  if(x %in% rownames(installed.packages()) == FALSE) {
    
    message(paste(x, "is required but not installed. Installing now"))
    Sys.sleep(1)
    install.packages(x)
    library(x)
    
  } else{ 
    
    library(x, character.only=TRUE)}
  
}

# install and load libraries required for these functions
install_if("dplyr")
install_if("assertthat")
install_if("gtools")

#'
#' @title n_unique()
#' 
#' @description Function that outputs the number of unique elements in a vector
#' 
#' @param x vector
#' 

n_unique <- function(x) { 
  
  # test if x is a numeric vector
  test_1 <- function(x) {
    
    is.vector(x)
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0("x is not a vector")
    
  }
  
  assertthat::assert_that(test_1(x = x))
  
  # return the number of unique elements
  return( length(unique(x)) ) 
  
  }

#' 
#' @title raw_cov()
#' 
#' @description Function to calculate raw covariance from two equal length,
#' numeric vectors
#' 
#' @param x numeric vector of length N
#' @param y numeric vector of length N
#' 

raw_cov <- function(x, y) {
  
  test_1 <- function(x, y) {
    
    assertthat::are_equal(length(x), length(y))
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(call$x, " & ", call$y, " do not have equal length")
    
  }
  
  assertthat::assert_that(test_1(x = x, y = y))
  
  # test if x and y are numeric vectors
  test_2 <- function(x, y) {
    
    (is.vector(x) & is.numeric(x)) & (is.vector(y) & is.numeric(y))
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0("either x or y are not a numeric vectors")
    
  }
  
  assertthat::assert_that(test_2(x = x, y = y))
  
  # calculate deviation from the mean of x and y
  c1 <- x - mean(x)
  c2 <- y - mean(y)
  
  # calculate the raw covariance
  c12 <- sum( (c1*c2) )/length(c1)
  
  return(c12)
  
}

#'
#' @title Isbell_2018_part()
#' 
#' @description Function to calculate Isbell et al.'s (2018, Ecology Letters) 
#' biodiversity effect partition using species mixture and monoculture data
#' across multiple times and places.
#' 
#' @param data data.frame in the format defined by Isbell et al. (2018, Ecology Letters):
#' column 1 - sample: variable specifying the unique place-time combination
#' column 2 - place: variable specifying the place
#' column 3 - time: variable specifying the time-point
#' column 4 - species: variable specifying the species name (all times and places must have all species names present)
#' column 5 - M: monoculture functioning
#' column 6 - Y: mixture function
#' @param RYe expected relative yields for the species across in all times and places:
#' numeric vector of length = (N species) and which sums to one
#' 
#' @symbol Mi - monoculture of each species
#' @symbol Yoi - observed yield of each species in mixture
#' @symbol Yo - observed mixture yield - sum(Yoi)
#' @symbol RYei - expected relative yield of each species (usually 1/n species but can be anything)
#' @symbol RYoi - observed relative yield (Yoi/Mi) i.e. measures the extent to which species i overyields
#' @symbol dRY - RYoi - RYei
#' @symbol N - n spp in mixture
#' @symbol Poi - observed proportion of each species in mixture (i.e. RYoi/sum(RYoi))
#' 

Isbell_2018_part <- function(data, RYe) {
  
  # test if the input data is a data.frame
  test_1 <- function(x) {
    
    is.data.frame(x)
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(call$x, " is not a data.frame")
    
  }
  
  assertthat::assert_that(test_1(x = data))
  
  # test if all the required columns are present
  test_2 <- function(x) {
    
    all(names(x) %in% c("sample", "place", "time", "species", "M", "Y"))
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " is missing one of the following columns: sample, place, time, species, M, Y")
    
  }
  
  assertthat::assert_that(test_2(x = data))
  
  # test if RYe is a number
  test_3 <- function(x) {
    
   is.vector(x) & is.numeric(x)
    
  }
  
  assertthat::on_failure(test_3) <- function(call, env){
    
    paste0(call$x, " is not a numeric vector")
    
  }
  
  assertthat::assert_that(test_3(x = RYe))
  
  # test if the length of the RYe vector equals the number of species
  test_4 <- function(x, y) {
    
    assertthat::are_equal(length(x), y)
    
  }
  
  assertthat::on_failure(test_4) <- function(call, env){
    
    paste0("x and y do not have equal length")
    
  }
  
  assertthat::assert_that(test_4(x = RYe, y = n_unique(data[["species"]])) )
  
  # test if the RYe vector sums to 1
  test_5 <- function(x) {
    
    sum(x) > 0.99
    
  }
  
  assertthat::on_failure(test_5) <- function(call, env){
    
    paste0("RYe does not sum to one")
    
  }
  
  assertthat::assert_that(test_5(x = RYe) )
  
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
  
  # define observed relative yields: Prevent dividing by zero
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
  
  # calculate total number of species, times and places
  N <- n_sp*n_t*n_p
  
  # 1. Net biodiversity (E10)
  NBE <- sum(df$dRY*df$M)
  
  # 2. Total complementarity (E10)
  TC <- N * mean(df$dRY) * mean(df$M)
  
  # 3. Total selection effect (E2)
  TS <- N*raw_cov(df$dRY, df$M) 
  
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
  
  # 8. Spatio-temporal insurance (E9)
  ST_df <- 
    merge(
      merge(
        merge(df, sm_t, by = c("time", "species")), 
        sm_p, by = c("place", "species")), 
      sm_s, by = "species")
  
  # derive the terms in E9 separately      
  ST_T1 <- with(ST_df, (d.Poi - d.Poi.t - d.Poi.p + d.Poi.s + mean(d.Poi)) )
  ST_T2 <- with(ST_df, (M - M.t - M.p + M.s + mean(M)))
  ST <- sum(ST_T1*ST_T2)
  
  # 9. Total insurance effect (E9)
  IT <- sum( (df$d.Poi - mean(df$d.Poi))*(df$M - mean(df$M)) )
  
  # 10. Local complementarity and local selection
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
  
  # check the internal consistency
  
  # NBE = TC + TS
  assertthat::are_equal(round(NBE, 1), round((TC + TS), 1) )
  
  # NBE = LC + LS
  assertthat::are_equal(round(NBE, 1), round((LC + LS), 1) )
  
  # TC - LC = LS - TS
  assertthat::are_equal(round((TC-LC), 1), round((LS - TS), 1) )
  
  # TS = NO + IT
  assertthat::are_equal(round(TS, 1), round((NO + IT), 1) )
  
  # IT = AS + TI + SI + ST
  assertthat::are_equal(round(IT, 1), round((AS + TI + SI + ST), 1) )

  # prepare output
  
  list(Beff = data.frame(Beff = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
                         Value = c(NBE, TC, TS, NO, IT, AS, TI, SI, ST)) ,
       L.Beff = data.frame(L.Beff = c("LC", "LS"),
                           Value = c(LC, LS)),
       RYe = RYe)
  
}

#'
#' @title Isbell_2018_sampler()
#' 
#' @description Function to calculate Isbell et al.'s (2018, Ecology Letters) 
#' biodiversity effect partition using species mixture and monoculture data
#' across multiple times and places whilst incorporating uncertainty in the RYe values
#' using the Dirichlet distribution.
#' 
#' @param data data.frame in the format defined by Isbell et al. (2018, Ecology Letters):
#' column 1 - sample: variable specifying the unique place-time combination
#' column 2 - place: variable specifying the place
#' column 3 - time: variable specifying the time-point
#' column 4 - species: variable specifying the species name (all times and places must have all species names present)
#' column 5 - M: monoculture functioning
#' column 6 - Y: mixture function
#' @param RYe expected relative yields for the species across in all times and places:
#' numeric vector of length = (N species) and which sums to one
#' @param RYe_post TRUE/FALSE deciding to calculate using samples from the Dirichlet distribution for the
#' RYe values. If TRUE, a distribution of each effects is outputted.
#' @param N number of samples to draw from the Dirichlet distribution
#' @param alpha_par alpha parameter determining the skewness of the Dirichlet distribution
#' 
#' 

Isbell_2018_sampler <- function(data, RYe_post = FALSE, N = 100, alpha_par = 4, RYe) {
  
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
      RYe <- sapply(RYe, function(x) x)
      
      x <- Isbell_2018_part(data = data, RYe = RYe)
      
      Beff_post[[i]] <- x$Beff
      
    }
    
    Beff_post <- bind_rows(Beff_post, .id = "sample")
    
    return(Beff_post)
    
  }

}

### END
