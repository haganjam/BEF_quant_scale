
#'@title: PI() 
#'
#'@description: function to calculate percentile confidence/credible interval
#'from a vector of data. The function is taken from the rethinking package
#'McElreath (2020, https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r)
#'
#'@param samples vector of data to calculate the percentile interval on
#'@param prob size of the percentile interval (default = 0.89)
#' 

PCI <- function( samples , prob=0.89 ) {
  x <- sapply( prob , function(p) {
    a <- (1-p)/2
    quantile( samples , probs=c(a,1-a) )
  } )
  # now order inside-out in pairs
  n <- length(prob)
  result <- rep(0,n*2)
  for ( i in 1:n ) {
    low_idx <- n+1-i
    up_idx <- n+i
    # lower
    result[low_idx] <- x[1,i]
    # upper
    result[up_idx] <- x[2,i]
    # add names
    a <- (1-prob[i])/2
    names(result)[low_idx] <- paste(round(a*100,0),"%")
    names(result)[up_idx] <- paste(round((1-a)*100,0),"%")
  }
  return(result)
}
PI <- PCI

#'@title: var2() 
#'
#'@description: function to calculate variance that does not use (n-1) as
#'the denominator. The function is taken from the rethinking package
#'McElreath (2020, https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r)
#'
#'@param x vector of data to calculate the percentile interval on
#' 

var2 <- function( x , na.rm = TRUE ) {
  # use E(x^2) - E(x)^2 form
  mean(x^2) - mean(x)^2
}

#'@title: HPDI() 
#'
#'@description: function to calculate the high density posterior interval
#'McElreath (2020, https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r)
#'
#'@param x vector of data to calculate the percentile interval on
#' 

# highest posterior density interval, sensu Box and Tiao
# requires coda library
HPDI <- function( samples , prob=0.89 ) {
  require(coda)
  coerce.list <- c( "numeric" , "matrix" , "data.frame" , "integer" , "array" )
  if ( inherits(samples, coerce.list) ) {
    # single chain for single variable
    samples <- coda::as.mcmc( samples )
  }
  x <- sapply( prob , function(p) coda::HPDinterval( samples , prob=p ) )
  # now order inside-out in pairs
  n <- length(prob)
  result <- rep(0,n*2)
  for ( i in 1:n ) {
    low_idx <- n+1-i
    up_idx <- n+i
    # lower
    result[low_idx] <- x[1,i]
    # upper
    result[up_idx] <- x[2,i]
    # add names
    names(result)[low_idx] <- paste("|",prob[i])
    names(result)[up_idx] <- paste(prob[i],"|")
  }
  return(result)
}

#' Sample LKJ correlation matrices.
#'
#' This function was copied from Richard McElreath's rethinking package hosted
#' at https://github.com/rmcelreath/rethinking. In turn, he appears to have
#' copied it from Ben Bolker's rLKJ function from the emdbook package, although
#' I cannot find it there (else I would have imported it).
#'
#' @param n Number of matrices to sample.
#' @param K dimenstion of matrix to sample.
#' @param eta Distribution parameter
#' @return matrix
#'
#' @importFrom stats rbeta rnorm
#' @export
rlkjcorr <- function (n, K, eta = 1) {
  
  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  #if (K == 1) return(matrix(1, 1, 1))
  
  f <- function() {
    alpha <- eta + (K - 2)/2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
    R[1,1] <- 1
    R[1,2] <- r12
    R[2,2] <- sqrt(1 - r12^2)
    if(K > 2) for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)
      
      # Draw uniformally on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])
      
      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)
    }
    return(crossprod(R))
  }
  R <- replicate( n , f() )
  if ( dim(R)[3]==1 ) {
    R <- R[,,1]
  } else {
    # need to move 3rd dimension to front, so conforms to array structure that Stan uses
    R <- aperm(R,c(3,1,2))
  }
  return(R)
}

### END
