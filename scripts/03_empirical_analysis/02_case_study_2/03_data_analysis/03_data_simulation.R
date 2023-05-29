#'
#' @title: Model missing monocultures
#' 
#' @description: Data simulation to test the model formulation along with
#' prior predictive simulation to choose appropriate priors.
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# import helper functions
source("scripts/03_empirical_analysis/helper_functions.R")

# load relevant libraries
library(dplyr)
library(MASS)

# create a simulated dataset

# simulate all combinations of S, C and T with three replicates
df_sim <- expand.grid(REP = 1:3, 
                      S = 1:5,
                      C = 1:10,
                      T = (1:3)/3)

# remove the replicate column
df_sim <- df_sim[, -1]

# simulate Y, PC1 and PC2 variables from the uniform distribution
df_sim$PC1 <- runif(n = nrow(df_sim), 0, 1)
df_sim$PC2 <- runif(n = nrow(df_sim), 0, 1)
df_sim$Y <- runif(n = nrow(df_sim), 0, 1)

# check the simulated data
head(df_sim)
dim(df_sim)

# data simulation

# set the number of species and clusters
Ns <- 5
Nc <- 10

# we want two sets of parameters (one for log-normal, one for logist = LG)
dist <- c("X", "LG")

# loop over the dist objects to draw one set for the lognormal
# and one set for the logistic model
for(j in dist) {
  
  # get the baseline parameters
  abar_par <- "abar"
  b15bar_par <- c("b1bar", "b2bar", "b3bar", "b4bar", "b5bar")
  sigma_list <- paste0("sigma_", c("a", "b1", "b2", "b3", "b4", "b5"))
  corr_list <- paste0("Rho_", c("a", "b1", "b2", "b3", "b4", "b5"))
  cov_list <- paste0("Cov_", c("a", "b1", "b2", "b3", "b4", "b5"))
  lm_offsets <- c("a", paste0("b", 1:5))
  
  # get the parameter names and add hu to them
  par_names <- c("abar_par", "b15bar_par", "sigma_list", "corr_list", "cov_list", "lm_offsets")
  if(j == "LG") {
    
    for(i in 1:length(par_names)) {
      
      assign(par_names[i], paste0(eval(parse(text = par_names[i])), "_hu"))
      
    }
    
  }
  
  # sample the alpha bar parameter
  assign(abar_par, round(runif(n = Ns, 0, 5), 1))
  
  # sample the betabar parameters
  for (i in b15bar_par) {
    assign(x = i, round(rnorm(n = Ns, mean = 0, sd = 2), 1) )
  }
  
  # sample the sigma sigma vals
  for(i in sigma_list) {
    assign(x = i, value = round(rexp(n = Ns, rate = 1), 1) ) 
  }
  
  # sample the correlation matrices
  for(i in corr_list) {
    x <- round(rlkjcorr(n = 1, K = Ns, eta = 2), 1)
    mx <- matrix(data = x, nrow = Ns, ncol = Ns)
    assign(x = i, value = mx)
  }
  
  # obtain the covariance matrices as follows
  for(i in 1:length(cov_list)) {
    
    x <- eval(parse(text = sigma_list[i]))
    y <- eval(parse(text = corr_list[i]))
    
    assign(cov_list[i], diag(x) %*% y %*% diag(x))
    
  }
  
  # get the sigma value for the log-normal distribution
  sigma <- rexp(n = 1, rate = 1)
  
  # draw varying effects from the multivariate normal distribution
  
  # get a list of linear model parameters
  lm_mean_pars <- c(abar_par, b15bar_par)
  
  for(i in 1:length(lm_mean_pars)) {
    
    # get the multivariate mean values
    mu <- eval(parse(text = lm_mean_pars[i]))
    
    # get the covariance matrix
    var_cov <- eval(parse(text = cov_list[i]))
    
    # sample the effects from the multivariate normal
    z <- mvrnorm( Nc , mu , var_cov )
    
    # assign to an object
    assign(lm_offsets[i], z)
    
  }
  
  
}


# simulate observations from these parameters
preds <- vector(length = nrow(df_sim))
for(i in 1:nrow(df_sim)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_sim, 
         (abar[S[i]] + a[C[i], S[i]]) + 
           (b1bar[S[i]] + b1[C[i], S[i]] * T[i]) + 
           (b2bar[S[i]] + b2[C[i], S[i]] * Y[i]) + 
           (b3bar[S[i]] + b3[C[i], S[i]] * PC1[i]) + 
           (b4bar[S[i]] + b4[C[i], S[i]] * PC2[i]) + 
           (b5bar[S[i]] + b5[C[i], S[i]] * Y[i] * PC1[i]))

  # draw from the log-normal distribution
  x <- rlnorm(n = length(x), x, sigma)

  # get the probability of 0
  y <- 
    with(df_sim, 
         (abar_hu[S[i]] + a_hu[C[i], S[i]]) + 
           (b1bar_hu[S[i]] + b1_hu[C[i], S[i]] * T[i]) + 
           (b2bar_hu[S[i]] + b2_hu[C[i], S[i]] * Y[i]) + 
           (b3bar_hu[S[i]] + b3_hu[C[i], S[i]] * PC1[i]) + 
           (b4bar_hu[S[i]] + b4_hu[C[i], S[i]] * PC2[i]) + 
           (b5bar_hu[S[i]] + b5_hu[C[i], S[i]] * Y[i] * PC1[i]) )
  y <- plogis(y)
  y <- 1-y
  
  # draw from the binomial distribution
  y <- rbinom(n = length(y), size = 1, prob = y)
  
  preds[i] <- (x*y)
  
}

# add the simulated values to the 











cov_ab <- sigma_a*sigma_b*rho
Sigma <- matrix( c(sigma_a^2,cov_ab,cov_ab,sigma_b^2) , ncol=2 )



# sample the priors

# set the number of samples to get from the priors
N <- 100

# log-normal model

# linear model priors
ln_lm <- c("abar", "b1bar", "b2bar", "b3bar", "b4bar", "b5bar")
for(i in ln_lm) {
  assign(x = i, value = rnorm(n = N, mean = 0, sd = 2))
}

# sigma priors
ln_sigma <- paste0("sigma_", c("a", "b1", "b2", "b3", "b4", "b5"))
for(i in ln_sigma) {
  assign(x = i, value = rexp(n = N, rate = 1) ) 
}

# cholesky factors of the lkj prior

# sample from lkj prior
x <- rlkjcorr(n = 1, K = 5, eta = 2)

# convert to a 5 x 5 matrix
mx <- matrix(data = x, nrow = 5, ncol = 5)

# get the Cholesky factor of the matrix
cx <- chol(x)

# convert back to a correlation matrix
crossprod(cx) 

# cholesky decomposition
# now for something positive semi-definite





L_Rho_a ~ lkj_corr_cholesky( 2 );
L_Rho_b1 ~ lkj_corr_cholesky( 2 );
L_Rho_b2 ~ lkj_corr_cholesky( 2 );
L_Rho_b3 ~ lkj_corr_cholesky( 2 );
L_Rho_b4 ~ lkj_corr_cholesky( 2 );
L_Rho_b5 ~ lkj_corr_cholesky( 2 );
// linear model priors: binomial model
abar_hu ~ normal( 0 , 2 );
b1bar_hu ~ normal( 0 , 2 );
b2bar_hu ~ normal( 0 , 2 );
b3bar_hu ~ normal( 0 , 2 );
b4bar_hu ~ normal( 0 , 2 );
b5bar_hu ~ normal( 0 , 2 );
sigma_a_hu ~ exponential( 1 );
sigma_b1_hu ~ exponential( 1 );
sigma_b2_hu ~ exponential( 1 );
sigma_b3_hu ~ exponential( 1 );
sigma_b4_hu ~ exponential( 1 );
sigma_b5_hu ~ exponential( 1 );
L_Rho_a_hu ~ lkj_corr_cholesky( 2 );
L_Rho_b1_hu ~ lkj_corr_cholesky( 2 );
L_Rho_b2_hu ~ lkj_corr_cholesky( 2 );
L_Rho_b3_hu ~ lkj_corr_cholesky( 2 );
L_Rho_b4_hu ~ lkj_corr_cholesky( 2 );
L_Rho_b5_hu ~ lkj_corr_cholesky( 2 );



# use the posterior predictive distribution
ppd <- TRUE

pred_mod <- vector("list", length = nrow(df_pred))
for(i in 1:nrow(df_pred)) {
  
  # get the mean prediction from the lognormal model on the natural scale
  x <- 
    with(df_pred, 
         (abar[-id, S[i]] + a[, , S[i]][-id, C[i]]) + 
           (b1bar[-id, S[i]] + b1[, , S[i]][-id, C[i]] * T[i]) + 
           (b2bar[-id, S[i]] + b2[, , S[i]][-id, C[i]] * Y[i]) + 
           (b3bar[-id, S[i]] + b3[, , S[i]][-id, C[i]] * PC1[i]) + 
           (b4bar[-id, S[i]] + b4[, , S[i]][-id, C[i]] * PC2[i]) + 
           (b5bar[-id, S[i]] + b5[, , S[i]][-id, C[i]] * Y[i] * PC1[i]))
  if(ppd) {
    x <- rlnorm(n = length(x), x, sigma)
  } else {
    x <- exp(x + (0.5*(sigma^2)))
  }
  
  # get the probability of 0
  y <- 
    with(df_pred, 
         (abar_hu[-id, S[i]] + a_hu[, , S[i]][-id, C[i]]) + 
           (b1bar_hu[-id, S[i]] + b1_hu[, , S[i]][-id, C[i]] * T[i]) + 
           (b2bar_hu[-id, S[i]] + b2_hu[, , S[i]][-id, C[i]] * Y[i]) + 
           (b3bar_hu[-id, S[i]] + b3_hu[, , S[i]][-id, C[i]] * PC1[i]) + 
           (b4bar_hu[-id, S[i]] + b4_hu[, , S[i]][-id, C[i]] * PC2[i]) + 
           (b5bar_hu[-id, S[i]] + b5_hu[, , S[i]][-id, C[i]] * Y[i] * PC1[i]) )
  y <- plogis(y)
  y <- 1-y
  if(ppd) {
    y <- rbinom(n = length(y), size = 1, prob = y)
  }
  
  pred_mod[[i]] <- (x*y)
  
}

# bind into a matrix
pred_mod <- do.call("cbind", pred_mod)

# pull into a data.frame for plotting
C <- factor(data$cluster_id)
levels(C) <- paste0("Cluster ", LETTERS[1:10])
df_plot <- data.frame(M_obs = data$M,
                      Obs_pred = ifelse(is.na(data$M), "Predicted", "Observed"),
                      S = as.character(df_pred$S),
                      C = C,
                      T = as.character(data$time),
                      Y = data$Y,
                      PC1 = data$PC1,
                      PC2 = data$PC2,
                      M_pred_mu = apply(pred_mod, 2, mean),
                      M_pred_PIlow = apply(pred_mod, 2, HPDI, 0.90)[1,],
                      M_pred_PIhigh = apply(pred_mod, 2, HPDI, 0.90)[2,])


# check the min and max
summary(df_plot)

# factor OTU order: Barn Bryo Bumpi Hydro Seasq

p1 <- 
  ggplot(data = df_pred,
         mapping = aes(x = M_obs, y = M_pred_mu, colour = M_S)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  geom_errorbar(mapping = aes(x = M_obs, ymin = M_pred_PIlow, ymax = M_pred_PIhigh),
                width = 0, alpha = 0.5, size = 0.25, show.legend = FALSE) +
  geom_point(shape = 16, alpha = 0.75) +
  ylab("Predicted monoculture (g)") +
  xlab("Observed monoculture (g)") +
  facet_wrap(~M_C, scales = "free", nrow = 2, ncol = 5) +
  scale_y_continuous(limits = c(0, 23)) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_colour_manual(name = "OTU",
                      labels = c("Barn", "Bryo", "Asci", "Hydro", "Ciona"),
                      values = viridis(n = 5, begin = 0.1, end = 0.9, option = "C")) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p1)
