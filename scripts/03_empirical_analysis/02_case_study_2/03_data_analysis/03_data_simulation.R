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
source("scripts/Function_plotting_theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)

# prior predictive simulation

# create a simulated dataset

# simulate all combinations of S, C and T with three replicates
df_preds <- expand.grid(REP = 1:3, 
                      S = 1:5,
                      C = 1:10,
                      T = (1:3)/3)

# remove the replicate column
df_preds <- df_sim[, -1]

# simulate Y, PC1 and PC2 variables from the uniform distribution
df_preds$PC1 <- runif(n = nrow(df_sim), 0, 1)
df_preds$PC2 <- runif(n = nrow(df_sim), 0, 1)
df_preds$Y <- runif(n = nrow(df_sim), 0, 1)

# check the simulated data
head(df_preds)
dim(df_preds)

# wrap the data simulation into a function so we can easily vary priors
data_simulation <- function(data,
                            abar_sd = 1,
                            b15bar_sd = 0.3,
                            sigma_ln = 4,
                            sigma_mvn = 4,
                            sigma_rate = 4,
                            lkj_rate = 3) {
  
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
    lm_offsets <- c("a", paste0("b", 1:5))
    
    # get the parameter names and add hu to them
    par_names <- c("abar_par", "b15bar_par", "sigma_list", "corr_list", "lm_offsets")
    if(j == "LG") {
      
      for(i in 1:length(par_names)) {
        
        assign(par_names[i], paste0(eval(parse(text = par_names[i])), "_hu"))
        
      }
      
    }
    
    # sample the alpha bar parameter
    assign(abar_par, round(rnorm(n = Ns, mean = 0, sd = abar_sd), 2))
    
    # sample the betabar parameters
    for (i in b15bar_par) {
      assign(x = i, round(rnorm(n = Ns, mean = 0, sd = b15bar_sd), 2) )
    }
    
    # sample the sigma sigma vals
    for(i in sigma_list) {
      assign(x = i, value = round(rexp(n = Ns, rate = sigma_mvn), 2) ) 
    }
    
    # get the log-normal sigma value
    sigma <- rexp(n = 1, rate = sigma_ln)
    
    # sample the correlation matrices
    for(i in corr_list) {
      
      # sample from the lkj distribution and convert to matrix
      x <- round(rlkjcorr(n = 1, K = Ns, eta = lkj_rate), 2)
      mx <- matrix(data = x, nrow = Ns, ncol = Ns)
      
      assign(x = i, value = mx)
    }
    
    # get the offsets
    for (i in 1:length(lm_offsets)) {
      
      # get a matrix of z-scores
      z <- matrix(rnorm(n = Ns*Nc, mean = 0, sd = 1), nrow = Ns, ncol = Nc)
      
      # get the sigmas
      sig <- eval(parse(text = sigma_list[i]))
      
      # get the correlation matrix
      R <- eval(parse(text = corr_list[i]))
      
      # obtain the Cholesky factor from this matrix
      L <- chol(R)
      
      # get the offsets from the z-scores
      v <- t(diag(sig) %*% L %*% z )
      
      assign(lm_offsets[i], v)
      
    }
    
  }
  
  # simulate observations from these parameters
  preds <- vector(length = nrow(data))
  for(i in 1:nrow(data)) {
    
    # get the mean prediction from the lognormal model on the natural scale
    x <- 
      with(data, 
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
      with(data, 
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
  
  # add the simulated values to the df_sim data
  data[["M"]] <- preds
  
  return(data)
  
}

# test the data simulation function
data_simulation(data = df_sim, abar_sd = 3, b15bar_sd = 0.3, 
                sigma_ln = 4, sigma_mvn = 3, lkj_rate = 3)

# try these priors and plot the resulting curves
sim_list <- 
  
  lapply(1:200, function(x) {
  
  data_simulation(df_sim,
                  abar_sd = 1, b15bar_sd = 0.25, 
                  sigma_ln = 4, sigma_mvn = 4, 
                  lkj_rate = 3)
  
})

# bind into a data.frame
sim_list <- bind_rows(sim_list, .id = "rep")

# check the summary statistics
sim_list %>%
  group_by(S) %>%
  summarise(mean_M = mean(M),
            min_M = min(M),
            max_M = max(M),
            n_zero = sum(M>0)/n())

# plot the data
p1 <- 
  ggplot(data = sim_list,
         mapping = aes(y = M, colour = as.character(S))) +
  geom_density() +
  facet_wrap(~C, scales = "free", nrow = 2, ncol = 5) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))
plot(p1)


# simulate a dataset to fit the model

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
  lm_offsets <- c("a", paste0("b", 1:5))
  
  # get the parameter names and add hu to them
  par_names <- c("abar_par", "b15bar_par", "sigma_list", "corr_list", "lm_offsets")
  if(j == "LG") {
    
    for(i in 1:length(par_names)) {
      
      assign(par_names[i], paste0(eval(parse(text = par_names[i])), "_hu"))
      
    }
    
  }
  
  # sample the alpha bar parameter
  assign(abar_par, round(rnorm(n = Ns, mean = 0, sd = 1), 2))
  
  # sample the betabar parameters
  for (i in b15bar_par) {
    assign(x = i, round(rnorm(n = Ns, mean = 0, sd = 0.3), 2) )
  }
  
  # sample the sigma sigma vals
  for(i in sigma_list) {
    assign(x = i, value = round(rexp(n = Ns, rate = 3), 2) ) 
  }
  
  # get the log-normal sigma value
  sigma <- rexp(n = 1, rate = 4)
  
  # sample the correlation matrices
  for(i in corr_list) {
    
    # sample from the lkj distribution and convert to matrix
    x <- round(rlkjcorr(n = 1, K = Ns, eta = 3), 2)
    mx <- matrix(data = x, nrow = Ns, ncol = Ns)
    
    assign(x = i, value = mx)
  }
  
  # get the offsets
  for (i in 1:length(lm_offsets)) {
    
    # get a matrix of z-scores
    z <- matrix(rnorm(n = Ns*Nc, mean = 0, sd = 1), nrow = Ns, ncol = Nc)
    
    # get the sigmas
    sig <- eval(parse(text = sigma_list[i]))
    
    # get the correlation matrix
    R <- eval(parse(text = corr_list[i]))
    
    # obtain the Cholesky factor from this matrix
    L <- chol(R)
    
    # get the offsets from the z-scores
    v <- t(diag(sig) %*% L %*% z )
    
    assign(lm_offsets[i], v)
    
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

# add the simulated values to the df_sim data
df_sim[["M"]] <- preds

# set-up the data list
# make a data.list with the training data
dat <- 
  list(N = length(df_sim$M),
       S_N = length(unique(df_sim$S)),
       C_N = length(unique(df_sim$C)),
       M = df_sim$M,
       Y = df_sim$Y,
       PC1 = df_sim$PC1,
       PC2 = df_sim$PC2,
       C = df_sim$C,
       T = df_sim$T,
       S = df_sim$S)

# compile the model
m_sim <- rstan::stan_model("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/03_model1.stan",
                        verbose = TRUE)

# sample the stan model
m_sim_fit <- rstan::sampling(m_sim, data = dat, 
                             iter = 1000, chains = 4, algorithm = c("NUTS"),
                             control = list(adapt_delta = 0.95, max_treedepth = 12),
                             seed = 54856, cores = 4)

print(m_sim_fit)

# check how similar the parameter estimates are

# check the model parameters
pars <- m_sim_fit@model_pars
pars <- pars[!(grepl("Rho", pars) | pars == "Z" | pars == "mu" | pars == "hu" | pars == "log_lik" | pars == "lp__")]

# parameter row names
pars_rows <- gsub(pattern = "\\[(.*?)\\]", replace = "", row.names(diag$summary))

# extract the diagnostic parameters
diag <- as.data.frame(rstan::summary(m_sim_fit)$summary)

# check simulated parameter
n <- 9

pars[n]
eval(parse(text = pars[n]))
pars_rows[pars[n] == pars_rows]
diag[pars[n] == pars_rows, ]$mean

cor(eval(parse(text = pars[n])), diag[pars[n] == pars_rows, ]$mean)






