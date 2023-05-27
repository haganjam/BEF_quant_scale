//
// LogNormal-Hurdle: Model 5
//
functions {
  /* hurdle lognormal log-PDF of a single response
   * taken from brms()
   * Args:
   *   y: the response value
   *   mu: mean parameter of the lognormal distribution
   *   sigma: sd parameter of the lognormal distribution
   *   hu: hurdle probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_lognormal_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             lognormal_lpdf(y | mu, sigma);
    }
  }
  /* hurdle lognormal log-PDF of a single response
   * logit parameterization of the hurdle part
   * taken from brms()
   * Args:
   *   y: the response value
   *   mu: mean parameter of the lognormal distribution
   *   sigma: sd parameter of the lognormal distribution
   *   hu: linear predictor for the hurdle part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_lognormal_logit_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | hu);
    } else {
      return bernoulli_logit_lpmf(0 | hu) +
             lognormal_lpdf(y | mu, sigma);
    }
  }
  // hurdle lognormal log-CCDF and log-CDF functions
  real hurdle_lognormal_lccdf(real y, real mu, real sigma, real hu) {
    return bernoulli_lpmf(0 | hu) + lognormal_lccdf(y | mu, sigma);
  }
  real hurdle_lognormal_lcdf(real y, real mu, real sigma, real hu) {
    return log1m_exp(hurdle_lognormal_lccdf(y | mu, sigma, hu));
  }
}
data{
     int<lower=1> N;
     int<lower=1> S_N;
     int<lower=1> C_N;
     vector[N] M;
    array[N] int C;
    array[N] int S;
}
parameters{
     // standard deviation: log-normal model
     real<lower=0> sigma;
     // standard normal deviations
     matrix[S_N,C_N] Z;
     // parameters: log-normal model
     cholesky_factor_corr[S_N] L_Rho_a;
     vector<lower=0>[S_N] sigma_a;
     vector[S_N] abar;
     // parameters: binomial model
     cholesky_factor_corr[S_N] L_Rho_a_hu;
     vector<lower=0>[S_N] sigma_a_hu;
     vector[S_N] abar_hu;
}
transformed parameters{
     // transformed parameters: log-normal model
     matrix[C_N,S_N] a;
    a = (diag_pre_multiply(sigma_a, L_Rho_a) * Z)';
     // transformed parameters: binomial model
     matrix[C_N,S_N] a_hu;
    a_hu = (diag_pre_multiply(sigma_a_hu, L_Rho_a_hu) * Z)';
}
model{
    // vector of means: log-normal linear model 
     vector[N] mu;
    // vector of means: binomial linear model 
     vector[N] hu;
    // standard deviation of the log-normal distribution
    sigma ~ exponential( 1 );
    // linear model priors: log-normal
    abar ~ normal( 0 , 2 );
    sigma_a ~ exponential( 1 );
    L_Rho_a ~ lkj_corr_cholesky( 2 );
    // linear model priors: binomial model
    abar_hu ~ normal( 0 , 2 );
    sigma_a_hu ~ exponential( 1 );
    L_Rho_a_hu ~ lkj_corr_cholesky( 2 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = (abar[S[i]] + a[C[i], S[i]]);
        hu[i] = (abar_hu[S[i]] + a_hu[C[i], S[i]]);
        target += hurdle_lognormal_logit_lpdf(M[i] | mu[i], sigma, hu[i]);
    }
}
generated quantities{
     // log-likelihood vector
     vector[N] log_lik;
     // mu vector: lognormal model
     vector[N] mu;
     // hu vector: binomial model
     vector[N] hu;
     // correlation matrices: lognormal model
     matrix[S_N,S_N] Rho_a;
    Rho_a = multiply_lower_tri_self_transpose(L_Rho_a);
     // correlation matrices: binomial model
     matrix[S_N,S_N] Rho_a_hu;
    Rho_a_hu = multiply_lower_tri_self_transpose(L_Rho_a_hu);
    for ( i in 1:N ) {
        mu[i] = (abar[S[i]] + a[C[i], S[i]]);
        hu[i] = (abar_hu[S[i]] + a_hu[C[i], S[i]]);
    }
    for ( i in 1:N ) log_lik[i] = hurdle_lognormal_logit_lpdf(M[i] | mu[i], sigma, hu[i]);
}
