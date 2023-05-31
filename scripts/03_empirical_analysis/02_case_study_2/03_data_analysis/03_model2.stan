//
// LogNormal-Hurdle: Model 2
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
     vector[N] Y;
     vector[N] PC1;
    array[N] int C;
    array[N] int S;
}
parameters{
     // standard deviation: log-normal model
     real<lower=0> sigma;
     // standard normal deviations: log-normal model
     matrix[2, S_N] Z;
     // parameters: log-normal model
     cholesky_factor_corr[2] L_Rho_v;
     vector<lower=0>[2] sigma_v;
     vector[2] vbar;
     // parameters: binomial model
     real a_hu;
     real b1_hu;
}
transformed parameters{
     // transformed parameters: log-normal model
     matrix[S_N, 2] v;
     vector[S_N] a;
     vector[S_N] b1;
     v = (diag_pre_multiply(sigma_v, L_Rho_v) * Z)';
     a = vbar[1] + v[,1];
     b1 = vbar[2] + v[,2];
}
model{
    // vector of means: log-normal linear model 
     vector[N] mu;
    // vector of means: binomial linear model 
     vector[N] hu;
    // standard deviation of the log-normal distribution
    sigma ~ exponential( 4 );
    // linear model priors: log-normal
    vbar ~ normal(0, 2);
    sigma_v ~ exponential( 2 );
    L_Rho_v ~ lkj_corr_cholesky( 2 );
    // linear model priors: binomial
    a_hu ~ normal(0, 2);
    b1_hu ~ normal(0, 2);
    // standard normal vectors
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = a[S[i]] + b1[S[i]] * Y[i];
        hu[i] =  a_hu + b1_hu * Y[i];
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
     matrix[2,2] Rho_v;
    Rho_v = multiply_lower_tri_self_transpose(L_Rho_v);
    for ( i in 1:N ) {
        mu[i] = a[S[i]] + b1[S[i]] * Y[i];
        hu[i] =  a_hu + b1_hu * Y[i];
    }
    for ( i in 1:N ) log_lik[i] = hurdle_lognormal_logit_lpdf(M[i] | mu[i], sigma, hu[i]);
}
