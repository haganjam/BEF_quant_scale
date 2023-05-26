//
// LogNormal hurdle model
//

functions {
  /* hurdle lognormal log-PDF of a single response
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
     vector[N] M;
    array[N] int C;
    array[N] int S;
}
parameters{
  matrix[5,10] Z;
  cholesky_factor_corr[5] L_Rho_a;
  vector<lower=0>[5] sigma_a;
  vector[5] abar;
  cholesky_factor_corr[5] L_Rho_a_hu;
  vector<lower=0>[5] sigma_a_hu;
  vector[5] abar_hu;
  real<lower=0> sigma;
}
transformed parameters{
     matrix[10,5] a;
     matrix[10,5] a_hu;
    a = (diag_pre_multiply(sigma_a, L_Rho_a) * Z)';
    a_hu = (diag_pre_multiply(sigma_a_hu, L_Rho_a_hu) * Z)';
}
model {
  vector[N] mu;
  vector[N] hu;
    sigma ~ exponential( 1 );
    abar ~ normal( 0 , 2 );
    sigma_a ~ exponential( 1 );
    L_Rho_a ~ lkj_corr_cholesky( 2 );
    abar_hu ~ normal( 0 , 2 );
    sigma_a_hu ~ exponential( 1 );
    L_Rho_a_hu ~ lkj_corr_cholesky( 2 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for (n in 1:N) {
      mu[n] = (abar[S[n]] + a[C[n], S[n]]);
      hu[n] = (abar_hu[S[n]] + a_hu[C[n], S[n]]);
      target += hurdle_lognormal_logit_lpdf(M[n] | mu[n], sigma, hu[n]);
    }
}
generated quantities{
    vector[N] log_lik;
    vector[N] mu;
    vector[N] hu;
    matrix[5,5] Rho_a;
    matrix[5,5] Rho_a_hu;
    Rho_a = multiply_lower_tri_self_transpose(L_Rho_a);
    Rho_a_hu = multiply_lower_tri_self_transpose(L_Rho_a_hu);
    for ( n in 1:N ) {
        mu[n] = (abar[S[n]] + a[C[n], S[n]]);
        hu[n] = (abar_hu[S[n]] + a_hu[C[n], S[n]]);
    }
    for ( n in 1:N ) log_lik[n] = hurdle_lognormal_logit_lpdf(M[n] | mu[n], sigma, hu[n]);
}