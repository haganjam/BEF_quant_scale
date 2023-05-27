//
// LogNormal-Hurdle: Model 4
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
     vector[N] T;
    array[N] int C;
    array[N] int S;
}
parameters{
     // standard deviation: log-normal model
     real<lower=0> sigma;
     // standard normal deviations
     matrix[S_N,C_N] Z;
     // parameters: log-normal model
     cholesky_factor_corr[S_N] L_Rho_b1;
     cholesky_factor_corr[S_N] L_Rho_a;
     vector<lower=0>[S_N] sigma_b1;
     vector<lower=0>[S_N] sigma_a;
     vector[S_N] b1bar;
     vector[S_N] abar;
     // parameters: binomial model
     cholesky_factor_corr[S_N] L_Rho_b1_hu;
     cholesky_factor_corr[S_N] L_Rho_a_hu;
     vector<lower=0>[S_N] sigma_b1_hu;
     vector<lower=0>[S_N] sigma_a_hu;
     vector[S_N] b1bar_hu;
     vector[S_N] abar_hu;
}
transformed parameters{
     // transformed parameters: log-normal model
     matrix[C_N,S_N] a;
     matrix[C_N,S_N] b1;
    b1 = (diag_pre_multiply(sigma_b1, L_Rho_b1) * Z)';
    a = (diag_pre_multiply(sigma_a, L_Rho_a) * Z)';
     // transformed parameters: binomial model
     matrix[C_N,S_N] a_hu;
     matrix[C_N,S_N] b1_hu;
    b1_hu = (diag_pre_multiply(sigma_b1_hu, L_Rho_b1_hu) * Z)';
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
    b1bar ~ normal( 0 , 2 );
    sigma_a ~ exponential( 1 );
    sigma_b1 ~ exponential( 1 );
    L_Rho_a ~ lkj_corr_cholesky( 2 );
    L_Rho_b1 ~ lkj_corr_cholesky( 2 );
    // linear model priors: binomial model
    abar_hu ~ normal( 0 , 2 );
    b1bar_hu ~ normal( 0 , 2 );
    sigma_a_hu ~ exponential( 1 );
    sigma_b1_hu ~ exponential( 1 );
    L_Rho_a_hu ~ lkj_corr_cholesky( 2 );
    L_Rho_b1_hu ~ lkj_corr_cholesky( 2 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = (abar[S[i]] + a[C[i], S[i]]) + (b1bar[S[i]] + b1[C[i], S[i]] * T[i]);
        hu[i] = (abar_hu[S[i]] + a_hu[C[i], S[i]]) + (b1bar_hu[S[i]] + b1_hu[C[i], S[i]] * T[i]);
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
     matrix[S_N,S_N] Rho_b1;
    Rho_b1 = multiply_lower_tri_self_transpose(L_Rho_b1);
    Rho_a = multiply_lower_tri_self_transpose(L_Rho_a);
     // correlation matrices: binomial model
     matrix[S_N,S_N] Rho_a_hu;
     matrix[S_N,S_N] Rho_b1_hu;
    Rho_b1_hu = multiply_lower_tri_self_transpose(L_Rho_b1_hu);
    Rho_a_hu = multiply_lower_tri_self_transpose(L_Rho_a_hu);
    for ( i in 1:N ) {
        mu[i] = (abar[S[i]] + a[C[i], S[i]]) + (b1bar[S[i]] + b1[C[i], S[i]] * T[i]);
        hu[i] = (abar_hu[S[i]] + a_hu[C[i], S[i]]) + (b1bar_hu[S[i]] + b1_hu[C[i], S[i]] * T[i]);
    }
    for ( i in 1:N ) log_lik[i] = hurdle_lognormal_logit_lpdf(M[i] | mu[i], sigma, hu[i]);
}
