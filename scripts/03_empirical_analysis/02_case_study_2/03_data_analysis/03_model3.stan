//
// LogNormal-Hurdle: Model 3
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
     matrix[S_N,C_N] Za;
     vector[S_N] Zb1;
     vector[S_N] Zb2;
     vector[S_N] Zb3;
     // parameters: log-normal model
     cholesky_factor_corr[S_N] L_Rho_a;
     vector<lower=0>[S_N] sigma_a;
     vector[S_N] abar;
     real<lower=0> sigma_b1;
     real b1bar;
     real<lower=0> sigma_b2;
     real b2bar;
     // parameters: binomial model
     vector[S_N] Za_hu;
     vector[S_N] Zb1_hu;
     vector[S_N] Zb2_hu;
     real<lower=0> sigma_a_hu;
     real abar_hu;
     real<lower=0> sigma_b1_hu;
     real b1bar_hu;
}
transformed parameters{
     // transformed parameters: log-normal model
     matrix[C_N,S_N] a;
     vector[S_N] b1;
     vector[S_N] b2;
     a = (diag_pre_multiply(sigma_a, L_Rho_a) * Za)';
     b1 = b1bar + (Zb1*sigma_b1);
     b2 = b2bar + (Zb2*sigma_b2);
     // transformed parameters: binomial model
     vector[S_N] a_hu;
     vector[S_N] b1_hu;
     a_hu = abar_hu + (Za_hu*sigma_a_hu);
     b1_hu = b1bar_hu + (Zb1_hu*sigma_b1_hu);
}
model{
    // vector of means: log-normal linear model 
     vector[N] mu;
    // vector of means: binomial linear model 
     vector[N] hu;
    // standard deviation of the log-normal distribution
    sigma ~ exponential( 4 );
    // linear model priors: log-normal
    abar ~ uniform( -4 , 4 );
    sigma_a ~ exponential( 2 );
    L_Rho_a ~ lkj_corr_cholesky( 2 );
    b1bar ~ normal(0, 1);
    sigma_b1 ~ exponential( 2 );
    b2bar ~ normal(0, 1);
    sigma_b2 ~ exponential( 2 );
    // linear model priors: binomial
    abar_hu ~ normal( 0 , 1 );
    sigma_a_hu ~ exponential( 2 );
    b1bar_hu ~ normal(0, 1);
    sigma_b1_hu ~ exponential( 2 );
    // standard normal vectors
    to_vector( Za ) ~ normal( 0 , 1 );
    to_vector( Zb1 ) ~ normal( 0 , 1 );
    to_vector( Zb2 ) ~ normal( 0 , 1 );
    to_vector( Za_hu ) ~ normal( 0 , 1 );
    to_vector( Zb1_hu ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = (abar[S[i]] + a[C[i], S[i]]) + b1[S[i]] * Y[i] + b2[S[i]] * PC1[i];
        hu[i] =  a_hu[S[i]] + b1_hu[S[i]] * Y[i];
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
    for ( i in 1:N ) {
        mu[i] = (abar[S[i]] + a[C[i], S[i]]) + b1[S[i]] * Y[i] + b2[S[i]] * PC1[i];
        hu[i] =  a_hu[S[i]] + b1_hu[S[i]] * Y[i];
    }
    for ( i in 1:N ) log_lik[i] = hurdle_lognormal_logit_lpdf(M[i] | mu[i], sigma, hu[i]);
}
