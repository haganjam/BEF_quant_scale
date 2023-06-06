//
// LogNormal-Hurdle: Model 3
//
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
    sigma ~ exponential( 5 );
    // linear model priors: log-normal
    vbar ~ normal(0, 1);
    sigma_v ~ exponential( 3 );
    L_Rho_v ~ lkj_corr_cholesky( 2 );
    // linear model priors: binomial
    a_hu ~ normal(0, 1.5);
    b1_hu ~ normal(0, 1.5);
    // standard normal vectors
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = a[S[i]] + b1[S[i]] * Y[i];
        hu[i] =  a_hu + b1_hu * Y[i];
      if (M[i] == 0) {
        target += bernoulli_lpmf(1 | inv_logit(hu[i]) );
      } else {
        target += bernoulli_lpmf(0 | inv_logit(hu[i]) ) +
                  lognormal_lpdf(M[i] | mu[i], sigma);
       }
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
      if (M[i] == 0) {
        log_lik[i] = bernoulli_lpmf(1 | inv_logit(hu[i]) );
      } else {
        log_lik[i] = bernoulli_lpmf(0 | inv_logit(hu[i]) ) +
                  lognormal_lpdf(M[i] | mu[i], sigma);
       }
    }
}
