//
// LogNormal-Hurdle: Model 1
//
data{
     int<lower=1> N;
     int<lower=1> S_N;
     vector[N] M;
     vector[N] Y;
     vector[N] PC1;
    array[N] int S;
}
parameters{
     // standard deviation: log-normal model
     real<lower=0> sigma;
     // standard normal deviations: log-normal model
     matrix[3, S_N] V;
     matrix[2, S_N] Z;
     // parameters: log-normal model
     cholesky_factor_corr[3] L_Rho_v;
     vector<lower=0>[3] sigma_v;
     vector[3] vbar;
     // parameters: binomial model
     cholesky_factor_corr[2] L_Rho_z;
     vector<lower=0>[2] sigma_z;
     vector[2] zbar;
}
transformed parameters{
     // transformed parameters: log-normal model
     matrix[S_N, 3] v;
     vector[S_N] a;
     vector[S_N] b1;
     vector[S_N] b2;
     v = (diag_pre_multiply(sigma_v, L_Rho_v) * V)';
     a = vbar[1] + v[,1];
     b1 = vbar[2] + v[,2];
     b2 = vbar[3] + v[,3];
     // transformed parameters: binomial model
     matrix[S_N, 2] z;
     vector[S_N] a_hu;
     vector[S_N] b1_hu;
     z = (diag_pre_multiply(sigma_z, L_Rho_z) * Z)';
     a_hu = zbar[1] + z[,1];
     b1_hu = zbar[2] + z[,2];
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
    zbar ~ normal(0, 1.5);
    sigma_z ~ exponential( 3 );
    L_Rho_z ~ lkj_corr_cholesky( 2 );
    // standard normal vectors
    to_vector( V ) ~ normal( 0 , 1 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = a[S[i]] + b1[S[i]] * Y[i] + b2[S[i]] * PC1[i];
        hu[i] =  a_hu[S[i]] + b1_hu[S[i]] * Y[i];
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
     matrix[3,3] Rho_v;
     Rho_v = multiply_lower_tri_self_transpose(L_Rho_v);
     // correlation matrices: binomial model
     matrix[2,2] Rho_z;
     Rho_z = multiply_lower_tri_self_transpose(L_Rho_z);
    for ( i in 1:N ) {
        mu[i] = a[S[i]] + b1[S[i]] * Y[i] + b2[S[i]] * PC1[i];
        hu[i] =  a_hu[S[i]] + b1_hu[S[i]] * Y[i];
      if (M[i] == 0) {
        log_lik[i] = bernoulli_lpmf(1 | inv_logit(hu[i]) );
      } else {
        log_lik[i] = bernoulli_lpmf(0 | inv_logit(hu[i]) ) +
                     lognormal_lpdf(M[i] | mu[i], sigma);
       }
    }
}
