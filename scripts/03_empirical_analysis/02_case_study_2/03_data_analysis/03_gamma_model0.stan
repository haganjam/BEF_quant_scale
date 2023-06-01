//
// LogNormal-Hurdle: Model 1
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
     // gamma model parameters
     // scale parameter
     real<lower=0> phi;
     // a parameters
     matrix[S_N, C_N] Za;
     cholesky_factor_corr[S_N] L_Rho_z;
     vector<lower=0>[S_N] sigma_z;
     row_vector[S_N] zbar;
     // b parameters
     vector[S_N] Zb1;
     real<lower=0> sigma_b1;
     real b1bar;
     // parameters: binomial model
     matrix[2, S_N] V;
     cholesky_factor_corr[2] L_Rho_v;
     vector<lower=0>[2] sigma_v;
     vector[2] vbar;
}
transformed parameters{
     // transformed parameters
     // a parameters
     matrix[C_N, S_N] z;
     matrix[C_N, S_N] zbar_mat;
     matrix[C_N, S_N] a;
     z = (diag_pre_multiply(sigma_z, L_Rho_z) * Za)';
     zbar_mat = rep_matrix(zbar, C_N);
     a = z + zbar_mat;
     // b parameters
     vector[S_N] b1;
     b1 = b1bar + (Zb1*sigma_b1);
     // transformed parameters: binomial model
     matrix[S_N, 2] v;
     vector[S_N] a_hu;
     vector[S_N] b1_hu;
     v = (diag_pre_multiply(sigma_v, L_Rho_v) * V)';
     a_hu = vbar[1] + v[,1];
     b1_hu = vbar[2] + v[,2];
}
model{
    // vector of means: gamma model
     vector[N] mu;
    // vector of means: binomial linear model 
     vector[N] hu;
    // standard deviation of the log-normal distribution
    phi ~ exponential( 4 );
    // linear model priors: gamma model
    // a parameters
    zbar ~ normal(0, 3);
    sigma_z ~ exponential( 2 );
    L_Rho_z ~ lkj_corr_cholesky( 1 );
    // b parameters
    b1bar ~ normal(0, 3);
    sigma_b1 ~ exponential(2);
    // linear model priors: binomial model
    vbar ~ normal(0, 2);
    sigma_v ~ exponential( 2 );
    L_Rho_v ~ lkj_corr_cholesky( 2 );
    // standard normal vectors
    to_vector( Za ) ~ normal( 0 , 1 );
    to_vector( Zb1 ) ~ normal( 0 , 1 );
    to_vector( V ) ~ normal(0, 1);
    for ( i in 1:N ) {
        mu[i] = a[C[i],S[i]] + b1[S[i]] * Y[i];
        hu[i] =  a_hu[S[i]] + b1_hu[S[i]] * Y[i];
      if (M[i] == 0) {
        target += bernoulli_lpmf(1 | inv_logit(hu[i]) );
      } else {
        target += bernoulli_lpmf(0 | inv_logit(hu[i]) ) +
                  gamma_lpdf(M[i] | exp(mu[i])*exp(mu[i])/phi, exp(mu[i])/phi);
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
     matrix[S_N,S_N] Rho_a;
    Rho_a = multiply_lower_tri_self_transpose(L_Rho_z);
     matrix[2,2] Rho_hu;
    Rho_hu = multiply_lower_tri_self_transpose(L_Rho_v);
    for ( i in 1:N ) {
        mu[i] = a[C[i],S[i]] + b1[S[i]] * Y[i];
        hu[i] =  a_hu[S[i]] + b1_hu[S[i]] * Y[i];
      if (M[i] == 0) {
        log_lik[i] = bernoulli_lpmf(1 | inv_logit(hu[i]) );
      } else {
        log_lik[i] = bernoulli_lpmf(0 | inv_logit(hu[i]) ) +
                    gamma_lpdf(M[i] | exp(mu[i])*exp(mu[i])/phi, exp(mu[i])/phi);
       }
    }
}
