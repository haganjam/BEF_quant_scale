//
// Gamma hurdle model
//
data{
     int<lower=1> N;
     int<lower=1> S_N;
     int<lower=1> C_N;
     vector[N] M;
     vector[N] Y;
    array[N] int C;
    array[N] int S;
}
parameters{
     // standard deviation: log-normal model
     real<lower=0> phi;
     // standard normal deviations: log-normal model
     vector[S_N] Za;
     // parameters: log-normal model
     real<lower=0> sigma_a;
     real abar;
     real b1;
     // parameters: binomial
     real a_hu;
     real b1_hu;
}
transformed parameters{
     // transformed parameters: log-normal model
     vector[S_N] a;
     a = abar + (Za*sigma_a);
}
model{
    // vector of means: log-normal linear model 
     vector[N] mu;
    // vector of means: binomial linear model 
     vector[N] hu;
    // standard deviation of the log-normal distribution
    phi ~ exponential( 2 );
    // linear model priors: log-normal
    abar ~ normal(0, 3);
    sigma_a ~ exponential( 3 );
    b1 ~ normal(0, 2);
    // linear model priors: binomial model
    a_hu ~ normal(0, 2);
    b1_hu ~ normal(0, 2);
    // standard normal vectors
    to_vector( Za ) ~ normal( 0 , 1 );
    // likelihood: Gamma hurdle
    for (i in 1:N) {
        mu[i] = a[S[i]] + b1 * Y[i];
        hu[i] = a_hu + b1_hu * Y[i];
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
    for ( i in 1:N ) {
        mu[i] = a[S[i]] + b1 * Y[i];
        hu[i] = a_hu + b1_hu * Y[i];
      if (M[i] == 0) {
        log_lik[i] = bernoulli_lpmf(1 | inv_logit(hu[i]) );
      } else {
        log_lik[i] = bernoulli_lpmf(0 | inv_logit(hu[i]) ) +
                    gamma_lpdf(M[i] | exp(mu[i])*exp(mu[i])/phi, exp(mu[i])/phi);
       }
    }
}
