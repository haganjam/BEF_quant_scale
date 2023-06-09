//
// LogNormal-Hurdle: Null model
//
data{
     int<lower=1> N;
     vector[N] M;
}
parameters{
     // standard deviation: log-normal model
     real<lower=0> sigma;
     // parameters: lognormal
     real a;
     // parameters: binomial
     real a_hu;
}
model{
    // vector of means: log-normal linear model 
     vector[N] mu;
    // vector of means: binomial linear model 
     vector[N] hu;
    // standard deviation of the log-normal distribution
    sigma ~ exponential( 5 );
    // linear model priors: log-normal
    a ~ normal(0, 1.5);
    // linear model priors: binomial model
    a_hu ~ normal(0, 1.5);
    for ( i in 1:N ) {
        mu[i] = a;
        hu[i] = a_hu;
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
    for ( i in 1:N ) {
        mu[i] = a;
        hu[i] = a_hu;
      if (M[i] == 0) {
        log_lik[i] = bernoulli_lpmf(1 | inv_logit(hu[i]) );
      } else {
        log_lik[i] = bernoulli_lpmf(0 | inv_logit(hu[i]) ) +
                  lognormal_lpdf(M[i] | mu[i], sigma);
       }
    }
}
