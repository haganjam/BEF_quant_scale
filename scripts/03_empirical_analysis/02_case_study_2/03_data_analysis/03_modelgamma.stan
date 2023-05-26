//
// Model 5
//

data{
     vector[287] M;
    array[287] int C;
    array[287] int S;
}
parameters{
     matrix[5,10] Z;
     cholesky_factor_corr[5] L_Rho_a;
     vector<lower=0>[5] sigma_a;
     vector[5] abar;
     real<lower=0.001,upper=100> shape;
}
transformed parameters{
     matrix[10,5] a;
    a = (diag_pre_multiply(sigma_a, L_Rho_a) * Z)';
}
model{
    abar ~ normal( 0 , 2 );
    sigma_a ~ exponential( 1 );
    L_Rho_a ~ lkj_corr_cholesky( 2 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:287 ) {
        M[i] ~ gamma(shape, ( shape/exp((abar[S[i]] + a[C[i], S[i]])) ));
    }
}
generated quantities{
     vector[287] log_lik;
     matrix[5,5] Rho_a;
     Rho_a = multiply_lower_tri_self_transpose(L_Rho_a);

     for ( i in 1:287 ) { 
       log_lik[i] = gamma_lpdf( M[i] | shape, ( shape/exp((abar[S[i]] + a[C[i], S[i]])) ) );
     }

}

