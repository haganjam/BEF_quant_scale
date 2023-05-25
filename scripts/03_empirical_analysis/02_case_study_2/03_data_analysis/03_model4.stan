//
// Model 4
//

data{
     vector[287] M;
     vector[287] T;
    array[287] int C;
    array[287] int S;
}
parameters{
     matrix[5,10] Z;
     cholesky_factor_corr[5] L_Rho_b1;
     cholesky_factor_corr[5] L_Rho_a;
     vector<lower=0>[5] sigma_b1;
     vector<lower=0>[5] sigma_a;
     vector[5] b1bar;
     vector[5] abar;
     real<lower=0> sigma;
}
transformed parameters{
     matrix[10,5] a;
     matrix[10,5] b1;
    b1 = (diag_pre_multiply(sigma_b1, L_Rho_b1) * Z)';
    a = (diag_pre_multiply(sigma_a, L_Rho_a) * Z)';
}
model{
     vector[287] u;
    sigma ~ exponential( 1 );
    abar ~ normal( 0 , 2 );
    b1bar ~ normal( 0 , 2 );
    sigma_a ~ exponential( 1 );
    sigma_b1 ~ exponential( 1 );
    L_Rho_a ~ lkj_corr_cholesky( 2 );
    L_Rho_b1 ~ lkj_corr_cholesky( 2 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:287 ) {
        u[i] = (abar[S[i]] + a[C[i], S[i]]) + (b1bar[S[i]] + b1[C[i], S[i]] * T[i]);
    }
    M ~ lognormal( u , sigma );
}
generated quantities{
    vector[287] log_lik;
     vector[287] u;
     matrix[5,5] Rho_a;
     matrix[5,5] Rho_b1;
    Rho_b1 = multiply_lower_tri_self_transpose(L_Rho_b1);
    Rho_a = multiply_lower_tri_self_transpose(L_Rho_a);
    for ( i in 1:287 ) {
        u[i] = (abar[S[i]] + a[C[i], S[i]]) + (b1bar[S[i]] + b1[C[i], S[i]] * T[i]);
    }
    for ( i in 1:287 ) log_lik[i] = lognormal_lpdf( M[i] | u[i] , sigma );
}
