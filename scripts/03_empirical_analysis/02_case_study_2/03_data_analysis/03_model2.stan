//
// Model 2
//

data{
     vector[287] M;
     vector[287] PC2;
     vector[287] PC1;
     vector[287] Y;
     vector[287] T;
    array[287] int C;
    array[287] int S;
}
parameters{
     matrix[5,10] Z;
     cholesky_factor_corr[5] L_Rho_b4;
     cholesky_factor_corr[5] L_Rho_b3;
     cholesky_factor_corr[5] L_Rho_b2;
     cholesky_factor_corr[5] L_Rho_b1;
     cholesky_factor_corr[5] L_Rho_a;
     vector<lower=0>[5] sigma_b4;
     vector<lower=0>[5] sigma_b3;
     vector<lower=0>[5] sigma_b2;
     vector<lower=0>[5] sigma_b1;
     vector<lower=0>[5] sigma_a;
     vector[5] b4bar;
     vector[5] b3bar;
     vector[5] b2bar;
     vector[5] b1bar;
     vector[5] abar;
     real<lower=0> sigma;
}
transformed parameters{
     matrix[10,5] a;
     matrix[10,5] b1;
     matrix[10,5] b2;
     matrix[10,5] b3;
     matrix[10,5] b4;
    b4 = (diag_pre_multiply(sigma_b4, L_Rho_b4) * Z)';
    b3 = (diag_pre_multiply(sigma_b3, L_Rho_b3) * Z)';
    b2 = (diag_pre_multiply(sigma_b2, L_Rho_b2) * Z)';
    b1 = (diag_pre_multiply(sigma_b1, L_Rho_b1) * Z)';
    a = (diag_pre_multiply(sigma_a, L_Rho_a) * Z)';
}
model{
     vector[287] u;
    sigma ~ exponential( 1 );
    abar ~ normal( 0 , 2 );
    b1bar ~ normal( 0 , 2 );
    b2bar ~ normal( 0 , 2 );
    b3bar ~ normal( 0 , 2 );
    b4bar ~ normal( 0 , 2 );
    sigma_a ~ exponential( 1 );
    sigma_b1 ~ exponential( 1 );
    sigma_b2 ~ exponential( 1 );
    sigma_b3 ~ exponential( 1 );
    sigma_b4 ~ exponential( 1 );
    L_Rho_a ~ lkj_corr_cholesky( 2 );
    L_Rho_b1 ~ lkj_corr_cholesky( 2 );
    L_Rho_b2 ~ lkj_corr_cholesky( 2 );
    L_Rho_b3 ~ lkj_corr_cholesky( 2 );
    L_Rho_b4 ~ lkj_corr_cholesky( 2 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:287 ) {
        u[i] = (abar[S[i]] + a[C[i], S[i]]) + (b1bar[S[i]] + b1[C[i], S[i]] * T[i]) + (b2bar[S[i]] + b2[C[i], S[i]] * Y[i]) + (b3bar[S[i]] + b3[C[i], S[i]] * PC1[i]) + (b4bar[S[i]] + b4[C[i], S[i]] * PC2[i]);
    }
    M ~ lognormal( u , sigma );
}
generated quantities{
    vector[287] log_lik;
     vector[287] u;
     matrix[5,5] Rho_a;
     matrix[5,5] Rho_b1;
     matrix[5,5] Rho_b2;
     matrix[5,5] Rho_b3;
     matrix[5,5] Rho_b4;
    Rho_b4 = multiply_lower_tri_self_transpose(L_Rho_b4);
    Rho_b3 = multiply_lower_tri_self_transpose(L_Rho_b3);
    Rho_b2 = multiply_lower_tri_self_transpose(L_Rho_b2);
    Rho_b1 = multiply_lower_tri_self_transpose(L_Rho_b1);
    Rho_a = multiply_lower_tri_self_transpose(L_Rho_a);
    for ( i in 1:287 ) {
        u[i] = (abar[S[i]] + a[C[i], S[i]]) + (b1bar[S[i]] + b1[C[i], S[i]] * T[i]) + (b2bar[S[i]] + b2[C[i], S[i]] * Y[i]) + (b3bar[S[i]] + b3[C[i], S[i]] * PC1[i]) + (b4bar[S[i]] + b4[C[i], S[i]] * PC2[i]);
    }
    for ( i in 1:287 ) log_lik[i] = lognormal_lpdf( M[i] | u[i] , sigma );
}
