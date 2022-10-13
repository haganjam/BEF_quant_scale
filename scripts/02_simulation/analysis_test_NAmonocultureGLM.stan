
data {
    int M[150];
    vector[150] E;
    vector[150] Y;
    int S[150];
}
parameters {
    vector[5] aS;
    vector[5] b_yS;
    vector[5] b_eS;
    vector[5] b_yeS;
}
model {
    vector[150] lambda;
    b_yeS ~ normal( 0 , 1 );
    b_eS ~ normal( 0 , 1 );
    b_yS ~ normal( 0 , 1 );
    aS ~ normal( 2 , 2 );
    for ( i in 1:150 ) {
        lambda[i] = aS[S[i]] + b_yS[S[i]] * Y[i] + b_eS[S[i]] * E[i] + b_yeS[S[i]] * E[i] * Y[i];
        lambda[i] = exp(lambda[i]);
    }
    M ~ poisson( lambda );
}
