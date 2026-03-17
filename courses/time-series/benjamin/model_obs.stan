data {
  int TT; 
  vector[TT] Y; 
}

parameters {
  real<lower=0> r; 
  real<lower=0> alpha; 
  real<lower=0> s_obs;
  real<lower=0> U1; 
}

model {
  vector[TT] U; // prediction variable

  // weak priors
  r ~ exponential(1);
  alpha ~ exponential(1);
  s_obs ~ normal(0, 1);
  U1 ~ normal(100, 10);

  // predictions
  U[1] = U1;
  for(t in 2:TT){
    U[t] = U[t-1]*exp(r*(1-alpha*U[t-1]));
  }

  // observation error
  for(t in 1:TT){
    Y[t] ~ lognormal(log(U[t]), s_obs);
  }
}
