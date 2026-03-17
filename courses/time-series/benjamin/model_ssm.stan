data {
  int TT;
  vector[TT] Y; 
}

parameters {
  real<lower=0> r; 
  real<lower=0> alpha; 
  real<lower=0> s_proc;
  real<lower=0> s_obs;
  vector[TT] Z; 
}

model {
  vector[TT] U; // prediction variable

  // weak priors
  r ~ exponential(1);
  alpha ~ exponential(1);
  s_proc ~ normal(0, 1);
  s_obs ~ normal(0, 1);
  Z[1] ~ normal(100, 10);

  // predictions and process error
  for(t in 2:TT){
    U[t] = Z[t-1]*exp(r*(1-alpha*Z[t-1]));
    Z[t] ~ lognormal(log(U[t]), s_proc);
  }

  // observation error
  for(t in 1:TT){
    Y[t] ~ lognormal(log(Z[t]), s_obs);
  }
}
