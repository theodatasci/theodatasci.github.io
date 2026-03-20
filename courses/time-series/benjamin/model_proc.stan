data {
  int TT;  
  vector[TT] Y;  
}

parameters {
  real<lower=0> r; 
  real<lower=0> alpha; 
  real<lower=0> s_proc;
}

model {
  vector[TT] U; // prediction variable

  // weak priors
  r ~ exponential(1);
  alpha ~ exponential(1);
  s_proc ~ normal(0, 1);

  // predictions
  for(t in 2:TT){
    U[t] = Y[t-1]*exp(r*(1-alpha*Y[t-1]));
  }

  // process error
  for(t in 2:TT){
    Y[t] ~ lognormal(log(U[t]), s_proc);
  }
}
