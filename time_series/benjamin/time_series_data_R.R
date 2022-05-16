rm(list=ls())
library("rstan")
library("coda")

set.seed(100) # random seed

# modeling ---------------------------------------------------------------------

TT = 14 # total time: 14 days 
Z = numeric(TT) # process time series

s_proc = 0.05 # process error standard deviation for environmental noise

r = 0.8 # growth rate 
K = 1000 # carrying capacity
alpha = 1/K # competition

# process equation: Ricker model
Z[1]=100
for(i in 1:(TT-1)){
  Z[i+1] = Z[i]*exp(r*(1-alpha*Z[i]))      # deterministic prediction
  Z[i+1] = Z[i+1]*exp(rnorm(1, 0, s_proc)) # environmental noise
  Z[i+1] = rpois(1, Z[i+1])                # demographic noise
}

# observation: count abundances in fraction p of total volume
p = 0.1
Y = rbinom(TT,round(Z),rep(p,TT))/p

plot(1:TT, Z,
     pch = 19, col="red", ty = "o",
     xlab = "time", ylab = "abundance",
     ylim=c(0,1.2*K), las=1)
points(1:TT, Y,
       pch=15, col="blue", ty="o", lty = 3)
legend("bottom",
       legend = c("Observation", "Process"),
       pch = c(15, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n")

# Fitting obs error only -------------------------------------------------------

code_obs = "
data {
  int TT; 
  int Y[TT]; 
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
"

data = list(Y=Y, TT=TT)

inits = list(r=0.5,
             alpha=1e-3,
             s_obs=0.1,
             U1=data$Y[1])

fit_obs = stan(model_code = code_obs,
               data = data,
               chains = 3, 
               cores = 3,
               init = rep(list(inits),3),
               iter = 4000)

print(fit_obs)

# plotting with coda package
samples = As.mcmc.list(fit_obs)
plot(samples[, 1:4])

# Fitting proc error only ------------------------------------------------------

code_proc = "
data {
  int TT;  
  int Y[TT];  
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
"

data = list(Y=Y, TT=TT)

inits = list(r=0.5,
             alpha=1e-3,
             s_proc=0.1)

fit_proc = stan(model_code = code_proc,
                data = data,
                chains = 3, 
                cores = 3,
                init = rep(list(inits),3),
                iter = 4000)

print(fit_proc)

# plotting with coda package
samples = As.mcmc.list(fit_proc)
plot(samples[, 1:3])

# Fitting a state-space model --------------------------------------------------

code_ssm = "
data {
  int TT;
  int Y[TT]; 
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
"

data = list(Y=Y, TT=TT)

inits = list(r=0.5,
             alpha=1e-3,
             s_proc=0.1,
             s_obs=0.1,
             Z=data$Y)

fit_ssm =  stan(model_code = code_ssm,
                data = data,
                chains = 3, 
                cores = 3,
                init = rep(list(inits),3),
                iter = 4000,
                control = list(adapt_delta=0.99))

print(fit_ssm)

# plotting with coda package
samples = As.mcmc.list(fit_ssm)
plot(samples[, 1:4])

# extract some quantiles for times series of latent states
Z_med = summary(fit_ssm)$summary[5:18, "50%"]
Z_lbd = summary(fit_ssm)$summary[5:18, "2.5%"]
Z_ubd = summary(fit_ssm)$summary[5:18, "97.5%"]

plot(1:TT, Z, type="n",
     xlab = "time", ylab = "abundance",
     ylim=c(0,1.2*K), las=1)
polygon(c(1:TT,TT:1),c(Z_lbd,rev(Z_ubd)), col="lightgrey" , border=NA)
lines(1:TT, Z, pch = 19, col="red", ty = "o")
lines(1:TT, Z_med)
points(1:TT, Y,
       pch=15, col="blue", ty="o", lty = 3)
legend("bottom",
       legend = c("Observation", "Process","Latent states"),
       pch = c(15, 19, NA),
       col = c("blue", "red", "black"),
       lty = c(3, 1, 1),
       horiz=TRUE, bty="n")
