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
Y = rbinom(TT,Z,rep(p,TT))/p

plot(1:TT, Z,
     pch = 19, col="red", ty = "o",
     xlab = "time", ylab = "abundance",
     ylim=c(0,1.2*K), las=1)
points(1:TT, Y,
       pch=15, col="blue", ty="o", lty = 3)
legend("bottomright",
       legend = c("Observation", "Process"),
       pch = c(15, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n")

# Fitting obs error only -------------------------------------------------------

data = list(Y=Y, TT=TT)

inits = list(r=0.5,
             alpha=1e-3,
             s_obs=0.1,
             U1=data$Y[1])

# fit_obs = stan(file = "model_obs.stan",
#                data = data,
#                chains = 3,
#                cores = 3,
#                init = rep(list(inits),3),
#                iter = 4000)
# save(fit_obs, file="fit_obs.RData")
load("fit_obs.RData")

print(fit_obs)

# plotting with coda package
samples = As.mcmc.list(fit_obs)
plot(samples[, 1:4])

## predictions -----------------------------------------------------------------

samples = as.matrix(fit_obs) # extract samples as matrix
# object for predictions: 1 row per sample, 1 column per observation 
U.pred = matrix(NA, nrow=1000, ncol=TT) 

for(j in 1:1000){ # loop over samples
  # extract parameters from posterior
  U.pred[j,1] = samples[j,"U1"] 
  r.pred = samples[j,"r"]
  a.pred = samples[j,"alpha"]
  for(i in 1:(TT-1)){ # loop through time
    # compute prediction
    U.pred[j,i+1] = U.pred[j,i]*exp(r.pred*(1-a.pred*U.pred[j,i]))      
  }
}

# extact quantiles
U.pred.q = apply(U.pred, 2, function(x) quantile(x, probs=c(0.025, 0.50, 0.975)))

plot(1:TT, Y, type="n",
     xlab = "time", ylab = "abundance",
     ylim=c(0,1.2*K), las=1)

# 95% quantiles:
polygon(c(1:TT,TT:1),c(U.pred.q[1,],rev(U.pred.q[3,])), 
        col=adjustcolor("red", alpha.f=0.25) , border=NA)
# median regression line
lines(1:TT, U.pred.q[2,], col="red", lwd=2)
# data
points(1:TT, Y, pch=15, col="blue", ty="o", lty = 3)

legend("bottomright",
       legend = c("Observation", "Prediction"),
       col = c("blue", "red"),
       lty = c(1, 1),
       horiz=TRUE, bty="n")


# Fitting proc error only ------------------------------------------------------

data = list(Y=Y, TT=TT)

inits = list(r=0.5,
             alpha=1e-3,
             s_proc=0.1)

# fit_proc = stan(file = "model_proc.stan",
#                 data = data,
#                 chains = 3,
#                 cores = 3,
#                 init = rep(list(inits),3),
#                 iter = 4000)
# save(fit_proc, file="fit_proc.RData")
load("fit_proc.RData")

print(fit_proc)

# plotting with coda package
samples = As.mcmc.list(fit_proc)
plot(samples[, 1:3])

## predictions -----------------------------------------------------------------

samples = as.matrix(fit_proc) # extract samples as matrix

# reshape data into Y_t and Y_t+1 of same length
Y.this = Y[1:(TT-1)]
Y.next = Y[2:TT]

# objects for predictions: 1 row per sample
U.this = seq(min(Y.this), max(Y.this), length.out=100)
U.next = matrix(NA, nrow=1000, ncol=100)

for(j in 1:1000){ # loop over samples
  # extract parameters from posterior
  r.pred = samples[j,"r"]
  a.pred = samples[j,"alpha"]
  # compute all one-step ahead predictions
  U.next[j, ] = U.this*exp(r.pred*(1-a.pred*U.this))
}

# extract quantiles
U.next.q = apply(U.next, 2, function(x) quantile(x, probs=c(0.025, 0.50, 0.975)))

plot(Y.this, Y.next, type="n", xlab="Y_t", ylab="Y_t+1")

# 95% quantiles:
polygon(c(U.this,rev(U.this)), c(U.next.q[1,],rev(U.next.q[3,])), 
        col=adjustcolor("red", alpha.f=0.25) , border=NA)
# median regression line:
lines(U.this, U.next.q[2,], col="red", lwd=2)
# data:
points(Y.this, Y.next, pch = 15, col="blue")

legend("bottomright",
       legend = c("Observation", "Prediction"),
       col = c("blue", "red"),
       lty = c(1, 1),
       horiz=TRUE, bty="n")

abline(0,1, col="darkgray", lty=2)

## alternative predictions -----------------------------------------------------
## log growth is linear model

samples = as.matrix(fit_proc) # extract samples as matrix

# reshape data into Y_t and Y_t+1 of same length
Y.this = Y[1:(TT-1)]
Y.next = Y[2:TT]

# objects for predictions: 1 row per sample
U.this = seq(min(Y), max(Y), length.out=100)
U.grow = matrix(NA, nrow=1000, ncol=100)

for(j in 1:1000){ # loop over samples
  # extract parameters from posterior
  r.pred = samples[j,"r"]
  a.pred = samples[j,"alpha"]
  # compute all one-step ahead predictions
  U.grow[j, ] = r.pred*(1-a.pred*U.this)
}

# extract quantiles
U.grow.q = apply(U.grow, 2, function(x) quantile(x, probs=c(0.025, 0.50, 0.975)))

plot(Y.this, log(Y.next/Y.this), type="n", xlab="Y_t", ylab="log(Y_t+1)-log(Y_t)")

# 95% quantiles:
polygon(c(U.this,rev(U.this)), c(U.grow.q[1,],rev(U.grow.q[3,])), 
        col=adjustcolor("red", alpha.f=0.25) , border=NA)
# median regression line:
lines(U.this, U.grow.q[2,], col="red", lwd=2)
# data:
points(Y.this, log(Y.next/Y.this), pch = 15, col="blue")

legend("bottomleft",
       legend = c("Observation", "Prediction"),
       col = c("blue", "red"),
       lty = c(1, 1),
       horiz=TRUE, bty="n")

abline(0,0, col="darkgray", lty=2)


# Fitting a state-space model --------------------------------------------------

data = list(Y=Y, TT=TT)

inits = list(r=0.5,
             alpha=1e-3,
             s_proc=0.1,
             s_obs=0.1,
             Z=data$Y)

# fit_ssm =  stan(file = "model_ssm.stan",
#                 data = data,
#                 chains = 3,
#                 cores = 3,
#                 init = rep(list(inits),3),
#                 iter = 4000,
#                 control = list(adapt_delta=0.99))
# save(fit_ssm, file="fit_ssm.RData")
load("fit_ssm.RData")

print(fit_ssm)

# plotting with coda package
samples = As.mcmc.list(fit_ssm)
plot(samples[, 1:4])

## predictions -----------------------------------------------------------------

# extract some quantiles for times series of latent states
Z_med = summary(fit_ssm)$summary[5:18, "50%"]
Z_lbd = summary(fit_ssm)$summary[5:18, "2.5%"]
Z_ubd = summary(fit_ssm)$summary[5:18, "97.5%"]

plot(1:TT, Z, type="n",
     xlab = "time", ylab = "abundance",
     ylim=c(0,1.2*K), las=1)
polygon(c(1:TT,TT:1),c(Z_lbd,rev(Z_ubd)), col=adjustcolor("red", alpha.f=0.25) , border=NA)
lines(1:TT, Z_med, col="red", lwd=2)
points(1:TT, Y,
       pch=15, col="blue", ty="o", lty = 3)
legend("bottomright",
       legend = c("Observation", "Latent states"),
       col = c("blue", "red"),
       lty = c(1, 1),
       horiz=TRUE, bty="n")

