

#Parameters
S=400  #number of species
meanA= .02  #average interaction
sigA=0.05 #std of interactions
sigK =.1 #std of K
Nmin = 0.00001 #Extinction threshold

#Carrying capacities
K=rnorm(S,1,sigK )
#Interactions
a=matrix(rnorm( S*S , 0,1   ),S,S)   #random with mean 0, var 1
A=meanA + sigA*a
diag(A)<-1

#Lotka volterra
eqs <- function(t,x,params){
  #x[x<Nmin]=0 #Avoid negative abundances
  dx = x *( K  - A%*%x ) 
  return(list(dx))
}

library(deSolve)

N0 = runif(S) #Random initial abundances
times=1:100 #Times at which the simulation is run

Nf=ode(N0, times=times, eqs ) #Run simulation

Neq=Nf[100,-1] #Abundances at final time step

#Abundance distribution at final time step
hist(Neq,freq=FALSE,breaks=20)

#Time series of abundances
plot(times,Nf[,2],type="l")
for (i in 3:S+1){
  lines(times,Nf[,i])
}

#Eigenvalues of matrix of survivors
plot(-eigen(A[Neq>0.001,Neq>0.001])$values)
