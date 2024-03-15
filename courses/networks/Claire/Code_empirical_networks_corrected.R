rm(list=ls())
library(devtools)
install.packages("MASS")
install.packages("plotrix")
library(MASS)
library(plotrix)

#################
# 14.03.2024
# Complexity-stability relationship in empirical food webs
###########################

###########################
#May's result in random communities
###########################
m<-matrix(rnorm(10^6),nrow=10^3)
diag(m)=-1
plot(eigen(m)$values,xlab="real part",ylab="imaginary part")
abline(h=0)
abline(v=0)

# Compute descriptors of complexity
Diversity = nrow(m) # Number of species in network i
Connectance = mean(m!=0) # Connectance
Connectance 
NNE_sd = sd(m) # Standard deviation of non-zero elements
# May's stability criterion (or complexity)
Complexity = NNE_sd*(sqrt(Diversity*Connectance))
Complexity
# Adding the stability criterion (the circle): use the function draw.circle() from plotrix package
draw.circle(mean(diag(m)),0,Complexity,lty=2)
max(as.numeric(eigen(m)$values)) # the maximum eigenvalue, the one driving stability

# Let's change the non-diagonal elements (the variance of interaction strengths):
m2<-matrix(rnorm(10^6,mean = 0, sd=0.1),nrow=10^3)
NNE_sd = sd(m2)
Diversity = nrow(m2) # Number of species in network i
Connectance = mean(m2!=0) # Connectance
Complexity = NNE_sd*(sqrt(Diversity*Connectance))
Complexity
plot(eigen(m2)$values,xlab="real part",ylab="imaginary part")
abline(h=0)
abline(v=0)
# Adding the stability criterion (the circle):
draw.circle(mean(diag(m2)),0,Complexity,lty=2)
max(as.numeric(eigen(m2)$values)) #the maximum eigenvalue, the one driving stability

# What if the matrix is smaller? E.g., 50 species instead of 1000? 
m<-matrix(rnorm(50^2),nrow=50)
plot(eigen(m)$values,xlab="real part",ylab="imaginary part")
abline(h=0)
abline(v=0)

# Descriptors of food-web complexity
Diversity = nrow(m) # Number of species in network i
Connectance =  mean(m!=0) # Connectance
Connectance 
NNE_sd = sd(m) # Standard deviation of non-zero elements
# May's stability criterion (or complexity)
Complexity = NNE_sd*(sqrt(Diversity*Connectance))
Complexity
# Adding the stability criterion (the circle):
draw.circle(mean(diag(m)),0,Complexity,lty=2)
max(as.numeric(eigen(m)$values))

###########################
# Data and package loading:
################
rm(list=ls())
# 1. Data frame with names of Ecopath models and habitat types (116 models):
load("Ecopath_models.Rdata")
ndat = nrow(Ecopath_models)
head(Ecopath_models)
# 2. List of species names for each model:
load("GroupName.Rdata")
# 3. List of species biomass for each model (tons/km^2):
load("B_vec.Rdata")
# 4. List of species production/biomass ratios (P/B) for each model (/year):
load("P_vec.Rdata")
# 5. List of species  consumption/biomass ratios (Q/B) for each model (/year):
load("Q_vec.Rdata")
# 6. List of DIET matrices (quantitative diet composition) for each model:
load("DIET.Rdata")
# Let's have a look at one food-web:
i = 61
Ecopath_models[i,]
GroupName[[i]]
Q_vec[[i]]
P_vec[[i]]
B_vec[[i]]
DIET[[i]]
apply(DIET[[i]],2,sum)
############################################################
# 1. Computing food-web complexity and stability
############################################################
A_mat = vector("list",length=ndat) # Empty list to store all interaction matrices
C_mat = vector("list",length=ndat) # Empty list to store all community matrices
Stab = numeric(length=ndat) # Local stability

# Descriptors of food-web complexity
Diversity = numeric(length=ndat) # Species richness
Connectance = numeric(length=ndat) # Connectance
NNE_sd = numeric(length=ndat) # Standard deviation of interaction strengths
Complexity = numeric(length=ndat) # May's stability criterion

for (i in 1 : ndat){
  
  D = DIET[[i]] # Diet matrix of network i
  B = as.numeric(B_vec[[i]]) # Species biomass
  P = as.numeric(P_vec[[i]]) # P/B
  Q = as.numeric(Q_vec[[i]]) # Q/B
  S = length(B) # Number of species in network i
  
  # We compute the outflows matrix (the biomass of each species that is consumed by another one)
  OUTFLOWS = -t(D%*%diag(Q)%*%diag(B)) # Negative effect of the consumer population on their resources
  
  # We compute the inflows matrix (the biomass consumed that is converted into predator biomass)
  e = P/Q # Resource conversion efficiency
  e[e=="Inf"] = 0 # if species i is not a consumer, then Q[i] = 0 and P/Q = Inf. We replace "Inf" by 0.
  INFLOWS = -t(OUTFLOWS*e) # Positive signs (i.e., add biomass to consumers)
  
  # We compute the per capita interaction matrix A, that is outflows and inflows matrices divided by predator and prey biomasses
  B2 = diag(B)%*%matrix(1,S,S)%*%diag(B) # B_rc x B_cr
  A = t((OUTFLOWS+INFLOWS)/B2)
  A_mat[[i]] = A # The per capita interaction matrix
  
  # We compute the community matrix (or Jacobian matrix):
  J = diag(B)%*%A
  C_mat[[i]] = J
  
  # We compute local stability
  Stab[i] = max(as.numeric(eigen(J)$values))

  # Species richness
  Diversity[i] = length(B) # Number of species in network i
  
  # Connectance
  Connectance[i] = mean(J!=0)
  
  # Standard deviation of non-zero off-diagonal elements
  vecJ = as.numeric(J)
  NNE_sd[i] = sd(vecJ[-c(which(vecJ==0))])
  
  # May's stability criterion (or complexity)
  Complexity[i] = NNE_sd[i]*(sqrt(Diversity[i]*Connectance[i]))
}

################################################
# 2. Analysing the relationship between complexity and stability in empirical food webs
###########################################################
# Plot the relationship between complexity and stability
plot(Complexity,Stab,las=1,
  xlab = expression(Complexity~sigma*sqrt(S*C)),
  ylab = expression(Instability~italic(Re(lambda*scriptstyle(max)))),
  pch=21)

# Test the significance of the relationship between complexity and stability
model = lm(Stab~Complexity)
summary(model) 

# Plotting the eigenvalues of Ecopath_models[61,]: Mid Atlantic Bight:
i=61
plot(eigen(C_mat[[i]])$values,xlab="real part",ylab="imaginary part")
abline(h=0)
abline(v=0)
Stab[i]
Complexity[i]
# Adding the stability criterion (the circle):
plot(eigen(C_mat[[i]])$values,xlab="real part",ylab="imaginary part",xlim=c(-25,25))
abline(h=0)
abline(v=0)
draw.circle(0,0,Complexity[i],lty=2)

## What if we fully randomize the food web? We pick up elements from a normal distribution with the same mean and sd
m<-matrix(rnorm((C_mat[[i]])^2,mean=mean(C_mat[[i]]),sd=sd(C_mat[[i]])),nrow=nrow(C_mat[[i]]))
plot(eigen(m)$values,xlab="real part",ylab="imaginary part",xlim=c(-25,25))
abline(h=0)
abline(v=0)
draw.circle(0,0,Complexity[i],lty=2)
# Complexity predicts pretty well the maximum eigenvalue (not perfect, but rather good)

#####################################################
# 3. Comparing the complexity-stability relationship of empirical and randomized food webs
#####################################################
# Let's focus on a small food-web first to look at useful function:
i = 6
J = C_mat[[i]]
S = nrow(J)
help(upper.tri) # useful function to select elements of a matrix
tJ = t(J) # useful too: tJ is the transpose matrix of J (rows become columns and columns become rows)
pairs = cbind(J[upper.tri(J)],tJ[upper.tri(tJ)]) # pairs of interactions (c_ij and c_ji)

help(sample) # useful function to do permutations
sample(pairs[,1])
sample(pairs[,2])

help(mvrnorm) # useful function to generate bivariate normal distribution (for example: with mean of positive elements, mean of negative elements and covariance matrix)
BD=mvrnorm(n = 100, c(mean(pairs[,1]),mean(pairs[,2])), Sigma=cov(pairs)) #Bivariate normal distribution with mean of positive elements, mean of negative elements and covariance matrix

######################################################################
# An example of randomization test:
# Changing the topological structure of the community matrix (H1)
Stab_H1 = numeric(length=ndat) # Store the stability of randomized matrices
rep =  numeric(length=100) # replications of the permutation

for (i in 1 : ndat){ # Loop over all food webs
  J = C_mat[[i]]
  diag(J)=0
  S = nrow(J)
  tJ = t(J)
  pairs = cbind(J[upper.tri(J)],tJ[upper.tri(tJ)]) # pairs of interactions
  
  for (oo in 1 :length(rep)){ # Replication loop
  vec = pairs[sample(1:nrow(pairs)),] # randomization of the position of pairs of interactions in the matrix   
  H1 = matrix(0,S,S)
  tH1 = t(H1)
  H1[upper.tri(H1)] = vec[,1]
  tH1[upper.tri((tH1))] = vec[,2]
  H1 = H1+t(tH1) # Community matrix with random topological structure
  rep[oo] =  max(as.numeric(eigen(H1)$values)) 
  }
  
  Stab_H1[i] = mean(rep) # Stability averaged over the replications
}
  
# Plot the relationship between complexity and stability
plot(Complexity,Stab_H1,las=1,
     xlab = expression(Complexity~sigma*sqrt(S*C)),
     ylab = expression(Instability~italic(Re(lambda*scriptstyle(max)))),
     pch=21,col=2)
points(Complexity,Stab,col=1) # Stability of non-randomized food webs

# Test the significance of the relationship between complexity and stability
model = lm(Stab_H1~Complexity)
summary(model) 

# Let's look at the distribution of eigenvalues:
plot(density(log10(Stab+1)),main="", xlab=expression(Stability~log[10](italic(Re(lambda*scriptstyle(max))))),cex.lab=1.4,cex.axis=1.4,col=1)
lines(density(log10(Stab_H1+1)),lwd=1.5,lty=1,col=2)

# Plot eigenvalues on the complex plane for one replicate:
i=61
plot(eigen(H1)$values,xlab="real part",ylab="imaginary part")
abline(h=0)
abline(v=0)
max(as.numeric(eigen(H1)$values))
Complexity[i]
# Adding the stability criterion (the circle):
plot(eigen(H1)$values,xlab="real part",ylab="imaginary part",xlim=c(-25,25))
abline(h=0)
abline(v=0)
draw.circle(0,0,Complexity[i],lty=2)

########################################################################

########################################################################
# Other randomization tests:
# Correlation between pairs of interactions c_ij and c_ji (H2)
Stab_H2 = numeric(length=ndat) # Store the stability of randomized matrices
rep =  numeric(length=100) # replications of the permutation

for (i in 1 : ndat){
  J = C_mat[[i]]
  S = nrow(J)
  tJ = t(J)
  pairs = cbind(J[upper.tri(J)],tJ[upper.tri(tJ)]) # pairs of interactions
  
  for (oo in 1 :length(rep)){ 
    pos = J[which(J>0)] # Isolate positive elements
    neg = J[which(J<0)] # Isolate negative elements
    n = sample(neg) # Randomization of positive elements
    p = sample(pos) # Randomization of negative elements
    H2 = sign(J)
    H2[which(H2<0)] = n
    H2[which(H2>0)] = p # Community matrix with no correlation between +/- elements
    rep[oo] =  max(as.numeric(eigen(H2)$values)) 
  }
  Stab_H2[i] = mean(rep)
}

# Plot the relationship between complexity and stability
plot(Complexity,Stab_H2,las=1,
     xlab = expression(Complexity~sigma*sqrt(S*C)),
     ylab = expression(Instability~italic(Re(lambda*scriptstyle(max)))),
     pch=21,col=3)
# points(Complexity,Stab_H1,col=2)
points(Complexity,Stab,col=1)

# Test the significance of the relationship between complexity and stability
model = lm(Stab_H2~Complexity)
summary(model) 

# Let's look a the distribution of eigenvalues:
plot(density(log10(Stab+1)),xlim=range(log10(Stab_H2+1),log10(Stab+1)),main="", xlab=expression(Stability~log[10](italic(Re(lambda*scriptstyle(max))))),cex.lab=1.4,cex.axis=1.4)
lines(density(log10(Stab_H1+1)),lwd=1.5,lty=1,col=2)
lines(density(log10(Stab_H2+1)),lwd=1.5,lty=1,col=3)

# Plot eigenvalues on the complex plane for one replicate:
i=61
J = C_mat[[i]]
S = nrow(J)
tJ = t(J)
pairs = cbind(J[upper.tri(J)],tJ[upper.tri(tJ)]) # pairs of interactions
pos = J[which(J>0)] # Isolate positive elements
neg = J[which(J<0)] # Isolate negative elements
n = sample(neg) # Randomization of positive elements
p = sample(pos) # Randomization of negative elements
H2 = sign(J)
H2[which(H2<0)] = n
H2[which(H2>0)] = p # Community matrix with no correlation between +/- elements

plot(eigen(H2)$values,xlab="real part",ylab="imaginary part")
abline(h=0)
abline(v=0)
max(as.numeric(eigen(H2)$values))
Complexity[i]
# Adding the stability criterion (the circle):
plot(eigen(H2)$values,xlab="real part",ylab="imaginary part",xlim=c(-25,25))
abline(h=0)
abline(v=0)
draw.circle(0,0,Complexity[i],lty=2)

#########################################################################
# Distribution of the elements of the community matrix   
Stab_H3 = numeric(length=ndat) # Store the stability of randomized matrices
rep =  numeric(length=100) # replications of the permutation

for (i in 1 : ndat){
  J = C_mat[[i]]
  S = nrow(J)
  tJ = t(J)
  pairs = cbind(J[upper.tri(J)],tJ[upper.tri(tJ)]) # pairs of interactions
  
  for (oo in 1 :length(rep)){ 
    cp = pairs[which(pairs[,1]!=0),] # Isolate non-zero elements, we want to keep the null elements at the same place 
    cp2=cp[which(cp[,1]>0),] # Isolate pairs where the positive element is in column 1
    cp3 = cbind(cp[which(cp[,1]<0),2],cp[which(cp[,1]<0),1]) 
    CP=rbind(cp2,cp3)
    if (nrow(cp2)==0)  CP=cp
    LN=mvrnorm(n = 5000, c(mean(CP[,1]),mean(CP[,2])), Sigma=cov(CP)) #Bivariate normal distribution with mean of positive elements, mean of negative elements and covariance matrix
    LN=LN[which(LN[,1]>0),]
    LN=LN[which(LN[,2]<0),]
    
    tt=pairs
    tt[which(pairs[,1]>0),1] = LN[c(1:length(which(pairs[,1]>0))),1]
    tt[which(pairs[,1]>0),2] = LN[c(1:length(which(pairs[,1]>0))),2]
    
    tt[which(pairs[,1]<0),1] = LN[c(1:length(which(pairs[,1]<0))),2]
    tt[which(pairs[,1]<0),2] = LN[c(1:length(which(pairs[,1]<0))),1]
    
    H3=matrix(0,S,S) 
    tH3=t(H3)
    H3[upper.tri(H3)]=tt[,1]
    tH3[upper.tri((tH3))]=tt[,2]
    H3=H3+t(tH3) # Community matrix with +/- elements picked from a bivariate normal distribution
    rep[oo] =  max(as.numeric(eigen(H3)$values)) 
  }
  Stab_H3[i] = mean(rep)
}

# Plot the relationship between complexity and stability
plot(Complexity,Stab_H3,las=1,
     xlab = expression(Complexity~sigma*sqrt(S*C)),
     ylab = expression(Instability~italic(Re(lambda*scriptstyle(max)))),
     pch=21,col=4)
points(Complexity,Stab,col=1)
#points(Complexity,Stab_H2,col=3)
#points(Complexity,Stab_H1,col=2)

# Test the significance of the relationship between complexity and stability
model = lm(Stab_H3~Complexity)
summary(model) 

# Let's look a the distribution of eigenvalues:
plot(density(log10(Stab+1)),xlim=range(log10(Stab_H3+1),log10(Stab+1)),main="", xlab=expression(Stability~log[10](italic(Re(lambda*scriptstyle(max))))),cex.lab=1.4,cex.axis=1.4)
lines(density(log10(Stab_H1+1)),lwd=1.5,lty=1,col=2)
lines(density(log10(Stab_H2+1)),lwd=1.5,lty=1,col=3)
lines(density(log10(Stab_H3+1)),lwd=1.5,lty=1,col=4)

# Plot eigenvalues on the complex plane for one replicate:
i=61
J = C_mat[[i]]
S = nrow(J)
tJ = t(J)
pairs = cbind(J[upper.tri(J)],tJ[upper.tri(tJ)]) # pairs of interactions

cp = pairs[which(pairs[,1]!=0),] # Isolate non-zero elements, we want to keep the null elements at the same place 
cp2=cp[which(cp[,1]>0),] # Isolate pairs where the positive element is in column 1
cp3 = cbind(cp[which(cp[,1]<0),2],cp[which(cp[,1]<0),1]) 
CP=rbind(cp2,cp3)
if (nrow(cp2)==0)  CP=cp
LN=mvrnorm(n = 5000, c(mean(CP[,1]),mean(CP[,2])), Sigma=cov(CP)) #Bivariate normal distribution with mean of positive elements, mean of negative elements and covariance matrix
LN=LN[which(LN[,1]>0),]
LN=LN[which(LN[,2]<0),]

tt=pairs
tt[which(pairs[,1]>0),1] = LN[c(1:length(which(pairs[,1]>0))),1]
tt[which(pairs[,1]>0),2] = LN[c(1:length(which(pairs[,1]>0))),2]

tt[which(pairs[,1]<0),1] = LN[c(1:length(which(pairs[,1]<0))),2]
tt[which(pairs[,1]<0),2] = LN[c(1:length(which(pairs[,1]<0))),1]

H3=matrix(0,S,S) 
tH3=t(H3)
H3[upper.tri(H3)]=tt[,1]
tH3[upper.tri((tH3))]=tt[,2]
H3=H3+t(tH3) # Community matrix with +/- elements picked from a bivariate normal distribution

plot(eigen(H3)$values,xlab="real part",ylab="imaginary part")
abline(h=0)
abline(v=0)
max(as.numeric(eigen(H3)$values))
Complexity[i]
# Adding the stability criterion (the circle):
plot(eigen(H3)$values,xlab="real part",ylab="imaginary part",xlim=c(-25,25))
abline(h=0)
abline(v=0)
draw.circle(0,0,Complexity[i],lty=2)