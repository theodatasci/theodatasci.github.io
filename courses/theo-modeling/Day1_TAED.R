
library(deSolve)


#---Lotka-Volterra model
#-----------------------
lv_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dN = r0*N - a*N*P
    dP = -m*P + e*a*N*P
    return(list(c(dN,dP)))
  })
}

#--- run the numerical integration
dyn <- data.frame(ode(y = c(N=10,P=1), times = seq(0,42,by=0.1), 
          func = lv_model, parms = c(r0=1,a=0.1,m=0.6,e=0.6), method = "rk4"))


#--- plot the result
par(mar=c(7,7,1,1))
plot(dyn$time,dyn$N,type="l",bty="l",lwd=3,xlab="time",ylab="density",col="orange",las=1,cex.lab=2,ylim=c(0,max(dyn[,2:3])), cex.axis=2,mgp=c(5,2,0))
lines(dyn$time,dyn$P,col="red",lwd=5)
legend("topright",bty="n",lwd=5,lty=1,col=c("orange","red"),legend=c("N","P"),cex=2)




#--- Logistic growth
#-------------------
logistic = function(r0,N0,t,K){K/(1+(K/N0-1)*exp(-r0*t))}

#--- parameters
rs = 0.5; k=30; N0=2; 
ts = seq(1,20,by=0.1); N = logistic(r0=rs,N0=N0,t=ts,K=k)

#--- plot
par(mar=c(5,5,5,5))
plot(ts,N,type="l",bty="l",lwd=7,ylim=c(0,40),las=1,cex.lab=3,cex.axis=3,xlab="time",xaxt="n",yaxt="n") #
lines(ts,logistic(r0=rs,N0=40,t=ts,K=k),lwd=7,lty=2)
lines(ts,logistic(r0=rs,N0=7,t=ts,K=k),lwd=7,lty=3)




#--- Logistic with 2 different initial conditions
#------------------------------------------------
rs = 0.5; k=30; N0=2; 
ts = seq(1,20,by=0.1); N = logistic(r0=rs,N0=N0,t=ts,K=k)
par(mar=c(5,5,5,5))
plot(ts,N,type="l",bty="l",lwd=7,ylim=c(0,40),las=1,cex.lab=3,cex.axis=3,xlab="time",xaxt="n",yaxt="n") #
lines(ts,logistic(r0=rs,N0=40,t=ts,K=k),lwd=7,lty=2)
lines(ts,logistic(r0=rs,N0=7,t=ts,K=k),lwd=7,lty=3)




#--- Lotka-Volterra 2-species competition 
#-----------------------------------------
LV_2sp= function(time,y,parms){
  with(as.list(c(y,parms)), {
    dN1 = r1*N1*(1-(N1+a_12*N2)/k1)
    dN2 = r2*N2*(1-(N2+a_21*N1)/k2)
    list(c(dN1,dN2))
  })
}

#--- parameters
pars2=c(r1=0.4,r2=0.3,a_12=0.2,a_21=0.1,k1=30,k2=20); times = seq(0,25,by=0.1)

#--- 2 dynamics with different initial conditions
dyn <- data.frame(ode(y = c(N1=1,N2=8), times = times, func = LV_2sp, parms = pars2, method = "rk4"))

dyn2 <- data.frame(ode(y = c(N1=2,N2=4), times = times, func = LV_2sp, parms = pars2, method = "rk4"))

#--- plot
par(mar=c(5,5,5,0))
plot(dyn$time,dyn$N1,type="l",bty="l",lwd=7,xlab="time",ylab="density",col="black",las=1,cex.lab=3,cex.axis=3,ylim=c(0,max(cbind(dyn[,2:3],dyn2[,2:3]))),xaxt="n",yaxt="n")
lines(dyn$time,dyn$N2,col="red",lwd=7)
lines(dyn2$time,dyn2$N1,col="black",lwd=7,lty=3)
lines(dyn2$time,dyn2$N2,col="red",lwd=7,lty=3)
legend("topleft",bty="n",lwd=5,lty=1,col=c("black","red"),legend=c("N1","N2"),cex=2)


#----------------------
#--- Analyse equilibria 
#----------------------

library(rootSolve)

#--- 2-species copetition function
EvalEq = function(t,y,pars) {
  with(as.list(c(y,pars)), {
    dN1 = N1*(r1+(a_11*N1+a_12*N2))
    dN2 = N2*(r2+(a_22*N2+a_21*N1))
    list(c(dN1,dN2),0)
  })
}

#--- parameter set
p = c(r1=0.1,r2=0.2,a_11=-1,a_22=-1,a_12=-0.5,a_21=0.4)

#--- evaluate equilibrium for this set, given certain initial densities
Eq = stode(y = c(N1=2,N2=3),func= EvalEq,parms=p,positive = TRUE)[[1]]
Eq


#--- Function to evaluate the jacobian
jac = function(t,y,pars) {
  with(as.list(c(y,pars)), {
    dN1 = N1*(r1+(a_11*N1+a_12*N2))
    dN2 = N2*(r2+(a_22*N2+a_21*N1))
    list(dN1,dN2)
  })
}

#--- evaluate the jacobian for this equilibrium 
J = jacobian.full(y=Eq,func=jac,parms=p) 

#--- calculate eigenvalues
eigens = eigen(J)$values 

#--- get the real part of the dominant eigenvalue
lmax=max(Re(eigens))      
lmax


#--- Plot equilibria and lambdamax while varying r1
r1s = seq(0.1,0.8,by=0.01);
tab = matrix(NA,ncol=4,nrow=length(r1s),dimnames=list(NULL,c("r1","N1","N2","lmax")))
i=1
for(i in 1:length(r1s)){
  p["r1"] = r1s[i]
  Eq = stode(y = c(N1=2,N2=3),func= EvalEq,parms=p,positive = TRUE)[[1]]
  J = jacobian.full(y=Eq,func=jac,parms=p) 
  tab[i,"r1"] = r1s[i]
  tab[i,"N1"] = Eq[1]
  tab[i,"N2"] = Eq[2]
  tab[i,"lmax"] = max(Re(eigen(J)$values))
}


layout(matrix(1:2,ncol=2))
par(mar=c(7,7,0,7))
plot(tab[,"r1"],tab[,"N1"],type="l",bty="l",ylim=range(tab[,2:3]),xlab="r1",ylab="Equilibrium",lwd=7,col="darkorange",las=1,cex.lab=4,cex.axis=4,xaxt="n",yaxt="n")
lines(tab[,"r1"],tab[,"N2"],lwd=4,col="navy")
legend("topleft",col=c("darkorange","navy"),lty=1,lwd=7,legend=c("N1","N2"),bty="n",cex=4)

plot(tab[,"r1"],tab[,"lmax"],type="l",bty="l",ylim=range(c(0,tab[,"lmax"])),xlab="r1",ylab="Lambda_max",lwd=7,las=1,cex.lab=4,cex.axis=4,xaxt="n",yaxt="n")
abline(h=0,lty=3,lwd=4)


#----------------
#--- Euler numerical integration for logistic growth
#---------------------------------------------------
logistic_deriv = function(N,r0,K){c(N=r0*N*(1-N/K))}             #--- Derivative

euler = function(N,r0,K,dt){return(N+dt*logistic_deriv(N,r0,K))} #--- Euler

dyn_euler = function(N0,tmin,tmax,dt,r0,K){                      #--- Iteration
  n=matrix(c(t=tmin,N=N0),1,2)
  for(t in 1:tmax) n=rbind(n,c(t=n[t,1]+dt,euler(N=n[t,2],r0=r0,K=K,dt=dt)))
  return(n)
}

dt = 2; tmax = 80  # time step and time end point

#--- run numerical integration for dt =2
dyn = dyn_euler(N0=1,tmin=0, tmax=tmax/dt,dt=dt,r0=0.1,K=10) 

#--- analytical solution for the logistic
logistic_analytic = function(ti,N0,r0,K)return(K/(1+(K/N0-1)*exp(-r0*ti)))

#--- dynamic with dt=5
dyn2 = dyn_euler(N0=1,tmin=0, tmax=tmax/5,dt=5,r0=0.1,K=10) 

#--- dynamic with dt=0.1
dyn3 = dyn_euler(N0=1,tmin=0, tmax=tmax/0.1,dt=0.1,r0=0.1,K=10) 

#--- plot results
par(mar=c(7,7,1,2))
plot(dyn[,"t"],dyn[,"N"],type="l",bty="l",xlab="time",ylab="density",las=1,col="orange",lwd=6,cex.lab=3)
tis=0:tmax
lines(dyn2[,"t"],dyn2[,"N"],col="red",lwd=6)
lines(dyn3[,"t"],dyn3[,"N"],col="purple",lwd=6)
lines(tis,logistic_analytic(tis,N0=1,r0=0.1,K=10),lty=2,lwd=6)
legend("topleft",legend=c("dt=0.1","dt=2","dt=5","solution"),col=c("purple","orange","red","black"),title="Euler",bty="n",lty=c(1,1,1,2),cex=2,lwd=3)


#--- Numerical integration of logistic with deSolve
#--------------------------------------------------
library(deSolve)
logistic_deriv2 = function(t,y,pars) {
  with(as.list(c(y,pars)), {
    dN = r0*N*(1-N/K)
    list(c(dN))
  })}
par(mar=c(7,7,1,2))
dyn=ode(y=c(N=1),times=seq(0,100,0.1),func=logistic_deriv2,parms=c(r0=0.1,K=10),method="euler")
plot(dyn[,"time"],dyn[,"N"],type="l",bty="l",xlab="time",ylab="density",las=1,lwd=6,cex.lab=3)


#--- Numerical integration of L-V model with deSolve
#---------------------------------------------------
library(deSolve)
lotka_Volterra_deriv = function(t,y,pars) {
  with(as.list(c(y,pars)), {
    dN = r0*N-a*N*P
    dP = e*a*N*P -m*P
    list(c(dN,dP))
  })}
par(mar=c(7,7,1,2))
dynLV=ode(y=c(N=2,P=1),times=seq(0,180,0.1),func=lotka_Volterra_deriv,
          parms=c(r0=0.2,a=0.3,e=0.5,m=0.1),method="rk4")
plot(dynLV[,"time"],dynLV[,"N"],type="l",bty="l",xlab="time",ylab="density",las=1,lwd=6,cex.lab=3,ylim=range(dynLV[,2:3]))
lines(dynLV[,"time"],dynLV[,"P"],lwd=6,col="red")
legend("topright",col = c("black","red"),legend=c("N","P"),cex=2,bty="n",lty=1,lwd=3)


#--- Rozensweig-McArthur model
#-----------------------------

LVK_deriv = function(t,y,pars) {
  with(as.list(c(y,pars)), {
    dN = r0*N*(1-N/K)-a*N*P/(1+h*N)
    dP = a*N*P/(1+h*N) -m*P
    list(c(dN,dP))
  })}

#--- with K=1
dynLVk1=ode(y=c(N=3,P=1),times=seq(0,180,0.1),func=LVK_deriv,parms=c(r0=0.2,K=1,a=0.3,e=0.5,m=0.1,h=0.5),method="rk4")

#--- with K=5
dynLVk5=ode(y=c(N=3,P=1),times=seq(0,180,0.1),func=LVK_deriv,parms=c(r0=0.2,K=5,a=0.3,e=0.5,m=0.1,h=0.5),method="rk4")

#--- plot
layout(matrix(1:2,ncol=2)); par(mar=c(7,7,2,2))
ylims = range(cbind(dynLVk1[,2:3],dynLVk5[,2:3]))
plot(dynLVk1[,"P"],dynLVk1[,"N"],type="l",bty="l",xlab="P",ylab="N",las=1,lwd=6,cex.lab=2,ylim=ylims)
plot(dynLVk5[,"P"],dynLVk5[,"N"],type="l",bty="l",xlab="P",ylab="N",las=1,lwd=6,cex.lab=2,ylim=ylims)


####----------------------------
#--- FITTING
#-----------

# REF: https://haddonm.github.io/URMQMF/model-parameter-estimation.html


#--- Least Sum Squared Residuals

library(MQMF)
#--- Models
vB = function(p, ages) return(p[1]*(1-exp(-p[2]*(ages-p[3]))))
Gz = function(p, ages) return(p[1]*exp(-p[2]*exp(p[3]*ages))) 
MM = function(p, ages) return((p[1]*ages)/(p[2] + ages^p[3]))

#--- Optimiaztion criterium SSQ
ssq <- function(funk,observed,...) {  
  predval <- funk(...) 
  return(sum((observed - predval)^2,na.rm=TRUE))  
}
#--- guess starting values
pars <- c("Linf"=27.0,"K"=0.15,"t0"=-2.0)

#--- run the optimization
bestvB <- nlm(f=ssq,funk=vB,observed=LatA$length,p=pars,ages=LatA$age,typsize=magnitude(pars))  

#--- compute the prediction
ages <- 1:max(LatA$age) # used in comparisons 
predvB <- vB(bestvB$estimate,ages)

#--- show optimization result
outfit(bestvB,backtran=FALSE,title="vB"); cat("\n")   


#--- fit Gomperz
pars <- c(26.0,0.7,-0.5) # Gompertz  
bestGz <- nlm(f=ssq,funk=Gz,observed=LatA$length,p=pars,  
              ages=LatA$age,typsize=magnitude(pars))  

#--- fit Michaelis-Menton 
pars <- c(23.2,1.0,1.0) # Michaelis-Menton - first start point  
bestMM1 <- nlm(f=ssq,funk=mm,observed=LatA$length,p=pars,  
               ages=LatA$age,typsize=magnitude(pars))  

#--- Predictions
predGz <- Gz(bestGz$estimate,ages) # using the outputs  
predmm <- mm(bestMM1$estimate,ages)

ymax <- getmax(LatA$length) #try ?getmax or getmax [no brackets]  
xmax <- getmax(LatA$age)  #there is also a getmin, not used here 

#--- plot
par(mar=c(7,7,1,1))
plot(LatA$age,LatA$length,type="p",pch=16, col=rgb(1,0,0,1/5),  
     cex=1.2,xlim=c(0,xmax),ylim=c(0,ymax),yaxs="i",xlab="Age",  
     ylab="Length (cm)",panel.first=grid(),las=1,cex.lab=2,cex.axis=2,mgp=c(4,1,0))  
lines(ages,predvB,lwd=5,col=4)       
lines(ages,predGz,lwd=5,col=1,lty=2) 
lines(ages,predmm,lwd=5,col=3,lty=3)  
legend("bottomright",cex=1.2,c("von Bertalanffy","Gompertz",  
                               "Michaelis-Menton"),col=c(4,1,3),lty=c(1,2,3),lwd=4,bty="n")  

#--- Compare normal and Log-normal residual distribution
#-------------------------------------------------------

pars <- c(27.25,0.15,-3.0)  

#--- fit with normal
bestvBN <- nlm(f=ssq,funk=vB,observed=LatA$length,p=pars,  
               ages=LatA$age,typsize=magnitude(pars),iterlim=1000)  

#--- function to calculate log residuals
ssqL <- function(funk,observed,...) {  
  predval <- funk(...)  
  return(sum((log(observed) - log(predval))^2,na.rm=TRUE))  
}

#--- fit with lognormal
bestvBLN <- nlm(f=ssqL,funk=vB,observed=LatA$length,p=pars,  
                ages=LatA$age,typsize=magnitude(pars),iterlim=1000)  

#--- predictions
predvBN <- vB(bestvBN$estimate,ages)   
predvBLN <- vB(bestvBLN$estimate,ages)   
ymax <- getmax(LatA$length)   
xmax <- getmax(LatA$age)  

#--- plot
parset()    
par(mar=c(7,7,1,1))
plot(LatA$age,LatA$length,type="p",pch=16, col=rgb(1,0,0,1/5),  
     cex=1.2,xlim=c(0,xmax),ylim=c(0,ymax),yaxs="i",xlab="Age",  
     ylab="Length (cm)",panel.first=grid(),cex.lab=3,cex.axis=3,mgp=c(4,1,0))  
lines(ages,predvBN,lwd=5,col=4,lty=2)   # add Normal dashed  
lines(ages,predvBLN,lwd=5,col=1)        # add Log-Normal solid  
legend("bottomright",c("Normal Errors","Log-Normal Errors"),  
       col=c(4,1),lty=c(2,1),lwd=4,bty="n",cex=1.2)  

#--- show fit summary results
outfit(bestvBN,backtran=FALSE,title="Normal errors"); cat("\n") 
outfit(bestvBLN,backtran=FALSE,title="Log-Normal errors") 


#----------------------------------------
#--- Fitting LSSQ with minpack.lm package
#----------------------------------------

#--- load data
dat = read.table("FITTING/ts_data.txt",h=T,sep="\t")
com = unique(dat$Community)
reps= unique(dat$Replicate)
xtime = unique(unlist(dat$"time_diff"))
sub = subset(dat,Community=="Tet")


#--- plot data for Tetrahymena
par(mar=c(7,8,2,2))
plot(dat$time,dat$Tet,type="n",bty="l",las=1,xlab="hours",ylab="density (x1000 ind.ml)", cex.axis=2,cex.lab=2, mgp=c(3,1,0),yaxt="n", main = "Tetrahymena")
axis(2,at=seq(0,10000,by=2000),labels =seq(0,10000,by=2000)/1000,las=1,cex.axis=2)
for(i in 1:length(reps)){
  sub = subset(dat,rep==reps[i])
  lines(sub$time,sub$Tet ,type="b",lty=2,lwd=3,cex=2)
}

#--- load package
library(minpack.lm)

#--- logistic growth
LogGrowth = function(time,y,parms){
  with(as.list(c(y,parms)), {
    dN = (r0 - alpha*N) * N
    list(c(dN))
  })
}

#--- function to calculate residuals
resids_LogG <- function(parms){	
  cinit=c(N=sub[1,"Tet"])
  output = ode(y=cinit, times=as.numeric(sub$time_diff), func=LogGrowth, parms=parms)
  return(output[,"N"]-c(sub[,"Tet"]))
}

#--- fit
fitval=nls.lm(par=c(r0=0.08,alpha=0.00001),fn=resids_LogG, 
              control=list(maxiter=1000), lower=c(r0=0, alpha=0))
fitval$par

#--- Predictions
yinit = c(N=mean(subset(sub,time_diff==0)[,"Tet"]))
fitted_output <- ode(y=yinit, times=c(0:max(xtime)), LogGrowth, fitval$par)

#--- Plot

par(mar=c(7,8,2,2))
plot(sub$time_diff,sub$Tet,type="n",bty="l",las=1,xlab="hours",ylab="density (x1000 ind.ml)", cex.axis=1.5,cex.lab=1.5, mgp=c(3,1,0),yaxt="n", main = "Tetrahymena")
axis(2,at=seq(0,10000,by=2000),labels =seq(0,10000,by=2000)/1000,las=1,cex.axis=1.5)
for(i in 1:length(reps)){
  subsub = subset(sub,Replicate==reps[i])
  lines(subsub$time_diff,subsub$Tet ,type="b",lty=2,lwd=3,cex=2)
}
lines(fitted_output[,2]~c(fitted_output[,1]),lwd=5)


#--- calculate confidence interval
p = 2 # number of parameters to fit
sqrt(diag(p*sum(resids_LogG(fitval$par)^2)/(length(xtime)-p)*solve(fitval$hessian*p)))* 1.96






#-----------------------------------
#--- MAXIMUM LIKELIHOOD
#-----------------------------------

#--- Fit with normally distributed residuals
#-------------------------------------------
pars <- c(Linf=27.0,K=0.15,t0=-3.0,sigma=2.5) # starting values  
ansvB <- nlm(f=negNLL,p=pars,funk=vB,observed=LatA$length,ages=ages,  
             typsize=magnitude(pars))  


#--- Fit with log-normal distribution of residuals
#--------------------------------------------------
data(tigers)

#--- Beveton-Holt model
lbh <- function(p,biom) return(log((p[1]*biom)/(p[2] + biom)))  
pars <- c("a"=25,"b"=4.5,"sigma"=0.4)   # includes a sigma  

#--- fit
best <- nlm(negNLL,pars,funk=lbh,observed=log(tigers$Recruit),biom=tigers$Spawn,typsize=magnitude(pars))  
outfit(best,backtran=FALSE,title="Beverton-Holt Recruitment")  


#--- prediction
predR <- exp(lbh(best$estimate,tigers$Spawn))   
#note exp(lbh(...)) is the median because no bias adjustment  
result <- cbind(tigers,predR,tigers$Recruit/predR)  


#--- plot
par(mar=c(7,7,1,1))
plot(tigers$Spawn,predR,xlab="Spawning Biomass",ylab="Recruitment",ylim=c(0,38),type="l",lwd=4,las=1,cex.lab=2,cex.axis=2,mgp=c(4,1,0))  
points(tigers$Spawn,tigers$Recruit,pch=16,cex=2,col=2) 

