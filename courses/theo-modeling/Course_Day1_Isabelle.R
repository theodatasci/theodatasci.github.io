
# clear workspace
rm(list=ls())

# load packages
library(deSolve)
library(phaseR)
library(rootSolve)
library(nleqslv)


############################################################################___
####.................................................
####----- WHAT TYPES OF MODELS ----------------------
####.................................................

###----- Deterministic versus Stochastic models -----
{
  #---- Model functions (Verhulst model, 2 versions) ----
  {
    #... Analytical solution of Verhulst's model
    verhulst_sol = function(N0,r0,K,t){
      N = K / (1 + (K / N0 - 1) * exp(-r0 * t))
      return(N)
    }
    
    #... Differential equation of Verhulst's model
    verhulst_ode = function(N,r0,K){
      dN = r0 * N * (1 - N / K)
      return(dN)
    }
    
    #... Runge-Kutta 4 algorithm for numerical integration of the Verhulst's model
    #.............................................................................
    #... n1: density at the previous time point; r0: growth rate; K: carrying capacity; h: time step
    rk4_verhulst = function(n1,r0,K,h){
      a1=verhulst_ode(N=n1,r0=r0,K=K)
      b1=verhulst_ode(N=n1+h*a1[1]/2,r0=r0,K=K)
      c1=verhulst_ode(N=n1+h*b1[1]/2,r0=r0,K=K)
      d1=verhulst_ode(N=n1+h*c1[1],r0=r0,K=K)
      return(c(n1+h*(a1[1]+2*b1[1]+2*c1[1]+d1[1])/6))
    }
    
    
    #... Function to calculate the dynamics for the deterministic version of the model (with the analytical solution)
    dyn_verhulst = function(n0,r0,K,tmax=250,h=0.1){
      tseq = seq(0,tmax,by = h)                               #... define the time sequence
      return(cbind(t=tseq,N=verhulst_sol(n0,r0,K,tseq)))
    }
    
    #... Function to calculate the dynamics for the stochastic version of the model (with the ODE)
    dyn_verhulst_stoch = function(n0,r0,K,tmax=250,h=0.1,a=0.25,sdstoch=0.5){
      n=c(n0)                                                     #... initialyse the vector of abundances
      tseq = seq(h,tmax,by = h)                                   #... define the time sequence   
      for(t in 1:length(tseq)){
        newn0 = rk4_verhulst(n[t],r0,K,h)                          #... calculate new density at time t with  ode
        newn = newn0 + a * rnorm(1,sd=sdstoch) * sqrt(newn0)       #... add stochastic demography term
        if(newn < 0)newn = 0                                       #... set to 0 if the new density is negative
        n=rbind(n,newn)                                            #... record the density
      }
      return(cbind("t"=c(0,tseq),"N"=n))
    }
    
  }
  
  
  #### Large populations ----
  {
    set.seed(1)          #... set a seed to reproduce the same figures
    
    #... Parameters
    Tmax = 500           #... maximum time
    h=0.1                #... time step (∆t)
    r0 = 0.1            #... population growth rate
    K = 2000             #... carrying capacity
    nrep_stoch = 10      #... number of stochastic replicates
    ast = 0.4            #... intensity of the stochastic noise
    sds = 1              #... standard deviation of the normal noise distribution          
    N0=30                #... initial density
    
    #... Stochastic dynamics (10 replicates)
    dynstoch = dyn_verhulst_stoch(N0,r0,K,Tmax,h,ast,sds)
    for(i in 2:nrep_stoch){
      dynnew = dyn_verhulst_stoch(N0,r0,K,Tmax,h,ast,sds)
      dynstoch = cbind(dynstoch,dynnew[,2])
    }
    meanstock = apply(dynstoch[,2:(nrep_stoch+1)],1,mean)    #... Mean of stochastic dynamics
    
    #... Deterministic dynamics
    dyndet = dyn_verhulst(n0=N0,r0,K,Tmax,h)
    
    
    #... Margins specifications
    par(mar = c(5,5,2,1))
    
    ylims = range(dynstoch[,2:(nrep_stoch+1)])  #... define Y-scale
    
    #... Plot one stochastic dynamics
    pdf(here::here("img","1_dynstoch_largepop1.pdf"), width = 5,height = 5)
    par(mar = c(5,5,2,1))
    plot(dynstoch[,1],dynstoch[,2],ylim=ylims,type="l",xlab="Time",ylab="Population density",lwd=1,col="chartreuse3",cex.lab=1.2,cex.axis=1.2,las=1,mgp=c(3.5,1,0))
    dev.off()
    
    #... Plot all stochastic dynamics
    pdf(here::here("img","1_dynstoch_largepop2.pdf"), width = 5,height = 5)
    par(mar = c(5,5,2,1))
    plot(dynstoch[,1],dynstoch[,2],ylim=ylims,type="l",xlab="Time",ylab="Population density",lwd=1,col="chartreuse3",cex.lab=1.2,cex.axis=1.2,las=1,mgp=c(3.5,1,0))
    for(i in 3:(nrep_stoch+1))lines(dynstoch[,1],dynstoch[,i],col="chartreuse3",lwd=1)
    dev.off()
    
    #... Plot all stochastic dynamics with mean
    pdf(here::here("img","1_dynstoch_largepop3.pdf"), width = 5,height = 5)
    par(mar = c(5,5,2,1))
    plot(dynstoch[,1],dynstoch[,2],ylim=ylims,type="l",xlab="Time",ylab="Population density",lwd=1,col="chartreuse3",cex.lab=1.2,cex.axis=1.2,las=1,mgp=c(3.5,1,0))
    for(i in 3:(nrep_stoch+1))lines(dynstoch[,1],dynstoch[,i],col="chartreuse3",lwd=1)
    lines(dynstoch[,1],meanstock,lwd=2,col="darkgreen")
    abline(h=K,lty=3)
    dev.off()
    
    
    #... Plot all stochastic dynamics with mean and deterministic dynamics
    pdf(here::here("img","1_dynstoch_largepop4.pdf"), width = 5,height = 5)
    par(mar = c(5,5,2,1))
    plot(dynstoch[,1],dynstoch[,2],ylim=ylims,type="l",xlab="Time",ylab="Population density",lwd=1,col="chartreuse3",cex.lab=1.2,cex.axis=1.2,las=1,mgp=c(3.5,1,0))
    for(i in 3:(nrep_stoch+1))lines(dynstoch[,1],dynstoch[,i],col="chartreuse3",lwd=1)
    lines(dynstoch[,1],meanstock,lwd=2,col="darkgreen")
    lines(dyndet[,1],dyndet[,2],lwd=2)
    abline(h=K,lty=3)
    dev.off()

    
  }
  
  
  #### Small populations ----
  {
    set.seed(2)          #... set a seed to reproduce the same figures
    
    K = 200             #... carrying capacity
    r0 = 0.04
    
    #... Stochastic dynamics (10 replicates)
    dynstoch = dyn_verhulst_stoch(N0,r0,K,Tmax,h,ast,sds)
    for(i in 2:nrep_stoch){
      dynnew = dyn_verhulst_stoch(N0,r0,K,Tmax,h,ast,sds)
      dynstoch = cbind(dynstoch,dynnew[,2])
    }
    meanstock = apply(dynstoch[,2:(nrep_stoch+1)],1,mean)    #... Mean of stochastic dynamics
    
    #... Deterministic dynamics
    dyndet = dyn_verhulst(N0,r0,K,Tmax,h)
    
    
    #... define Y-scale
    ylims = range(dynstoch[,2:(nrep_stoch+1)])  
    
    #... Plot one stochastic dynamics
    pdf(here::here("img","2_dynstoch_smallpop1.pdf"), width = 5,height = 5)
    par(mar = c(5,5,2,1))
    plot(dynstoch[,1],dynstoch[,2],ylim=ylims,type="l",xlab="Time",ylab="Population density",lwd=1,col="chartreuse3",cex.lab=1.2,cex.axis=1.2,las=1,mgp=c(3.5,1,0))
    dev.off()
    
    #... Plot all stochastic dynamics
    pdf(here::here("img","2_dynstoch_smallpop2.pdf"), width = 5,height = 5)
    par(mar = c(5,5,2,1))
    plot(dynstoch[,1],dynstoch[,2],ylim=ylims,type="l",xlab="Time",ylab="Population density",lwd=1,col="chartreuse3",cex.lab=1.2,cex.axis=1.2,las=1,mgp=c(3.5,1,0))
    for(i in 3:(nrep_stoch+1))lines(dynstoch[,1],dynstoch[,i],col="chartreuse3",lwd=1)
    dev.off()
    
    #... Plot all stochastic dynamics with mean
    pdf(here::here("img","2_dynstoch_smallpop3.pdf"), width = 5,height = 5)
    par(mar = c(5,5,2,1))
    plot(dynstoch[,1],dynstoch[,2],ylim=ylims,type="l",xlab="Time",ylab="Population density",lwd=1,col="chartreuse3",cex.lab=1.2,cex.axis=1.2,las=1,mgp=c(3.5,1,0))
    for(i in 3:(nrep_stoch+1))lines(dynstoch[,1],dynstoch[,i],col="chartreuse3",lwd=1)
    lines(dynstoch[,1],meanstock,lwd=2,col="darkgreen")
    abline(h=K,lty=3)
    dev.off()
    
    #... Plot all stochastic dynamics with mean and deterministic dynamics
    pdf(here::here("img","2_dynstoch_smallpop4.pdf"), width = 5,height = 5)
    par(mar = c(5,5,2,1))
    plot(dynstoch[,1],dynstoch[,2],ylim=ylims,type="l",xlab="Time",ylab="Population density",lwd=1,col="chartreuse3",cex.lab=1.2,cex.axis=1.2,las=1,mgp=c(3.5,1,0))
    for(i in 3:(nrep_stoch+1))lines(dynstoch[,1],dynstoch[,i],col="chartreuse3",lwd=1)
    lines(dynstoch[,1],meanstock,lwd=2,col="darkgreen")
    lines(dyndet[,1],dyndet[,2],lwd=2)
    abline(h=K,lty=3)
    dev.off()

  }
}






############################################################################___
####.................................................
####----- HOW TO BUILD A MODEL ----------------------
####.................................................

############################################################################___
###----- NUMERICAL INTEGRATION, EXAMPLES OF INTEGRATION ERRORS ----------------


###----- Logistic Growth -----

{
  ###---- Function for the analytical solution of the Verhulst model
  verhulst_sol = function(N0,r0,K,t){
    N = K / (1 + (K / N0 - 1) * exp(-r0 * t))
    return(N)
  }
  
  
  ###---- Differential equation of Verhulst's model (derivative)
  verhulst_ode = function(N,r0,K){
    dN = r0 * N * (1 - N / K)
    return(dN)
  }
  
  
  #--- Function to apply Euler algorithm for the numerical integration of the Verhulst model (get a matrix of the dynamic)
  dyn_verhulst_euler = function(n0,r0,K,tmax=250,h=0.1){
    tseq = seq(0,tmax,by = h)                               #... define the time sequence
    density = numeric(length(tseq))
    density[1] = n0
    for(i in 2:length(tseq)){
      density[i] = density[i-1] +  verhulst_ode(density[i-1],r0,K) * h
    }
    return(cbind(t=tseq,N=density))
  }
  
  
  ###---- Parameters
  tmax = 200                                     #--- maximal time
  h=0.1                                          #--- timestep
  tseq = seq(0,tmax,by = h)                      #--- time sequence
  hseq = c(0.1,1,10)                             #--- sequence to vary the timestep
  col_hseq = c("chartreuse3", "orange", "red")   #--- colors for the timesteps
  
  K = 2000                                       #--- Carrying capacity
  r0 = 0.1                                       #--- intrinsic growth rate
  
  Ninit = 1                                      #--- Initial population density
  
  
  ###---- Plot the analytical solution
  pdf(here::here("img","3_euler_Logistic.pdf"), width = 5,height = 5.5)
  plot(tseq,verhulst_sol(Ninit,r0,K,tseq),type="l",xlab="Time", ylab= "Population density",lwd=2,las=1, main="Example 1, Logistic growth")
  
  ###---- Add the lines corresponding to numerical integration using each timestep value
  for(i in 1:length(hseq)){
    mat = dyn_verhulst_euler(Ninit,r0,K,tmax=tmax,h=hseq[i])
    lines(mat[,"t"],mat[,"N"],lwd=2,col=col_hseq[i],lty=2)
  }
  legend("topleft",lty=c(1,2,2,2),lwd=rep(2,4),col=c("black",col_hseq),legend=c("solution",hseq),title=expression(paste(delta, "t")),bty="n")
  dev.off()
  
}


###---- Predator-prey Lotka-Volterra model -----

{
  ###---- Differential equation of Verhulst's model (derivative)
  LotkaVolterra_ode = function(N,P,r0,a,e,m){
    dN = r0 * N - a * N * P
    dP = e * a * N * P - m * P
    return(c(dN,dP))
  }
  
  
  
  ###---- Function to apply Euler's algorithm for the numerical integration of the Lotka-Volterra model (get a matrix of the dynamic)
  dyn_lotkaVolterra_euler = function(n0,p0,r0,a,e,m,tmax=250,h=0.1){
    tseq = seq(0,tmax,by = h)                               #... define the time sequence
    density_N = density_P = numeric(length(tseq))
    density_N[1] = n0
    density_P[1] = p0
    for(i in 2:length(tseq)){
      deriv = LotkaVolterra_ode(density_N[i-1],density_P[i-1],r0,a,e,m)
      density_N[i] = density_N[i-1] +  deriv[1] * h
      density_P[i] = density_P[i-1] +  deriv[2] * h
    }
    return(cbind(t=tseq,N=density_N,P=density_P))
  }
  
  
  ###---- Parameters
  tmax = 200                                     #--- maximal time
  hseq = c(0.01,0.1,1)                             #--- sequence to vary the timestep
  col_hseq = c("chartreuse3", "orange", "red")   #--- colors for the timesteps
  
  
  r0 = 0.1                                       #--- intrinsic growth rate
  a=0.5
  e=0.5
  m=0.2
  
  Ninit = 1                                      #--- Initial population density of the prey
  Pinit = 1                                      #--- Initial population density of the predator
  
  
  mat = dyn_lotkaVolterra_euler(Ninit,Pinit,r0,a,e,m,tmax=tmax,h=hseq[1])
  
  ###---- Plot the dynamic for the smaller timestep
  pdf(here::here("img","4_euler_LV.pdf"), width = 5,height = 5.5)
  plot(mat[,"t"],mat[,"N"],type="l",col=col_hseq[1],lty=2,xlab="Time", ylab= "Population density",lwd=2,las=1, main="Example 2, Lotka-Volterra",ylim=c(-1,4.5))
  lines(mat[,"t"],mat[,"P"],type="l",col=col_hseq[1],lwd=2)
  
  ###---- Add the lines corresponding to numerical integration using each timestep value
  for(i in 2:length(hseq)){
    mat = dyn_lotkaVolterra_euler(Ninit,Pinit,r0,a,e,m,tmax=tmax,h=hseq[i])
    lines(mat[,"t"],mat[,"N"],lwd=2,col=col_hseq[i],lty=2)
    lines(mat[,"t"],mat[,"P"],lwd=2,col=col_hseq[i],lty=1)
  }
  legend("topleft",lty=c(2,2,2),lwd=rep(2,4),col=col_hseq,legend=hseq,title=expression(paste(delta, "t")),bty="n")
  dev.off()
  
}



############################################################################___
####----- EXPLORE THE ROSENZWEIG AND MACARTHUR MODEL --------------------------

{
  ###---- Rosenzweig and MacArthur model, the function of the derivatives
  rma_model <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      dP = r0 * P * (1 - P / K) - a * P * H / (1 + a * h * P)
      dH =  eps * a * P * H / (1 + a * h * P) - m * H
      list(c(dP,dH))
    })
  }
  
  ###---- Set the parameter values
  intrinsic_growth_rate = 0.5   
  carrying_capacity = 10
  attack_rate=0.7
  conversion_efficiency=0.3
  mortality_rate=0.4
  handling_time=0.4
  
  ###---- Make a vector of parameters, with the names corresponding to the symbols in the function of the derivatives
  parameters = c(r0 = intrinsic_growth_rate, a = attack_rate, eps = conversion_efficiency, m=mortality_rate, h=handling_time,K=carrying_capacity)
  
  ###---- Set initial conditions for state variables
  N_init = c(P=4,H=2) 
  
  
  ###---- Set the timestep
  timestep =0.1   
  
  ###---- a vector of time steps: times at which to estimate the model
  tmax = 250
  simtime = seq(0,tmax,by = timestep) 
  
  
  ###---- Run the dynamics with the solver and put it in a dataframe
  dyn = data.frame(ode(y = N_init, times = simtime, func = rma_model, parms = parameters, method = "rk4"))
  head(dyn)
  
  ###---- Plot trajectories
  pdf(here::here("img","5_rma_dynamics.pdf"), width = 7,height = 4)
  {
    plot(dyn[, "time"], dyn[, "P"], col="chartreuse3",type = "l", ylim=range(dyn[,c("P","H")]), ylab = "Density", xlab = "Time", main = "Dynamics", las = 1,lwd=2)
    lines(dyn[, "time"], dyn[, "H"], type = "l", col = "darkorange3",lwd=2)
    legend("topright", c("P","H"), lty = 1, col = c("chartreuse3","darkorange3"), bty = "n",lwd=2)  
  }
  dev.off()
  
  ###---- look at the five last time steps
  tail(dyn)
  
  
  ###---- EXPLORE
  
  ###---- Vary initial conditions line 35 => Does long-term result change?
  
  ###---- Vary the timestep
  
  ###---- Vary the parameters => what strategy to explore the model ?
  
}







############################################################################___
####..................................................................
####----- HOW TO ANALYSE A MODEL -------------------------------------
####..................................................................


############################################################################___
####----- Isoclines and phase plane -------------------------------------------
{
  
  ####----- Isocline functions for the coexistence equilibrium
  nullclineP_Hstar = function(r0,K,a,h, P){r0/a*(1+P*a*h)*(1-P/K)}
  nullclineH_Pstar = function(eps,a,h,m){m/(a*(eps-m*h))}
  
  
  ####----- Example of stable point -------------------------------------
  {
    ####----- Parameter values
    r0 = 0.5   
    K = 10
    a=0.7
    eps=0.3
    m=0.4         # m=0.4
    h=0.4
    
    pars = c(r0 = r0, a = a, eps = eps, m=m, h=h,K=K) #--- parameter vector
    
    
    ####----- Set the sequence of P values (x-scale range)
    Ps = seq(0,K,by=0.1)
    
    
    ####----- calculate the isoclines with the functions
    nullclineP = nullclineP_Hstar(r0,K,a,h,Ps)
    nullclineH = nullclineH_Pstar(eps,a,h,m)
    
    
    ####----- Colors
    colo = c("chartreuse3","darkorange3")
    
    
    ####----- Plot of isoclines with trajectories and flow field
    pdf(here::here("img","6_isoclines_rma_3.pdf"), width = 6,height = 5)
    {
      xmax=max(K,max(nullclineH)); ymax=1.6
      plot(NA,xlim = c(0,xmax),ylim = c(0,ymax),xlab="P",ylab="H",las=1,lwd=2)
      abline(h=0,col=colo[2],lwd=2)
      abline(v=c(0,nullclineH),col=colo,lwd=2)
      lines(Ps[nullclineP>=0],nullclineP[nullclineP>=0],col=colo[1],lwd=2)
      legend("topright",lty=1,col=colo,legend = c("P nullclines","H nullclines"),box.col="white",bg="white",lwd=2)
      
      #--- Matrix of all equilibria
      Pstar = nullclineH_Pstar(eps,a,h,m)
      Hstar = nullclineP_Hstar(r0,K,a,h, Pstar)
      alleqs = matrix(c(0,0,K,0,Pstar,Hstar),ncol=2,byrow=T,dimnames=list(NULL,c("P","H")))
      
      #--- Plot equilibria
      points(alleqs[,"P"],alleqs[,"H"],lwd=2)
      
      
      # Run and trace the dynamics (to plot on the phase plane)
      simtime = seq(0,300,by = 0.1)     # simulation time
      yini = c(P=0.5,H=0.5)
      dyn = data.frame(ode(y = yini, times = simtime, func = rma_model, parms = pars, method = "rk4"))
      
      points(yini[1],yini[2],pch=16)                                              #... plot the initial point
      lines(dyn[,"P"],dyn[,"H"])                                                  #... plot the trajectory of the dynamics
      points(dyn[nrow(dyn),"P"],dyn[nrow(dyn),"H"],pch=4,cex=2,lwd=3,col="red")   #... plot the endpoint


      # Run and trace the dynamics (to plot on the phase plane)
      simtime = seq(0,300,by = 0.1)     # simulation time
      yini = c(P=9,H=1.2)
      dyn = data.frame(ode(y = yini, times = simtime, func = rma_model, parms = pars, method = "rk4"))

      points(yini[1],yini[2],pch=16,col="purple")                                              #... plot the initial point
      lines(dyn[,"P"],dyn[,"H"],lty=2,col="purple")                                                  #... plot the trajectory of the dynamics
      points(dyn[nrow(dyn),"P"],dyn[nrow(dyn),"H"],pch=4,cex=2,lwd=3,col="red")   #... plot the endpoint

    }
    dev.off()
    
    
    
    ####----- Plot of isoclines with trajectories using Phase R-package functions
    pdf(here::here("img","6_isoclines_rma_4.pdf"), width = 6,height = 5)
    {
      isocline = nullclines(rma_model,xlim = c(0,xmax) ,ylim = c(0,ymax) , parameters = pars, points = 400, col = colo,state.names = c("P","H"),add=F,add.legend = F,lwd=2)
      legend("topright",lty=1,col=colo,legend = c("P nullclines","H nullclines"),box.col="white",bg="white",lwd=2)
      grid()
      traj = trajectory(rma_model,y0=c(P=0.5,H=0.5),tlim = c(0,300) , parameters = pars, state.names = c("P","H"),add=T)
    }
    dev.off()
    
    
    
    ####----- Plot corresponding dynamics
    pdf(here::here("img","7_dynamics_stable.pdf"), width = 5,height = 3)
    {
      plot(dyn[, "time"], dyn[, "P"], col="chartreuse3",type = "l", ylim=range(dyn[,c("P","H")]), ylab = "Density", xlab = "Time",lwd=2)
      lines(dyn[, "time"], dyn[, "H"], type = "l", col = "darkorange3",lwd=2)
      #legend("topright", c("P","H"), lty = 1, col = c("chartreuse3","darkorange3"), bty = "n",lwd=2)  
    }
    dev.off()
  }
  
  ####----- Example of oscillatory dynamics  -------------------------------------
  {
    ####----- decrease predator mortality
    m=0.35         # m=0.4
    pars = c(r0 = r0, a = a, eps = eps, m=m, h=h,K=K) #--- parameter vector
    yini = c(P=0.5,H=0.5)
    
    
    ####----- Plot of isoclines with trajectories and flow field V2: using Phase R-package functions
    pdf(here::here("img","8_isoclines_rma_V1_oscillatory.pdf"), width = 6,height = 5)
    {
      ####----- calculate the isoclines with the functions
      nullclineP = nullclineP_Hstar(r0,K,a,h,Ps)
      nullclineH = nullclineH_Pstar(eps,a,h,m)
      
      xmax=max(K,max(nullclineH)); ymax=1.6
      plot(NA,xlim = c(0,xmax),ylim = c(0,ymax),xlab="P",ylab="H",las=1,lwd=2)
      abline(h=0,col=colo[2],lwd=2)
      abline(v=c(0,nullclineH),col=colo,lwd=2)
      lines(Ps[nullclineP>=0],nullclineP[nullclineP>=0],col=colo[1],lwd=2)
      legend("topright",lty=1,col=colo,legend = c("P nullclines","H nullclines"),box.col="white",bg="white",lwd=2)
      
      #--- Matrix of all equilibria
      Pstar = nullclineH_Pstar(eps,a,h,m)
      Hstar = nullclineP_Hstar(r0,K,a,h, Pstar)
      alleqs = matrix(c(0,0,K,0,Pstar,Hstar),ncol=2,byrow=T,dimnames=list(NULL,c("P","H")))
      
      #--- Plot equilibria
      points(alleqs[,"P"],alleqs[,"H"],lwd=2)
      
      
      # Run and trace the dynamics (to plot on the phase plane)
      simtime = seq(0,300,by = 0.1)     # simulation time
      yini = c(P=0.5,H=0.5)
      dyn = data.frame(ode(y = yini, times = simtime, func = rma_model, parms = pars, method = "rk4"))
      
      points(yini[1],yini[2],pch=16)                                              #... plot the initial point
      lines(dyn[,"P"],dyn[,"H"])                                                  #... plot the trajectory of the dynamics
      points(dyn[nrow(dyn),"P"],dyn[nrow(dyn),"H"],pch=4,cex=2,lwd=3,col="red")   #... plot the endpoint

    }
    dev.off()
    
    ####----- Plot of isoclines with trajectories and flow field V2: using Phase R-package functions
    pdf(here::here("img","8_isoclines_rma_V2_oscillatory.pdf"), width = 6,height = 5)
    {
      isocline = nullclines(rma_model,xlim = c(0,xmax) ,ylim = c(0,ymax) , parameters = pars, points = 400, col = colo,state.names = c("P","H"),add=F,add.legend = F,lwd=2)
      legend("topright",lty=1,col=colo,legend = c("P nullclines","H nullclines"),box.col="white",bg="white",lwd=2)
      grid()
      traj = trajectory(rma_model,y0=yini,tlim = c(0,300), parameters = pars, state.names = c("P","H"),add=T)
    }
    dev.off()
    
    
    pdf(here::here("img","9_dynamics_oscillatory.pdf"), width = 5,height = 3)
    {
      dyn = data.frame(ode(y = yini, times = simtime, func = rma_model, parms = pars, method = "rk4"))
      plot(dyn[, "time"], dyn[, "P"], col="chartreuse3",type = "l", ylim=range(dyn[,c("P","H")]), ylab = "Density", xlab = "Time",lwd=2)
      lines(dyn[, "time"], dyn[, "H"], type = "l", col = "darkorange3",lwd=2)
      #legend("topright", c("P","H"), lty = 1, col = c("chartreuse3","darkorange3"), bty = "n",lwd=2)  
    }
    dev.off()
    
  }
}


############################################################################___
####----- Local stability analysis --------------------------------------------

{
  # Initial values of the variables
  yini = c(P=1,H=1) 
  m=0.4         # m=0.4
  pars = c(r0 = r0, a = a, eps = eps, m=m, h=h,K=K) #--- parameter vector
  
  
  ####----- Function to evaluate the jacobian matrix
  jac_rma <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      dP = r0 * P * (1 - P / K) - a * P * H / (1 + a * h * P)
      dH =  eps * a * P * H / (1 + a * h * P) - m * H
      list(dP,dH)
    })
  }
  
  
  ####----- Calculate the numerical values of equilibria with the given parameter set
  Pstar = nullclineH_Pstar(eps,a,h,m)
  Hstar = nullclineP_Hstar(r0,K,a,h, Pstar)
  list_eq = list(eq1 = c(P=0,H=0), eq2 = c(P=K,H=0), eq3 = c(P=Pstar, H=Hstar))
  
  # The function stode try to find the stable coexistence one but sometimes fail
  #Eq = stode(y = yini,func=rma_model,parms=pars,positive = TRUE)[[1]]
  
  
  ####----- Calculate the eigenvalues for a given equilibrium
  J = jacobian.full(y=list_eq[[3]],func=jac_rma,parms=pars)   #--- Calculate the jacobian matrix 
  eigen(J)$values                                             #--- Calculate the eigenvalue
  max(Re(eigen(J)$values))                                    #--- Stability criteria: sign of the real part of the dominant eigenvalue
  
  
}


############################################################################___
####----- Bifurcation diagram - Extrema -------------------------------------
####----- Set the parameter values
{
  r0 = 0.5     # r0 = 0.1
  K = 20
  a=0.7
  eps=0.3
  m=0.4
  h=0.4
  pars = c(r0 = r0, a = a, eps = eps, m=m, h=h,K=K)
  Ks = seq(0.5,20,length=70)
  
  simtime = seq(0,2000,by = 0.1)     # simulation time
  yini = c("P"=1,"H"=1)              # initial
  
  ####----- Record the extrema at equilibrium
  eqs = data.frame(P_min=NULL,P_max=NULL,H_min=NULL,H_max=NULL)
  for(i in 1: length(Ks)){
    #---  Update Carrying capacity
    pars["K"] = Ks[i]
    
    #--- Run the dynamics
    dyn = data.frame(ode(y = yini, times = simtime, func = rma_model, parms = pars, method = "rk4"))
    
    #--- Record the min and max at the end of the dynamics
    dynend = dyn[(nrow(dyn)-1000):nrow(dyn),]
    eqs = rbind(eqs,data.frame(P_min=min(dynend[,"P"]),P_max=max(dynend[,"P"]),H_min=min(dynend[,"H"]),H_max=max(dynend[,"H"])))
  }
  
  ####----- Plot the bifurcation diagram
  pdf(here::here("img","10_bifurcation_extrema.pdf"), width = 5,height = 4)
  colos = rep(colo,each=2)
  plot(NA,xlim = range(Ks),ylim=range(eqs),xlab="K",ylab="Extremum densities",las=1)
  for(i in 1:4)points(Ks,eqs[,i],pch=16,col=colos[i])
  legend("topleft",pch=16,col=colo,legend=c("P","H"),bty="n")
  dev.off()
  
}

############################################################################___
####----- Bifurcation diagram - Stability criteria ----------------------------
{
  ####----- Record stability criteria for equilibrium 2 (plant alone) and 3 (coexistence)
  lambdamax2 = lambdamax3 = NULL 
  for(i in 1: length(Ks)){
    #--- Update Carrying capacity
    pars["K"] = Ks[i]
    
    #--- Calculate equilibium 2 and 3 for the given parameter set
    Eq2 = c(P=Ks[i],H=0)
    Pstar=nullclineH_Pstar(eps,a,h,m)
    Eq3 = c(P = Pstar, H= nullclineP_Hstar(r0,Ks[i],a,h, Pstar))
    
    #--- Evaluate the Jacobian for the given equilibrium and parameter set
    J2 = jacobian.full(y=Eq2,func=jac_rma,parms=pars)
    J3 = jacobian.full(y=Eq3,func=jac_rma,parms=pars)
    
    #--- Record dominant eigenvalue
    lambdamax2 = c(lambdamax2,max(Re(eigen(J2)$values)))
    lambdamax3 = c(lambdamax3,max(Re(eigen(J3)$values)))
  }
  
  ####----- Plot the stability criterium
  pdf(here::here("img","11_bifurcation_lambda.pdf"), width = 5,height = 4)
  plot(Ks,lambdamax3,xlab="K",type="n",ylab="real part of dominant eigenvalue",pch=16, las=1,ylim=range(c(lambdamax2,lambdamax3)))
  
  #--- Color in red when unstable
  col_Eq2 = col_Eq3 = rep("black",length(Ks))
  col_Eq2[which(lambdamax2>0)] = "red"
  col_Eq3[which(lambdamax3>0)] = "red"
  points(Ks,lambdamax2,col=col_Eq2)
  points(Ks,lambdamax3,col=col_Eq3,pch=16)
  abline(h=0,lty=3)
  legend("bottomright",pch=c(1,16),legend=c("Eq (2)","Eq (3)"),bty="n")
  dev.off()
  
  
  ####----- Record stability criteria for the steday state found by the stode function
  lambdamax = NULL 
  for(i in 1: length(Ks)){
    pars["K"] = Ks[i]
    Eq = stode(y = yini,func=rma_model,parms=pars,positive = TRUE)[[1]] #!!! Be careful with this function: sometimes doesn't provide the stable equilibrium
    J = jacobian.full(y=Eq,func=jac_rma,parms=pars)
    lambdamax = c(lambdamax,max(Re(eigen(J)$values)))
  }
  
  ####----- Plot the stability criterium
  pdf(here::here("img","12_bifurcation_lambda_stode.pdf"), width = 5,height = 4)
  plot(Ks,lambdamax,xlab="K",ylab="real part of dominant eigenvalue",pch=16, las=1,ylim=range(c(lambdamax2,lambdamax3)))
  abline(h=0,lty=3)
  dev.off()
}


############################################################################___
####----- Sensibility to initial conditions -----------------------------------
{
  
  #--- Lotka-Volterra competition model
  lv_compet_ode <- function(t, y,parms) {
    with(as.list(c(y, parms)), {
      return(list(c(r1*y[1]*(1-(y[1]+a12*y[2])/K1),
                    r2*y[2]*(1-(y[2]+a21*y[1])/K2))))
    })
  }
  
  pdf(here::here("img","13_bistability.pdf"), width = 8,height = 4)
  {
    # One parameter set
    pars=c(r1=0.1,r2=0.2,a12=1.5,a21=1.5,K1=15,K2=20)
    coloN = c("darkorange","cyan3")
    
    #--- Dynamics with the same set but 2 different initial conditions
    simtime = seq(0,600,by = 0.1)
    init1 = c(5,5)
    dyn1 = data.frame(ode(y = init1, times = simtime, func = lv_compet_ode, parms = pars, method = "rk4"))
    init2 = c(10,1)
    dyn2 = data.frame(ode(y = init2, times = simtime, func = lv_compet_ode, parms = pars, method = "rk4"))
    
    #--- Plot the dynamics
    layout(matrix(1:2,ncol=2))
    plot(dyn1[,1],dyn1[,2],xlab="times",ylim=range(dyn1[,2:3]),las=1,ylab="densities",main=paste("P0 =",init1[1],"H0 =",init1[2]),type="l",col=coloN[1],lwd=2)
    lines(dyn1[,1],dyn1[,3],col=coloN[2],lwd=2)
    legend("right",lwd=2,col=coloN,legend=c("N1","N2"),bty="n")
    
    plot(dyn2[,1],dyn2[,2],xlab="times",ylim=range(dyn2[,2:3]),las=1,ylab="densities",main=paste("P0 =",init2[1],"H0 =",init2[2]),type="l",col=coloN[1],lwd=2)
    lines(dyn2[,1],dyn2[,3],col=coloN[2],lwd=2)
  }
  dev.off()
  
  
  #................................................................................
  #---- Systematic search of all equilibria and stability with nleqslv package ----
  #................................................................................
  
  #--- Define a matrix of initial densities from uniform distributions
  yinimat = cbind(N1=runif(1000,min=0,max=20),N2=runif(1000,min=0,max=20))
  head(yinimat)
  
  #--- Define your Lotka-Volterra model
  lv_compet_model <- function(x, parms) {
    with(as.list(parms), {
      return(c(r1*x[1]*(1-(x[1]+a12*x[2])/K1),
               r2*x[2]*(1-(x[2]+a21*x[1])/K2)))
    })
  }
  
  
  #--- STEP 1: Search for all equilibria
  #.....................................
  pars=c(r1=0.1,r2=0.2,a12=1.5,a21=1.5,K1=15,K2=20)
  eqs = searchZeros(yinimat,fn=lv_compet_model,parms=pars,control=list(xtol=10e-11,ftol=10e-11,allowSingular=TRUE))$x
  eqs = ifelse(eqs<0.0000001,0,eqs)           # set to zero the near-zero values (precision issue)
  eqs
  
  #--- STEP 2: Find the stability of each equilibrium
  #..................................................
  
  #--- Define functions to determine stability
  {
    # formulation for jacobian
    jac_lv <- function(t, x, parms) {
      with(as.list(c(x, parms)), {
        list(r1*x[1]*(1-(x[1]+a12*x[2])/K1),
             r2*x[2]*(1-(x[2]+a21*x[1])/K2))
      })
    }
    
    # find the stability criteria
    stability_lv = function(Eq,pars){
      Jac = jacobian.full(y=Eq,func=jac_lv,parms=pars)
      lambdamax = max(Re(eigen(Jac)$values))
      if(lambdamax < 0) return(1)
      else return(0)
    }
  }
  
  
  #--- Find the stability for each equilibrium
  stab = numeric(nrow(eqs))
  for(i in 1:nrow(eqs))stab[i] = stability_lv(eqs[i,],pars)
  cbind(eqs,stab)
  
  # Particularly useful for intractable ode models where multiple stable states arise
  
}



############################################################################___
####---- 2-D parameter exploration for the Rosenzweig-MacArthur model a / h ---------

{
  
  final_state_types = c("none","P","P-H")
  col_coex = c("grey","chartreuse3","purple")
  
  #--- Function to calculate all equilibria from the analytical expressions
  rma_eqs=function(r0,K,a,h,eps,m){
    z = eps-h*m
    P=c(0,K,m/(a*z))
    H=c(0,0,(eps*r0*(a*K*z-m))/(a^2*K*z^2))
    return(cbind(P,H))
  }
  
  #--- Set the parameter values
  r0 = 0.5     
  K = 20
  a=0.7
  eps=0.3
  m=0.4
  h=0.4
  
  #--- Sequences of variation
  n=30
  a_s = seq(0.05,1.2,length=n)
  h_s = seq(0.05,1.2,length=n)

  
  
  #--- Record coexistence and stability while varying a and h
  res = c(h=NULL,a=NULL,coex=NULL,stab=NULL)
  
  for(a in a_s){
    for(h in h_s){
      eqs = rma_eqs(r0,K,a,h,eps,m)                               # calculate the equilibria
      pars = c(r0 = r0, a = a, eps = eps, m=m, h=h,K=K)           # vector of parameters
      stab_eq = numeric(3)
      for(i in 1:3)stab_eq[i]=ifelse(max(Re(eigen(jacobian.full(y=eqs[i,],func=jac_rma,parms=pars))$value))<0,1,0) # Calculate stability for each equilibrium
      if(sum(stab_eq)==1){                             #--- if there is only one stable equilibrium
        if(which(stab_eq==1)==2){                      #--- and this stable equilibrium is the one with only plants
          res = rbind(res,c(h=h,a=a,coex=2,stab=1))    #--- record coexistence and stability
        }else{
          res = rbind(res,c(h=h,a=a,coex=3,stab=1))    #--- else record it is the coexistence equilirbium
        }
      }else{
        res = rbind(res,c(h=h,a=a,coex=3,stab=0))      #--- if there is no stable equilibrium, it is the coexistence one
      }
    }
  } 
  
  ####---- Plot coexistence and stability in 2-2 parameter space
  pdf(here::here("img","14_dim2_traits.pdf"), width = 5,height = 5)
  {
    plot(res[,"h"],res[,"a"],pch=16,xlab="h",ylab="a",col=col_coex[res[,"coex"]],main="coexistence and stability",las=1)
    points(res[which(res[,"stab"]==0),"h"],res[which(res[,"stab"]==0),"a"],pch=16,cex=0.5,col="white")
    legend("topright",pch=16,col=col_coex[c(2,3,3)],legend=c("P","PH stable","PH unstable"),bg="white")
    legend("topright",pch=16,col=c(col_coex[c(2,3)],"white"),pt.cex=c(1,1,0.5),legend=c("P","PH stable","PH unstable"),bty="n")
  }
  dev.off()
  
  # ####---- Plot a vertical transect (h=0.4) of extrema to check stability
  # pars = c(r0 = r0, a = a, eps = eps, m=m, h=0.4,K=K)
  # simtime = seq(0,2000,by = 0.1)     # simulation time
  # yini = c("P"=1,"H"=1)              # initial
  # a_s = seq(0.05,1.2,length=100)
  # 
  # ####----- Record the extrema at equilibrium
  # eqs = data.frame(P_min=NULL,P_max=NULL,H_min=NULL,H_max=NULL)
  # for(a in a_s){
  #   pars["a"] = a
  #   dyn = data.frame(ode(y = yini, times = simtime, func = rma_model, parms = pars, method = "rk4"))
  #   dynend = dyn[(nrow(dyn)-1000):nrow(dyn),]
  #   eqs = rbind(eqs,data.frame(P_min=min(dynend[,"P"]),P_max=max(dynend[,"P"]),H_min=min(dynend[,"H"]),H_max=max(dynend[,"H"])))
  # }
  # 
  # ####----- Plot the bifurcation diagram
  # pdf(here::here("img","15_transect.pdf"), width = 5,height = 4)
  # {
  #   colos = rep(colo,each=2)
  #   plot(NA,xlim = range(a_s),ylim=range(eqs),xlab="a",ylab="Extremum densities",las=1)
  #   for(i in 1:4)points(a_s,eqs[,i],pch=16,col=colos[i])
  #   legend("right",pch=16,col=colo,legend=c("P","H"),bty="n")
  # }
  # dev.off()
  # 
  
}




############################################################################___
####---- Model comparison -----------------------------------------------------

####---- 2-D parameter exploration for the Rosenzweig-MacArthur model a / K ---------
{
  #--- Set the parameter values
  r0 = 0.5     
  K = 20
  a=0.7
  eps=0.3
  m=0.4
  h=0.4

  #--- Set the parameter values
  pars = c(r0 = 0.5, a = 0.7, eps = 0.3, m=0.4, h=0.4,K=20)
  
  #--- Sequences of variation
  n=30
  K_s = seq(0.05,20,length=n)
  a_s = seq(0.05,1.2,length=n)
  
  
  #--- Record coexistence and stability while varying a and h
  res = c(K=NULL,a=NULL,coex=NULL,stab=NULL)
  K = K_s[30]; a=a_s[30]
  for(a in a_s){
    for(K in K_s){
      eqs = rma_eqs(r0,K,a,h,eps,m)                               # calculate the equilibria
      pars = c(r0 = r0, a = a, eps = eps, m=m, h=h,K=K)           # vector of parameters
      stab_eq = numeric(3)
      for(i in 1:3)stab_eq[i]=ifelse(max(Re(eigen(jacobian.full(y=eqs[i,],func=jac_rma,parms=pars))$value))<0,1,0) # Calculate stability for each equilibrium
      if(sum(stab_eq)==1){                             #--- if there is only one stable equilibrium
        if(which(stab_eq==1)==2){                      #--- and this stable equilibrium is the one with only plants
          res = rbind(res,c(K=K,a=a,coex=2,stab=1))    #--- record coexistence and stability
        }else{
          res = rbind(res,c(K=K,a=a,coex=3,stab=1))    #--- else record it is the coexistence equilirbium
        }
      }else{
        res = rbind(res,c(K=K,a=a,coex=3,stab=0))      #--- if there is no stable equilibrium, it is the coexistence one
      }
    }
  }  
  
  ####---- Plot coexistence and stability in 2-2 parameter space
  pdf(here::here("img","16_Compare_2.pdf"), width = 5,height = 5)
  {
    plot(res[,"K"],res[,"a"],pch=16,xlab="K",ylab="a",col=col_coex[res[,"coex"]],main="Holling type II",las=1)
    points(res[which(res[,"stab"]==0),"K"],res[which(res[,"stab"]==0),"a"],pch=16,cex=0.5,col="white")
    legend("topright",pch=16,col=col_coex[c(2,3,3)],legend=c("P","PH stable","PH unstable"),bg="white")
    legend("topright",pch=16,col=c(col_coex[c(2,3)],"white"),pt.cex=c(1,1,0.5),legend=c("P","PH stable","PH unstable"),bty="n")
  }
  dev.off()

}

####---- 2-D parameter exploration for the linear equivalent model a / K ---------
{
   # Function to calculate equilibria from analytical formulas for the model with linear functional response
  rma_linear_eqs=function(r0,K,a,eps,m){
    P=c(0,K,m/(a*eps))
    H=c(0,0,(eps*a*K-m)/(a^2*K*eps))
    return(cbind(P,H))
  }
  
  
  # Function to evaluate the jacobian matrix of the model with linear functional response
  jac_linear <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      dP = r0*P*(1-P/K) - a*P*H
      dH =  eps*a*P*H - m*H
      list(dP,dH)
    })
  }
  
  
  
  #--- Set the parameter values
  pars = c(r0 = 0.5, a = 0.7, eps = 0.3, m=0.4,K=20)
  
  
  #--- Record coexistence and stability while varying a and h
  res = c(K=NULL,a=NULL,coex=NULL,stab=NULL)
  
  for(K in K_s){
    for(a in a_s){
      eqs = rma_linear_eqs(r0,K,a,eps,m)                               # calculate the equilibria
      pars = c(r0 = r0, a = a, eps = eps, m=m, h=h,K=K)           # vector of parameters
      stab_eq = numeric(3)
      for(i in 1:3)stab_eq[i]=ifelse(max(Re(eigen(jacobian.full(y=eqs[i,],func=jac_linear,parms=pars))$value))<0,1,0) # Calculate stability for each equilibrium
      if(sum(stab_eq)==1){
        if(which(stab_eq==1)==2){
          res = rbind(res,c(K=K,a=a,coex=2,stab=1))
        }else{
          res = rbind(res,c(K=K,a=a,coex=3,stab=1))
        }
      }else{
        res = rbind(res,c(K=K,a=a,coex=3,stab=0))
      }
    }
  } 
  
  ####---- Plot coexistence and stability in 2-2 parameter space
  pdf(here::here("img","16_Compare_1_linear.pdf"), width = 5,height = 5)
  {
    plot(res[,"K"],res[,"a"],pch=16,xlab="K",ylab="a",col=col_coex[res[,"coex"]],main="Holling type I",las=1)
    points(res[which(res[,"stab"]==0),"K"],res[which(res[,"stab"]==0),"a"],pch=16,cex=0.5,col="white")
  }
  dev.off()
  
}



