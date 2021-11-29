rm(list=ls())

library(rootSolve)
library(mvtnorm)


## Set up simulation model ------------------------------------------------------------------------

## tmax is the length of the simulation; f0 is the fecundity rate at low density; kB is the
## saturation parameter during the breeding season; s0 is the overwintering survival rate at low density;
## kW is the saturation parameter during the overwintering season; cor.ebij is the spatial synchrony of
## breeding season environment; cor.ewij is the spatial synchrony of the overwintering environment;
## cor.ebew is the correlation between breeding and overwintering environments;
## sd.e is the standard deviation of the environmental noises; N0 are initial population abundances.

simmod_main<-function(tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                      sd.e, dfrac=0, getBt=FALSE){

  
  if(exp(s0) > 1){
    stop("exp(s0) must be <= 1")
  }
  
  #create environmental noise time series
  sigma<-matrix(0, 4, 4)
  
  sigma[2,1] <- cor.ebij*(sd.e^2)
  sigma[3:4,1] <- cor.ebew*(sd.e^2)
  sigma[3:4,2] <- cor.ebew*(sd.e^2)
  sigma[4,3] <- cor.ewij*(sd.e^2)  
  sigma <- sigma+t(sigma)
  diag(sigma) <- rep(sd.e^2, 4)
  
  env <- rmvnorm(tmax,sigma=sigma)
  colnames(env) <- c("ebi","ebj","ewi","ewj")
  eb <- env[,1:2]
  ew <- env[,3:4]
  
  #if needed, create dispersal matrix
  if(dfrac > 0){
    dmat <- matrix(c(1-dfrac, dfrac, dfrac, 1-dfrac), 2, 2)
    
  }
  
  N0 <- rnorm(2, mean(c(kB, kW)), sd.e)
  
  Nt <- matrix(NA, tmax, 2)
  Nt[1,] <- N0
  
  if(getBt){
    Bt <- matrix(NA, tmax, 2)
  }
  
  mymin <- function(x){
    out <- rep(NA, length(x))
    for(ii in 1:length(x)){
      out[ii] <- min(c(x[ii]),1)
    }
    return(out)
  }
  
  for(tt in 2:tmax){
    fN <- exp(f0)*exp(-Nt[tt-1,]/kB)*exp(eb[tt,])
    sN <- mymin(exp(s0)*exp(-(Nt[tt-1,]*fN)/kW)*exp(ew[tt,]))
    Nt[tt,] <- Nt[tt-1,]*fN*sN
    
    if(dfrac > 0){
      Nt[tt,] <- colSums(dmat*Nt[tt,])
    }
    
    if(getBt){
      Bt[tt,] <- Nt[tt-1,]*fN
    }
  }
  if(getBt){
    out <- list(Nt=Nt, Bt=Bt, env=env)
  }
  else{
    out <- list(Nt=Nt, env=env)
  }
  return(out)
}



##-------------------------------------------------------------------------------------------------
## Set up function for analytical solution

## parameter definitions as above.


analytical.solution<-function(f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e){
  
  #find equilibrium value of N; this is when rate = 0
  rate <- function(N){
    exp(f0)*exp(-N/kB) * exp(s0)*exp(-(N*exp(f0)*exp(-N/kB))/kW) - 1
  }
  
  Eq <- uniroot.all(rate, c(1, max(c(kB,kW))*1.5))
  
  #test stability of equilibrium
  eig <- vector()
  for (i in 1:length(Eq)){
    eig[i] <- sign (gradient(rate, Eq[i]))
  }
  
  if(!any(eig==-1)){
    stop("No stable equilibria between N=1 and N=1.5*max(kB, kW)")
  }
  
  Neq <- Eq[eig==-1] #take the stable equilibrium if multiple
  
  
  # #compute derivatives at N
  g <- expression(Neq*exp(f0)*exp(-Neq/kB)*exp(eB)*exp(s0)*exp(-(Neq*exp(f0)*exp(-Neq/kB)*exp(eB))/kW)*exp(eW))
  dgdeb <- D(g, 'eB')
  dgdew <- D(g, 'eW')
  eB <- eW <- 0
  PB <- as.numeric(eval(dgdeb))
  PW <- as.numeric(eval(dgdew))
  
  #convert between correlation/sd and covariance/variance
  cov.ebij <- cor.ebij*(sd.e^2)
  cov.ewij <- cor.ewij*(sd.e^2)
  cov.ebew <- cor.ebew*(sd.e^2)
  var.e <- sd.e^2
  
  # this is equation 26 from Dan's math document
  return((PB^2*cov.ebij + PW^2*cov.ewij + 2*PB*PW*cov.ebew)/(PB^2*var.e + PW^2*var.e + 2*PB*PW*cov.ebew))
  
}


## define parameters and run ----------------------------------------------------------------------
tmax = 10000
burn = 200
f0 = 1.1
kB = 20
s0 = -0.01
kW = 10
cor.ebij = 0.7
cor.ewij = 0.4
cor.ebew = -0.3
sd.e = 0.01
dfrac = 0



linear.test<-simmod_main(tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e, dfrac, getBt=TRUE)
cor(linear.test$Nt[-c(1:burn),])

analytical.solution(f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)


