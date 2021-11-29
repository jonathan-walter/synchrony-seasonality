#library(parallel)
#library(lhs)
library(rootSolve)

## Set up simulation model ----------------------------

simmod_main<-function(tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                      sd.e, dfrac=0, N0=NULL, getBt=FALSE){
  
  library(mvtnorm)
  
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
  
  if(is.null(N0)){
    N0 <- rnorm(2, mean(kB, kW), sd.e)
  }
  
  Nt <- matrix(NA, tmax, 2)
  Nt[1,] <- N0
  
  if(getBt){
    Bt <- matrix(NA, tmax, 2)
  }
  
  for(tt in 2:tmax){
    fN <- exp(f0)*exp(-Nt[tt-1,]/kB)*exp(eb[tt,])
    sN <- min(exp(s0)*exp(-(Nt[tt-1,]*fN)/kW)*exp(ew[tt,]), 1)
    Nt[tt,] <- Nt[tt-1]*fN*sN
    
    if(dfrac > 0){
      Nt[tt,] <- colSums(dmat*Nt[tt,])
    }
    
    if(getBt){
      Bt[tt,] <- Nt[tt-1]*fN
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



## Set parameters ---------------------------------------------------------------------------------

## start with one example, then move to a search over the parameter space
tmax = 1000
burn = 100
f0 = 1.1
kB = 20
s0 = -0.01
kW = 40
cor.ebij = 0.5
cor.ewij = 0.5
cor.ebew = 0.2
sd.e = 0.01
dfrac = 0



## Run a test with above params and do basic exploration/plotting ---------------------------------

linear.test<-simmod_main(tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e, dfrac, getBt=TRUE)


cor(linear.test$Nt[-c(1:burn),])
cor(linear.test$Bt[-c(1:burn),], use="pairwise.complete.obs")
cor(linear.test$env[-c(1:burn),])


plot(linear.test$Nt[,1])
mean(linear.test$Nt[-c(1:burn),1])


## Set up deterministic skeleton, plot, attempt to find equilibria --------------------------------

## here, eB and eW are assumed to be zero, so excluded



rate <- function(N){
  (exp(f0)*exp(-N/kB)) * (exp(s0)*exp(-(N*exp(f0)*exp(-N/kB))/kW))-1
}


# #find derivative
# g <- expression( (N*exp(f0)*exp(-N/kB)) * (exp(s0)*exp(-(N*exp(f0)*exp(-N/kB))/kW)) )
# D(g,"N")
# 
# rate <- function(N){
#   return( (exp(f0) * exp(-N/kB) - N * exp(f0) * (exp(-N/kB) * (1/kB))) * 
#             (exp(s0) * exp(-(N * exp(f0) * exp(-N/kB))/kW)) - (N * exp(f0) * 
#             exp(-N/kB)) * (exp(s0) * (exp(-(N * exp(f0) * exp(-N/kB))/kW) * 
#            ((exp(f0) * exp(-N/kB) - N * exp(f0) * (exp(-N/kB) * (1/kB)))/kW))) )
# }
# 
plot(N, rate(N))

Eq <- uniroot.all(rate, c(0,55))
print(Eq)
# 
eig <- vector()
for (i in 1:length(Eq)){
  eig[i] <- sign (gradient(rate, Eq[i]))
}

print(eig) 
