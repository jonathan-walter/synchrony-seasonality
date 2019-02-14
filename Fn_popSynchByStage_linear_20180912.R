## Implement a model that allows for different levels of synchronizing influences at different stages in an annual cycle

## Version with linear dynamics

## this is valuable as a starting point because it's linear so without stage specificity we know from Moran (1953)
## that the population correlation will equal the environmental correlation.

popSynchByStage.linear<-function(beta, Sw, rho.b, rho.w, sd.b, sd.w, 
                                 tmax=100, nlocs=10, ret.Bt=FALSE, ret.eb=FALSE, ret.ew=FALSE)
{
  
  #initialize state variables
  Nt<-matrix(NA, nlocs, tmax)
  Nt[,1]<-rnorm(nlocs)
  Bt<-matrix(NA, nlocs, tmax)
  Wt<-matrix(NA, nlocs, tmax)
  
  #make time series of synchronous variation
  synchts<-function(rho, sd, nlocs, tmax){
    library(MASS)
    cormat<-matrix(rho*sd^2, nlocs, nlocs)
    diag(cormat)<-sd^2
    out<-t(mvrnorm(tmax, mu=rep(0,nlocs), Sigma=cormat, empirical=T))
    return(out)
  }
  
  eb<-synchts(rho.b, sd.b, nlocs, tmax)
  ew<-synchts(rho.w, sd.w, nlocs, tmax)
  #run model
  for(tt in 2:tmax){
    Bt[,tt] = beta*Nt[,(tt-1)] + eb[,tt]
    Wt[,tt] = (1-Sw)*(Bt[,tt]) + ew[,tt]
    Nt[,tt] = Bt[,tt]-Wt[,tt]
  }
  out<-list(Nt=Nt)
  if(ret.Bt){out[["Bt"]]<-Bt}
  if(ret.eb){out[["eb"]]<-eb}
  if(ret.ew){out[["ew"]]<-ew}
  
  return(out)
}
