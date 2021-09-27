synchmod_basicricker_multnoise<-function(r, K, rho.b, rho.w, sd.b, sd.w, Db=0, 
                                   tmax=100, nlocs=2, ret.eb=FALSE, ret.ew=FALSE)
{
  
  #initialize state variables
  Nt<-matrix(NA, nlocs, tmax)
  Nt[,1]<-runif(nlocs,1,K)


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
  
  #set up dispersal matrices
  Dmat.b<-matrix(Db/nlocs, nlocs, nlocs); diag(Dmat.b)<-1-Db
  
  #run model
  for(tt in 2:tmax){
    Nt[,tt] = Nt[,(tt-1)]*exp(r*(1-(Nt[,(tt-1)]/K)) + eb[,tt] + ew[,tt]) 
    if(!Db==0){Nt[,tt] = Nt[,tt] %*% Dmat.b} #dispersal
    Nt[,tt][Nt[,tt]<0]<-0 #make it so that population can't drop below 0
    if(sum(Nt[,tt])==0){break}
  }
  
  out<-list(Nt=Nt)
  if(ret.eb){out[["eb"]]<-eb}
  if(ret.ew){out[["ew"]]<-ew}
  return(out)
  
}