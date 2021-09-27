synchsim_main<-function(r, Kb, Kw, rho.b, rho.w, sd.b, sd.w, Db=0, Dw=0, tmax=100, nlocs=2,
                               ret.Bt=TRUE, ret.eb=FALSE, ret.ew=FALSE)
{
  
  #initialize state variables
  Nt<-matrix(NA, nlocs, tmax)
  Nt[,1]<-runif(nlocs,1,mean(c(Kb,Kw)))
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
  
  #set up dispersal matrices
  Dmat.b<-matrix(Db/nlocs, nlocs, nlocs); diag(Dmat.b)<-1-Db
  Dmat.w<-matrix(Dw/nlocs, nlocs, nlocs); diag(Dmat.w)<-1-Dw
  
  #run model
  for(tt in 2:tmax){
    Bt[,tt] = Nt[,(tt-1)]*exp((r + eb[,tt])*(1-(Nt[,(tt-1)]/Kb)))
    Bt[,tt][Bt[,tt]<0]<-0
    if(!Db==0){Bt[,tt] = Bt[,tt] %*% Dmat.b} #dispersal
    
    Wt[,tt] = max(0, Bt[,tt]/(Kw + ew[,tt])-1)*Bt[,tt] #Winter mortality
    Nt[,tt] = Bt[,tt]-Wt[,tt]
    if(Dw==0){
      Nt[,tt][Nt[,tt]<=0]<-NA #make it so that population can't drop below 0
    }
    else{
      Nt[,tt][Nt[,tt]<0]<-0 #make it so that population can't drop below 0
    }
    
    if(!Dw==0){Nt[,tt] = Nt[,tt] %*% Dmat.w} #Winter dispersal
    if(all(is.na(Nt[,tt])) | sum(Nt[,tt]==0, na.rm=T)){break}
  }
  
  out<-list(Nt=Nt)
  if(ret.Bt){out[["Bt"]]<-Bt}
  if(ret.eb){out[["eb"]]<-eb}
  if(ret.ew){out[["ew"]]<-ew}
  return(out)
  
}