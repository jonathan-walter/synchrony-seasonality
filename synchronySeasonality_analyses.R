## Theoretical analysis of impacts of seasonality in population processes 
## and environmental drivers on spatial synchrony and metapopulation variability

## setup ------------------------------------------------------------------------------------------

rm(list=ls())

setwd("~/GitHub/synchrony-seasonality")

library(rootSolve)
library(parallel)
library(lhs)
library(dplyr)
library(viridis)

## source and define helper functions -------------------------------------------------------------
source("simmod_main.R")
source("simmod_alt_nowinter.R")
source("simmod_alt_sameenv.R")

getNtcor<-function(inlist){ #extract synchrony in Nt from model output
  out<-rep(NA, length(inlist))
  for(ii in 1:length(inlist)){
    tmp<-inlist[[ii]]$Nt[-c(1:burn),]
    out[ii]<-cor(tmp)[2,1]
  }
  return(out)
}

getBtcor<-function(inlist){ #extract synchrony in Bt from model output
  out<-rep(NA, length(inlist))
  for(ii in 1:length(inlist)){
    tmp<-inlist[[ii]]$Bt[-c(1:burn),]
    out[ii]<-cor(tmp)[2,1]
  }
  return(out)
}

getCV <- function(inlist){
  out <- rep(NA, length(inlist))
  for(ii in 1:length(inlist)){
    tmp<-inlist[[ii]]$Nt[-c(1:burn),]
    out[ii] <- var(rowSums(tmp))/mean(rowSums(tmp))
  }
  return(out)
}

# compute analytical solution to linearized model
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
  g <- expression(Neq*exp(f0)*exp(-Neq/kB)*exp(eB)*
                    exp(s0)*exp(-(Neq*exp(f0)*exp(-Neq/kB)*exp(eB))/kW)*exp(eW))
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
  return((PB^2*cov.ebij + PW^2*cov.ewij + 2*PB*PW*cov.ebew)/
           (PB^2*var.e + PW^2*var.e + 2*PB*PW*cov.ebew))
  
}


## Analytical stuff (to be populated from LGS code) -----------------------------------------------


## Simulation -- sensitivity of synchrony to model parameters -------------------------------------

# nn<-200
# set.seed(666)
# hypercube<-optimumLHS(n=nn, k=9)
# saveRDS(hypercube, "hypercube.rds")
hypercube <- readRDS("hypercube.rds")

f0=qunif(hypercube[,1],0.3,2.45)
kB=qunif(hypercube[,2],10,200)
s0=qunif(hypercube[,3],-0.9, 0)*f0
kW=kB*qunif(hypercube[,4],0.1,1)
cor.ebij=qunif(hypercube[,5],0,1)
cor.ewij=qunif(hypercube[,6],0,1)
sd.e=hypercube[,7]
dfrac=qunif(hypercube[,8],0,0.2)
cor.ebew=qunif(hypercube[,9],-0.5,0.5)
cor.eij <- cor.ebij

tmax=10000
burn=1000

ncores=detectCores()-4

cl<-makeCluster(ncores)
main.sample<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                      sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)

nowint.sample<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                        sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)

sameenv.sample<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                         getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)

#helper function for getting cor(Nt)
getNtcor<-function(inlist){ #extract synchrony in Nt
  out<-rep(NA, length(inlist))
  for(ii in 1:length(inlist)){
    tmp<-inlist[[ii]]$Nt[-c(1:burn),]
    out[ii]<-cor(tmp)[2,1]
  }
  return(out)
}


for(ii in 1:length(main.sample)){
  
  tmp <- main.sample[[ii]]$Nt
  plot(tmp[-c(1:burn),1], type="l", main=paste(ii))
  lines(tmp[-c(1:burn),2], col="blue")
  
}


#make data frames for analysis  
sampres.main <- data.frame(rhoN=getNtcor(main.sample), 
                           #scale predictors to standardize effect sizes
                           f0=scale(f0), kB=scale(kB), s0=scale(s0), kW=scale(hypercube[,4]), cor.ebij=scale(cor.ebij),
                           cor.ewij=scale(cor.ewij), sd.e=scale(sd.e), dfrac=scale(dfrac))

sampres.nowint <- data.frame(rhoN=getNtcor(nowint.sample), 
                             #scale predictors to standardize effect sizes
                             f0=scale(f0), kB=scale(kB), cor.ebij=scale(cor.ebij),
                             cor.ewij=scale(cor.ewij), sd.e=scale(sd.e), dfrac=scale(dfrac))

sampres.sameenv <- data.frame(rhoN=getNtcor(sameenv.sample), 
                              #scale predictors to standardize effect sizes
                              f0=scale(f0), kB=scale(kB), s0=scale(s0), kW=scale(hypercube[,4]), cor.ebij=scale(cor.eij),
                              sd.e=scale(sd.e), dfrac=scale(dfrac))

#regression analysis
lm.main.rhoN<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + cor.ewij + cor.ebew + sd.e + dfrac + 
                   cor.ebij:cor.ewij + kB:kW + f0:cor.ebij + s0:cor.ewij, 
                 data=sampres.main, na.action=na.exclude)
summary(lm.main.rhoN)
#hist(resid(lm.main.rhoN))
modsumm.main<-summary(lm.main.rhoN)

lm.nowint.rhoN<-lm(rhoN ~ f0 + kB + cor.ebij + cor.ewij + sd.e + dfrac + cor.ebew +
                     cor.ebij:cor.ewij + f0:cor.ebij, 
                   data=sampres.nowint, na.action=na.exclude)
summary(lm.nowint.rhoN)
#hist(resid(lm.nowint.rhoN))
modsumm.nowint<-summary(lm.nowint.rhoN)

lm.sameenv.rhoN<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + sd.e + dfrac +
                      kB:kW + f0:cor.ebij + s0:cor.ebij, 
                    data=sampres.main, na.action=na.exclude)
summary(lm.sameenv.rhoN)
hist(resid(lm.sameenv.rhoN))
modsumm.sameenv<-summary(lm.sameenv.rhoN)

#compile results into a figure
effects.main <- modsumm.main$coefficients[-1,1]
effects.main <- as.data.frame(effects.main)
effects.main$param <- row.names(effects.main)

effects.sameenv <- modsumm.sameenv$coefficients[-1,1]
effects.sameenv <- as.data.frame(effects.sameenv)
effects.sameenv$param <- row.names(effects.sameenv)

effects.nowint <- modsumm.nowint$coefficients[-1,1]
effects.nowint <- as.data.frame(effects.nowint)
effects.nowint$param <- row.names(effects.nowint)

effects.comb <- left_join(effects.main, effects.sameenv)
effects.comb <- left_join(effects.comb, effects.nowint)

saveRDS(effects.comb, file="~/GitHub/synchrony-seasonality/sim_sensitivity.RDS")

barplot(t(as.matrix(effects.comb[,-2])), beside=T, names.arg=effects.comb$param,
        legend.text=c("main", "same environment", "no overwintering"),
        args.legend=list(x="topleft",bty="n"), las=2)

## TODO: make nice figure for manuscript


## Simulations -- vary environmental synchrony, parameter set 1 -----------------------------------

scentxt<-"set1"

#parameter set 1
rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 1.8
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=25)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=25)
cor.ebew = 0
sd.e = 0.1
dfrac = 0
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                       sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                         sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                          getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)


Ntcor.main<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos))
Ntcor.main<-aggregate(Ntcor.main, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main<-matrix(Ntcor.main$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos))
Ntcor.nowint<-aggregate(Ntcor.nowint, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint<-matrix(Ntcor.nowint$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos))
Ntcor.sameenv<-aggregate(Ntcor.sameenv, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor


CV.main<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, CV=getCV(main.varrhos))
CV.main<-aggregate(CV.main, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
CV.main<-matrix(CV.main$CV, nrow=length(rho), ncol=length(rho))

CV.nowint<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, CV=getCV(nowint.varrhos))
CV.nowint<-aggregate(CV.nowint, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
CV.nowint<-matrix(CV.nowint$CV, nrow=length(rho), ncol=length(rho))

CV.sameenv<-data.frame(cor.ebij=cor.ebij, CV=getCV(sameenv.varrhos))
CV.sameenv<-aggregate(CV.sameenv, by=list(cor.ebij), FUN=mean, na.rm=T)$CV


Btcor.main<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos))
Btcor.main<-aggregate(Btcor.main, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main<-matrix(Btcor.main$Btcor, nrow=length(rho), ncol=length(rho))


diff.NtBtcor.main <- Ntcor.main - Btcor.main

pal<-colorRampPalette(colors=c("red","white","blue"))


## Make figure comparing syncrhony in different seasons
png(paste0("~/GitHub/synchrony-seasonality/fig4_compare_seasons.png"), units="in", res=300, width=6.1, height=2.3)

par(mfrow=c(1,3), mar=c(1.8,1.8,1.5,0), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(1.1,1,0,1))
image(rho, rho, Ntcor.main, xlab="", ylab="", asp=1, main="cor(Nt)",
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main, add=T)
text(0.05,0.95,"a)")

image(rho, rho, Btcor.main, xlab="", ylab="", asp=1, main="cor(Bt)",
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main, add=T)
text(0.05,0.95,"b)")

image(rho, rho, diff.NtBtcor.main, xlab="", ylab="", asp=1, main="cor(Nt)-cor(Bt)",
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main, add=T)
text(0.05,0.95,"c)")

mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.8)
mtext("Overwintering synchrony",2,outer=T,cex=0.8)

dev.off()



#make figure comparing main model to alternate

png(paste0("~/GitHub/synchrony-seasonality/fig3_compare_alternates_",scentxt,".png"), 
    units="in", res=300, width=6.1, height=4.7)

par(mfcol=c(2,3), mar=c(2,2,3.1,1), oma=c(1.5,1.5,0,0),mgp=c(2,0.6,0), tcl=-0.4)

#main model
image(rho, rho, Ntcor.main, xlab="", ylab="", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main, add=T)
mtext("Main model",3,line=0.1,cex=0.7)
text(0.05,0.95,"a)")

image(rho, rho, CV.main, xlab="", ylab="", asp=1,
      main="", col=viridis(50))
contour(rho, rho, CV.main, add=T)
mtext("Main model",3,line=0.1, cex=0.7)
text(0.05,0.95,"d)")

#alternate - same environment
plot(rho, Ntcor.sameenv, pch=16, xlab="", ylab="")
mtext("Alt: same environment",3,line=0.1,cex=0.7)
mtext("Spatial Synchrony",3,line=1.5, cex=0.9)
text(0.05,0.95,"b)")

plot(rho, CV.sameenv, pch=16, xlab="", ylab="")
mtext("Alt: same environment",3,line=0.1,cex=0.7)
mtext("Metapopulation CV",3,line=1.5, cex=0.9)
text(0.05,0.95,"b)")

#alternate - no overwintering
image(rho, rho, Ntcor.nowint, xlab="", ylab="", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.7)
text(0.05,0.95,"a)")

image(rho, rho, CV.nowint, xlab="", ylab="", asp=1,
      main="", col=viridis(50))
contour(rho, rho, CV.nowint, add=T)
mtext("Alt: no overwintering",3,line=0.1, cex=0.7)
text(0.05,0.95,"d)")

mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.9,line=0.1)
mtext("Spatial synchrony of overwintering season environment",2,outer=T,cex=0.9,line=0.1)

dev.off()


