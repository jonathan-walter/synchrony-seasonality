rm(list=ls())

library(rootSolve)
library(mvtnorm)

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

simmod_det<-function(tmax, f0, kB, s0, kW, getBt=FALSE){
  
  if(exp(s0) > 1){
    stop("exp(s0) must be <= 1")
  }
  
  N0 <- 1
  
  Nt <- rep(NA, tmax)
  Nt[1] <- N0
  
  if(getBt){
    Bt <- rep(NA, tmax)
  }
  
  for(tt in 2:tmax){
    fN <- exp(f0)*exp(-Nt[tt-1]/kB)
    sN <- min(exp(s0)*exp(-(Nt[tt-1]*fN)/kW))
    Nt[tt] <- Nt[tt-1]*fN*sN
    
    if(getBt){
      Bt[tt] <- Nt[tt-1,]*fN
    }
  }
  if(getBt){
    out <- list(Nt=Nt, Bt=Bt)
  }
  else{
    out <- list(Nt=Nt)
  }
  return(out)
}

#helper function for checking regime
regime.check <- function(inlist){
  out <- rep(NA, length(inlist))
  for(ii in 1:length(out)){
    tmp <- inlist[[ii]]
    if(tmp[length(tmp)] <= 0.1){
      out[ii] <- "check"
    }
    else if(any(tmp > tmp[length(tmp)])){
      out[ii] <- "overcompensatory"
    }
    else{
      out[ii] <- "undercompensatory"
    }
  }
  return(out)
}

#### define parameters across runs ####
tmax = 2000
burn = 1000

#### 1. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .8
cor.ewij = .2
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.a <- matrix(NA, nrow=length(f0), ncol=length(prop))
regimes <- matrix(NA, nrow=length(f0), ncol=length(prop))
#det.a <- array(NA, c(length(f0), length(prop), tmax))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.a[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
    temp <- simmod_det(tmax, f0[xx], kB, s0, kW)
    #det.a[xx,yy,] <- simmod_det(tmax, f0[xx], kB, s0, kW)
    regimes[xx,yy] <- regime.check(temp)
  }
}
regimes_clean <- ifelse(regimes == "overcompensatory", 1, 0)

quartz(height=6, width=6)
pal<-colorRampPalette(colors=c("red","white","blue"))
par(mfrow=c(3,2), mar=c(1,1.5,2.5,1), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(3,3,3,1))
# sim 1
# analytical
image(f0, prop, results.a, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(f0, prop, results.a, add=T)
text(0.35,.08,"a)", xpd=NA)
text(1.3,-1.1, paste(expression("cor.ebij="), cor.ebij, 
                     expression("and cor.ewij="), cor.ewij), xpd=NA)

image(f0, prop, regimes_clean, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)

#### 2. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .2
cor.ewij = .8
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.b <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.b[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}


image(f0, prop, results.b, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(f0, prop, results.b, add=T)
text(0.35,.08,"b)", xpd=NA)
text(1.3,-1.1, paste(expression("cor.ebij="), cor.ebij, 
                     expression("and cor.ewij="), cor.ewij), xpd=NA)

#### 3. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .2
cor.ewij = .1
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.c <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.c[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}

image(f0, prop, results.c, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(f0, prop, results.c, add=T)
text(0.35,.08,"c)", xpd=NA)
text(1.3,-1.1, paste(expression("cor.ebij="), cor.ebij, 
                     expression("and cor.ewij="), cor.ewij), xpd=NA)


#### 4. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .1
cor.ewij = .2
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.d <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.d[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}

image(f0, prop, results.d, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(f0, prop, results.d, add=T)
text(0.35,.08,"d)", xpd=NA)
text(1.3,-1.1, paste(expression("cor.ebij="), cor.ebij, 
                     expression("and cor.ewij="), cor.ewij), xpd=NA)


#### 5. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .6
cor.ewij = .1
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.e <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.e[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}

image(f0, prop, results.e, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(f0, prop, results.e, add=T)
text(0.35,.08,"e)", xpd=NA)
text(1.3,-1.1, paste(expression("cor.ebij="), cor.ebij, 
                     expression("and cor.ewij="), cor.ewij), xpd=NA)


#### 6. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .1
cor.ewij = .6
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.f <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.f[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}

image(f0, prop, results.f, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(f0, prop, results.f, add=T)
text(0.35,.08,"f)", xpd=NA)
text(1.3,-1.1, paste(expression("cor.ebij="), cor.ebij, 
                     expression("and cor.ewij="), cor.ewij), xpd=NA)
