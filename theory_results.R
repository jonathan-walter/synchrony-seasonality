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

#### define parameters across runs ####
tmax = 2000
burn = 1000

#### eB eW heat map ####

cor.ebij = seq(0, 1, .05)
cor.ewij = seq(0, 1, .05)

f0 = 2.0
kB = 100
s0 = -0.15
kW = 85
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.a <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))

for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    results.a[xx,yy] <- analytical.solution(f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e)
  }
}

#### f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .8
cor.ewij = .1
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


### Plotting

quartz(height=6, width=3)
pal<-colorRampPalette(colors=c("red","white","blue"))
par(mfrow=c(3,1), mar=c(1,.5,.5,.5), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(3,3,3,1))
# sim 1
# analytical
image(cor.ebij, cor.ewij, results.a, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(cor.ebij, cor.ewij, results.a, add=T)
text(0.0,1.08,"a)", xpd=NA)


# sim 2
# analytical
image(f0, prop, results.b, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", xaxt="n", cex=1.25)
contour(f0, prop, results.b, add=T)
text(0.0,1.08,"b)", xpd=NA)


# sim 3
# analytical
image(cor.ebij, cor.ewij, results.c, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(cor.ebij, cor.ewij, results.c, add=T)
text(0.0,1.08,"e)", xpd=NA)


mtext(expression(paste("Spatial synchrony of breeding season environment (", epsilon[b], ")")), 
      1, outer=T,cex=0.8, line=1.2)
mtext(expression(paste("Spatial synchrony of overwintering season environment (", epsilon[w], ")")),
      2,outer=T,cex=0.8, line=1.2)
mtext(expression(paste("Analytical")), 3, outer=T, cex=1, line=.3, adj=.2)
mtext(expression(paste("Simuation")), 3, outer=T,cex=1, line=.5, adj=.85)


