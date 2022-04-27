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

#### Plotting set up ####

#quartz(height=6, width=4)
pdf(file="fig_analytical.pdf", width=3, height=6)
pal<-colorRampPalette(colors=c("red","white","blue"))
par(mfrow=c(3,1), mar=c(3,1.5,1.5,.5), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(2,3,1,1))

# text(0.0,1.08,"a)", xpd=NA)

#### eB eW heat map ####
f0 = 1.6
kB = 100
s0 = 0
kW = 50
cor.ebew = 0.2
sd.e = 0.01
cor.ebij = seq(0, 1, .05)
cor.ewij = seq(0, 1, .05)
dfrac = 0

results.a <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))

for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    results.a[xx,yy] <- analytical.solution(f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e)
  }
}

image(cor.ebij, cor.ewij, results.a, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(cor.ebij, cor.ewij, results.a, add=T)
text(0.02,1.08,"a)", xpd=NA)
mtext(expression(paste("Synchrony of breeding season environ. (", epsilon[B], ")")), 
      1, outer=F,cex=0.75, line=1.7)
mtext(expression(paste("Synchrony of overwintering")),
      2,outer=F,cex=0.75, line=2.5)
mtext(expression(paste("season environment (", epsilon[W], ")")),
      2,outer=F,cex=0.75, line=1.4)

#### f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .8
cor.ewij = .2
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
#text(1.5, .08, paste(expression("cor.ebij="), cor.ebij, 
#                     expression("& cor.ewij="), cor.ewij), xpd=NA)

mtext(expression(paste("Growth rate (", italic('f')[0], ")")), 1, outer=F,cex=0.75, line=1.7)
mtext(expression(paste("Relative survival (", italic(hat('s'))[0], ")")), 2, outer=F,cex=0.75, line=1.7)

#### cor(eB,eW) figure ####

f0_1 <- 1.2
prop_1 <- -.4
s0_1 <- f0_1*prop_1

f0_2 <- 1.8
prop_2 <- 0
s0_2 <- f0_2*prop_2

f0_3 <- 2.4
prop_3 <- -.1
s0_3 <- f0_3*prop_3

f0_4 <- .6
prop_4 <- -.7
s0_4 <- f0_4*prop_4


kB = 100
kW = 85
cor.ebij = .3
cor.ewij = .1
cor.ebew = seq(-.5,.5,.05)
sd.e = 0.05
dfrac = 0

results.c <- rep(NA, length=length(cor.ebew))
results.d <- rep(NA, length=length(cor.ebew))
results.e <- rep(NA, length=length(cor.ebew))
results.f <- rep(NA, length=length(cor.ebew))

for(xx in 1:length(cor.ebew)) {
    results.c[xx] <- analytical.solution(f0_1, kB, s0_1, kW, cor.ebij, cor.ewij, cor.ebew[xx], sd.e)
    results.d[xx] <- analytical.solution(f0_2, kB, s0_2, kW, cor.ebij, cor.ewij, cor.ebew[xx], sd.e)
    results.e[xx] <- analytical.solution(f0_3, kB, s0_3, kW, cor.ebij, cor.ewij, cor.ebew[xx], sd.e)
    results.f[xx] <- analytical.solution(f0_4, kB, s0_4, kW, cor.ebij, cor.ewij, cor.ebew[xx], sd.e)
    
}

my_min <- min(c(results.c, results.d, results.e, results.f))
my_max <- max(c(results.c, results.d, results.e, results.f))

plot(cor.ebew, results.c, ylim=c(my_min-0.05, my_max+0.1), col="black", lwd=2, type="l",
     xlab="", ylab="")
lines(cor.ebew, results.d,  col="blue", lwd=2)
lines(cor.ebew, results.e,  col="red", lwd=2)
lines(cor.ebew, results.f,  col="darkgreen", lwd=2)
text(-.5,.7,"c)", xpd=NA)
#text(0.35,.08,"b)", xpd=NA)

mtext(expression(paste("Cross season synchrony, cor(", epsilon[B], ",", epsilon[W], ")")), 1, outer=F,cex=0.75, line=1.7)
mtext(expression(paste("Population synchrony")), 2, outer=F,cex=0.75, line=1.7)

dev.off()
