## investigating covariance between two environmental drivers in two locations
## using PRISM seasonal average temperature

rm(list=ls())

library(raster)
library(spdep)
library(rgdal)
library(ncf)
library(stringr)

set.seed(11)

states<-readOGR("/Users/jonathanwalter/Documents/Research/DATA/Basemaps/statesp020.shp")
states<-states[!states$STATE %in% c("Alaska","U.S. Virgin Islands","Puerto Rico","Hawaii"),]

samplepts<-spsample(states, n=1000, type="random")
plot(states)
points(samplepts)

distKm<-gcdist(samplepts@coords[,1], samplepts@coords[,2])
hist(distKm)

years<-1989:2009
months=1:12
basepath="/Users/jonathanwalter/Documents/Research/DATA/PRISM_4km2_gridded/tmean/"
basename="PRISM_tmean_stable_4kmM2_"
ext="_bil.bil"


filelist<-list()

for(yy in years){
  for(mm in months){
    filelist<-c(filelist, paste0(basepath,basename,yy,str_pad(mm,2,"left","0"),ext))
  }
}

tavg_stack<-raster::stack(filelist)

tavg_sampled<-tavg_stack[samplepts,]


years2=1990:2009
summer_tavg<-matrix(nrow=length(samplepts), ncol=length(years2))
winter_tavg<-matrix(nrow=length(samplepts), ncol=length(years2))

years.long<-rep(years, each=12)
months.long<-rep(months, times=21)

for(ii in 1:length(samplepts)){
  for(jj in 1:length(years2)){
    
    summer_tavg[ii,jj]<-mean(tavg_sampled[ii,years.long==years2[jj] & months.long %in% c(6,7,8)])
    winter_tavg[ii,jj]<-mean(c(tavg_sampled[ii,years.long==years2[jj] & months.long %in% c(1,2)],
                               tavg_sampled[ii,years.long==(years2[jj]-1)] & months.long==12))
    
  }
}


cor_b1w2<-matrix(NA, length(samplepts), length(samplepts))

for(ii in 1:length(samplepts)){
  for(jj in 1:ii){
    #correlate breeding season with next winter
    cor_b1w2[ii,jj]<-cor(summer_tavg[ii,1:19], winter_tavg[jj,2:20])
    #cor_b1w2[ii,jj]<-cor(summer_tavg[ii,1:19], winter_tavg[jj,2:20])
  }
}

hist(cor_b1w2)
summary(c(cor_b1w2))

quantile(cor_b1w2[lower.tri(cor_b1w2)], c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=T)

quantile(diag(cor_b1w2),c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=T)
#plot(distKm, cor_b1w2)


# cor_w1b2<-matrix(NA, length(samplepts), length(samplepts))
# 
# for(ii in 2:length(samplepts)){
#   for(jj in 1:ii){
#     cor_w1b2[ii,jj]<-cor(summer_tavg[jj,1:19], winter_tavg[ii,2:20])
#   }
# }
# 
# hist(cor_w1b2)
# plot(distKm, cor_w1b2)

good <- complete.cases(summer_tavg)

summer_tavg <- summer_tavg[good,]
winter_tavg <- winter_tavg[good,]

# sync_eB <- Sncf(x=samplepts@coords[good,1], y=samplepts@coords[good,2], z=summer_tavg, latlon=TRUE, resamp=100)
# sync_eW <- Sncf(x=samplepts@coords[good,1], y=samplepts@coords[good,2], z=winter_tavg, latlon=TRUE, resamp=100)
# sync_eBeW <- Sncf(x=samplepts@coords[good,1], y=samplepts@coords[good,2], z=summer_tavg[,1:19], w=winter_tavg[,2:20], latlon=TRUE, resamp=100)


# pdf("~/GitHub/synchrony-seasonality/cross_season_synch_roughFigs.pdf", onefile=T)
# 
# plot(sync_eB, main="synchrony in breeding season")
# plot(sync_eW, main="synchrony in overwintering season")
# plot(sync_eBeW, main="cross-season synchrony")
# hist(cor_b1w2, main="cross-season synchrony")
# abline(v=median(cor_b1w2, na.rm=T), lwd=2, lty=2)
# hist(diag(cor_b1w2))
# abline(v=median(diag(cor_b1w2), na.rm=T), lwd=2, lty=2)
# #hist(diag(cor_w1b2))
# 
# pal <- colorRampPalette(colors=c("red","lightgrey","blue"))
# pal <- pal(100)
# 
# rampscale <- function(x){
#   
#   out <- x - min(x, na.rm=T)
#   out <- out / max(out, na.rm=T)
#   out <- round(out*100)
#   return(out)
# }
# 
# 
# layout(matrix(1:2,ncol=2,byrow=T), widths=c(0.9,0.1))
# 
# par(mar=c(1,1,1,1))
# plot(samplepts, pch=16, col=pal[rampscale(diag(cor_b1w2))])
# mtext("cor(bi, wi)")
# image(t(matrix(1:100)), col=pal, xaxt="n", yaxt="n")
# axis(2, at=seq(0,1,length.out=5), las=2,
#      labels=round(seq(min(diag(cor_b1w2), na.rm=T), max(diag(cor_b1w2), na.rm=T), length.out=5), digits=2))
# 
# dev.off()
# 
# 
# 
# ## nicer figures for manuscript
# 
# png("~/GitHub/synchrony-seasonality/fig_sup_crossCorHist.png", units="in", res=300, width=4.5, height=4.5)
# 
# par(mar=c(4.1,4.1,1.1,1.1))
# hist(cor_b1w2, xlab="Cross-variable synchrony\n(Pearson correlation)", main="")
# abline(v=median(cor_b1w2, na.rm=T), lwd=2, lty=2)
# 
# dev.off()
# 
# 
# png("~/GitHub/synchrony-seasonality/fig_sup_distdecay.png", units="in", res=300, width=4.5, height=4.5)
# 
# par(mar=c(4.1,4.1,1.1,1.1))
# plot(distKm, cor_b1w2, pch=16, col="lightgrey", cex=0.2, ylim=c(-1,1), xlab="Distance (Km)", 
#      ylab="Cross-variable synchrony")
# lines(sync_eBeW$real$predicted$x, sync_eBeW$real$predicted$y, lwd=2)
# lines(sync_eBeW$real$predicted$x, sync_eBeW$boot$boot.summary$predicted$y[2,], lwd=1, lty=2)
# lines(sync_eBeW$real$predicted$x, sync_eBeW$boot$boot.summary$predicted$y[10,], lwd=1, lty=2)
# 
# dev.off()
# 
