## investigating covariance between two environmental drivers in two locations
## using PRISM seasonal average temperature

rm(list=ls())

library(raster)
library(spdep)
library(rgdal)
library(ncf)
library(stringr)

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


cov_b1w2<-matrix(NA, length(samplepts), length(samplepts))

for(ii in 2:length(samplepts)){
  for(jj in 1:ii){
    cov_b1w2[ii,jj]<-cov(summer_tavg[ii,1:19], winter_tavg[jj,2:20])
  }
}

hist(cov_b1w2)


cov_w1b2<-matrix(NA, length(samplepts), length(samplepts))

for(ii in 2:length(samplepts)){
  for(jj in 1:ii){
    cov_w1b2[ii,jj]<-cov(summer_tavg[jj,1:19], winter_tavg[ii,2:20])
  }
}

hist(cov_w1b2)
