###################### 1. GCM data ######################

# https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.derived.surface.html
# fileurl="ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc"
# download.file(fileurl,destfile="./data/air_mon_mean.nc")

# library(ncdf4)
# ncin=nc_open(fileloc)
# library(RNetCDF)
# ncin=open.nc(fileloc)

library(rhdf5)
fileloc="./data/air_mon_mean.nc"
h5ls(fileloc)
# group name       otype dclass            dim
# 0     /  air H5I_DATASET  FLOAT 144 x 73 x 830
# 1     /  lat H5I_DATASET  FLOAT             73
# 2     /  lon H5I_DATASET  FLOAT            144
# 3     / time H5I_DATASET  FLOAT            830

LON=h5read(fileloc,"lon")
LAT=h5read(fileloc,"lat")
X=h5read(fileloc,"air") # lon*lat*time unit degC
NLON=length(LON)
NLAT=length(LAT)

LAT=rev(LAT)
LON=c(LON[(NLON/2+1):NLON]-360,LON[1:(NLON/2)])
X=X[, length(LAT):1, ]
X=X[c((NLON/2+1):NLON,1:(NLON/2)),,]
X=X[,,1:(dim(X)[3]-2)] # 1948-2016

# remove monthly mean
X_mon_mean=X[,,1:12]
for(i in 0:11){
  X0=X[,,(1:dim(X)[3])%%12==i]
  X_mon_mean[,,i+1]=apply(X0,c(1,2),mean) # [12,1,2,3,...,11]
}

for(i in 1:dim(X)[3]){
  X[,,i]=X[,,i]-X_mon_mean[,,(i%%12)+1]
}

############### 1.1 Get subsample ################
idsellon=2*(1:(NLON/2))-1
idsellat=2*(1:(NLAT/2))-1
LON=LON[idsellon]
LAT=LAT[idsellat]
NLON=length(LON)
NLAT=length(LAT)
X=X[idsellon,idsellat,]
save(LON,LAT,NLON,NLAT,X,file="./data/air_mon_mean_mon_mean_removed_sub.RData")

############### 1.2 EDA ################
# Plots
X1=X[,,1]
#X1=apply(X,c(1,2),mean)
image(LON,LAT,X1)
library(maptools)
data(wrld_simpl)
plot(wrld_simpl,add=TRUE)


#

x1=X[which(LON==-25),which(LAT==37.5),]
x2=X[which(LON==-17.5),which(LAT==65),]


x1=X[which(LON==-120),which(LAT==0),]
x2=X[which(LON==-70),which(LAT==45),]
plot(x1,type="l",col="blue")
points(x2,type="l",col="red")


plot(density(x1),col="blue")
lines(density(x2),col="red")

plot(x1,x2)
cor(x1,x2)

# plot(X_mon_mean[100,30,],col="blue",ylim=c(min(X_mon_mean),max(X_mon_mean)))
# for(i in 100:144){
#   for(j in 30:36){
#     lines(X_mon_mean[i,j,],col=rgb((i+j)/250,0,0))
#   }
# }
# 
# 
# 
# plot(density(X[1,1,]),col="white",ylim=c(0,1))
# for(i in 1:144){
#   lines(density(X[i,30,]),col=rgb((i)/144,0,0))
# }


#
D=matrix(1,144*73,144*73)
for(i1 in 1:length(LON)){
  for(j1 in 1:length(LAT)){
    for(i2 in 1:length(LON)){
      for(j2 in 1:length(LAT)){
        i=73*(i1-1)+j1
        j=73*(i2-1)+j2
        x=X[i1,j1,]
        y=X[i2,j2,]
        if(i>j){
          D[i,j]=D[j,i]=cor(x,y)
        }
        
      }
    }
  }
}
  

save(D,file="air_mon_mean_cor.RData")


###################### 2. OBS data ######################
############### 2.1 Daily to Monthly with CDO ###############
# cannot run in windows, since cdo is a linux command tool
for(i in 1979:2016){
  cat("CDO: year",i,"...\n")
  cmd=paste0("cdo -monavg ","tmmx_",i,".nc monavg_tmmx_",i,".nc")
  system(cmd)
  cmd=paste0("cdo -monavg ","tmmn_",i,".nc monavg_tmmn_",i,".nc")
  system(cmd)
}

############### 2.2 Load Monthly data ###############
library(RNetCDF)
library(abind)
# lon: 1386
# lat: 585
# day: 12
# air_temperature: (lat, lon, day)
x=NULL
for(i in 1979:2016){
  fileloc=paste0("D:/works/obs_data_mon_avg/monavg_tmmx_",i,".nc")
  nc=open.nc(fileloc) # print.nc(nc)
  x=abind(x,aperm(var.get.nc(nc, "air_temperature"),c(2,1,3)),along=3)
  
}
lon=var.get.nc(nc, "lon")
lat=var.get.nc(nc, "lat")
lat=rev(lat)
nlon=length(lon)
nlat=length(lat)







