##### graph_geographical_data
# plot girded phenology data
# SOS, EOS, LGS, average and trend

library(ncdf4)
library(maptools)
library(raster)
library(rgdal)

lst_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_color.csv",header = F)/255))
ano_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_ano.csv",header = F)/255))
ano_discrete_ramp<-c(ano_color_ramp[c(1,26,51,76,101)],"grey80","grey80",ano_color_ramp[c(140,165,190,215,240)])
ano_discrete_ramp_leg<-rep(ano_discrete_ramp,each=20)

longlat <-  CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
ae<-"+proj=aeqd +lat_0=90 +lon_0=-0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

plot_lat<-function(lat){
  norths<-rep(lat,181)
  easts<- -90:90*2
  
  lati <- cbind(easts, norths)	
  lati <- rbind(lati, lati[1,])	# creates a matrix with two colums
  proj_lat <- SpatialPolygons(list(Polygons(list(Polygon(lati)), 1)))	# create a polgon from this matrix
  proj4string(proj_lat) <- longlat
  replat<-spTransform(proj_lat,ae)
  return(replat)
}

plot_long<-function(lon){
  # if (lon==-180){
  #   norths<- 10:16*5
  #   easts<- rep(lon,7)
  # }else{
  # 
  # }
  norths<- 6:16*5
  easts<- rep(lon,11)
  
  long <- cbind(easts, norths)	
  long <- rbind(long, long[1,])	# creates a matrix with two colums
  proj_lon <- SpatialPolygons(list(Polygons(list(Polygon(long)), 1)))	# create a polgon from this matrix
  proj4string(proj_lon) <- longlat
  replon<-spTransform(proj_lon,ae)
  return(replon)
}


nc2ae<-function(dat){
  latlongdat<-raster(apply(dat,1,rev))
  extent(latlongdat)<-c(-180,180,30,90)
  projection(latlongdat)<-longlat
  extae <- projectExtent(latlongdat, ae)
  aedat<-projectRaster(latlongdat,extae,method = "bilinear")
  return(aedat)
}

plotlatlong<-function(){
  
  for (i in 4:8){
    prjlat<-plot_lat(i*10)
    lines(prjlat,lty=8,col='grey50',lwd=0.5)
  }
  
  for (i in -6:5){
    prolon<-plot_long(i*30)
    lines(prolon,lty=8,col='grey50',lwd=0.5)
  }
  
  lines(repcoa,col='black',lwd=0.2)
}

setrange<-function(dat,a,b){
  dat[dat<a]<-a
  dat[dat>b]<-b
  return(dat)
}


border<-plot_lat(30)
#lat60<-plot_lat(60)

### coastlines
coastline<-shapefile("/Users/yzhang/Data/GIS_data/global/ne_110m_coastline.shp")
cropcoast<-crop(coastline,extent(-180,180,30,90))
repcoa<-spTransform(cropcoast,ae)




indicator<-c("VI",'SIF')
k=1
indi<-indicator[k]

setwd("/Users/yzhang/Project/SIF_phenology/analysis/")
ncin<-nc_open(paste("./",indi,"_SOS_30N_fixed_stat.nc",sep=""))
sossif<-ncvar_get(nc = ncin,varid = "MEAN")
nc_close(ncin)

ncin<-nc_open(paste("./",indi,"_EOS_30N_fixed_stat.nc",sep=""))
eossif<-ncvar_get(nc = ncin,varid = "MEAN")
nc_close(ncin)

ncin<-nc_open(paste("./",indi,"_LGS_30N_fixed_stat.nc",sep=""))
lgssif<-ncvar_get(nc = ncin,varid = "MEAN")
nc_close(ncin)

k=2
indi<-indicator[k]

ncin<-nc_open(paste("./",indi,"_SOS_30N_fixed_stat.nc",sep=""))
sossif2<-ncvar_get(nc = ncin,varid = "MEAN")
nc_close(ncin)

ncin<-nc_open(paste("./",indi,"_EOS_30N_fixed_stat.nc",sep=""))
eossif2<-ncvar_get(nc = ncin,varid = "MEAN")
nc_close(ncin)

ncin<-nc_open(paste("./",indi,"_LGS_30N_fixed_stat.nc",sep=""))
lgssif2<-ncvar_get(nc = ncin,varid = "MEAN")
nc_close(ncin)



aesos_diff<-nc2ae(sossif2-sossif)-6/365
aeeos_diff<-nc2ae(eossif2-eossif)-6/365
aelgs_diff<-nc2ae(lgssif2-lgssif)


pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/pheno_map_fix_diff.pdf",sep=""),width=11/3*1.2,height=11)

#############################
par(fig=c(0,0.9,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
image(setrange(aesos_diff,-0.1643836,0.1643836),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-0.1643836,0.1643836))
text(0, 5500000,expression("SOS"[SIF]-"SOS"[NDVI]),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.66,1),new=T)
plot(aesos, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-0.1643836,0.1643836),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-60, 60, 10)/365,
                    labels=c('-60','','-40','','-20','','0','','20','','40','','60'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))



par(fig=c(0,0.9,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aeeos_diff,-0.2465753,0.2465753),add=T,col=ano_discrete_ramp,
      axes=F,zlim=c(-0.2465753,0.2465753))
text(0, 5500000,expression("EOS"[SIF]-"EOS"[NDVI]),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.33,0.66),new=T)
plot(aesos, legend.only=TRUE, col=ano_discrete_ramp_leg,horizontal=F,zlim=c(-0.2465753,0.2465753),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-90, 90, 15)/365,
                    labels=c('-90','','-60','','-30','','0','','30','','60','','90'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))


par(fig=c(0,0.9,0.00,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aelgs_diff,-0.3287671,0.3287671),add=T,col=ano_discrete_ramp,axes=F,zlim=c(-0.3287671,0.3287671))
text(0, 5500000,expression("LGS"[SIF]-"LGS"[NDVI]),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"c",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0.5,0.96,0,0.33),new=T)
plot(aelgs_t, legend.only=TRUE, 
     col=ano_discrete_ramp_leg,horizontal=F,zlim=c(-0.3287671,0.3287671),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-120, 120, 20)/365,
                    labels=c('-120','','-80','','-40','','0','','40','','80','','120'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))
#box()
#axis(1)
par(fig=c(0,1,0,1),new=T)
plot(NA,axes=F,xlim=c(-1,1),ylim=c(-1.5,1.5),xaxs="i",yaxs='i')

text(0.94, 1.02, labels = expression('Day'), xpd = NA, srt = -90,cex=1.1) 
text(0.94, 0, labels = expression('Day'), xpd = NA, srt = -90,cex=1.1)    
text(0.94, -1.03, labels = expression('Day'), xpd = NA, srt = -90,cex=1.1)
#text(12000000, 0, labels = expression('Day y'^-1), xpd = NA, srt = -90,cex=1.1,col='red')

dev.off()

