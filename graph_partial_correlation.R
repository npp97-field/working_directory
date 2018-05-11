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




setwd("/Users/yzhang/Project/SIF_phenology/analysis/")
ncin<-nc_open(paste("./p_corr/pcor_sos_csif_p2s.nc",sep=""))
sos_p2s<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef"))
nc_close(ncin)

ncin<-nc_open(paste("./p_corr/pcor_sos_eos.nc",sep=""))
sos_eos<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef"))
nc_close(ncin)

ncin<-nc_open(paste("./p_corr/pcor_sos_csif_s2p.nc",sep=""))
sos_s2p<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef"))
nc_close(ncin)


pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/pcor_sos_eos_sif.pdf",sep=""),width=11/3*1.2,height=11)

#############################
par(fig=c(0,0.9,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
image(setrange(sos_p2s,-1,1),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-1,1))
text(0, 5500000,expression("SOS~SIF"[Peak_to_sene]),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.66,1),new=T)
plot(sos_p2s, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1, 1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))



par(fig=c(0,0.9,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_eos,-1,1),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-1,1))
text(0, 5500000,expression("SOS~EOS"),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.33,0.66),new=T)
plot(sos_eos, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1,1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))


par(fig=c(0,0.9,0.00,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_s2p,-1,1),add=T,col=rev(ano_discrete_ramp),axes=F,zlim=c(-1,1))
text(0, 5500000,expression("SOS~SIF"[Start_to_peak]),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"c",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0.5,0.96,0,0.33),new=T)
plot(sos_s2p, legend.only=TRUE, 
     col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1,1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))
#box()
#axis(1)
par(fig=c(0,1,0,1),new=T)
plot(NA,axes=F,xlim=c(-1,1),ylim=c(-1.5,1.5),xaxs="i",yaxs='i')

text(0.94, 1.02, labels = expression('partial r'), xpd = NA, srt = -90,cex=1.1) 
text(0.94, 0, labels = expression('partial r'), xpd = NA, srt = -90,cex=1.1)    
text(0.94, -1.03, labels = expression('partial r'), xpd = NA, srt = -90,cex=1.1)
#text(12000000, 0, labels = expression('Day y'^-1), xpd = NA, srt = -90,cex=1.1,col='red')

dev.off()
