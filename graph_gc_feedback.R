####### plot the figure for the feedback

library(ncdf4)
library(maptools)
library(raster)
library(rgdal)
library(R.matlab)

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

Asign<-readMat("/Users/yzhang/Project/SIF_GC/analysis/ERA_north/A_sign.mat")
Sig<-readMat("/Users/yzhang/Project/SIF_GC/analysis/ERA_north/significance5_mat.mat")
gc<-readMat("/Users/yzhang/Project/SIF_GC/analysis/ERA_north/g_causality.mat")
landmask_f<-nc_open("/Users/yzhang/Project/SIF_phenology/data/global_land_mask.nc")
landmask<-ncvar_get(landmask_f)[,241:360]

PAR_SIF<-Asign$A.sign[,,1,3]*Sig$significance5[,,1,3]*gc$g.causality[,,1,3]*landmask
par_sif_ae<-nc2ae(PAR_SIF)
SIF_PAR<-Asign$A.sign[,,3,1]*Sig$significance5[,,3,1]*gc$g.causality[,,3,1]*landmask
sif_par_ae<-nc2ae(SIF_PAR)
loop_ae<-nc2ae(PAR_SIF*SIF_PAR)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/SIF_PAR_GC_era.pdf",sep=""),width=11/3*1.2,height=11)

#############################
par(fig=c(0,0.9,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),oma=c(0,0,0,2))
plot(border)
image(setrange(par_sif_ae,-0.3,0.3),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-0.3,0.3))
text(0, 5500000,expression("PAR"%->%"SIF"),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.66,1),new=T)
plot(par_sif_ae, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-0.3,0.3),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-0.3, 0.3, 0.05),
                    c('-0.3','','-0.2','','-0.1','','0','','0.1','','0.2','','0.3'),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

par(fig=c(0,0.9,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sif_par_ae,-0.3,0.3),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-0.3,0.3))
text(0, 5500000,expression("SIF"%->%"PAR"),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,1/3,2/3),new=T)
plot(sif_par_ae, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-0.3,0.3),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-0.3, 0.3, 0.05),
                    c('-0.3','','-0.2','','-0.1','','0','','0.1','','0.2','','0.3'),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

par(fig=c(0,0.9,0,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(loop_ae,-0.06,0.06),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-0.06,0.06))
text(0, 5500000,expression("PAR"%->%"SIF"%->%"PAR"),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"c",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0,1/3),new=T)
plot(loop_ae, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-0.06,0.06),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-0.06, 0.06, 0.01),
                    c('-0.06','','-0.04','','-0.02','','0','','0.02','','0.04','','0.06'),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))
dev.off()
