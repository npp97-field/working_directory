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

Asign<-readMat("/Users/yzhang/Project/SIF_GC/analysis/rs_north/A_sign.mat")
Sig<-readMat("/Users/yzhang/Project/SIF_GC/analysis/rs_north/significance5_mat.mat")
sub_mon<-readMat("/Users/yzhang/Project/SIF_GC/analysis/rs_north/sub_mon.mat")
mon_seasonal<-readMat("/Users/yzhang/Project/SIF_GC/analysis/rs_north/mon_seasonal.mat")
seasonal_annual<-readMat("/Users/yzhang/Project/SIF_GC/analysis/rs_north/seasonal_annual.mat")
inter_annual<-readMat("/Users/yzhang/Project/SIF_GC/analysis/rs_north/inter_annual.mat")
landmask_f<-nc_open("/Users/yzhang/Project/SIF_phenology/data/global_land_mask.nc")
landmask<-ncvar_get(landmask_f)[,241:360]

var_SIF<-list()
###PAR
var_SIF[[1]]<-Asign$A.sign[,,1,3]*Sig$significance5[,,1,3]*sub_mon$sub.mon[,,1,3]*landmask
var_SIF[[2]]<-Asign$A.sign[,,1,3]*Sig$significance5[,,1,3]*mon_seasonal$mon.seasonal[,,1,3]*landmask
var_SIF[[3]]<-Asign$A.sign[,,1,3]*Sig$significance5[,,1,3]*seasonal_annual$seasonal.annual[,,1,3]*landmask
var_SIF[[4]]<-Asign$A.sign[,,1,3]*Sig$significance5[,,1,3]*inter_annual$inter.annual[,,1,3]*landmask
#PREC
var_SIF[[5]]<-Asign$A.sign[,,1,2]*Sig$significance5[,,1,2]*sub_mon$sub.mon[,,1,2]*landmask
var_SIF[[6]]<-Asign$A.sign[,,1,2]*Sig$significance5[,,1,2]*mon_seasonal$mon.seasonal[,,1,2]*landmask
var_SIF[[7]]<-Asign$A.sign[,,1,2]*Sig$significance5[,,1,2]*seasonal_annual$seasonal.annual[,,1,2]*landmask
var_SIF[[8]]<-Asign$A.sign[,,1,2]*Sig$significance5[,,1,2]*inter_annual$inter.annual[,,1,2]*landmask
###temp
var_SIF[[9]]<-Asign$A.sign[,,1,4]*Sig$significance5[,,1,4]*sub_mon$sub.mon[,,1,4]*landmask
var_SIF[[10]]<-Asign$A.sign[,,1,4]*Sig$significance5[,,1,4]*mon_seasonal$mon.seasonal[,,1,4]*landmask
var_SIF[[11]]<-Asign$A.sign[,,1,4]*Sig$significance5[,,1,4]*seasonal_annual$seasonal.annual[,,1,4]*landmask
var_SIF[[12]]<-Asign$A.sign[,,1,4]*Sig$significance5[,,1,4]*inter_annual$inter.annual[,,1,4]*landmask


SIF_var<-list()
### PAR
SIF_var[[1]]<-Asign$A.sign[,,3,1]*Sig$significance5[,,3,1]*sub_mon$sub.mon[,,3,1]*landmask
SIF_var[[2]]<-Asign$A.sign[,,3,1]*Sig$significance5[,,3,1]*mon_seasonal$mon.seasonal[,,3,1]*landmask
SIF_var[[3]]<-Asign$A.sign[,,3,1]*Sig$significance5[,,3,1]*seasonal_annual$seasonal.annual[,,3,1]*landmask
SIF_var[[4]]<-Asign$A.sign[,,3,1]*Sig$significance5[,,3,1]*inter_annual$inter.annual[,,3,1]*landmask
### PAR
SIF_var[[5]]<-Asign$A.sign[,,2,1]*Sig$significance5[,,2,1]*sub_mon$sub.mon[,,2,1]*landmask
SIF_var[[6]]<-Asign$A.sign[,,2,1]*Sig$significance5[,,2,1]*mon_seasonal$mon.seasonal[,,2,1]*landmask
SIF_var[[7]]<-Asign$A.sign[,,2,1]*Sig$significance5[,,2,1]*seasonal_annual$seasonal.annual[,,2,1]*landmask
SIF_var[[8]]<-Asign$A.sign[,,2,1]*Sig$significance5[,,2,1]*inter_annual$inter.annual[,,2,1]*landmask
### PAR
SIF_var[[9]]<-Asign$A.sign[,,4,1]*Sig$significance5[,,4,1]*sub_mon$sub.mon[,,4,1]*landmask
SIF_var[[10]]<-Asign$A.sign[,,4,1]*Sig$significance5[,,4,1]*mon_seasonal$mon.seasonal[,,4,1]*landmask
SIF_var[[11]]<-Asign$A.sign[,,4,1]*Sig$significance5[,,4,1]*seasonal_annual$seasonal.annual[,,4,1]*landmask
SIF_var[[12]]<-Asign$A.sign[,,4,1]*Sig$significance5[,,4,1]*inter_annual$inter.annual[,,4,1]*landmask
###
sif_var_ae<-list()
var_sif_ae<-list()
for (i in 1:12){
  sif_var_ae[[i]]<-nc2ae(SIF_var[[i]])
  var_sif_ae[[i]]<-nc2ae(var_SIF[[i]])
}

x_title<-c("Submonthly",'Subseasonal',"Seasonal","Interannual")
pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/SIF_PAR_GC_freqency.pdf",sep=""),width=4*4.3,height=4*6.1)

#############################
par(mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),oma=c(1,3,3,6))
for (i in 1:4){
  par(fig=c(i/4-0.25,i/4,5/6,1),new=T)
  plot(border)
  image(setrange(var_sif_ae[[i]],-0.3,0.3),add=T,col=rev(ano_discrete_ramp),
        axes=F,zlim=c(-0.3,0.3))
  mtext(x_title[i],3,line=1,cex=2.5)
  plotlatlong()
  if (i==1){
    mtext(expression("PAR"%->%"SIF"),2,line=1,cex=2.5)
  }
  
  par(fig=c(i/4-0.25,i/4,2/3,5/6),new=T)
  plot(border)
  image(setrange(sif_var_ae[[i]],-0.3,0.3),add=T,col=rev(ano_discrete_ramp),
        axes=F,zlim=c(-0.3,0.3))
  plotlatlong()
  if (i==1){
    mtext(expression("SIF"%->%"PAR"),2,line=1,cex=2.5)
  }
  #############
  par(fig=c(i/4-0.25,i/4,1/2,2/3),new=T)
  plot(border)
  image(setrange(var_sif_ae[[i+4]],-0.3,0.3),add=T,col=rev(ano_discrete_ramp),
        axes=F,zlim=c(-0.3,0.3))
  plotlatlong()
  if (i==1){
    mtext(expression("Precip."%->%"SIF"),2,line=1,cex=2.5)
  }
  
  par(fig=c(i/4-0.25,i/4,1/3,1/2),new=T)
  plot(border)
  image(setrange(sif_var_ae[[i+4]],-0.3,0.3),add=T,col=rev(ano_discrete_ramp),
        axes=F,zlim=c(-0.3,0.3))
  plotlatlong()
  if (i==1){
    mtext(expression("SIF"%->%"Precip."),2,line=1,cex=2.5)
  }
  ################
  par(fig=c(i/4-0.25,i/4,1/6,1/3),new=T)
  plot(border)
  image(setrange(var_sif_ae[[i+8]],-0.3,0.3),add=T,col=rev(ano_discrete_ramp),
        axes=F,zlim=c(-0.3,0.3))
  plotlatlong()
  if (i==1){
    mtext(expression("Temp."%->%"SIF"),2,line=1,cex=2.5)
  }
  
  par(fig=c(i/4-0.25,i/4,0,1/6),new=T)
  plot(border)
  image(setrange(sif_var_ae[[i+8]],-0.3,0.3),add=T,col=rev(ano_discrete_ramp),
        axes=F,zlim=c(-0.3,0.3))
  plotlatlong()
  if (i==1){
    mtext(expression("SIF"%->%"Temp."),2,line=1,cex=2.5)
  }
}

par(fig=c(0.5,0.98,1/3,2/3),oma = c(1,1,1,1),new=T)
plot(sif_var_ae[[1]], legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-0.3,0.3),
     legend.width=3.3, legend.shrink=0.75,
     axis.args=list(at=seq(-0.3, 0.3, 0.05),
                    lab=c('-0.3','','-0.2','','-0.1','','0','','0.1','','0.2','','0.3'),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=2.5))
dev.off()
