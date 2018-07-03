##### graph_geographical_data
# plot girded phenology data
# SOS, EOS, LGS, average and trend

library(ncdf4)
library(maptools)
library(raster)
library(rgdal)

lst_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_color.csv",header = F)/255))
ano_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_ano.csv",header = F)/255))
#ano_discrete_ramp<-c(ano_color_ramp[c(1,26,51,76,101)],"grey80","grey80",ano_color_ramp[c(140,165,190,215,240)])
#ano_discrete_ramp_leg<-rep(ano_discrete_ramp,each=20)
jet.colorslow <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(10)
jet.colors<-rep(jet.colorslow[6:10],each=10)


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


cropland<-shapefile("/Users/yzhang/Project/SIF_phenology/data/crop_north.shp")
projection(cropland)<-longlat
repcrop<-spTransform(cropland,ae)


setwd("/Users/yzhang/Project/SIF_phenology/analysis/")
ncin<-nc_open(paste("./multivariate_regression/dsos_with_trend.nc",sep=""))
sos_reduce<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_Rsquare"))
nc_close(ncin)

ncin<-nc_open(paste("./multivariate_regression/dsos_with_eos.nc",sep=""))
sos_full<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_Rsquare"))
nc_close(ncin)

ncin<-nc_open(paste("./multivariate_regression/dsos_diff.nc",sep=""))
sos_diff<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef"))
nc_close(ncin)


ncin<-nc_open(paste("./multivariate_regression/deos_with_trend.nc",sep=""))
eos_reduce<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_Rsquare"))
nc_close(ncin)

ncin<-nc_open(paste("./multivariate_regression/deos_with_sos.nc",sep=""))
eos_full<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_Rsquare"))
nc_close(ncin)

ncin<-nc_open(paste("./multivariate_regression/deos_diff.nc",sep=""))
eos_diff<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef"))
nc_close(ncin)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/dmultivariate_regression.pdf",sep=""),width=11/3*2.2,height=11)

#############################
par(fig=c(0,0.46,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
image(setrange(sos_reduce,0,1),add=T,col=jet.colors,
      axes=F,zlim=c(0,1))
text(0, 5500000,expression("SOS~T,P,R"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)



par(fig=c(0,0.46,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_full,0,1),add=T,col=jet.colors,
      axes=F,zlim=c(0,1))
text(0, 5500000,expression("SOS~T,P,R,EOS_pre"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)


par(fig=c(0,0.46,0.00,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_diff,0,1),add=T,col=jet.colors,axes=F,zlim=c(0,1))
text(0, 5500000,expression("R"^2~"difference"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"c",cex=1.8,font=2,padj=-7,las=2)


######################
########################

par(fig=c(0.46,0.92,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(eos_reduce,0,1),add=T,col=jet.colors,
      axes=F,zlim=c(0,1))
text(0, 5500000,expression("EOS~T,P,R"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"d",cex=1.8,font=2,padj=-7,las=2)


par(fig=c(0.46,0.92,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(eos_full,0,1),add=T,col=jet.colors,
      axes=F,zlim=c(0,1))
text(0, 5500000,expression("EOS~T,P,R,SOS"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"e",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0.46,0.92,0.00,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(eos_diff,0,1),add=T,col=jet.colors,axes=F,zlim=c(0,1))
text(0, 5500000,expression("R"^2~"difference"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"f",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0.5,0.98,0.33,0.66),new=T)
plot(eos_diff, legend.only=TRUE, 
     col=jet.colors,horizontal=F,zlim=c(0,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(0,1, 1/5),
                    seq(0,1, 1/5), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

dev.off()