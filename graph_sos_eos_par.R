#### graph correlation between three climate variables

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
ncin<-nc_open(paste("./correlation_clear/cor_eos_p2e_par.nc",sep=""))
#ncin<-nc_open(paste("./climate_intercor/cor_pre_end_par_pre_prec.nc",sep=""))
#ncin<-nc_open(paste("./climate_intercor/cor_pre_start_par_pre_prec.nc",sep=""))
eos_par<-nc2ae(ncvar_get(nc = ncin,varid = "cor_coef"))
nc_close(ncin)

ncin<-nc_open(paste("./correlation_clear/cor_sos_p2e_par.nc",sep=""))
#ncin<-nc_open(paste("./climate_intercor/cor_pre_end_par_pre_temp.nc",sep=""))
#ncin<-nc_open(paste("./climate_intercor/cor_pre_start_par_pre_temp.nc",sep=""))
sos_par<-nc2ae(ncvar_get(nc = ncin,varid = "cor_coef"))
nc_close(ncin)

ncin<-nc_open(paste("./correlation_clear/cor_p2e_sif_par.nc",sep=""))
#ncin<-nc_open(paste("./climate_intercor/cor_pre_end_prec_pre_temp.nc",sep=""))
#ncin<-nc_open(paste("./climate_intercor/cor_pre_start_prec_pre_temp.nc",sep=""))
sif_par<-nc2ae(ncvar_get(nc = ncin,varid = "cor_coef"))
nc_close(ncin)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/sos_par_eos.pdf",sep=""),width=11/3*1.2,height=11)

#############################
par(fig=c(0,0.9,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
image(setrange(sos_par,-1,1),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-1,1))
text(0, 5500000,expression("p2e_PAR~SOS"),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.66,1),new=T)
plot(sos_par, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1, 1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))



par(fig=c(0,0.9,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(eos_par,-1,1),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-1,1))
text(0, 5500000,expression("p2e_PAR~EOS"),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.33,0.66),new=T)
plot(eos_par, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1,1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))


par(fig=c(0,0.9,0.00,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sif_par,-1,1),add=T,col=rev(ano_discrete_ramp),axes=F,zlim=c(-1,1))
text(0, 5500000,expression("p2e_PAR~p2e_SIF"),cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"c",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0,0.33),new=T)
plot(sif_par, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1,1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

dev.off()

