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


setwd("/Users/yzhang/Project/SIF_phenology/")
# #ncin<-nc_open(paste("./climate_variability/clear_daily_SOS_30N_pre_start_par_stat.nc",sep=""))
# ncin<-nc_open(paste("./data/monthly_temp.nc",sep=""))
# sd_var<-ncvar_get(nc = ncin,varid = "sd_temp")
# nc_close(ncin)
# 
# pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/month_temp_var.pdf",sep=""),width=11/3*4.3,height=11)
# mon<-c("Jan",'Feb',"Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# #############################
# for (i in 1:12){
#   row = ceiling(i/4)
#   col = i-row*4+4
#   par(fig=c((col-1)/4,col/4,1-row/3,1-row/3+1/3),mar=c(0.4,0.4,0.4,0.4),oma=c(0,0,0,4),mgp=c(3,0.3,0),new=T)
#   plot(border)
#   # image(setrange(par_var,0,10),add=T,col=rev(lst_color_ramp),
#   #       axes=F,zlim=c(0,10))
#   image(setrange(nc2ae(sd_var[,,i]),0,4),add=T,col=rev(lst_color_ramp),
#         axes=F,zlim=c(0,4))
#   text(0, 5500000,mon[i],cex=1.3)
#   plotlatlong()
# 
# }
# par(fig=c(0.5,0.97,0,1),mar=c(0,0,0,0),oma=c(0,0,0,0),mgp=c(3,0.3,0))
# plot(par_var, legend.only=TRUE, col=rev(lst_color_ramp),horizontal=F,zlim=c(0,4),#zlim=c(0,10),#for std
#      legend.width=1.3, legend.shrink=0.75,
#      axis.args=list(#at=seq(-1, 1, 1/6),
#        #c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'),
#        mgp=c(3,0.2,0),tck=0.3,
#        cex.axis=1))
# 
# dev.off()
# 
# 



#ncin<-nc_open(paste("./climate_variability/clear_daily_SOS_30N_pre_start_par_stat.nc",sep=""))
ncin<-nc_open(paste("./analysis/all_daily_SOS_30N_fixed_stat.nc",sep=""))
sd_SOS<-ncvar_get(nc = ncin,varid = "SD")
nc_close(ncin)

ncin<-nc_open(paste("./analysis/all_daily_EOS_30N_fixed_stat.nc",sep=""))
sd_EOS<-ncvar_get(nc = ncin,varid = "SD")
nc_close(ncin)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/SOS_EOS_var.pdf",sep=""),width=11/3*2.1,height=11/3)
#############################


par(fig=c(0,0.5,0,1),mar=c(0.4,0.4,0.4,0.4),oma=c(0,0,0,4),mgp=c(3,0.3,0))
plot(border)
image(setrange(nc2ae(sd_SOS)*365,0,15),add=T,col=rev(lst_color_ramp),
      axes=F,zlim=c(0,15))
text(0, 5500000,"SD SOS",cex=1.3)
plotlatlong()




par(fig=c(0.5,1,0,1),mar=c(0.4,0.4,0.4,0.4),oma=c(0,0,0,4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(nc2ae(sd_EOS)*365,0,15),add=T,col=rev(lst_color_ramp),
      axes=F,zlim=c(0,15))
text(0, 5500000,"SD EOS",cex=1.3)
plotlatlong()


par(fig=c(0.5,0.97,0,1),mar=c(0,0,0,0),oma=c(0,0,0,0),mgp=c(3,0.3,0))
plot(par_var, legend.only=TRUE, col=rev(lst_color_ramp),horizontal=F,zlim=c(0,15),#zlim=c(0,10),#for std
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(#at=seq(-1, 1, 1/6),
       #c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'),
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=1))

dev.off()



ncin<-nc_open(paste("./analysis/all_daily_SOS_30N_fixed_stat.nc",sep=""))
sd_SOS<-ncvar_get(nc = ncin,varid = "SD")
nc_close(ncin)

ncin<-nc_open(paste("./analysis/all_daily_EOS_30N_fixed_stat.nc",sep=""))
sd_EOS<-ncvar_get(nc = ncin,varid = "SD")
nc_close(ncin)
