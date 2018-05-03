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
 
  for (i in 3:8){
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
ncin<-nc_open("./SOS_30N_fixed_stat.nc")
sossif<-ncvar_get(nc = ncin,varid = "MEAN")
sostrend<-ncvar_get(nc = ncin,varid = "TREND")
nc_close(ncin)

ncin<-nc_open("./EOS_30N_fixed_stat.nc")
eossif<-ncvar_get(nc = ncin,varid = "MEAN")
eostrend<-ncvar_get(nc = ncin,varid = "TREND")
nc_close(ncin)

ncin<-nc_open("./LGS_30N_fixed_stat.nc")
lgssif<-ncvar_get(nc = ncin,varid = "MEAN")
lgstrend<-ncvar_get(nc = ncin,varid = "TREND")
nc_close(ncin)


aesos<-nc2ae(sossif)
aesos_t<-nc2ae(sostrend)
aeeos<-nc2ae(eossif)
aeeos_t<-nc2ae(eostrend)
aelgs<-nc2ae(lgssif)
aelgs_t<-nc2ae(lgstrend)


pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/pheno_map.pdf",width=11/3*2.4,height=11)
par(fig=c(0,0.45,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
image(setrange(aesos+2.5/365,0.2192,0.5479),add=T,col=lst_color_ramp,axes=F,zlim=c(0.2192,0.5479))
text(0, 5500000,"SOS",cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0,0.48,0.66,1),new=T)
plot(aesos, legend.only=TRUE, col=lst_color_ramp,horizontal=F,zlim=c(0.2191,0.5480),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(80, 200, 20)/365,
                    labels=seq(80, 200, 20), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))


par(fig=c(0,0.45,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aeeos+2.5/365,0.5479,0.8768),add=T,col=rev(lst_color_ramp),axes=F,zlim=c(0.5479,0.8768))
text(0, 5500000,"EOS",cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0,0.48,0.33,0.66),new=T)
plot(aesos, legend.only=TRUE, col=rev(lst_color_ramp),horizontal=F,zlim=c(0.5479,0.8768),
     legend.width=1.3, legend.shrink=0.75,label="SOS DOY",
     axis.args=list(at=seq(200, 320, 20)/365,
                    labels=seq(200, 320, 20), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))#,
     #legend.args=list(text='DOY', side=4, font=1, line=2, cex=1.1))


par(fig=c(0,0.45,0.00,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aelgs,0.2191,0.5480),add=T,col=rev(lst_color_ramp),axes=F,zlim=c(0.2191,0.5480))
text(0, 5500000,"LGS",cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"c",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0,0.48,0,0.33),new=T)
plot(aesos, legend.only=TRUE, col=rev(lst_color_ramp),horizontal=F,zlim=c(0.2191,0.5480),
     legend.width=1.3, legend.shrink=0.75,label="SOS DOY",
     axis.args=list(at=seq(80, 200, 20)/365,
                    labels=seq(80, 200, 20), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))#,
#legend.args=list(text='DOY', side=4, font=1, line=2, cex=1.1))


#############################
par(fig=c(0.5,0.95,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aesos_t,-0.004109589,0.004109589),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-0.004109589,0.004109589))
text(0, 5500000,"SOS",cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"d",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.98,0.66,1),new=T)
plot(aesos, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-0.004109589,0.004109589),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1.5, 1.5, 0.25)/365,
                    labels=c('-1.5','','-1','','-0.5','','0','','0.5','','1','','1.5'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))



par(fig=c(0.5,0.95,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aeeos_t,-0.004109589,0.004109589),add=T,col=ano_discrete_ramp,
      axes=F,zlim=c(-0.004109589,0.004109589))
text(0, 5500000,"EOS",cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"e",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.98,0.33,0.66),new=T)
plot(aesos, legend.only=TRUE, col=ano_discrete_ramp_leg,horizontal=F,zlim=c(-0.004109589,0.004109589),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1.5, 1.5, 0.25)/365,
                    labels=c('-1.5','','-1','','-0.5','','0','','0.5','','1','','1.5'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))


par(fig=c(0.5,0.95,0.00,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aelgs_t,-0.00822,0.00822),add=T,col=ano_discrete_ramp,axes=F,zlim=c(-0.00822,0.00822))
text(0, 5500000,"LGS",cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"f",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0.5,0.98,0,0.33),new=T)
plot(aelgs_t, legend.only=TRUE, 
     col=ano_discrete_ramp_leg,horizontal=F,zlim=c(-0.00822,0.00822),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-3, 3, 0.5)/365,
                    labels=c('-3','','-2','','-1','','0','','1','','2','','3'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))
#box()
#axis(1)
par(fig=c(0,1,0,1),new=T)
plot(NA,axes=F,xlim=c(-1,1),ylim=c(-1.5,1.5),xaxs="i",yaxs='i')
text(-0.03, 1.02, labels = 'DOY', xpd = T, srt = -90,cex=1.1)     
text(-0.03, 0, labels = 'DOY', xpd = NA, srt = -90,cex=1.1)
text(-0.03,-1.03, labels = 'DOY', xpd = NA, srt = -90,cex=1.1)
text(0.97, 1.02, labels = expression('Day y'^-1), xpd = NA, srt = -90,cex=1.1) 
text(0.97, 0, labels = expression('Day y'^-1), xpd = NA, srt = -90,cex=1.1)    
text(0.97, -1.03, labels = expression('Day y'^-1), xpd = NA, srt = -90,cex=1.1)
#text(12000000, 0, labels = expression('Day y'^-1), xpd = NA, srt = -90,cex=1.1,col='red')    


dev.off()

