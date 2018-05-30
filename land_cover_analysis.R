###### get the land cover data
library(raster)
library(ncdf4)
library(plotKML)
# ncin<-nc_open('/Users/yzhang/Project/SIF_phenology/data/mcd12c1_landcover1_majority.nc')
# lccmg<-ncvar_get(ncin,varid='landcover1')
# ras_lc <- raster(lccmg)
# 
# ## aggregate by a factor of 3, with "modal"
# m <- aggregate(ras_lc, fact = 10, fun = modal, na.rm = TRUE)
# m_mat<-getValues(m)
# dim(m_mat)<-c(360,720)
# aoi<-t(m_mat[241:360,])
# outfile<-"/Users/yzhang/Project/SIF_phenology/data/North_mcd12c1_landcover1_majority.nc"
# export_nc(aoi,outfile,'lc1')
# 
export_nc<-function(dat,outfile,varname){
  latmin<- 30
  latmax<- 90
  latd<- 0.5
  lonmin<- -180
  lonmax<- 180
  lond<- 0.5

  lat<- seq(latmin+latd/2,latmax-latd/2,latd)
  long<-seq(lonmin+lond/2,lonmax-lond/2,lond)

  dimlat<-ncdim_def('latitude','deg',lat)
  dimlong<-ncdim_def('longitude','deg',long)
  ncvar<-ncvar_def(varname,'NA',list(dimlong,dimlat),-9999,longname=varname,prec='float',compression=9)

  if (file.exists(outfile)){
    file.remove(outfile)
  }
  ncout<-nc_create(outfile,list(ncvar))
  ncvar_put(ncout,varid=ncvar,dat)
  nc_close(ncout)
}

to_raster<-function(nc){
  ras<-raster(apply(nc,1,rev))
  extent(ras)<-c(-180,180,30,90)
  return(ras)
}

wet_r<-to_raster(wet)
plot(wet_r)
plot(coastline,add=T)
barren_r<-to_raster(barren)
plot(barren_r)



setwd("/Users/yzhang/Project/SIF_phenology/")

ncin<-nc_open("./data/North_mcd12c1_landcover1_majority.nc")
lc<-ncvar_get(ncin,"lc1")
nc_close(ncin)
#### 0 wat 1 ENF 2 EBF 3 DNF 4 DBF 5 MF 
#### 6 CSH 7 OSH 8 WSA
#### 9 SAV 10 GRA 11 WET
#### 12 CRO 14 CNV

lc_col<-worldgrids_pal$IGBP[2:15]
lc_type<-c("ENF","EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO","URB","CNV")

barren<-!(lc==0|lc==15|lc==16)
barren[barren==0]<-NA
export_nc(barren,"./data/North_barren_mask.nc","barren")

lc<-lc*barren

ncin<-nc_open("./data/ai_cru.nc")
ai<-ncvar_get(ncin,'ai')/1000
nc_close(ncin)

ncin<-nc_open("./analysis/p_corr/pcor_sos_csif_p2s.nc")
pcor_sos_p2s<-ncvar_get(ncin,"pcor_coef")
nc_close(ncin)

ncin<-nc_open("./analysis/p_corr/pcor_sos_eos.nc")
pcor_sos_eos<-ncvar_get(ncin,"pcor_coef")
nc_close(ncin)

cor_lc_ai<-cbind(as.vector(pcor_sos_p2s),as.vector(pcor_sos_eos),as.vector(lc),as.vector(ai))
good_dat<-cor_lc_ai[complete.cases(cor_lc_ai), ]

#plot(good_dat[,1],log(good_dat[,4]),col=lc_col[good_dat[,3]],cex=0.3,pch=15)
#tapply(good_dat[,1], good_dat[,3], mean)

#######################################

longlat <-  CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
ae<-"+proj=aeqd +lat_0=90 +lon_0=-0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline<-shapefile("/Users/yzhang/Data/GIS_data/global/ne_110m_coastline.shp")
cropcoast<-crop(coastline,extent(-180,180,30,90))
repcoa<-spTransform(cropcoast,ae)



###### plot functions
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
  aedat<-projectRaster(latlongdat,extae,method = "ngb")
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
border<-plot_lat(30)
lc_rep<-nc2ae(lc)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/landcover.pdf",sep=""),width=5,height=4)

par(fig=c(0,1,0,1),mar=c(0.4,0.4,0.4,4.4),mgp=c(3,0.3,0))
plot(border)
image(lc_rep,add=T,col=lc_col,axes=F,zlim=c(1,14))
plot(repcoa,add=T)
#text(0, 5500000,"SOS",cex=1.3)
plotlatlong()
#mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0,1,0,1),mar=c(0.4,0.4,0.4,0.4),new=T)
plot(lc_rep, legend.only=TRUE, col=rep(rev(lc_col),each=20),horizontal=F,zlim=c(1,14*20),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=(1:14)*20-10,
                    labels=rev(lc_type), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

dev.off()


# ######################
# ### only plot the forest ecosystems
# forests<-good_dat[good_dat[,3]==1|good_dat[,3]==2|good_dat[,3]==3|good_dat[,3]==4|good_dat[,3]==5,]
# plot(forests[,1],log(forests[,4]),col=lc_col[forests[,3]],cex=0.3,pch=15)
# 
# non_forest<-good_dat[good_dat[,3]==6|good_dat[,3]==7|good_dat[,3]==8|good_dat[,3]==9|good_dat[,3]==10|good_dat[,3]==11,]
# plot(non_forest[,1],log(non_forest[,4]),col=lc_col[non_forest[,3]],cex=0.3,pch=15)
# 
# crop<-good_dat[good_dat[,3]==12|good_dat[,3]==14,]
# plot(crop[,1],log(crop[,4]),col=lc_col[crop[,3]],cex=0.3,pch=15)