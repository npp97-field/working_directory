#### freeze thaw analysis

library(raster)
library(ncdf4)
library()

freeze<-as.matrix(read.table("/Users/yzhang/Data/freeze_thaw/freez_a.dat"))
coastline<-shapefile("/Users/yzhang/Data/GIS_data/global/ne_110m_coastline.shp")

freeze_raster<-raster(freeze)
extent(freeze_raster)<-c(-180,180,-90,90)
image(freeze_raster)
plot(coastline,add=T)

thaw<-as.matrix(read.table("/Users/yzhang/Data/freeze_thaw/thaw_a.dat"))
thaw_raster<-raster(thaw)
extent(thaw_raster)<-c(-180,180,-90,90)
image(thaw_raster,zlim=c(0,6000))
plot(coastline,add=T)

