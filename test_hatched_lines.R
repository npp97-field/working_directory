library(maptools)
library(raster)
library(rgdal)
library(ncdf4)

longlat <-  CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
ae<-"+proj=aeqd +lat_0=90 +lon_0=-0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

coastline<-shapefile("/Users/yzhang/Data/GIS_data/global/ne_110m_coastline.shp")
cropcoast<-crop(coastline,extent(-180,180,30,90))
repcoa<-spTransform(cropcoast,ae)

# landcover_f<-nc_open("/Users/yzhang/Project/SIF_phenology/data/North_mcd12c1_landcover1_majority.nc")
# landcover<-ncvar_get(landcover_f,'lc1')
# lc<-raster("/Users/yzhang/Project/SIF_phenology/data/north_landcover.tif")
# crop = lc==12|lc==14
# writeRaster(crop,"/Users/yzhang/Project/SIF_phenology/data/north_crop.tif")
cropland<-shapefile("/Users/yzhang/Project/SIF_phenology/data/filtered_result.shp")
projection(cropland)<-longlat
repcrop<-spTransform(cropland,ae)
plot(repcrop[c(1:288,290:291),],density=110,border=NA)
plot(repcoa,add=T)
