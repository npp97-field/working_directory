#### the snow thaw and SOS, snow fall and eos
library(raster)
library(rgdal)
library(ncdf4)
setwd("/Users/yzhang/Project/SIF_phenology/")

load("./analysis/snow_data.RData")

### coastlines
coastline<-shapefile("/Users/yzhang/Data/GIS_data/global/ne_110m_coastline.shp")
cropcoast<-crop(coastline,extent(-180,180,30,90))
repcoa<-spTransform(cropcoast,ae)

setwd("/Users/yzhang/Project/SIF_phenology/analysis/clear_daily_phenology/")
ncin<-nc_open(paste("./clear_daily_SOS_30N_fixed_stat.nc",sep=""))
sossif<-ncvar_get(nc = ncin,varid = "MEAN")
nc_close(ncin)

aesos<-nc2ae(sossif)
aesnow<-nc2ae(average_snow_end)
pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/pheno_map_fix_",indi,".pdf",sep=""),width=11/3*2.4,height=11/3)
par(fig=c(0,0.45,0,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
image(setrange(aesos+(14.5-6*k)/365,0.2192,0.5479),add=T,col=lst_color_ramp,axes=F,zlim=c(0.2192,0.5479))
text(0, 5500000,"SOS",cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0,0.48,0,1),new=T)
plot(aesos, legend.only=TRUE, col=lst_color_ramp,horizontal=F,zlim=c(0.2191,0.5480),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(80, 200, 20)/365,
                    labels=seq(80, 200, 20), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

par(fig=c(0.5,0.95,0,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aesnow/365,0.2192,0.5479),add=T,col=lst_color_ramp,axes=F,zlim=c(0.2192,0.5479))
text(0, 5500000,"Snow thaw",cex=1.3)
plotlatlong()
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.98,0,1),new=T)
plot(aesos, legend.only=TRUE, col=lst_color_ramp,horizontal=F,zlim=c(0.2191,0.5480),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(80, 200, 20)/365,
                    labels=seq(80, 200, 20), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))



