##### graph_geographical_data
# plot girded phenology data
# SOS, EOS, LGS, average and trend
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")

k=3
indicator<-c("VI",'all_daily',"clear_daily")
indi<-indicator[k]

border<-plot_lat(30)
#lat60<-plot_lat(60)

### coastlines
coastline<-shapefile("/Users/yzhang/Data/GIS_data/global/ne_110m_coastline.shp")
cropcoast<-crop(coastline,extent(-180,180,30,90))
repcoa<-spTransform(cropcoast,ae)

setwd("/Users/yzhang/Project/SIF_phenology/analysis/clear_daily_phenology/")
ncin<-nc_open(paste("./",indi,"_EOS_30N_fixed_stat.nc",sep=""))
eossif<-ncvar_get(nc = ncin,varid = "MEAN")
eostrend<-ncvar_get(nc = ncin,varid = "TREND")
nc_close(ncin)

aeeos<-nc2ae(eossif)
aeeos_t<-nc2ae(eostrend)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/pheno_map_fix_",indi,".pdf",sep=""),width=11/3*2.4,height=11/3)

par(fig=c(0,0.45,0,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),oma=c(0,0,0,1))
plot(border)
image(setrange(aeeos+(14.5-6*k)/365,0.5479,0.8768),add=T,col=rev(lst_color_ramp),axes=F,zlim=c(0.5479,0.8768))
text(0, 5500000,"EOP",cex=1.3)
text(9700000, 0, labels = 'DOY', xpd = NA, srt = -90,cex=1.1)   
plotlatlong()
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0,0.48,0,1),new=T)
plot(aeeos, legend.only=TRUE, col=rev(lst_color_ramp),horizontal=F,zlim=c(0.5479,0.8768),
     legend.width=1.3, legend.shrink=0.75,label="SOS DOY",
     axis.args=list(at=seq(200, 320, 20)/365,
                    labels=seq(200, 320, 20), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))#,

par(fig=c(0.5,0.95,0,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(aeeos_t,-0.004109589,0.004109589),add=T,col=ano_discrete_ramp,
      axes=F,zlim=c(-0.004109589,0.004109589))
text(0, 5500000,"EOP",cex=1.3)
text(9700000, 0, labels = expression('Day y'^-1), xpd = NA, srt = -90,cex=1.1) 
plotlatlong()
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.98,0,1),new=T)
plot(aeeos, legend.only=TRUE, col=ano_discrete_ramp_leg,horizontal=F,zlim=c(-0.004109589,0.004109589),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1.5, 1.5, 0.25)/365,
                    labels=c('-1.5','','-1','','-0.5','','0','','0.5','','1','','1.5'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

dev.off()

