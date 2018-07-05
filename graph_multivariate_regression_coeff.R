##### graph_geographical_data
# plot girded phenology data
# SOS, EOS, LGS, average and trend
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")

lst_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_color.csv",header = F)/255))
ano_discrete_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_ano.csv",header = F)/255))
#ano_discrete_ramp<-c(jet.colors[c(1,26,51,76,101)],"grey80","grey80",jet.colors[c(140,165,190,215,240)])
#ano_discrete_ramp_leg<-rep(ano_discrete_ramp,each=20)
jet.colorslow <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100)
jet.colors<-jet.colorslow[51:100]

cropland<-shapefile("/Users/yzhang/Project/SIF_phenology/data/crop_north.shp")
projection(cropland)<-longlat
repcrop<-spTransform(cropland,ae)


setwd("/Users/yzhang/Project/SIF_phenology/analysis/")

# ## for SOS
# ncin<-nc_open(paste("./multivariate_regression/dsos_with_trend.nc",sep=""))
# dsos_reduce_temp<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef")[,,2])
# nc_close(ncin)
# 
# ncin<-nc_open(paste("./multivariate_regression/dsos_with_eos.nc",sep=""))
# dsos_full_temp<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef")[,,2])
# nc_close(ncin)
# 
# ncin<-nc_open(paste("./multivariate_regression/sos_with_trend.nc",sep=""))
# sos_reduce_temp<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef")[,,2])
# nc_close(ncin)
# 
# ncin<-nc_open(paste("./multivariate_regression/sos_with_eos.nc",sep=""))
# sos_full_temp<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef")[,,2])
# nc_close(ncin)


##for eos
ncin<-nc_open(paste("./multivariate_regression/deos_with_trend.nc",sep=""))
dsos_reduce_temp<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef")[,,2])
nc_close(ncin)

ncin<-nc_open(paste("./multivariate_regression/deos_with_sos.nc",sep=""))
dsos_full_temp<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef")[,,2])
nc_close(ncin)

ncin<-nc_open(paste("./multivariate_regression/eos_with_trend.nc",sep=""))
sos_reduce_temp<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef")[,,2])
nc_close(ncin)

ncin<-nc_open(paste("./multivariate_regression/eos_with_sos.nc",sep=""))
sos_full_temp<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef")[,,2])
nc_close(ncin)


pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/dmultivariate_regression_coef_eos.pdf",sep=""),width=11/3*2.2,height=11)

#############################
par(fig=c(0,0.46,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
image(setrange(sos_reduce_temp,0,3),add=T,col=jet.colors,
      axes=F,zlim=c(0,3))
text(0, 5500000,expression(gamma["T"]^"EOS"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0,0.46,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(dsos_reduce_temp,0,3),add=T,col=jet.colors,
                  axes=F,zlim=c(0,3))
text(0, 5500000,expression(gamma["dT"]^"dEOS"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0,0.46,0,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_reduce_temp-dsos_reduce_temp,-0.8,0.8),add=T,col=ano_color_ramp,
      axes=F,zlim=c(-0.8,0.8))
text(0, 5500000,expression(paste(Delta,gamma["T"]^"EOS"),sep=""),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"c",cex=1.8,font=2,padj=-7,las=2)

######################
########################

par(fig=c(0.46,0.92,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_full_temp,0,3),add=T,col=jet.colors,
      axes=F,zlim=c(0,3))
text(0, 5500000,expression(gamma["T"]^"EOS,SOS"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"d",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0.46,0.92,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(dsos_full_temp,0,3),add=T,col=jet.colors,
      axes=F,zlim=c(0,3))
text(0, 5500000,expression(gamma["dT"]^"dEOS,SOS"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"e",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0.46,0.92,0,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_full_temp-dsos_full_temp,-0.8,0.8),add=T,col=ano_color_ramp,
      axes=F,zlim=c(-0.8,0.8))
text(0, 5500000,expression(paste(Delta,gamma["T"]^"EOS,SOS"),sep=""),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"f",cex=1.8,font=2,padj=-7,las=2)



par(fig=c(0.5,0.98,0.5,0.83),new=T)
plot(sos_full_temp, legend.only=TRUE, 
     col=jet.colors,horizontal=F,zlim=c(0,3),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(0,3, 1),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

par(fig=c(0.5,0.98,0,0.33),new=T)
plot(sos_full_temp, legend.only=TRUE, 
     col=ano_discrete_ramp,horizontal=F,zlim=c(-0.8,0.8),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-0.8,0.8,0.2),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))


dev.off()
