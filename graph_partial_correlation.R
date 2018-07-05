##### graph_geographical_data
# plot girded phenology data
# SOS, EOS, LGS, average and trend

source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")

lst_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_color.csv",header = F)/255))
ano_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_ano.csv",header = F)/255))
ano_discrete_ramp<-c(ano_color_ramp[c(1,26,51,76,101)],"grey80","grey80",ano_color_ramp[c(140,165,190,215,240)])
ano_discrete_ramp_leg<-rep(ano_discrete_ramp,each=20)

cropland<-shapefile("/Users/yzhang/Project/SIF_phenology/data/filtered_result.shp")
projection(cropland)<-longlat
repcrop<-spTransform(cropland,ae)


setwd("/Users/yzhang/Project/SIF_phenology/analysis/")
# ncin<-nc_open(paste("./correlation_clear/pcor_sos_csif_p2s.nc",sep=""))
# sos_p2s<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef"))
# nc_close(ncin)
# 
# ncin<-nc_open(paste("./correlation_clear/pcor_sos_eos.nc",sep=""))
# sos_eos<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef"))
# nc_close(ncin)
# 
# ncin<-nc_open(paste("./correlation_clear/pcor_sos_csif_s2p.nc",sep=""))
# sos_s2p<-nc2ae(ncvar_get(nc = ncin,varid = "pcor_coef"))
# nc_close(ncin)


ncin<-nc_open(paste("./correlation_clear/cor_sos_csif_p2s.nc",sep=""))
sos_p2s<-nc2ae(ncvar_get(nc = ncin,varid = "cor_coef"))
nc_close(ncin)

ncin<-nc_open(paste("./correlation_clear/cor_sos_eos.nc",sep=""))
sos_eos<-nc2ae(ncvar_get(nc = ncin,varid = "cor_coef"))
nc_close(ncin)

ncin<-nc_open(paste("./correlation_clear/cor_sos_csif_s2p.nc",sep=""))
sos_s2p<-nc2ae(ncvar_get(nc = ncin,varid = "cor_coef"))
nc_close(ncin)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/cor_sos_eos_sif_clear.pdf",sep=""),width=11/3*1.2,height=11)

#############################
par(fig=c(0,0.9,2/3,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
image(setrange(sos_p2s,-1,1),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-1,1))
text(0, 5500000,expression("SOS~SIF"[Peak_to_sene]),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.66,1),new=T)
plot(sos_p2s, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1, 1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))



par(fig=c(0,0.9,1/3,2/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_eos,-1,1),add=T,col=rev(ano_discrete_ramp),
      axes=F,zlim=c(-1,1))
text(0, 5500000,expression("SOS~EOS"),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"b",cex=1.8,font=2,padj=-7,las=2)
par(fig=c(0.5,0.96,0.33,0.66),new=T)
plot(sos_eos, legend.only=TRUE, col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1,1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))


par(fig=c(0,0.9,0.00,1/3),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(sos_s2p,-1,1),add=T,col=rev(ano_discrete_ramp),axes=F,zlim=c(-1,1))
text(0, 5500000,expression("SOS~SIF"[Start_to_peak]),cex=1.3)
plotlatlong()
plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
mtext(side=2,line=-1.5,"c",cex=1.8,font=2,padj=-7,las=2)

par(fig=c(0.5,0.96,0,0.33),new=T)
plot(sos_s2p, legend.only=TRUE, 
     col=rev(ano_discrete_ramp_leg),horizontal=F,zlim=c(-1,1),
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(at=seq(-1,1, 1/6),
                    c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'), mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))
#box()
#axis(1)
par(fig=c(0,1,0,1),new=T)
plot(NA,axes=F,xlim=c(-1,1),ylim=c(-1.5,1.5),xaxs="i",yaxs='i')

# text(0.94, 1.02, labels = expression('partial r'), xpd = NA, srt = -90,cex=1.1) 
# text(0.94, 0, labels = expression('partial r'), xpd = NA, srt = -90,cex=1.1)    
# text(0.94, -1.03, labels = expression('partial r'), xpd = NA, srt = -90,cex=1.1)
#text(12000000, 0, labels = expression('Day y'^-1), xpd = NA, srt = -90,cex=1.1,col='red')

dev.off()

