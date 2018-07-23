###### graph figure 5 for the future change of the temperature limited and precipitation limited EOP
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")

setwd("/Users/yzhang/Project/SIF_phenology/")

ncin<-nc_open("./analysis/future_t.nc")
currentT<-ncvar_get(ncin,varid="current_tas")
futureT<-ncvar_get(ncin,varid="rcp85_tas")
nc_close(ncin)
ncin<-nc_open("./analysis/future_AI.nc")
currentAI<-ncvar_get(ncin,varid="currentai")
futureAI<-ncvar_get(ncin,varid="rcp85ai")
nc_close(ncin)

######### start the graph here.

CTae<-nc2ae(currentT)
FTae<-nc2ae(futureT)
CAIae<-nc2ae(currentAI)
FAIae<-nc2ae(futureAI)

get_diff<-function(val1,val2){
  out<-(val1+1)*(val2+3)
  lim<-out
  lim[lim==3]<-1
  lim[lim==6]<-2
  lim[lim==4]<-3
  lim[lim==8]<-4
  return(lim)
}

t_pred<-get_diff(CTae>2.5,FTae>2.5)
ai_pred<-get_diff(CAIae<0.6,FAIae<0.6)


dis_col_gmt<-rep(ano_col_gmt[c(2,5,11,14)],each=15)
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Figure_future_AI_TEMP.pdf",height=4.5,width=10)
par(oma=c(0,0,0,8),mgp=c(3,0.3,0))
par(fig=c(0,0.5,0,1),mar=c(0.4,0.4,0.4,0.4))
plot(border)
image(setrange(ai_pred,1,4),add=T,col=rev(dis_col_gmt),
      axes=F,zlim=c(1,4))
text(0, 5500000,"Aridity Index",cex=1)
plotlatlong()
mtext(side=2,line=-1.5,'a',cex=1.3,font=2,padj=-9.5,las=2)

par(fig=c(0.5,1,0,1),mar=c(0.4,0.4,0.4,0.4),new=T)
plot(border)
image(setrange(t_pred,1,4),add=T,col=rev(dis_col_gmt),
      axes=F,zlim=c(1,4))
text(0, 5500000,"Mean annual Temperature",cex=1)
text(11300000,0, labels = "Change of limitation on EOP", xpd = NA, srt = -90,cex=1.1)   
plotlatlong()
mtext(side=2,line=-1.5,'b',cex=1.3,font=2,padj=-9.5,las=2)


par(fig=c(0.5,0.9,0.1,0.9),oma=c(0,0,0,0),new=T)
plot(ai_pred, legend.only=TRUE, col=rev(dis_col_gmt),horizontal=F,zlim=c(0.5,4.5),
     legend.width=1.5, legend.shrink=0.75,
     axis.args=list(at=1:4,label=c("T",expression("P"%->%"T"),expression("T"%->%"P"),"P"),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))

dev.off()

