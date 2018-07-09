#### graph and analyze the EOS climate relationship with aridity index.
## 
library("LSD")
## use the heatscatter function
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")
setwd("/Users/yzhang/Project/SIF_phenology/")
#ncin<-nc_open("./data/ai_cru_P_PET.nc")
ncin<-nc_open("./data/ai_cru_P_PET_grow.nc")
ai_data<-ncvar_get(ncin,"ai")
nc_close(ncin)

ncin<-nc_open("./data/North_mcd12c1_landcover1_majority.nc")
lc_data<-ncvar_get(ncin,"lc1")
nc_close(ncin)

ncin<-nc_open("./analysis/clear_daily_phenology/clear_daily_EOS_30N_fixed_stat.nc")
eos_data<-ncvar_get(ncin,'MEAN')*365
nc_close(ncin)

ncin<-nc_open("./analysis/correlation_clear_rs/max_correlation_eos_tday.nc")
tday<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
ncin<-nc_open("./analysis/correlation_clear_rs/max_significance_eos_tday.nc")
tday_sig<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
tday[tday_sig>0.05]<-NA

ncin<-nc_open("./analysis/correlation_clear_rs/max_correlation_eos_prec.nc")
prec<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
ncin<-nc_open("./analysis/correlation_clear_rs/max_significance_eos_prec.nc")
prec_sig<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
prec[prec_sig>0.05]<-NA

ncin<-nc_open("./analysis/correlation_clear_rs/max_correlation_eos_par.nc")
par<-ncvar_get(ncin,"max_cor")
nc_close(ncin)

# ##### get the correlation for all points, ignore land cover tyeps
# lc_type<-c("ENF","EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO","URB","CNV")
# id<-c(1,4,5,7,8,9,10,12)
# bind_dat<-cbind(as.vector(ai_data),as.vector(tday),as.vector(lc_data))
# pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/ai_climate_by_veg_grow.pdf",width=8,height=8)
# par(mfrow=c(3,3),mar=c(3,3,3,3),oma=c(3,3,1,1))
# for (i in 1:8){
#   lc_sub<-bind_dat[bind_dat[,3]==id[i],]
#   heatscatter(lc_sub[,1],lc_sub[,2],xlim=c(0,1.5),ylim=c(-1,1),main="")
#   mtext(side=3,line=0,lc_type[id[i]],cex=1.3,font=2)
#   abline(v=0.5)
# }
# mtext(side=2,line=0,"Correlation",outer=T,cex=1.6)
# mtext(side=1,line=0,"Aridity Index",outer=T,cex=1.6)
# dev.off()

##### analysis 2:
####    get the mean eos, EOS-Tday R, and AI and graph the EOS ai relationship
bufferline<-function(v=v){
  for (i in 1:10){
    abline(v=v,lwd=i*i*i*0.01,col=adjustcolor("grey30",alpha.f = 0.7-i*0.06))
  }
}
# use the coarse vegetation map
lc_data1<- lc_data
lc_data1[lc_data>=1&lc_data<=5]<-1  #forest
lc_data1[lc_data>=6&lc_data<=8]<-2  # woody land
lc_data1[lc_data>=9&lc_data<=11]<-3 # grass wet sav
lc_data1[lc_data==12|lc_data==14]<-4 #cropland
lc_data1[lc_data1>4|lc_data<1]<-NA
lc_ras<-raster(lc_data)

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Pre_tday_control.pdf",height=11.2,width=8)
par(oma=c(0,1,1,1),mgp=c(3,0.3,0))
par(fig=c(0,1,0.7,1),mar=c(3,4,1,8))
all<-cbind(as.vector(lc_data1),as.vector(ai_data),as.vector(tday),as.vector(eos_data))
all_na<-all[complete.cases(all),]
col_ano<-adjustcolor(ano_col_gmt[round((all_na[,3]+1)*8)],alpha.f = 0.6)
plot(all_na[,2],all_na[,4],col=col_ano,cex=0.3,pch=15,xlim=c(0,1.8),ylim=c(0,365),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
mtext(side=1,line=1.8,"Aridity Index")
mtext(side=2,line=2.2,"End of Photosynthesis (DOY)")
mtext(side=2,line=4,'a',cex=1,font=2,padj=-10,las=2)
axis(1,at=c(0,0.05,0.2,0.5,0.65,1,1.8),tck = -0.02)
axis(2,at=c((0:7)*50),las=2,tck = -0.02)
axis(2,at=c((0:36)*10),label=rep("",37),tck = -0.01)
text(2.3,187, labels = expression(italic(r)[paste(T[day],",EOP",sep="")]), xpd = NA, srt = -90,cex=1.1)   
text(1.98,50, labels = "Density", xpd = NA, srt = -90,cex=1.1)
arrows(x0 = 0.55,x1 = 0.25, y0 = 340,y1=340,length = 0.1)
arrows(x0 = 0.65,x1 = 0.95, y0 = 340,y1=340,length = 0.1)
text(0.13,340,"dry")
text(1.07,340,"wet")
text(0.95,155,expression("RED:"),col="red",pos=2)
text(0.9,150,expression("higher T"[day] %=>%"delayed EOP"),pos=4)
text(0.95,125,expression("BLUE:"),col="blue",pos=2)
text(0.9,120,expression("higher T"[day] %=>%"advanced EOP"),pos=4)
####
negativeall<-all_na[all_na[,3]<0,]
positiveall<-all_na[all_na[,3]>0,]
negativeall[negativeall[,2]>1.8,2]=NA
positiveall[positiveall[,2]>1.8,2]=NA
par(new=T)
histneg<-hist(negativeall[,2],breaks=30,plot=F)
barplot(histneg$density,ylim=c(0,7),width=0.05,xlim=c(0,1.8),axes=F,space=0,xaxs='i',yaxs="i",
     col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border="white")
histpos<-hist(positiveall[,2],breaks=30,plot=F)
barplot(histpos$density,ylim=c(0,7),width=0.05,xlim=c(0,1.8),axes=F,space=0,xaxs='i',yaxs="i",
        col=adjustcolor(ano_col_gmt[15],alpha.f=0.7),border="white",add=T)
bufferline(v=0.6)
axis(4,at=(0:4)/2,las=2,tck = -0.02)
box()
par(fig=c(0.2,1,0.7,1),mar=c(3,4,1,8),new=T)
plot(lc_ras, legend.only=TRUE, col=ano_col_gmt,horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))

##### prec
all<-cbind(as.vector(lc_data1),as.vector(ai_data),as.vector(prec),as.vector(eos_data))
all_na<-all[complete.cases(all),]
col_ano<-adjustcolor(ano_col_gmt[round((all_na[,3]+1)*8)],alpha.f = 0.6)
par(fig=c(0,1,0.4,0.7),mar=c(3,4,1,8),new=T)
plot(all_na[,2],all_na[,4],col=col_ano,cex=0.3,pch=15,xlim=c(0,1.8),ylim=c(0,365),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
mtext(side=1,line=1.8,"Aridity Index")
mtext(side=2,line=2.2,"End of Photosynthesis (DOY)")
mtext(side=2,line=4,'b',cex=1,font=2,padj=-10,las=2)
axis(1,at=c(0,0.05,0.2,0.5,0.65,1,1.8),tck = -0.02)
axis(2,at=c((0:7)*50),las=2,tck = -0.02)
axis(2,at=c((0:36)*10),label=rep("",37),tck = -0.01)
text(2.3,187, labels = expression(italic(r)["Prec.,EOP"]), xpd = NA, srt = -90,cex=1.1)   
text(1.98,50, labels = "Density", xpd = NA, srt = -90,cex=1.1)  
arrows(x0 = 0.55,x1 = 0.25, y0 = 340,y1=340,length = 0.1)
arrows(x0 = 0.65,x1 = 0.95, y0 = 340,y1=340,length = 0.1)
text(0.13,340,"dry")
text(1.07,340,"wet")
text(0.95,155,expression("RED:"),col="red",pos=2)
text(0.9,150,expression("higher Prec."%=>%"delayed EOP"),pos=4)
text(0.95,125,expression("BLUE:"),col="blue",pos=2)
text(0.9,120,expression("higher Prec." %=>%"advanced EOP"),pos=4)
##
negativeall<-all_na[all_na[,3]<0,]
positiveall<-all_na[all_na[,3]>0,]
negativeall[negativeall[,2]>1.8,2]=NA
positiveall[positiveall[,2]>1.8,2]=NA
par(new=T)
histneg<-hist(negativeall[,2],breaks=30,plot=F)
barplot(histneg$density,ylim=c(0,7),width=0.05,xlim=c(0,1.8),axes=F,space=0,xaxs='i',yaxs="i",
        col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border="white")
histpos<-hist(positiveall[,2],breaks=30,plot=F)
barplot(histpos$density,ylim=c(0,7),width=0.05,xlim=c(0,1.8),axes=F,space=0,xaxs='i',yaxs="i",
        col=adjustcolor(ano_col_gmt[15],alpha.f=0.7),border="white",add=T)
axis(4,at=(0:4)/2,las=2,tck = -0.02)
box()
bufferline(v=0.6)
par(fig=c(0.2,1,0.4,0.7),mar=c(3,4,1,8),new=T)
plot(lc_ras, legend.only=TRUE, col=ano_col_gmt,horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))
#######
par(fig=c(0,1,0.1,0.4),mar=c(3,4,1,8),new=T)
biome_col<-c("darkgreen","red","coral","gold2")
plot(NA,xlim=c(0,1.8),ylim=c(0,3),xaxs="i",yaxs="i",xlab="",ylab="",axes=F)
mtext(side=1,line=1.8,"Aridity Index")
mtext(side=2,line=2.2,"Density")
mtext(side=2,line=4,'c',cex=1,font=2,padj=-10,las=2)
axis(1,at=c(0,0.05,0.2,0.5,0.65,1,1.8),tck = -0.02)
axis(2,at=c(0:6/2),tck = -0.02,las=2)
arrows(x0 = 0.55,x1 = 0.25, y0 = 2.794521,y1=2.794521,length = 0.1)
arrows(x0 = 0.65,x1 = 0.95, y0 = 2.794521,y1=2.794521,length = 0.1)
text(0.13,2.794521,"dry")
text(1.07,2.794521,"wet")
for (i in 1:4){
  x<-density(all_na[all_na[,1]==i,2])
  lines(x$x,x$y,col=biome_col[i],lwd=4)
}
bufferline(v=0.6)
box()
legend(1.3,2.8,c("Forest","Woodland","Grassland","Cropland"),
       lty=rep(1,4),lwd=rep(4,4),col=biome_col,bty = 'n')

par(fig=c(0,1,0,0.1),mar=c(1,4,0,8),new=T)
plot(NA,axes=F,xlab="",ylab="",xlim=c(0,1.8),ylim=c(0,1),xaxs="i",yaxs='i')
ai_range<-c(0,0.05,0.2,0.5,0.65,1.8)
ai_color<-c('red',"firebrick1","darkorange","cyan4","blue")
raw_col<-c(0.8,0,0,
           0.9,0.5,0,
           1,1,0,
           0.3,1,0,
           0,0.6,0)
dim(raw_col)<-c(3,5)
ai_col<-rgb(t(raw_col))
for (i in 1:5){
  lines(c(ai_range[i],ai_range[i+1]),c(0.6,0.6),col=ai_col[i],lend=1,lwd=5)
}
text(0.025,0.4,"hyperarid",xpd = NA)
text(0.125,0.8,"arid")
text(0.35,0.4,"semi-arid")
text(0.575,0.8,"dry subhumid")
text(1.225,0.4,"humid")
dev.off()
### need to add some descriptions


