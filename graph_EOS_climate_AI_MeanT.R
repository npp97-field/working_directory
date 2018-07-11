#### graph and analyze the EOS climate relationship with aridity index.
## 
library("LSD")
## use the heatscatter function
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")
setwd("/Users/yzhang/Project/SIF_phenology/")
ncin<-nc_open("./data/ai_cru_P_PET.nc")
#ncin<-nc_open("./data/ai_cru_P_PET_grow.nc")
ai_data<-ncvar_get(ncin,"ai")
nc_close(ncin)
ncin<-nc_open("./data/mean_annual_temp.nc")
T_data<-ncvar_get(ncin,"mean_annual_temp")
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
bufferline<-function(v=NA,h=NA){
  for (i in 1:10){
    abline(v=v,lwd=i*i*i*0.01,col=adjustcolor("grey30",alpha.f = 0.7-i*0.06))
  }
  for (i in 1:10){
    abline(h=h,lwd=i*i*i*0.01,col=adjustcolor("grey30",alpha.f = 0.7-i*0.06))
  }
}
# use the coarse vegetation map
lc_data1<- lc_data
lc_data1[lc_data>=1&lc_data<=5]<-1  #forest
lc_data1[lc_data>=6&lc_data<=8]<-2  # woody land
lc_data1[lc_data>=9&lc_data<=11]<-3 # grass wet sav
lc_data1[lc_data==12|lc_data==14]<-NA #cropland
lc_data1[lc_data1>4|lc_data<1]<-NA
lc_ras<-raster(lc_data)

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Pre_tday_control_AI_TEMP.pdf",height=12.5,width=5)
par(oma=c(1,1,0,1),mgp=c(3,0.3,0))
par(fig=c(0,0.7,0.7,0.93),mar=c(3,4,0,0))
all<-cbind(as.vector(lc_data1),as.vector(ai_data),as.vector(T_data),as.vector(tday),as.vector(eos_data))
all_na<-all[complete.cases(all),]
col_ano<-adjustcolor(ano_col_gmt[round((all_na[,4]+1)*8)],alpha.f = 0.6)
plot(all_na[,2],all_na[,3],col=col_ano,cex=0.3,pch=15,ylim=c(-25,25),xlim=c(0,1.8),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
mtext(side=1,line=1.4,"Aridity Index")
mtext(side=2,line=2.2,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
mtext(side=2,line=4,'a',cex=1,font=2,padj=-11,las=2)
text(2.7,0, labels = expression(italic(r)[paste(T[day],",EOP",sep="")]), xpd = NA, srt = -90,cex=1.1)   
axis(1,at=c(0,0.05,0.2,0.5,0.65,1,1.8),tck = -0.02)
axis(2,at=(-5:5)*5,tck = -0.02,las=2)
box()

par(fig=c(0,0.7,0.93,1),mar=c(0,4,1,0),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,2]> 1.8|negativeall[,2]< 0,2]=NA
positiveall[positiveall[,2]> 1.8|positiveall[,2]< 0,2]=NA
histneg<-hist(negativeall[,2],breaks=0:30*0.06,plot=F)
plot(NA,xlim=c(0,1.8),ylim=c(0,3),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
rect(xleft=histneg$mids-0.023, ybottom=0, xright=histneg$mids+0.023, ytop=histneg$density,
     col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border=NA)
histpos<-hist(positiveall[,2],breaks=0:30*0.06,plot=F)
rect(xleft=histpos$mids-0.023, ybottom=0, xright=histpos$mids+0.023, ytop=histpos$density,
     col=adjustcolor(ano_col_gmt[15],alpha.f=0.7),border=NA)
bufferline(v=0.6)
arrows(x0 = 0.55,x1 = 0.2, y0 = 2.5,y1=2.5,length = 0.1)
arrows(x0 = 0.65,x1 = 1, y0 = 2.5,y1=2.5,length = 0.1)
text(0.08,2.5,"dry",xpd = NA)
text(1.12,2.5,"wet",xpd = NA)
text(1.65,1.1,expression("RED:"),col="red",pos=2,xpd = NA,cex=0.78)
text(1.53,1.0,expression("higher T"[day] %=>%"delayed EOP"),pos=4,xpd = NA,cex=0.78)
text(1.65,1.85,expression("BLUE:"),col="blue",pos=2,xpd = NA,cex=0.78)
text(1.53,1.75,expression("higher T"[day] %=>%"advanced EOP"),pos=4,xpd = NA,cex=0.78)
####
par(fig=c(0.7,0.9,0.7,0.93),mar=c(3,0,1,2),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,3]>25|negativeall[,3]< -25,3]=NA
positiveall[positiveall[,3]>25|positiveall[,3]< -25,3]=NA
histneg<-hist(negativeall[,3],breaks=-15:15*(5/3),plot=F)
plot(NA,xlim=c(0,0.1),ylim=c(-25,25),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
rect(xleft=0, ybottom=histneg$mids-0.766666, xright=histneg$density, ytop=histneg$mids+0.766666,
     col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border=NA)
histpos<-hist(positiveall[,3],breaks=-15:15*(5/3),plot=F)
rect(xleft=0, ybottom=histpos$mids-0.766666, xright=histpos$density, ytop=histpos$mids+0.766666,
     col=adjustcolor(ano_col_gmt[15],alpha.f=0.7),border=NA)
bufferline(h=2.5)
arrows(x0 = 0.08333333,x1 = 0.08333333, y0 = 0,y1=-10,length = 0.1)
arrows(x0 = 0.08333333,x1 = 0.08333333, y0 = 5,y1=15,length = 0.1)
text(0.08333333,-14.33,"cold",xpd = NA,srt = -90)
text(0.08333333,21.33,"warm",xpd = NA,srt = -90)

par(fig=c(0.2,1,0.7,0.93),mar=c(3,4,1,3),new=T)
plot(lc_ras, legend.only=TRUE, col=ano_col_gmt,horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))

##### prec
par(oma=c(0,1,1,1),mgp=c(3,0.3,0))
par(fig=c(0,0.7,0.4,0.63),mar=c(3,4,0,0),new=T)
all<-cbind(as.vector(lc_data1),as.vector(ai_data),as.vector(T_data),as.vector(prec),as.vector(eos_data))
all_na<-all[complete.cases(all),]
col_ano<-adjustcolor(ano_col_gmt[round((all_na[,4]+1)*8)],alpha.f = 0.6)
plot(all_na[,2],all_na[,3],col=col_ano,cex=0.3,pch=15,ylim=c(-25,25),xlim=c(0,1.8),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
mtext(side=1,line=1.4,"Aridity Index")
mtext(side=2,line=2.2,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
text(2.7,0, labels = expression(italic(r)[paste("Prec.",",EOP",sep="")]), xpd = NA, srt = -90,cex=1.1)   
mtext(side=2,line=4,'b',cex=1,font=2,padj=-11,las=2)
axis(1,at=c(0,0.05,0.2,0.5,0.65,1,1.8),tck = -0.02)
axis(2,at=(-5:5)*5,tck = -0.02,las=2)
box()
par(fig=c(0,0.7,0.63,0.7),mar=c(0,4,1,0),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,2]> 1.8|negativeall[,2]< 0,2]=NA
positiveall[positiveall[,2]> 1.8|positiveall[,2]< 0,2]=NA
histneg<-hist(negativeall[,2],breaks=0:30*0.06,plot=F)
plot(NA,xlim=c(0,1.8),ylim=c(0,3),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
rect(xleft=histneg$mids-0.023, ybottom=0, xright=histneg$mids+0.023, ytop=histneg$density,
     col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border=NA)
histpos<-hist(positiveall[,2],breaks=0:30*0.06,plot=F)
rect(xleft=histpos$mids-0.023, ybottom=0, xright=histpos$mids+0.023, ytop=histpos$density,
     col=adjustcolor(ano_col_gmt[15],alpha.f=0.7),border=NA)
bufferline(v=0.6)
arrows(x0 = 0.55,x1 = 0.2, y0 = 2.5,y1=2.5,length = 0.1)
arrows(x0 = 0.65,x1 = 1, y0 = 2.5,y1=2.5,length = 0.1)
text(0.08,2.5,"dry",xpd = NA)
text(1.12,2.5,"wet",xpd = NA)
text(1.65,1.1,expression("RED:"),col="red",pos=2,xpd = NA,cex=0.78)
text(1.53,1.0,expression("higher Prec."%=>%"delayed EOP"),pos=4,xpd = NA,cex=0.78)
text(1.65,1.85,expression("BLUE:"),col="blue",pos=2,xpd = NA,cex=0.78)
text(1.53,1.75,expression("higher Prec."%=>%"advanced EOP"),pos=4,xpd = NA,cex=0.78)

####
par(fig=c(0.7,0.9,0.4,0.63),mar=c(3,0,1,2),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,3]>25|negativeall[,3]< -25,3]=NA
positiveall[positiveall[,3]>25|positiveall[,3]< -25,3]=NA
histneg<-hist(negativeall[,3],breaks=-15:15*(5/3),plot=F)
plot(NA,xlim=c(0,0.1),ylim=c(-25,25),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
rect(xleft=0, ybottom=histneg$mids-0.766666, xright=histneg$density, ytop=histneg$mids+0.766666,
     col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border=NA)
histpos<-hist(positiveall[,3],breaks=-15:15*(5/3),plot=F)
rect(xleft=0, ybottom=histpos$mids-0.766666, xright=histpos$density, ytop=histpos$mids+0.766666,
     col=adjustcolor(ano_col_gmt[15],alpha.f=0.7),border=NA)
bufferline(h=2.5)
arrows(x0 = 0.08333333,x1 = 0.08333333, y0 = 0,y1=-10,length = 0.1)
arrows(x0 = 0.08333333,x1 = 0.08333333, y0 = 5,y1=15,length = 0.1)
text(0.08333333,-14.33,"cold",xpd = NA,srt = -90)
text(0.08333333,21.33,"warm",xpd = NA,srt = -90)

par(fig=c(0.2,1,0.4,0.63),mar=c(3,4,1,3),new=T)
plot(lc_ras, legend.only=TRUE, col=ano_col_gmt,horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))

#######

biome_col<-adjustcolor(c("darkgreen","blue","red"),alpha.f = 0.6)

par(fig=c(0,0.7,0.1,0.33),mar=c(3,4,0,0),new=T)
all<-cbind(as.vector(lc_data1),as.vector(ai_data),as.vector(T_data),as.vector(tday),as.vector(eos_data))
all_na<-all[complete.cases(all),]
col_ano<-biome_col[all_na[,1]]
plot(all_na[,2],all_na[,3],col=col_ano,cex=0.3,pch=15,ylim=c(-25,25),xlim=c(0,1.8),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
mtext(side=1,line=1.4,"Aridity Index")
mtext(side=2,line=2.2,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
mtext(side=2,line=4,'c',cex=1,font=2,padj=-11,las=2)
axis(1,at=c(0,0.05,0.2,0.5,0.65,1,1.8),tck = -0.02)
axis(2,at=(-5:5)*5,tck = -0.02,las=2)
box()

par(fig=c(0,0.7,0.33,0.4),mar=c(0,4,1,0),new=T)
forest<-all_na[all_na[,1]==1,]
woodland<-all_na[all_na[,1]==2,]
grassland<-all_na[all_na[,1]==3,]

forest[forest[,2]> 1.8|forest[,2]< 0,2]=NA
woodland[woodland[,2]> 1.8|woodland[,2]< 0,2]=NA
grassland[grassland[,2]> 1.8|grassland[,2]< 0,2]=NA

histfor<-hist(forest[,2],breaks=0:30*0.06,plot=F)
plot(NA,xlim=c(0,1.8),ylim=c(0,3),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
rect(xleft=histfor$mids-0.023, ybottom=0, xright=histfor$mids+0.023, ytop=histfor$density,
     col=adjustcolor(biome_col[1],alpha.f=0.7),border=NA)
histwoo<-hist(woodland[,2],breaks=0:30*0.06,plot=F)
rect(xleft=histwoo$mids-0.023, ybottom=0, xright=histwoo$mids+0.023, ytop=histwoo$density,
     col=adjustcolor(biome_col[2],alpha.f=0.7),border=NA)
histgra<-hist(grassland[,2],breaks=0:30*0.06,plot=F)
rect(xleft=histgra$mids-0.023, ybottom=0, xright=histgra$mids+0.023, ytop=histgra$density,
     col=adjustcolor(biome_col[3],alpha.f=0.7),border=NA)
bufferline(v=0.6)
arrows(x0 = 0.55,x1 = 0.2, y0 = 2.5,y1=2.5,length = 0.1)
arrows(x0 = 0.65,x1 = 1, y0 = 2.5,y1=2.5,length = 0.1)
text(0.08,2.5,"dry",xpd = NA)
text(1.1,2.5,"wet",xpd = NA)
####
par(fig=c(0.7,0.9,0.1,0.33),mar=c(3,0,1,2),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,3]>25|negativeall[,3]< -25,3]=NA
positiveall[positiveall[,3]>25|positiveall[,3]< -25,3]=NA
histneg<-hist(negativeall[,3],breaks=-15:15*(5/3),plot=F)
plot(NA,xlim=c(0,0.12),ylim=c(-25,25),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
forest<-all_na[all_na[,1]==1,]
woodland<-all_na[all_na[,1]==2,]
grassland<-all_na[all_na[,1]==3,]

forest[forest[,3]> 25|forest[,3]< -25,3]=NA
woodland[woodland[,3]> 25|woodland[,3]< -25,3]=NA
grassland[grassland[,3]> 25|grassland[,3]< -25,3]=NA

histfor<-hist(forest[,3],breaks=-15:15*5/3,plot=F)
histwoo<-hist(woodland[,3],breaks=-15:15*5/3,plot=F)
histgra<-hist(grassland[,3],breaks=-15:15*5/3,plot=F)
rect(xleft=0, ybottom=histfor$mids-0.766666, xright=histfor$density, ytop=histfor$mids+0.766666,
     col=adjustcolor(biome_col[1],alpha.f=0.7),border=NA)
rect(xleft=0, ybottom=histwoo$mids-0.766666, xright=histwoo$density, ytop=histwoo$mids+0.766666,
     col=adjustcolor(biome_col[2],alpha.f=0.7),border=NA)
rect(xleft=0, ybottom=histgra$mids-0.766666, xright=histgra$density, ytop=histgra$mids+0.766666,
     col=adjustcolor(biome_col[3],alpha.f=0.7),border=NA)

bufferline(h=2.5)
arrows(x0 = 0.1,x1 = 0.1, y0 = 0,y1=-10,length = 0.1)
arrows(x0 = 0.1,x1 = 0.1, y0 = 5,y1=15,length = 0.1)
text(0.1,-14.33,"cold",xpd = NA,srt = -90)
text(0.1,21.33,"warm",xpd = NA,srt = -90)

#text(31,73, labels = "Density", xpd = NA, srt = -90,cex=1.1)

# text(-18.5,185,expression("RED:"),col="red",pos=2)
# text(-20,180,expression("higher T"[day] %=>%"delayed EOP"),pos=4)
# text(-18.5,155,expression("BLUE:"),col="blue",pos=2)
# text(-20,150,expression("higher T"[day] %=>%"advanced EOP"),pos=4)
# ####

par(fig=c(0.2,1,0.1,0.33),mar=c(3,4,1,3),new=T)
plot(lc_ras, legend.only=TRUE, col=rep(c("darkgreen","blue","red"),each=20),horizontal=F,zlim=c(0.5,3.5),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=1:3,label=c("Forest","Woodland","Grassland"),
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))
####### aridity 
par(fig=c(0,0.7,0,0.1),mar=c(0,4,0,0),new=T)
plot(NA,axes=F,xlab="",ylab="",xlim=c(0,1.8),ylim=c(0,1),xaxs="i",yaxs='i')
ai_range<-c(0,0.05,0.2,0.5,0.65,1.8)
raw_col<-c(0.8,0,0,
           0.9,0.5,0,
           1,1,0,
           0.3,1,0,
           0,0.6,0)
dim(raw_col)<-c(3,5)
ai_col<-rgb(t(raw_col))
for (i in 1:5){
  lines(c(ai_range[i],ai_range[i+1]),c(0.9,0.9),col=ai_col[i],lend=1,lwd=5)
}
text(0.025-0.07,0.85,"hyperarid",xpd = NA,srt=-90,pos=4)
text(0.125-0.07,0.85,"arid",xpd = NA,srt=-90,pos=4)
text(0.35-0.07,0.85,"semi-arid",xpd = NA,srt=-90,pos=4)
text(0.575-0.07,0.85,"dry subhumid",xpd = NA,srt=-90,pos=4)
text(1.225-0.07,0.85,"humid",xpd = NA,srt=-90,pos=4)
dev.off()
### need to add some descriptions


