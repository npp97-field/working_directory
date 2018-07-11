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
lc_data1[lc_data==12|lc_data==14]<-NA #cropland
lc_data1[lc_data1>4|lc_data<1]<-NA
lc_ras<-raster(lc_data)

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Pre_tday_control_temperature.pdf",height=11.2,width=8)
par(oma=c(0,1,1,1),mgp=c(3,0.3,0))
par(fig=c(0,1,0.7,1),mar=c(3,4,1,8))
all<-cbind(as.vector(lc_data1),as.vector(T_data),as.vector(tday),as.vector(eos_data))
all_na<-all[complete.cases(all),]
col_ano<-adjustcolor(ano_col_gmt[round((all_na[,3]+1)*8)],alpha.f = 0.6)
plot(all_na[,2],all_na[,4],col=col_ano,cex=0.3,pch=15,xlim=c(-25,25),ylim=c(0,365),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
### add regression for <2.5 deg C
cold<-all_na[all_na[,2]<2.5,]
cor.test(cold[,2],cold[,4])
reg<-lm(cold[,4]~cold[,2])   #r 0.611,    slope =1.036
xmin<-min(cold[,2])
xmax<-2.5
coef_reg<-coefficients(reg)
lines(c(xmin,xmax),c(coef_reg[1]+coef_reg[2]*xmin,coef_reg[1]+coef_reg[2]*xmax),lwd=4,lty=2,
      col=adjustcolor("grey3",alpha.f = 0.8))
text(-20,270,pos=4,expression(paste(italic(r),"=0.61 ",italic(p),'<0.001',sep="")))

mtext(side=1,line=1.8,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
mtext(side=2,line=2.2,"End of Photosynthesis (DOY)")
mtext(side=2,line=4,'a',cex=1,font=2,padj=-10,las=2)
axis(1,at=(-5:5)*5,tck = -0.02)
axis(2,at=c((0:7)*50),las=2,tck = -0.02)
axis(2,at=c((0:36)*10),label=rep("",37),tck = -0.01)
text(39,187, labels = expression(italic(r)[paste(T[day],",EOP",sep="")]), xpd = NA, srt = -90,cex=1.1)   
text(31,73, labels = "Density", xpd = NA, srt = -90,cex=1.1)
arrows(x0 = 0,x1 = -10, y0 = 340,y1=340,length = 0.1)
arrows(x0 = 5,x1 = 15, y0 = 340,y1=340,length = 0.1)
text(-12,340,"cold")
text(17,340,"hot")
text(-18.5,185,expression("RED:"),col="red",pos=2)
text(-20,180,expression("higher T"[day] %=>%"delayed EOP"),pos=4)
text(-18.5,155,expression("BLUE:"),col="blue",pos=2)
text(-20,150,expression("higher T"[day] %=>%"advanced EOP"),pos=4)
####
negativeall<-all_na[all_na[,3]<0,]
positiveall<-all_na[all_na[,3]>0,]
negativeall[negativeall[,2]> 25|negativeall[,2]< -25,2]=NA
positiveall[positiveall[,2]>25|positiveall[,2]< -25,2]=NA
par(new=T)
histneg<-hist(negativeall[,2],breaks=-15:15*5/3,plot=F)
# barplot(histneg$density,ylim=c(0,0.2),width=1,xlim=c(0,50),axes=F,space=0,xaxs='i',yaxs="i",
#         col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border="white")
plot(NA,xlim=c(-25,25),ylim=c(0,0.2),axes=F,space=0,xaxs='i',yaxs="i",xlab="",ylab="")
rect(xleft=histneg$mids-0.66667, ybottom=0, xright=histneg$mids+0.66667, ytop=histneg$density,
     col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border="white")
histpos<-hist(positiveall[,2],breaks=-15:15*5/3,plot=F)
rect(xleft=histpos$mids-0.66667, ybottom=0, xright=histpos$mids+0.66667, ytop=histpos$density,
     col=adjustcolor(ano_col_gmt[15],alpha.f=0.7),border="white")
bufferline(v=2.5)
axis(4,at=(0:4)/50,las=2,tck = -0.02)
box()
par(fig=c(0.2,1,0.7,1),mar=c(3,4,1,3),new=T)
plot(lc_ras, legend.only=TRUE, col=ano_col_gmt,horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))

##### prec
all<-cbind(as.vector(lc_data1),as.vector(T_data),as.vector(prec),as.vector(eos_data))
all_na<-all[complete.cases(all),]
col_ano<-adjustcolor(ano_col_gmt[round((all_na[,3]+1)*8)],alpha.f = 0.6)
par(fig=c(0,1,0.4,0.7),mar=c(3,4,1,8),new=T)
plot(all_na[,2],all_na[,4],col=col_ano,cex=0.3,pch=15,xlim=c(-25,25),ylim=c(0,365),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
mtext(side=1,line=1.8,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
mtext(side=2,line=2.2,"End of Photosynthesis (DOY)")
mtext(side=2,line=4,'b',cex=1,font=2,padj=-10,las=2)
axis(1,at=-5:5*5,tck = -0.02)
axis(2,at=c((0:7)*50),las=2,tck = -0.02)
axis(2,at=c((0:36)*10),label=rep("",37),tck = -0.01)
text(39,187, labels = expression(italic(r)["Prec.,EOP"]), xpd = NA, srt = -90,cex=1.1)   
text(31,73, labels = "Density", xpd = NA, srt = -90,cex=1.1)  
arrows(x0 = 0,x1 = -10, y0 = 340,y1=340,length = 0.1)
arrows(x0 = 5,x1 = 15, y0 = 340,y1=340,length = 0.1)
text(-12,340,"cold")
text(17,340,"hot")
text(-18.5,185,expression("RED:"),col="red",pos=2)
text(-20,180,expression("higher Prec."%=>%"delayed EOP"),pos=4)
text(-18.5,155,expression("BLUE:"),col="blue",pos=2)
text(-20,150,expression("higher Prec." %=>%"advanced EOP"),pos=4)
##
negativeall<-all_na[all_na[,3]<0,]
positiveall<-all_na[all_na[,3]>0,]
negativeall[negativeall[,2]> 25|negativeall[,2]< -25,2]=NA
positiveall[positiveall[,2]>25|positiveall[,2]< -25,2]=NA
par(new=T)
histneg<-hist(negativeall[,2],breaks=-15:15*5/3,plot=F)
plot(NA,xlim=c(-25,25),ylim=c(0,0.2),axes=F,space=0,xaxs='i',yaxs="i",xlab="",ylab="")
rect(xleft=histneg$mids-0.66667, ybottom=0, xright=histneg$mids+0.66667, ytop=histneg$density,
     col=adjustcolor(ano_col_gmt[2],alpha.f=0.7),border="white")
histpos<-hist(positiveall[,2],breaks=-15:15*5/3,plot=F)
rect(xleft=histpos$mids-0.66667, ybottom=0, xright=histpos$mids+0.66667, ytop=histpos$density,
     col=adjustcolor(ano_col_gmt[15],alpha.f=0.7),border="white")
axis(4,at=(0:4)/50,las=2,tck = -0.02)
box()
bufferline(v=2.5)
par(fig=c(0.2,1,0.4,0.7),mar=c(3,4,1,8),new=T)
plot(lc_ras, legend.only=TRUE, col=ano_col_gmt,horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))
#######
par(fig=c(0,1,0.1,0.4),mar=c(3,4,1,8),new=T)
biome_col<-c("darkgreen","red","coral","gold2")
plot(NA,xlim=c(-25,25),ylim=c(0,0.12),xaxs="i",yaxs="i",xlab="",ylab="",axes=F)
mtext(side=1,line=1.8,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
mtext(side=2,line=2.2,"Density")
mtext(side=2,line=4,'c',cex=1,font=2,padj=-10,las=2)
axis(1,at=(-5:5)*5,tck = -0.02)
axis(2,at=c(0:4*0.03),tck = -0.02,las=2)
arrows(x0 = 0,x1 = -10, y0 = 0.1117808,y1=0.1117808,length = 0.1)
arrows(x0 = 5,x1 = 15, y0 = 0.1117808,y1=0.1117808,length = 0.1)
text(-12,0.1117808,"cold")
text(17,0.1117808,"hot")
for (i in 1:3){
  x<-density(all_na[all_na[,1]==i,2])
  lines(x$x,x$y,col=biome_col[i],lwd=4)
}
bufferline(v=2.5)
box()
legend(-23,0.11,c("Forest","Woodland","Grassland"),#,"Cropland"),
       lty=rep(1,3),lwd=rep(4,3),col=biome_col[1:3],bty = 'n')

dev.off()
### need to add some descriptions


