##########
## this code plot the temperature dependence for the EOS and SOS
# to see if a different pattern between the date and temperature exist.
# this only applies for the temperature limited regions, i.e., based on the 
# AI > 0.6
# A density plot will be generated,the x axis is the interannual anomaly
# for each pixel, and the x axis is the anomaly of EOS and SOS for each 
# pixel.
# the correlation is first calculate with 1 month lags, but also calculated for other lags.
#

library(rgdal)
library(ncdf4)
library(raster)
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")
####  read the data for each pixel.
# including 1 month preseason temp for both SOS and EOS. and the SOS and EOS for each pixel.
pheno_files<-list.files("./pheno_hd_fixed_threshold_clear/",full.names = T)
climate_files<-list.files("./pheno_hd_fixed_threshold_climate/clear_era/",full.names = T)

# data matrix
pre_T_sos_1mon<-array(NA,dim=c(720,120,16))
pre_T_eos_1mon<-array(NA,dim=c(720,120,16))
sos<-array(NA,dim=c(720,120,16))
eos<-array(NA,dim=c(720,120,16))

### read in the data
for (i in 1:length(pheno_files)){
  ncin<-nc_open(pheno_files[i])
  sos[,,i]<-ncvar_get(ncin, varid = "SOS")
  eos[,,i]<-ncvar_get(ncin, varid = "EOS")
  ncin<-nc_open(climate_files[i])
  pre_T_sos_1mon[,,i]<-ncvar_get(ncin, varid = "pre_start1_temp")
  pre_T_eos_1mon[,,i]<-ncvar_get(ncin, varid = "pre_end1_temp")
  nc_close(ncin)
}

#### calcualte the anomaly
## this function read in the ts of TEMP, SOS/EOS and return their anomalies
ano_for_all<-function(vec){
  if (sum(is.na(vec))>10){
    return(rep(NA,length(vec)))
  }
  mean_vec<-mean(vec,na.rm=T)
  return(vec-mean_vec)
}

ano_sos<-apply(sos,c(1,2),ano_for_all)
ano_eos<-apply(eos,c(1,2),ano_for_all)
ano_pre_T_sos<-apply(pre_T_sos_1mon,c(1,2),ano_for_all)
ano_pre_T_eos<-apply(pre_T_eos_1mon,c(1,2),ano_for_all)

#### mask out the barren and ocean areas. use the aridity threshold and get the 
#### temperature limited regions
mask_out<-function(dat,lmask,aridity_mask){
  #### the data is a 720 by 120 by 16 matrix
  aridity_mask[aridity_mask<0.6]<-NA
  aridity_mask[!is.na(aridity_mask)]<-1
  totalmask<-rep(lmask*aridity_mask,each=16)
  dim(totalmask)<-c(16,720,120)
  return(totalmask*dat)
}

ncin<-nc_open("./analysis/North_barren_mask.nc")
landmask<-ncvar_get(ncin,"barren")
ncin<-nc_open("./analysis/climate/ai_cru_P_PET.nc")
ai<-ncvar_get(ncin,"ai")
ncin<-nc_open("./analysis/North_mcd12c1_landcover1_majority.nc")
lc<-ncvar_get(ncin,"lc1")
m_sos_ano<-mask_out(ano_sos,landmask,ai)
m_eos_ano<-mask_out(ano_eos,landmask,ai)
m_ano_pre_T_sos<-mask_out(ano_pre_T_sos,landmask,ai)
m_ano_pre_T_eos<-mask_out(ano_pre_T_eos,landmask,ai)

dim(m_sos_ano)<-c(16,86400)
dim(m_eos_ano)<-c(16,86400)
dim(m_ano_pre_T_sos)<-c(16,86400)
dim(m_ano_pre_T_eos)<-c(16,86400)

out_dat<-t(rbind(m_sos_ano,m_ano_pre_T_sos,m_eos_ano,m_ano_pre_T_eos,as.vector(lc)))
out_complete<-out_dat[complete.cases(out_dat),]   #### something wrong here, need to check


numpix<-rep(NA,14)

sos_stat_out<-array(NA, dim=c(16,14,4))
#names(sos_stat_out)<-c("mean_date","sd_date","mean_temp","sd_temp")
eos_stat_out<-array(NA, dim=c(16,14,4))

for (i in 1:16){  # i for the 16 years
  for (j in 1:14){   # j for the 16 land cover types
    sos_stat_out[i,j,1]<-mean(out_complete[out_complete[,65]==j,i])*365
    sos_stat_out[i,j,2]<-sd(out_complete[out_complete[,65]==j,i])*365
    sos_stat_out[i,j,3]<-mean(out_complete[out_complete[,65]==j,i+16])
    sos_stat_out[i,j,4]<-sd(out_complete[out_complete[,65]==j,i+16])
    
    eos_stat_out[i,j,1]<-mean(out_complete[out_complete[,65]==j,i+32])*365
    eos_stat_out[i,j,2]<-sd(out_complete[out_complete[,65]==j,i+32])*365
    eos_stat_out[i,j,3]<-mean(out_complete[out_complete[,65]==j,i+48])
    eos_stat_out[i,j,4]<-sd(out_complete[out_complete[,65]==j,i+48])
  }
}
for (j in 1:14){
  numpix[j]<-sum(out_complete[,65]==j)
}

save(numpix,sos_stat_out,eos_stat_out,file="./analysis/ano_pheno_climate.RData")
## get all values. and plot them into a density map.





######################################################################
#load("/Users/yzhang/Project/SIF_phenology/analysis/ano_pheno_climate.RData")
#setwd("/Users/yzhang/habanero/PROJECT/SIF_phenology/")
setwd("/Users/yzhang/PROJECT/SIF_phenology/")
load("./analysis/ano_pheno_climate.RData")
# sos_stat_out[,,1]<-sos_stat_out[,,1]*365
# sos_stat_out[,,2]<-sos_stat_out[,,2]*365
# eos_stat_out[,,1]<-eos_stat_out[,,1]*365
# eos_stat_out[,,2]<-eos_stat_out[,,2]*365
# library(LSD)
# par(mfcol=c(1,2))
# heatscatter(as.vector(out_complete[,1:16])*365,as.vector(out_complete[,17:32]),xlim=c(-7,7))
# heatscatter(as.vector(out_complete[,33:48])*365,as.vector(out_complete[,49:64]),xlim=c(-7,7))
# library(plotKML)
# data(worldgrids_pal)
# lc_col<-worldgrids_pal$IGBP[2:15]
lc_col<-c("palegreen4","green4","green2","yellowgreen","lightseagreen",
               'darkmagenta',"coral3","tan4","orange2","darkolivegreen",'royalblue4',"darkorange3")
lc_used<-c(1,4,5,7,8,9,10,11)
lc_type<-c("ENF","EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO","URB","CNV")
s_slope<-c()
e_slope<-c()

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Fig/temp_sensitivity.pdf",width=8,height=4)
par(fig=c(0,0.5,0,1),mar=c(3.5,3.5,1,1),oma=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(NA,xlim=c(-3,3),ylim=c(-12,12),axes=F,xaxs='i',yaxs='i',xlab="",ylab="")
mtext(side=1,line=2,expression(paste("Preseason T anomaly (",degree,"C)",sep="")))
mtext(side=2, line=2,"SOP anomaly (days)")
box()
axis(1,tck=-0.01)
axis(2,las=2,tck=-0.01)
for (j in 1:length(lc_used)){
  points(sos_stat_out[,lc_used[j],3],sos_stat_out[,lc_used[j],1],col=lc_col[lc_used[j]],cex=0.4)
  reg<-lm(sos_stat_out[,lc_used[j],1]~sos_stat_out[,lc_used[j],3])
  s_slope[j]<-reg$coefficients[2]
  range<-c(min(sos_stat_out[,lc_used[j],3],na.rm=T),max(sos_stat_out[,lc_used[j],3],na.rm=T))
  lines(range,range*reg$coefficients[2]+reg$coefficients[1],col=lc_col[lc_used[j]])
}

legend("bottomleft",lc_type[lc_used],col=lc_col[lc_used],pch=rep(1,9),bty="n",cex=0.7)
reg_all<-lm(as.vector(sos_stat_out[,lc_used,1])~as.vector(sos_stat_out[,lc_used,3]))
corr<-cor.test(as.vector(sos_stat_out[,lc_used,1]),as.vector(sos_stat_out[,lc_used,3]))
text(3,10,pos=2,substitute(paste(italic(gamma),"=",slo%+-%std),
                           list(slo=round(mean(s_slope),2),std=round(sd(s_slope),2))))
text(3,7,pos=2,substitute(paste(italic(r),"=",slo,"   P<0.001"),
                           list(slo=round(corr$estimate,2))))

par(fig=c(0.5,1,0,1),new=T)
plot(NA,xlim=c(-3,3),ylim=c(-12,12),axes=F,xaxs='i',yaxs='i',xlab="",ylab="")
mtext(side=1,line=2,expression(paste("Preseason T anomaly (",degree,"C)",sep="")))
mtext(side=2, line=2,"EOP anomaly (days)")
box()
axis(1,tck=-0.01)
axis(2,las=2,tck=-0.01)
for (j in 1:length(lc_used)){
  points(eos_stat_out[,lc_used[j],3],eos_stat_out[,lc_used[j],1],col=lc_col[lc_used[j]],cex=0.4)
  reg<-lm(eos_stat_out[,lc_used[j],1]~eos_stat_out[,lc_used[j],3])
  e_slope[j]<-reg$coefficients[2]
  range<-c(min(eos_stat_out[,lc_used[j],3],na.rm=T),max(eos_stat_out[,lc_used[j],3],na.rm=T))
  lines(range,range*reg$coefficients[2]+reg$coefficients[1],col=lc_col[lc_used[j]])
}

reg_all<-lm(as.vector(eos_stat_out[,lc_used,1])~as.vector(eos_stat_out[,lc_used,3]))
corr<-cor.test(as.vector(eos_stat_out[,lc_used,1]),as.vector(eos_stat_out[,lc_used,3]))
text(3,10,pos=2,substitute(paste(italic(gamma),"=",slo%+-%std),
                           list(slo=round(mean(e_slope),2),std=round(sd(e_slope),2))))
text(3,7,pos=2,substitute(paste(italic(r),"=",slo,"   P<0.001"),
                          list(slo=round(corr$estimate,2))))

dev.off()

