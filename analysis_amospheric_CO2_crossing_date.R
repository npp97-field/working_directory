##### analyze the atmospheric CO2 data
co2_brw<-read.table("/Users/yzhang/Data/CO2_data/co2_brw_surface-insitu_1_ccgg_DailyData.txt",sep=' ')
names(co2_brw)<-c('site_code','year','month','day','hour','minute','second','value','value_unc','nvalue',
                  'latitude','longitude','altitude','elevation','intake_height','instrument','qcflag')
co2_brw<-co2_brw[co2_brw$year>=1980&co2_brw$year<=2016,]
date_brw<-strptime(paste(co2_brw$year,'-',co2_brw$month,'-',co2_brw$day,sep=""),format="%Y-%m-%d")
#plot(date_brw,co2_brw$value)
co2_brw$value[co2_brw$value<0]<-NA


co2_mlo<-read.table('/Users/yzhang/Data/CO2_data/co2_mlo_surface-insitu_1_ccgg_DailyData.txt',sep=" ")
names(co2_mlo)<-c('site_code','year','month','day','hour','minute','second','value','value_unc','nvalue',
                  'latitude','longitude','altitude','elevation','intake_height','instrument','qcflag')
co2_mlo<-co2_mlo[co2_mlo$year>=2000&co2_mlo$year<=2016,]
date_mlo<-strptime(paste(co2_mlo$year,'-',co2_mlo$month,'-',co2_mlo$day,sep=""),format="%Y-%m-%d")
#points(date_mlo,co2_mlo$value,col='red')
co2_mlo$value[co2_mlo$value<0]<-NA
plot(date_brw,co2_brw$value,type='l')
lines(date_mlo,co2_mlo$value,col='red')



library(TSA)
ts_brw<-ts(data = co2_brw$value,start = c(1980,1),frequency=365.2425)
har_brw<-harmonic(ts_brw,4)
xt<-(1:length(co2_brw$value))/365.2425

fit_brw<-lm(co2_brw$value~har_brw+poly(xt,4))

plot(date_brw,co2_brw$value,type="l")
lines(date_brw[!is.na(co2_brw$value)],fit_brw$fitted.values,col='red')
plot(date_brw[!is.na(co2_brw$value)],fit_brw$residuals,type="l")
resid<-fit_brw$residuals
library(latticeExtra)
library(zoo)
smoothed_resi<-simpleSmoothTs(resid,c=3,sides=2,width=30)
lines(date_brw[!is.na(co2_brw$value)],smoothed_resi,col='red')

###add the smoothed residual back to fitted to reconstruct co2
smoothed_brw<-smoothed_resi+fit_brw$fitted.values
plot(date_brw[!is.na(co2_brw$value)],smoothed_brw)

fit_coef<-fit_brw$coefficients
deseasonal<-fit_coef[1]+poly(xt,4)%*%fit_coef[10:13]
plot(date_brw,deseasonal)

seasonal = smoothed_brw-deseasonal[!is.na(co2_brw$value)]
plot(seasonal)
co2_brw$seasonal[!is.na(co2_brw$value)]<-seasonal
plot(date_brw,co2_brw$seasonal,type="l")
co2_brw$no_gap_seasonal<-na.fill(co2_brw$seasonal,fill="extend")


szc<-c()
fzc<-c()
for (y in 1980:2016){
  dat<-co2_brw$no_gap_seasonal[co2_brw$year==y]
  eop_start<-length(dat)-which.max(rev(dat)<0)+1
  sop_start<-which.max(dat<0)-1
  fzc[y-1979]<-eop_start+dat[eop_start]/(dat[eop_start]-dat[eop_start+1])
  szc[y-1979]<-sop_start+dat[sop_start]/(dat[sop_start]-dat[sop_start+1])
}

##########compare with the average EOP date during 2000-2016
library(ncdf4)
library(raster)
phen_files<-list.files("/Users/yzhang/Project/SIF_phenology/analysis/pheno_hd_fixed_threshold_clear/",full.names = T)
weightedeos<-list()
weightedsos<-list()
land_area<-raster('/Users/yzhang/Data/land/0.5deg.area.tif')
landval<-getValues(land_area)
dim(landval)<-c(720,360)
northland_area<-landval[241:360]
average_eos<-c()
average_sos<-c()
for (i in 1:16){
  ncin<-nc_open(phen_files[i])
  if(i ==1){
    thresh<-ncvar_get(ncin,varid="THESHOLD")
  }
  weightedeos[[i]]<-ncvar_get(ncin,varid='EOS')*thresh*northland_area
  weightedsos[[i]]<-ncvar_get(ncin,varid='SOS')*thresh*northland_area
  average_eos[i]<-mean(weightedeos[[i]],na.rm=T)/mean(thresh*northland_area,na.rm=T)*365
  average_sos[i]<-mean(weightedsos[[i]],na.rm=T)/mean(thresh*northland_area,na.rm=T)*365
}


phen_files<-list.files("/Users/yzhang/Project/SIF_phenology/analysis/pheno_fluxcom/ANN/",full.names = T)
weightedeos_flux<-list()
weightedsos_flux<-list()
average_eos_flux<-c()
average_sos_flux<-c()
for (i in 1:13){
  ncin<-nc_open(phen_files[i])
  if(i ==1){
    thresh<-ncvar_get(ncin,varid="THESHOLD")
  }
  weightedeos_flux[[i]]<-ncvar_get(ncin,varid='EOS')*thresh*northland_area
  weightedsos_flux[[i]]<-ncvar_get(ncin,varid='SOS')*thresh*northland_area
  average_eos_flux[i]<-mean(weightedeos_flux[[i]],na.rm=T)/mean(thresh*northland_area,na.rm=T)*365
  average_sos_flux[i]<-mean(weightedsos_flux[[i]],na.rm=T)/mean(thresh*northland_area,na.rm=T)*365
}

plot(average_eos*365,ylim=c(200,340))
plot(szc[1:23],col='red',type="l")

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/correlation_between_SOP_EOP_and_crosing_date.pdf",width=8,height=6)
par(fig=c(0,1,0.5,1),mgp=c(3,0.3,0),mar=c(4,4,0.4,0.4))
plot(2001:2016,average_sos-mean(average_sos),xlim=c(2001,2016),
     ylim=c(-15,15),ylab="",xlab="",axes=F,col="blue",lwd=2)
points(2001:2013,average_sos_flux-mean(average_sos_flux),col="darkgreen",lwd=2)
lines(2001:2016,szc[22:37]-mean(szc[22:37]),col="red",lwd=2)
axis(1)
axis(2,las=2)
mtext(side=1,line=2,"Year")
mtext(side=2,line=2.3,"DOY")
box()
text(2013,13,"r=0.48, p=0.059",col="blue")
text(2013,10,"r=-0.02, p=0.95",col="darkgreen")
legend("topleft",c("average start of photosynthesis (CSIF)","average start of photosynthesis (FLUXCOM)","spring zero crossing date"),
       lty=c(NA,NA,1),pch=c(1,1,NA),col=c("blue","darkgreen","red"),bty="n")
cor.test(szc[22:37],average_sos)
cor.test(szc[22:34],average_sos_flux)

par(fig=c(0,1,0,0.5),mgp=c(3,0.3,0),mar=c(4,4,0.4,0.4),new=T)
plot(2001:2016,average_eos-mean(average_eos),xlim=c(2001,2016),
     ylim=c(-25,25),ylab="",xlab="",axes=F,col="blue",lwd=2)
points(2001:2013,average_eos_flux-mean(average_eos_flux),col="darkgreen",lwd=2)
lines(2001:2016,fzc[22:37]-mean(fzc[22:37]),col="red",lwd=2)
axis(1)
axis(2,las=2)
mtext(side=1,line=2,"Year")
mtext(side=2,line=2.3,"DOY")
box()
cor.test(fzc[22:37],average_eos)

text(2010,21,"r=-0.15, p=0.59",col="blue")
text(2010,17,"r=-0.30, p=0.31",col="darkgreen")
dev.off()


