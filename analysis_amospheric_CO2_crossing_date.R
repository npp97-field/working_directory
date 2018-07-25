##### analyze the atmospheric CO2 data
co2_brw<-read.table("/Users/yzhang/Data/CO2_data/co2_brw_surface-insitu_1_ccgg_DailyData.txt",sep=' ')
names(co2_brw)<-c('site_code','year','month','day','hour','minute','second','value','value_unc','nvalue',
                  'latitude','longitude','altitude','elevation','intake_height','instrument','qcflag')
co2_brw<-co2_brw[co2_brw$year>=2000&co2_brw$year<=2016,]
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
ts_brw<-ts(data = co2_brw$value,start = c(2000,1),frequency=365.2425)
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
for (y in 2000:2016){
  dat<-co2_brw$no_gap_seasonal[co2_brw$year==y]
  pos_start<-length(dat)-which.max(rev(dat)<0)+1
  szc[y-1999]<-pos_start+dat[pos_start]/(dat[pos_start]-dat[pos_start+1])
}

##########compare with the average EOP date during 2000-2016
library(ncdf4)
library(raster)
phen_files<-list.files("/Users/yzhang/Project/SIF_phenology/analysis/pheno_hd_fixed_threshold_clear/",full.names = T)
weightedeos<-list()
land_area<-raster('/Users/yzhang/Data/land/0.5deg.area.tif')
landval<-getValues(land_area)
dim(landval)<-c(720,360)
northland_area<-landval[,241:360]
average_eos<-c()
for (i in 1:16){
  ncin<-nc_open(phen_files[i])
  if(i ==1){
    thresh<-ncvar_get(ncin,varid="THESHOLD")
  }
  weightedeos[[i]]<-ncvar_get(ncin,varid='EOS')*thresh*northland_area
  average_eos[i]<-mean(weightedeos[[i]],na.rm=T)/mean(thresh*northland_area,na.rm=T)*365
}

plot(average_eos*365,ylim=c(200,340))
lines(szc[2:17],col='red')

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/correlation_between_EOS_and_crosing_date.pdf",width=5,height=4)
plot(average_eos*365,szc[2:17],xlab='average EOP',ylab="zero crossing at fall")
text(263,328,"r=-0.16")
dev.off()
