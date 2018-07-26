##### analyze the atmospheric CO2 data
co2_brw<-read.table("/Users/yzhang/Data/CO2_data/co2_brw_surface-insitu_1_ccgg_DailyData.txt",sep=' ')
names(co2_brw)<-c('site_code','year','month','day','hour','minute','second','value','value_unc','nvalue',
                  'latitude','longitude','altitude','elevation','intake_height','instrument','qcflag')
co2_brw<-co2_brw[co2_brw$year>=1975&co2_brw$year<=2016,]
date_brw<-strptime(paste(co2_brw$year,'-',co2_brw$month,'-',co2_brw$day,sep=""),format="%Y-%m-%d")
#plot(date_brw,co2_brw$value)
co2_brw$value[co2_brw$value<0]<-NA


co2_mlo<-read.table('/Users/yzhang/Data/CO2_data/co2_mlo_surface-insitu_1_ccgg_DailyData.txt',sep=" ")
names(co2_mlo)<-c('site_code','year','month','day','hour','minute','second','value','value_unc','nvalue',
                  'latitude','longitude','altitude','elevation','intake_height','instrument','qcflag')
co2_mlo<-co2_mlo[co2_mlo$year>=1975&co2_mlo$year<=2016,]
date_mlo<-strptime(paste(co2_mlo$year,'-',co2_mlo$month,'-',co2_mlo$day,sep=""),format="%Y-%m-%d")
#points(date_mlo,co2_mlo$value,col='red')
co2_mlo$value[co2_mlo$value<0]<-NA
plot(date_brw,co2_brw$value,type='l')
lines(date_mlo,co2_mlo$value,col='red')



library(TSA)


get_critical<-function(co2_brw){

  ts_brw<-ts(data = co2_brw$value,start = c(1975,1),frequency=365.2425)
  har_brw<-harmonic(ts_brw,4)
  xt<-(1:length(co2_brw$value))/365.2425
  
  #fit the mean seasonal + interannual trend
  fit_brw<-lm(co2_brw$value~har_brw+poly(xt,4))
  
  plot(date_brw,co2_brw$value,type="l")
  lines(date_brw[!is.na(co2_brw$value)],fit_brw$fitted.values,col='red')
  
  co2_brw$residuals[!is.na(co2_brw$value)]<-fit_brw$residuals
  residual_std<-sd(co2_brw$residuals,na.rm=T)
  threshold<-quantile(co2_brw$residuals,c(0.01,0.99),na.rm=T)
  res<-co2_brw$residuals
  res[res<threshold[1]|res>threshold[2]]<-NA
  co2_brw$residuals_no_outlier<-res
  plot(date_brw,co2_brw$residuals,type="l")
  lines(date_brw,co2_brw$residuals_no_outlier,type="l",col='red')
  
  resid<-co2_brw$residuals_no_outlier
  
  ####### this is the function for the gaussian filter with a sigma value
  gaussian_filter<-function(ts,fwhm,cutoff=3){
    sigma<-fwhm/2.3548
    size_filter<-round(fwhm/2.3548*cutoff)
    #generate the kernal
    kernal<-rep(NA,size_filter*2-1)
    for (i in 1:(size_filter*2-1)){
      kernal[i]<-exp((-(i-size_filter)^2)/(2*sigma^2))/(sigma*sqrt(2*pi))
    }
    #### apply the kernal
    ts_filtered<-ts
    for (i in size_filter:(length(ts)-size_filter+1)){
      ts_filtered[i]<-sum(ts[(i-size_filter+1):(i+size_filter-1)]*kernal,na.rm=T)/
        (sum((!is.na(ts[(i-size_filter+1):(i+size_filter-1)]))*kernal))
    }
    return(ts_filtered)
  }

  # smoothed_resi<-simpleSmoothTs(resid,c=3,sides=2,width=30)
  # smoothed_resi<-smth.gaussian(resid,window = 30)
  smoothed_resi<-gaussian_filter(ts = c(resid,resid[(length(resid)-364):length(resid)]),fwhm = 30)[1:length(resid)]
  lines(date_brw,smoothed_resi,col='blue')
  
  fit_coef<-fit_brw$coefficients
  ###add the smoothed residual back to fitted to reconstruct co2
  fitted<-fit_coef[1]+poly(xt,4)%*%fit_coef[10:13]+har_brw%*%fit_coef[2:9]
  
  smoothed_brw<-smoothed_resi+fitted
  deseasonal<-gaussian_filter(c(smoothed_brw,rep(smoothed_brw[(length(smoothed_brw)-364):length(smoothed_brw)],2)), 
                              fwhm=390)[1:length(smoothed_brw)]

  plot(date_brw,smoothed_brw,type='l')
  lines(date_brw,co2_brw$value,col='red')
  lines(date_brw,deseasonal,col="orange")
  plot(date_brw,smoothed_brw-deseasonal)
  fit_coef<-fit_brw$coefficients
  
  co2_brw$deseasonal<-deseasonal
  seasonal = smoothed_brw-deseasonal
  plot(date_brw,seasonal)
  co2_brw$seasonal<-seasonal

  szc<-rep(NA,37)
  fzc<-rep(NA,37)
  though<-rep(NA,37)
  for (y in 1980:2016){
    dat<-co2_brw$seasonal[co2_brw$year==y|co2_brw$year==y+1][100:500]
    eop_start<-length(dat)-which.max(rev(dat)<0)+1
    sop_start<-which.max(dat<0)-1
    fzc[y-1979]<-eop_start+dat[eop_start]/(dat[eop_start]-dat[eop_start+1])+100
    szc[y-1979]<-sop_start+dat[sop_start]/(dat[sop_start]-dat[sop_start+1])+100
    though[y-1979]<-min(dat,na.rm=T)
  }
  return(list(co2_brw[co2_brw$year>=1980,],szc,fzc,though))
}

brw<-get_critical(co2_brw)
mlo<-get_critical(co2_mlo)

plot(brw[[1]]$deseasonal,mlo[[1]]$deseasonal)
abline(0,1)
plot(mlo[[1]]$deseasonal)
lines(brw[[1]]$deseasonal)

plot(brw[[4]],ylim=c(-15,0))
lines(mlo[[4]])
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


phen_files<-list.files("/Users/yzhang/Project/SIF_phenology/analysis/pheno_fluxcom/ANN/",full.names = T)[2:35]
weightedeos_flux<-list()
weightedsos_flux<-list()
average_eos_flux<-c()
average_sos_flux<-c()
for (i in 1:34){
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



co2_brw<-brw[[1]]
szc<-brw[[2]]
fzc<-brw[[3]]
though<-brw[[4]]

library(pracma)
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/correlation_between_SOP_EOP_and_crosing_date.pdf",width=8,height=6)
par(fig=c(0,1,0.5,1),mgp=c(3,0.3,0),mar=c(4,4,0.4,0.4))
ano_sos<-average_sos-mean(average_sos,na.rm=T)
ano_flux<-average_sos_flux-mean(average_sos_flux,na.rm=T)
ano_szc<-szc-mean(szc,na.rm=T)

years<-1980:2016
plot(years[22:37],ano_sos,xlim=c(1980,2016),
     ylim=c(-10,10),ylab="",xlab="",axes=F,col="blue",type="l",lwd=2,xaxs="i",yaxs="i")
lines(years[1:34],ano_flux,col="darkgreen",lwd=2)
lines(years,szc-mean(szc,na.rm=T),col="red",lwd=2)
axis(1,at=1980:2016,label=rep("",37),tck=-0.01)
axis(1,at=0:7*5+1980,tck=-0.03)
axis(2,las=2,tck=-0.03)
mtext(side=1,line=2,"Year")
mtext(side=2,line=2.3,"DOY")
box()
text(2002,13,"r=0.48, p=0.059",col="blue")
text(2002,10,"r=-0.02, p=0.95",col="darkgreen")
legend("topright",c("Average SOP (CSIF)","Average SOP (FLUXCOM)","spring zero crossing date"),
       lty=c(1,1,1),lwd=c(2,2,2),col=c("blue","darkgreen","red"),bty="n")
cor.test(szc[22:37],average_sos)
cor.test(szc[1:34],average_sos_flux)
##### regression slopes
reg1<-lm(ano_szc~years)
lines(years,reg1$fitted.values,col="red",lty=2)
reg4<-lm(ano_szc[22:37]~years[22:37])
lines(years[22:37],reg4$fitted.values,col="red",lty=2)
reg2<-lm(ano_sos~years[22:37])
lines(years[22:37],reg2$fitted.values,col="blue",lty=2)
reg3<-lm(ano_flux~years[1:34])
lines(years[1:34],reg3$fitted.values,col="darkgreen",lty=2)
##### significance test


par(fig=c(0,1,0,0.5),mgp=c(3,0.3,0),mar=c(4,4,0.4,0.4),new=T)
ano_eos<-average_eos-mean(average_eos,na.rm=T)
ano_eos_flux<-average_eos_flux-mean(average_eos_flux,na.rm=T)
ano_fzc<-fzc-mean(fzc,na.rm=T)
plot(2001:2016,ano_eos,xlim=c(1980,2016),
     ylim=c(-15,15),ylab="",xlab="",axes=F,col="blue",lwd=2,type='l',xaxs="i",yaxs="i")
lines(1980:2013,ano_eos_flux,col="darkgreen",lwd=2)
lines(1980:2016,ano_fzc,col="red",lwd=2)
axis(1,at=1980:2016,label=rep("",37),tck=-0.01)
axis(1,at=0:7*5+1980,tck=-0.03)
axis(2,las=2,tck=-0.03)
mtext(side=1,line=2,"Year")
mtext(side=2,line=2.3,"DOY")
box()
cor.test(fzc[22:37],average_eos)
reg1<-lm(ano_fzc~years)
lines(years,reg1$fitted.values,col="red",lty=2)
reg2<-lm(ano_eos~years[22:37])
lines(years[22:37],reg2$fitted.values,col="blue",lty=2)
reg3<-lm(ano_eos_flux~years[1:34])
lines(years[1:34],reg3$fitted.values,col="darkgreen",lty=2)
text(2002,21,"r=-0.15, p=0.59",col="blue")
text(2002,16,"r=-0.30, p=0.31",col="darkgreen")
dev.off()


plot_reg<-function(years,ano,col_used){
  reg<-lm(ano~years)
  lines(years,reg$fitted.values,col=col_used,lty=2)
}

jet.color <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(37)
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/BRW.pdf",height=8,width=8)
par(fig=c(0,1,0.5,1),mar=c(4,4,1,4))
plot(NA,xlim=c(1,365),ylim=c(-15,10),axes=F,xaxs="i",yaxs="i",xlab="",ylab="")
axis(1)
axis(2,las=2)
box()
mtext(side=1,line=2.3,"DOY")
mtext(side=2,line=2.3,"Detrended CO2 (ppm)")
abline(h=0,lty=2)
for (y in 1980:2016){
  year_data<-co2_brw[co2_brw$year==y,]
  lines(1:length(year_data$seasonal),year_data$seasonal,col=jet.color[y-1979])  
}
temp<-raster(array(1:10,dim=c(2,5)))
par(fig=c(0.5,1,0.5,1),mar=c(4,4,1,4))
plot(temp, legend.only=TRUE, col=rev(jet.color),horizontal=F,zlim=c(1980,2016),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=1980+(0:7)*5,
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))
ano_though<-though-mean(though)
par(fig=c(0,1,0,0.5),new=T)
plot(years,ano_szc,ylim=c(-15,15),lwd=2,col="red",type='l',axes=F,xaxs="i",yaxs="i",xlab="",ylab="")
axis(1)
axis(2,las=2)
box()
mtext(side=1,line=2.3,"DOY")
mtext(side=2,line=2.3,"Zero crossing date (day)")
plot_reg(years,ano_szc,"red")
lines(years,ano_fzc,lwd=2,col="blue")
plot_reg(years,ano_fzc,"blue")
par(new=T)
plot(years,ano_though,lwd=2,col='darkgreen',type='l',ylim=c(-5,5),axes=F,xlab="",ylab="")
plot_reg(years,ano_though,"darkgreen")
axis(4,las=2,col="darkgreen")
mtext(side=4,line=2.3,"Amplitude (ppm)",col="darkgreen")
dev.off()

