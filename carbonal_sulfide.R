#### this is for the carbonal sulfide data analysis

ocs_files<-list.files("/Users/yzhang/Data/OCS/",full.names = T)
brw_ocs<-read.table(ocs_files[2],skip=1,header = T,stringsAsFactors = F)

brw_ocs<-brw_ocs[as.numeric(brw_ocs$yyyymmdd)>20000000,]
obsdate<-strptime(paste(brw_ocs$yyyymmdd,brw_ocs$hhmmss,sep="-"),"%Y%m%d-%H%M")

ocs<-as.numeric(as.character(brw_ocs$OCS_))

plot(obsdate,ocs)


all_year<-as.data.frame(1:round(365.2425*19))
all_year$dat<-NA
all_year$dat[round((as.numeric(brw_ocs$dec_date)-2000)*365.2425)]<-as.numeric(brw_ocs$OCS_)

library(TSA)
ts_brw<-ts(data = all_year$dat,start = c(2000,1),frequency=365.2425)
har_brw<-harmonic(ts_brw,4)
xt<-(1:length(all_year$dat))/365.2425
reg<-lm(all_year$dat~har_brw+poly(xt,4))
fit_coef<-reg$coefficients
fitted<-fit_coef[1]+poly(xt,4)%*%fit_coef[10:13]+har_brw%*%fit_coef[2:9]

plot(2000+1:length(all_year$dat)/365.2425,all_year$dat,ylab="OCS",xlab="year")
lines(2000+1:length(all_year$dat)/365.2425,fitted)

deseasonal<-fit_coef[1]+poly(xt,4)%*%fit_coef[10:13]
residual<-all_year$dat-fitted
plot(residual)
thresh<-quantile(residual,c(0.01,0.99),na.rm=T)
residual[residual<thresh[1]|residual>thresh[2]]<-NA
points(residual,col='red')

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
smoothed<-gaussian_filter(residual,fwhm = 30)
lines(smoothed)

reconstruct<-smoothed+fitted
plot(reconstruct)
interannual<-gaussian_filter(reconstruct,370)
lines(interannual)

seasonal<-reconstruct-interannual
par(fig=c(0,1,0,1),mar=c(4,4,1,4))
plot(NA,xlim=c(1,365),ylim=c(-100,100),axes=F,xaxs="i",yaxs="i",xlab="",ylab="")
jet.color <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(19)
axis(1)
axis(2,las=2)
box()
mtext(side=1,line=2.3,"DOY")
mtext(side=2,line=2.3,"Detrended OCS (pptv)")
abline(h=0,lty=2)
for (i in 1:19){
  year_data<-seasonal[round(365.2425*i-364.2425):round(365.2425*i)]
  lines(1:length(year_data),year_data,col=jet.color[i])  
}

library(raster)
temp<-raster(array(1:10,dim=c(2,5)))
par(fig=c(0.5,1,0,1),mar=c(4,4,1,4))
plot(temp, legend.only=TRUE, col=jet.color,horizontal=F,zlim=c(2000,2018),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=2000+(0:7)*5,
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))


