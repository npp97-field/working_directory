##### map phenology using the polynomial maximum method for VI
###
library(ncdf4)
library(zoo)
library(parallel)
get_two_threshold<-function(ts){
  ### get the mean seasonal cycle for each pixel.
  ts[ts<= 0]<-NA
  if (sum(is.na(ts))>20*14){
    return(c(NA,NA))
  }
  ts_fill<-na.fill(ts,fill = "extend")
  dim(ts_fill)<-c(23,14)
  meants<-apply(ts_fill,1,mean,na.rm=T)
  ### calcualte the NDVI ratio
  ndvi_ratio<-(meants[2:23]-meants[1:22])/meants[1:22]
  #get max (positive) for sos, and min (negative) for eos
  return(meants[c(which.max(ndvi_ratio),which.min(ndvi_ratio)+1)])
}

extract_phenology<-function(ts_thresh){
  if (is.na(sum(ts_thresh[24:25]))){
    return(c(NA,NA))
  }
  sos_thresh<-ts_thresh[24]
  eos_thresh<-ts_thresh[25]
  #### fit the ts with a 6 order polynomial
  doy<-0.5:22.5
  if (sum(!is.na(ts_thresh[1:23]))<=6){
    return(c(NA,NA))
  }
  ts_thresh[1:23]<-na.fill(ts_thresh[1:23],fill="extend")
  ts_fit<-lm(ts_thresh[1:23]~poly(doy,6,raw=T))
  b<-coef(ts_fit)
  predicted<-b[1] + poly(1:365/16,6,raw=T)%*%b[-1]
  ### get the date
  peak<-max(predicted,na.rm=T)
  if (is.infinite(peak)){
    return(c(NA,NA))
  }
  increase<-c(rep(1,which(predicted==peak)[1]),rep(0,365-which(predicted==peak)[1]))
  decrease<-c(rep(0,which(predicted==peak)[1]-1),rep(1,365-which(predicted==peak)[1]+1))
  if (peak<sos_thresh){
    sos_acc<-NA
  }else{
    sos_int<-max(which(predicted<sos_thresh&increase))
    sos_acc<-sos_int+(sos_thresh-predicted[sos_int])/(predicted[sos_int+1]-predicted[sos_int])  
  }
  if (peak<eos_thresh){
    eos_acc<-NA
  }else{
    eos_int<-min(which(predicted<eos_thresh&decrease))-1
    eos_acc<-eos_int+(predicted[eos_int]-eos_thresh)/(predicted[eos_int]-predicted[eos_int+1])
  }
  return(c(sos_acc,eos_acc))
}
# for (j in 1:86400){
#   extract_phenology(year_with_thresh[j,])
# }

export_nc<-function(outdata,outfile){
  latmin<- 30
  latmax<- 90
  latd<- 0.5
  lonmin<- -180
  lonmax<- 180
  lond<- 0.5
  
  lat<- seq(latmin+latd/2,latmax-latd/2,latd)
  long<-seq(lonmin+lond/2,lonmax-lond/2,lond)
  
  dimlat<-ncdim_def('latitude','deg',lat)
  dimlong<-ncdim_def('longitude','deg',long)
  ncsos<-ncvar_def('SOS','NA',list(dimlong,dimlat),-9999,longname="start of season",prec='float',compression=9)
  nceos<-ncvar_def('EOS','',list(dimlong,dimlat),-9999,longname="end of season",prec='float',compression=9)
  ncsos_t<-ncvar_def('THESHOLD SOS','',list(dimlong,dimlat),-9999,longname="threshold of SOS",prec='float',compression=9)
  nceos_t<-ncvar_def('THESHOLD EOS','',list(dimlong,dimlat),-9999,longname="threshold of EOS",prec='float',compression=9)
  
  if(file.exists(outfile)){
    file.remove(outfile)
  }
  
  ncout<-nc_create(outfile,list(ncsos,nceos,ncsos_t,nceos_t))
  
  dim(outdata)<-c(4,720,120)
  ncvar_put(ncout,varid=ncsos,outdata[1,,])
  ncvar_put(ncout,varid=nceos,outdata[2,,])
  ncvar_put(ncout,varid=ncsos_t,outdata[3,,])
  ncvar_put(ncout,varid=nceos_t,outdata[4,,])
  nc_close(ncout)
}


setwd("/rigel/glab/users/zy2309/")
####read the 14 year data
NDVI<-array(NA,dim=c(720,120,23*14))
vi_files<-list.files("./DATA/MOD13C1_HD/",full.names = T)

cl<-makeCluster(getOption("cl.cores",10))
clusterEvalQ(cl, {
  library(zoo)
})
for (i in 1:length(vi_files)){
  ncf<-nc_open(vi_files[i])
  NDVI_temp<-ncvar_get(ncf,varid = "NDVI")
  NDVI[,,((i-1)*23+1):(i*23)]<-NDVI_temp[,241:360,]
}
dim(NDVI)<-c(86400,322)
NDVI[is.nan(NDVI)]<-NA
#NDVIthresh<-apply(NDVI,1,get_two_threshold)
NDVIthresh<-parRapply(cl,NDVI,get_two_threshold)
dim(NDVIthresh)<-c(2,86400)
# fit the time series of each year using a six order polynomial
for (i in 1:length(vi_files)){
  year_ndvi<-NDVI[,((i-1)*23+1):(i*23)]
  year_with_thresh<-cbind(year_ndvi,t(NDVIthresh))
  #sos_eos<-apply(year_with_thresh,1,extract_phenology)
  sos_eos<-parRapply(cl,year_with_thresh,extract_phenology)
  dim(sos_eos)<-c(2,86400)
  ### export to nc file
  outfile<-paste("./PROJECT/SIF_phenology/pheno_VI_hd_PM/VI_SOS_EOS_PM_",i+2002,".nc",sep="")
  outdata<-rbind(sos_eos,NDVIthresh)
  export_nc(outdata,outfile)
}

stopCluster(cl)

