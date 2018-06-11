####map phenology using the server.
#library(phenopix)
library(zoo)
library(ncdf4)
library(parallel)
###
setwd("/rigel/glab/users/zy2309/")
year_files<-list.files("./DATA/clear_inst_SIF_4day_HD_BISE/",full.names = T)
##read the theshold from each year
#phenology_files<-list.files("./PROJECT/SIF_phenology/pheno_hd_var_threshold/",pattern="30N",full.names = T)
#plot(1:255,pch=16,cex=2,col=lst_color_ramp)

years<-2001:2016
#years=2003

retrieve_pheno_fix<-function(xt){
  thresh_val<-xt[93]
  extract_thresh<-function(var_ts,thresh_val){
    #var_ts<-coredata(in_var_ts)
    n <- length(var_ts)
    avg <- mean(var_ts, na.rm=TRUE)
    
    peak <- max(var_ts, na.rm=TRUE)
    trough <- max(min(var_ts, na.rm=TRUE),0)
    ampl <- peak - trough
    if (peak<=0|peak<thresh_val){
      return(c(-9999,-9999,-9999,-9999))
    }
    pos<-which(var_ts==peak)
    increase<-c(rep(1,pos),rep(0,n-pos))
    decrease<-c(rep(0,pos-1),rep(1,n-pos+1))
    
    
    sos_acc<- -9999
    eos_acc<- -9999
    #thresh_val<- -9999
    
    #thresh_val<-trough+thresh*ampl
    sos_int<-max(which(var_ts<thresh_val&increase))
    sos_acc<-sos_int+(thresh_val-var_ts[sos_int])/(var_ts[sos_int+1]-var_ts[sos_int])  
    eos_int<-min(which(var_ts<thresh_val&decrease))-1
    eos_acc<-eos_int+(var_ts[eos_int]-thresh_val)/(var_ts[eos_int]-var_ts[eos_int+1])
    
    return(c((sos_acc-1)/n,which(var_ts==peak)/n,(eos_acc-1)/n,thresh_val))
  }
  
  phen<-extract_thresh(xt[1:92],thresh_val)
  return(as.vector(phen))
}


export_nc<-function(pheno_dat,outfile){
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
  ncsos<-ncvar_def('SOS','',list(dimlong,dimlat),-9999,longname="start of season",prec='float',compression=9)
  ncpos<-ncvar_def('POS','',list(dimlong,dimlat),-9999,longname="peak of season",prec='float',compression=9)
  nceos<-ncvar_def('EOS','',list(dimlong,dimlat),-9999,longname="end of season",prec='float',compression=9)
  ncthr<-ncvar_def('THESHOLD','',list(dimlong,dimlat),-9999,longname="threshold of CSIF",prec='float',compression=9)
  
  ncout<-nc_create(outfile,list(ncsos,ncpos,nceos,ncthr))
  ncvar_put(ncout,varid=ncsos,pheno_dat[1,,])
  ncvar_put(ncout,varid=ncpos,pheno_dat[2,,])
  ncvar_put(ncout,varid=nceos,pheno_dat[3,,])
  ncvar_put(ncout,varid=ncthr,pheno_dat[4,,])
  nc_close(ncout)
}

### get the threshold 
### to be consistent with the site level-method. modified on May 25th 2018,us multi-year_average seasonal value
get_threshold<-function(var_ts,thresh=0.3){
  n<-length(var_ts)
  dim(var_ts)<-c(92,n/92)
  season_mean<-apply(var_ts,1,mean,na.rm=T)
  
  peak <- max(season_mean, na.rm=TRUE)
  trough<-max(min(season_mean, na.rm=TRUE),0)
  if (trough<0.1*peak){
    trough<-0
  }
  ampl <- peak - trough
  if (peak<=0){
    return(-9999)
  }
  thresh_val<- -9999
  thresh_val<-trough+thresh*ampl
  return(thresh_val)
}

####smooth SIF for each year
smooth_SIF<-function(var_ts){
  var_ts[var_ts< -900]<-0
  var_ts[is.nan(var_ts)]<-0
  var_ts[var_ts< -3|var_ts>3]<-NA
  if (sum(is.na(var_ts))>20){
    return(rep(-9999,length(var_ts)))
  }
  x_fill<-na.fill(var_ts,fill = "extend")
  sp_fit<-smooth.spline(x_fill,df=9)
  return(sp_fit$y)
}

#yeardata<-array(NA,dim=c(7200*1200,92))
cl<-makeCluster(getOption("cl.cores",10))
clusterEvalQ(cl, {
  library(zoo)
})

ncbar<-nc_open('/rigel/glab/users/zy2309/PROJECT/SIF_phenology/analysis/climate/North_barren_mask.nc')
barren<-ncvar_get(ncbar,"barren")
nc_close(ncbar)
barrenmask<-rep(barren,each =4)
dim(barrenmask)<-c(4,720,120)


# first get annual smoothed VI and get threshold
smoothed_CSIF<-array(NA,dim=c(86400,92*length(year_files)))
for(i in 1:length(year_files)){
  year<-substr(basename(year_files[i]),16,19)
  ncf<-nc_open(year_files[i])
  ncdat<-ncvar_get(ncf,varid = "clear_daily_sif")
  nc_close(ncf)
  yeardata<-ncdat[,241:360,]
  dim(yeardata)<-c(86400,92)
  year_smoothed<-parRapply(cl,yeardata,smooth_SIF)
  dim(year_smoothed)<-c(92,86400)
  smoothed_CSIF[,(i*92-91):(i*92)]<-t(year_smoothed)
}
multiyear_threshold<-parRapply(cl,smoothed_CSIF,get_threshold)


for(i in 1:length(year_files)){
  year<-substr(basename(year_files[i]),16,19)
  yeardata_with_thresh<-cbind(smoothed_CSIF[,(i*92-91):(i*92)],multiyear_threshold)
  sif_pheno<-parRapply(cl,yeardata_with_thresh,retrieve_pheno_fix)
  #sif_pheno<-apply(yeardata_with_thresh,1,retrieve_pheno_fix)
  #sif_res<-unlist(sif_pheno)
  dim(sif_pheno)<-c(4,720,120)
  sif_pheno<-sif_pheno*barrenmask
  sif_pheno[is.na(sif_pheno)]<--9999
  outfile<-paste("./PROJECT/SIF_phenology/pheno_hd_fixed_threshold_clear/clear_daily_SOS_POS_EOS_threshold_",year,".nc",sep="")
  if(file.exists(outfile))
    file.remove(outfile)
  export_nc(sif_pheno,outfile)
}
stopCluster(cl)

