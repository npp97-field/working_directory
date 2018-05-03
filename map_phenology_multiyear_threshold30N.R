####map phenology using the server.
#library(phenopix)
library(zoo)
library(ncdf4)
library(parallel)
###
setwd("/rigel/glab/users/zy2309/")
year_files<-list.files("./DATA/all_daily_SIF_4day_HD/",full.names = T)
##read the theshold from each year
phenology_files<-list.files("./PROJECT/SIF_phenology/pheno_hd_var_threshold/",pattern="30N",full.names = T)
#plot(1:255,pch=16,cex=2,col=lst_color_ramp)

years<-2003:2016
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
      return(c(-9999,-9999,-9999))
    }
    
    increase<-c(rep(1,which(var_ts==peak)),rep(0,n-which(var_ts==peak)))
    decrease<-c(rep(0,which(var_ts==peak)-1),rep(1,n-which(var_ts==peak)+1))
    
    
    sos_acc<- -9999
    eos_acc<- -9999
    #thresh_val<- -9999
    
    #thresh_val<-trough+thresh*ampl
    sos_int<-max(which(var_ts<thresh_val&increase))
    sos_acc<-sos_int+(thresh_val-var_ts[sos_int])/(var_ts[sos_int+1]-var_ts[sos_int])  
    eos_int<-min(which(var_ts<thresh_val&decrease))-1
    eos_acc<-eos_int+(var_ts[eos_int]-thresh_val)/(var_ts[eos_int]-var_ts[eos_int+1])
    
    return(c((sos_acc-1)/n,(eos_acc-1)/n,thresh_val))
  }
  x<-xt[1:92]
  x[is.nan(x)]<-NA
  x[x< -5]<-NA
  if (sum(is.na(x))>20|is.na(thresh_val)){
    return(c(-9999,-9999,-9999))
  }
  x_fill<-na.fill(x,fill = "extend")
  sp_fit<-smooth.spline(x_fill,df=9)
  phen<-extract_thresh(sp_fit$y,thresh_val)
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
  ncsos<-ncvar_def('SOS','NA',list(dimlong,dimlat),-9999,longname="start of season",prec='float',compression=9)
  nceos<-ncvar_def('EOS','',list(dimlong,dimlat),-9999,longname="end of season",prec='float',compression=9)
  ncthr<-ncvar_def('THESHOLD','',list(dimlong,dimlat),-9999,longname="threshold of CSIF",prec='float',compression=9)
  
  ncout<-nc_create(outfile,list(ncsos,nceos,ncthr))
  ncvar_put(ncout,varid=ncsos,pheno_dat[1,,])
  ncvar_put(ncout,varid=nceos,pheno_dat[2,,])
  ncvar_put(ncout,varid=ncthr,pheno_dat[3,,])
  nc_close(ncout)
}

### get the threshold
annual_threshold<-array(NA, dim=c(86400,14))
for(i in 1:length(phenology_files)){
  year<-substr(basename(phenology_files[i]),20,23)
  ncf<-nc_open(phenology_files[i])
  ncdat<-ncvar_get(ncf,varid = "THESHOLD")
  nc_close(ncf)
  #yeardata<-ncdat[,241:360,]
  ## get the threshold 
  annual_threshold[,i]<-ncdat
}
annual_threshold[annual_threshold< -100]<-NA
multiyear_threshold<-apply(annual_threshold,1,mean,na.rm=T)

#yeardata<-array(NA,dim=c(7200*1200,92))
cl<-makeCluster(getOption("cl.cores",6))
clusterEvalQ(cl, {
  library(zoo)
})

for(i in 1:length(year_files)){
  year<-substr(basename(year_files[i]),20,23)
  ncf<-nc_open(year_files[i])
  ncdat<-ncvar_get(ncf,varid = "all_daily_sif")
  nc_close(ncf)
  yeardata<-ncdat[,241:360,]
  dim(yeardata)<-c(86400,92)
  yeardata_with_thresh<-cbind(yeardata,multiyear_threshold)
  sif_pheno<-parRapply(cl,yeardata_with_thresh,retrieve_pheno_fix)
  #sif_pheno<-apply(yeardata_with_thresh,1,retrieve_pheno_fix)
  #sif_res<-unlist(sif_pheno)
  dim(sif_pheno)<-c(3,720,120)
  outfile<-paste("./PROJECT/SIF_phenology/pheno_hd_fixed_threshold/SOS_EOS_threshold_",year,".nc",sep="")
  if(file.exists(outfile))
    file.remove(outfile)
  export_nc(sif_pheno,outfile)
}
stopCluster(cl)

