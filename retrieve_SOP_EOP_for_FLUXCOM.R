#### retrieve SOP EOP for fluxcom
#library(phenopix)
library(zoo)
library(ncdf4)
library(parallel)
library(abind)
###

models<-c("ANN","MARS","RF")

get_threshold<-function(var_ts,thresh=0.3){
  n<-length(var_ts)
  if (sum(is.na(var_ts))>1000){
    return(NA)
  }
  dim(var_ts)<-c(365,13)
  season_mean<-apply(var_ts,1,mean,na.rm=T)
  
  peak <- max(season_mean, na.rm=TRUE)
  trough<-max(min(season_mean, na.rm=TRUE),0)
  if (trough<0.1*peak){
    trough<-0
  }
  ampl <- peak - trough
  if (peak<=0){
    return(NA)
  }
  thresh_val<- NA
  thresh_val<-trough+thresh*ampl
  return(thresh_val)
}

retrieve_pheno_fix<-function(xt){
  thresh_val<-xt[366]  
  extract_thresh<-function(var_ts,thresh_val){
    #var_ts<-coredata(in_var_ts)
    n <- length(var_ts)
    avg <- mean(var_ts, na.rm=TRUE)
    
    peak <- max(var_ts, na.rm=TRUE)
    trough <- max(min(var_ts, na.rm=TRUE),0)
    ampl <- peak - trough
    if (peak<=0|peak<thresh_val){
      return(rep(NA,4))
    }
    pos<-which(var_ts==peak)[1]
    
    increase<-c(rep(1,pos),rep(0,n-pos))
    decrease<-c(rep(0,pos-1),rep(1,n-pos+1))
    
    sos_acc<- NA
    eos_acc<- NA
    
    sos_int<-max(which(var_ts<thresh_val&increase))
    sos_acc<-sos_int+(thresh_val-var_ts[sos_int])/(var_ts[sos_int+1]-var_ts[sos_int])  
    eos_int<-min(which(var_ts<thresh_val&decrease))-1
    eos_acc<-eos_int+(var_ts[eos_int]-thresh_val)/(var_ts[eos_int]-var_ts[eos_int+1])
    
    return(c((sos_acc-1)/n,pos/n,(eos_acc-1)/n,thresh_val))
  }
  if (is.na(thresh_val)){
    return(rep(NA,4))
  }
  phen<-extract_thresh(xt[1:365],thresh_val)
  return(as.vector(phen))
}

smooth_GPP<-function(var_ts){
  if (sum(is.na(var_ts))==365){
    return(rep(NA,365))
  }
  sp_fit<-smooth.spline(var_ts,df=9)
  return(sp_fit$y)
}


export_nc<-function(pheno_dat,outfile,xdim,ydim){
  ncsos<-ncvar_def('SOS','',list(xdim,ydim),-9999,longname="start of season",prec='float',compression=9)
  ncpos<-ncvar_def('POS','',list(xdim,ydim),-9999,longname="peak of season",prec='float',compression=9)
  nceos<-ncvar_def('EOS','',list(xdim,ydim),-9999,longname="end of season",prec='float',compression=9)
  ncthr<-ncvar_def('THESHOLD','',list(xdim,ydim),-9999,longname="threshold of CSIF",prec='float',compression=9)
  
  ncout<-nc_create(outfile,list(ncsos,ncpos,nceos,ncthr))
  ncvar_put(ncout,varid=ncsos,pheno_dat[1,,])
  ncvar_put(ncout,varid=ncpos,pheno_dat[2,,])
  ncvar_put(ncout,varid=nceos,pheno_dat[3,,])
  ncvar_put(ncout,varid=ncthr,pheno_dat[4,,])
  nc_close(ncout)
}


setwd("/rigel/glab/users/zy2309/")
cl<-makeCluster(getOption("cl.cores",10))
clusterEvalQ(cl, {
  library(zoo)
})

gpp_day_f<-list.files("./DATA/FLUXCOM/daily/",pattern="GPP_HB",full.names = T)
gpp_night_f<-list.files("./DATA/FLUXCOM/daily/",pattern="GPP.daily",full.names = T)
retrieve_sopeop<-function(gppfile,mod){
  ###read in the 20 year's GPP data from 198501-200412
  modelf<-gppfile[grep(mod,gppfile)]
  #get the multiyear average first.
  #
  num_f<-length(modelf)
  smoothed_north_gpp<-array(NA,dim=c(365*13,720*120))
  for (i in (num_f-12):num_f){
    ncin<-nc_open(modelf[i])
    gpp<-ncvar_get(ncin,varid = ncin$var[[1]])
    xdim<-ncin$dim[["lon"]]
    yval<-ncin$dim[["lat"]]$val[1:120]
    nc_close(ncin)
    smoothed_north_gpp[((i-num_f+12)*365+1):((i-num_f+13)*365),]<-parApply(cl,gpp[,1:120,1:365],c(1,2),smooth_GPP)
  }
  ### get the threshold
  multiyear_threshold<-parCapply(cl,smoothed_north_gpp,get_threshold)
  ydim<-ncdim_def("lat","",vals = yval)
  ##  retrieve phenology for each year!
  sop<-array(NA,dim=c(720,120,13))
  eop<-array(NA,dim=c(720,120,13))
  pop<-array(NA,dim=c(720,120,13))
  threshold<-array(NA,dim=c(720,120,13))
  dir.create(paste("./PROJECT/SIF_phenology/pheno_fluxcom/",mod,"/",sep=""),recursive = T)
  for (y in 1:13){
    pheno<-parCapply(cl,rbind(smoothed_north_gpp[(y*365-364):(y*365),],multiyear_threshold),retrieve_pheno_fix)
    dim(pheno)<-c(4,720,120)
    file.remove(paste("./PROJECT/SIF_phenology/pheno_fluxcom/",mod,"/pheno_",mod,"_",y+2000,".nc",sep=""))
    export_nc(pheno,paste("./PROJECT/SIF_phenology/pheno_fluxcom/",mod,"/pheno_",mod,"_",y+2000,".nc",sep=""),xdim,ydim)
    sop[,,y]<-pheno[1,,]
    pop[,,y]<-pheno[2,,]
    eop[,,y]<-pheno[3,,]
    threshold[,,y]<-pheno[4,,]
  }
  
  # get the mean average SOP and EOP
  mean_sop<-apply(sop,c(1,2),mean,na.rm=T)
  mean_pop<-apply(pop,c(1,2),mean,na.rm=T)
  mean_eop<-apply(eop,c(1,2),mean,na.rm=T)
  mean_threshold<-apply(threshold,c(1,2),mean,na.rm=T)
  outdat<-abind(mean_sop,mean_pop,mean_eop,mean_threshold,along=3)
  ## output to nc files
  file.remove(paste("./PROJECT/SIF_phenology/pheno_fluxcom/",mod,"/pheno_",mod,"_mean.nc",sep=""))
  export_nc(pheno,paste("./PROJECT/SIF_phenology/pheno_fluxcom/",mod,"/pheno_",mod,"_mean.nc",sep=""),xdim,ydim)
}

for (i in 1:3){
  retrieve_sopeop(gppfile = gpp_night_f,models[i])
}


