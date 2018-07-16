####map phenology using the server.
#library(phenopix)
library(zoo)
library(ncdf4)
library(parallel)
library(abind)
###
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")

models<-c("CanESM2","CMCC-CESM","GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR",
          "MIROC-ESM","MPI-ESM2-LR","NorESM1-ME")
for (i in 1:length(models)){
  dir.create(paste("./pheno_CMIP5/",models[i],sep=""))
}
#yeardata<-array(NA,dim=c(7200*1200,92))
cl<-makeCluster(getOption("cl.cores",10))
clusterEvalQ(cl, {
  library(zoo)
})

### get the threshold 
### to be consistent with the site level-method. modified on May 25th 2018,us multi-year_average seasonal value
get_threshold<-function(var_ts,thresh=0.3){
  n<-length(var_ts)
  if (sum(is.na(var_ts))>120){
    return(NA)
  }
  dim(var_ts)<-c(12,20)
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
  thresh_val<-xt[13]  
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
  phen<-extract_thresh(xt[1:12],thresh_val)
  return(as.vector(phen))
}

smooth_GPP<-function(var_ts){
  #var_ts[var_ts< -900]<-0
  #var_ts[is.nan(var_ts)]<-0
  #var_ts[var_ts< -3|var_ts>100]<-NA
  #x_fill<-na.fill(var_ts,fill = "extend")
  if (sum(is.na(var_ts))==12){
    return(rep(NA,12))
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


#### write a function for the phenology retrieval for each models.
gpp_all_f<-list.files("./CMIP5/historical/",pattern="gpp_",full.names = T)
retrieve_sopeop<-function(mod){
  ###read in the 20 year's GPP data from 198501-200412
  ncin<-nc_open(gpp_all_f[grep(mod,gpp_all_f)])
  gpp_mod<-ncvar_get(ncin,varid="gpp")
  xydim<-dim(gpp_mod)[1:2]
  xdim<-ncin$dim[["lon"]]
  ydim<-ncin$dim[["lat"]]
  nc_close(ncin)
  ## reform it to a 2 dimensional array
  dim(gpp_mod)<-c(xydim[1]*xydim[2],240)
  smoothed_gpp<-array(NA,dim=c(xydim[1]*xydim[2],240))
  ## smooth the GPP for each year
  for (i in 1:20){
    temp<-parRapply(cl,gpp_mod[,(i*12-11):(i*12)],smooth_GPP)
    dim(temp)<-c(12,xydim[1]*xydim[2])
    smoothed_gpp[,(i*12-11):(i*12)]<-t(temp)
  }
  ### get the 30% threshold
  multiyear_threshold<-parRapply(cl,smoothed_gpp,get_threshold)
  
  ##  retrieve phenology for each year!
  sop<-array(NA,dim=c(xydim,20))
  eop<-array(NA,dim=c(xydim,20))
  pop<-array(NA,dim=c(xydim,20))
  threshold<-array(NA,dim=c(xydim,20))
  dir.create(paste("./pheno_CMIP5/",mod,"/pheno/",sep=""))
  for (y in 1:20){
    pheno<-parRapply(cl,cbind(smoothed_gpp[,(y*12-11):(y*12)],multiyear_threshold),retrieve_pheno_fix)
    dim(pheno)<-c(4,xydim)
    export_nc(pheno,paste("./pheno_CMIP5/",mod,"/pheno/pheno_",y+1984,".nc",sep=""),xdim,ydim)
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
  export_nc(pheno,paste("./pheno_CMIP5/",mod,"/pheno_mean.nc",sep=""),xdim,ydim)
}

# first get annual smoothed VI and get threshold
for (i in 1:8){
  retrieve_sopeop(models[i])
}
stopCluster(cl)

