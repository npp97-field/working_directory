##### get (1) average, (2) standard deviation (3) trend
##### of SOS EOS, and threshold 
library(trend)
library(ncdf4)

### trend calculation
mk_test<-function(x){
  x[is.nan(x)]<-NA
  if(sum(is.na(x))>0){
    return(NA)
  }
  return(mk.test(x)$p.value)
}

###sen's slope
senslope<-function(x){
  x[is.nan(x)]<-NA
  if(sum(is.na(x))>0){
    return(NA)
  }
  return(sens.slope(x)$estimate)
}


calculate_stat<-function(dat,outfile){
  mean_dat<-apply(dat,1,mean,na.rm=T)
  std_dat<-apply(dat,1,sd,na.rm=T)
  sen_trend<-apply(dat,1,senslope)
  sig_trend<-apply(dat,1,mk_test)
  sif_analysis<-cbind(mean_dat,std_dat,sen_trend,sig_trend)
  dim(sif_analysis)<-c(720,120,4)
  export_nc(sif_analysis,outfile)
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
  ncmean<-ncvar_def('MEAN','NA',list(dimlong,dimlat),-9999,longname="mean",prec='float',compression=9)
  ncsd<-ncvar_def('SD','NA',list(dimlong,dimlat),-9999,longname="standard deviation",prec='float',compression=9)
  nctrend<-ncvar_def('TREND','NA',list(dimlong,dimlat),-9999,longname="sen_slope_trend",prec='float',compression=9)
  ncsig<-ncvar_def('SIG','NA',list(dimlong,dimlat),-9999,longname="trend significance",prec='float',compression=9)
  
  ncout<-nc_create(outfile,list(ncmean,ncsd,nctrend,ncsig))
  ncvar_put(ncout,varid=ncmean,pheno_dat[,,1])
  ncvar_put(ncout,varid=ncsd,pheno_dat[,,2])
  ncvar_put(ncout,varid=nctrend,pheno_dat[,,3])
  ncvar_put(ncout,varid=ncsig,pheno_dat[,,4])
  nc_close(ncout)
}

########################

setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")
indicator<-c('var','fixed')

stat_dataset<-function(indi){
  pheno_files<-list.files(paste("./pheno_hd_",indi,"_threshold/",sep=""),pattern=".nc",full.names = T)
  sos<-array(NA,dim=c(86400,14))
  eos<-array(NA,dim=c(86400,14))
  lgs<-array(NA,dim=c(86400,14))
  thresh<-array(NA,dim=c(86400,14))
  
  for (i in 1:length(pheno_files)){
    ncf<-nc_open(pheno_files[i])
    ncsos<-ncvar_get(ncf,varid = "SOS")
    nceos<-ncvar_get(ncf,varid = "EOS")
    ncthresh<-ncvar_get(ncf,varid = "THESHOLD")
    nc_close(ncf)
    sos[,i]<-ncsos
    eos[,i]<-nceos
    thresh[,i]<-ncthresh
  }
  lgs<-eos-sos
  calculate_stat(sos,paste("./analysis/SOS_30N_",indi,"_stat.nc",sep=""))
  calculate_stat(eos,paste("./analysis/EOS_30N_",indi,"_stat.nc",sep=""))
  calculate_stat(lgs,paste("./analysis/LGS_30N_",indi,"_stat.nc",sep=""))
  calculate_stat(thresh,paste("./analysis/THRESH_30N_",indi,"_stat.nc",sep=""))
}

#### variable threshold for each year

stat_dataset(indicator[1])

#### fixed threshold for each year

stat_dataset(indicator[2])








