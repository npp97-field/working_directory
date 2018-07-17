##### calculate CMIP5 pre season climate 
# 
### calcuate the peak to senescence mean SIF or NDVI for each year
library(ncdf4)
## readin the CMIP5 phenology
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology")

models<-c("CanESM2","CMCC-CESM","GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR",
          "MIROC-ESM","MPI-ESM2-LR","NorESM1-ME")

cmip_pheno<-list.files("./pheno_CMIP5/",full.names = T,recursive = T,pattern = "pheno_mean.nc")

cmip_pr<-list.files("./CMIP5/historical/",pattern="pr_",full.names = T)
cmip_tas<-list.files("./CMIP5/historical/",pattern="tas_",full.names = T)
cmip_tasmax<-list.files("./CMIP5/historical/",pattern="tasmax_",full.names = T)


#### interpolate monthly obs to daily
spline_interpolate<-function(x){
  nyears<-length(x)/12
  if (sum(is.na(x))>100){
    return(rep(NA,365*nyears))
  }
  x_in<-(1:length(x)-0.5)/12
  x_out<-(1:(365*nyears)-0.5)/365
  spl_temp<-spline(x=x_in,y=x,xout=x_out)$y
  return(spl_temp)
}
#### get preseason climate variables
get_pre_var<-function(var){
  eos<-var[366]
  if (is.na(eos)){
    return(rep(NA,4))
  }
  end<-round(eos*365)
  pre0<-mean(var[max(1,end-15):end],na.rm=T)
  pre1<-mean(var[max(1,end-30):end],na.rm=T)
  pre2<-mean(var[max(1,end-60):end],na.rm=T)
  pre3<-mean(var[max(1,end-90):end],na.rm=T)
  return(c(pre0,pre1,pre2,pre3))
}

model_pre_var<-function(var,eos,xydim,i){
  var_y<-var[,,(i*12-11):(i*12)]
  dim(var_y)<-c(xydim[1]*xydim[2],12)
  var_d<-apply(var_y,1,spline_interpolate)
  #####  get the pre eos with 15, 30, 60, 90days
  comdat<-rbind(var_d,as.vector(eos))
  dat<-apply(comdat,2,get_pre_var)
  dim(dat)<-c(4,xydim[1],xydim[2])
  return(dat)
}

export_nc<-function(prdat,tasdat,tasmaxdat,outfile,xdim,ydim){
  zdim<-ncdim_def("lags","",vals = 0:3,longname = "months of lags")
  ncpr_pre<-ncvar_def('pre_pr','mm day-1',list(zdim,xdim,ydim),-9999,
                      longname="pre eos precipitation 0-3 mon",prec='float',compression=9)
  nctas_pre<-ncvar_def('pre_tas','deg C',list(zdim,xdim,ydim),-9999,
                       longname="pre eos air temperature 0-3 mon",prec='float',compression=9)
  nctasmax_pre<-ncvar_def('pre_tasmax','deg C',list(zdim,xdim,ydim),-9999,
                          longname="pre eos max air temperature 0-3 mon",prec='float',compression=9)
  
  ncout<-nc_create(outfile,list(ncpr_pre,nctas_pre,nctasmax_pre))
  ncvar_put(ncout,varid=ncpr_pre,prdat)
  ncvar_put(ncout,varid=nctas_pre,tasdat)
  ncvar_put(ncout,varid=nctasmax_pre,tasmaxdat)
  nc_close(ncout)
}


get_preseason_climate<-function(mod){
  ncin<-nc_open(cmip_pheno[grep(mod,cmip_pheno)])
  #sop<-ncvar_get(ncin,varid="SOS")
  eos<-ncvar_get(ncin,varid="EOS")
  xydim<-dim(eos)[1:2]
  xdim<-ncin$dim[["lon"]]
  ydim<-ncin$dim[["lat"]]
  nc_close(ncin)
  ncin<-nc_open(cmip_pr[grep(mod,cmip_pr)])
  pr<-ncvar_get(ncin,varid="pr")
  nc_close(ncin)
  ncin<-nc_open(cmip_tas[grep(mod,cmip_tas)])
  tas<-ncvar_get(ncin,varid="tas")
  nc_close(ncin)
  tasmax_f<-cmip_tasmax[grep(mod,cmip_tasmax)]
  if(length(tasmax_f)==0){
    tasmax<-tas
  }else{
    ncin<-nc_open(tasmax_f)
    tasmax<-ncvar_get(ncin,varid="tasmax")
    nc_close(ncin)
  }
  dir.create(paste("./pheno_CMIP5_climate/",mod,sep=""))
  #### for each year, the climate variables is interpolated to daily values.
  for (i in 1:20){
    prdat<-model_pre_var(pr,eos,xydim,i)
    tasdat<-model_pre_var(tas,eos,xydim,i)
    tasmaxdat<-model_pre_var(tasmax,eos,xydim,i)
    outfilename<-paste("./pheno_CMIP5_climate/",mod,"/",i+1984,"_pre_eos_climate.nc",sep="")
    export_nc(prdat,tasdat,tasmaxdat,outfilename,xdim,ydim)
  }
}

#####calcualte for the 8 models,
for (i in 1:8){
  get_preseason_climate(models[i])
}

#########-----------------------------------------------------
###--------- this is for the correlation analysis.
library(ncdf4)
#library(abind)
## readin the CMIP5 phenology
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology")

models <-c("CanESM2","CMCC-CESM","GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR",
           "MIROC-ESM","MPI-ESM2-LR","NorESM1-ME")

calculate_correlation<-function(ts_var){
  eos<-ts_var[1:20*5-4]
  if (sum(is.na(eos))>5){
    return(rep(NA,8))
  }
  coef<-c()
  pv<-c()
  for (i in 1:4){
    var_val<-ts_var[1:20*5-4+i]
    if (sum(is.na(var_val))>5){
      return(rep(NA,8))
    }
    cor_t<-cor.test(eos,var_val)
    coef[i]<-cor_t$estimate
    pv[i]<-cor_t$p.value
  }
  return(c(coef,pv))
}


get_eos_var_correlation<-function(varname,cmip_pheno,cmip_cli){
  for (i in 1:20){
    ncin<-nc_open(cmip_pheno[i])
    eos<-ncvar_get(ncin,var="EOS")
    nc_close(ncin)
    
    ncin<-nc_open(cmip_cli[i])
    pre_var<-ncvar_get(ncin,var=varname)
    nc_close(ncin)
    
    xydim<-dim(eos)
    dim(eos)<-c(xydim[1]*xydim[2])
    dim(pre_var)<-c(4,xydim[1]*xydim[2])
    if (i==1){
      dat_all<-rbind(eos,pre_var)
    }else{
      temp<-rbind(eos,pre_var)
      dat_all<-rbind(dat_all,temp)
    }
  }
  cor_var<-apply(dat_all,2,calculate_correlation)
  dim(cor_var)<-c(8,xydim)
  return(cor_var)
}

### function for export
export_nc<-function(pr_rel,tas_rel,tasmax_rel,outfile,xdim,ydim){
  zdim<-ncdim_def("lags","",vals = 0:3,longname = "months of lags")
  ncpr_cor<-ncvar_def('r_pr','',list(zdim,xdim,ydim),-9999,
                      longname="correlation between EOS and pre eos precipitation 0-3 mon",
                      prec='float',compression=9)
  nctas_cor<-ncvar_def('r_tas','',list(zdim,xdim,ydim),-9999,
                       longname="correlation between EOS and pre eos air temperature 0-3 mon",
                       prec='float',compression=9)
  nctasmax_cor<-ncvar_def('r_tasmax','',list(zdim,xdim,ydim),-9999,
                          longname="correlation between EOS and pre eos max air temperature 0-3 mon",
                          prec='float',compression=9)
  
  ncpr_pv<-ncvar_def('pv_pr','',list(zdim,xdim,ydim),-9999,
                      longname="p-value between EOS and pre eos precipitation 0-3 mon",
                      prec='float',compression=9)
  nctas_pv<-ncvar_def('pv_tas','',list(zdim,xdim,ydim),-9999,
                       longname="p-value between EOS and pre eos air temperature 0-3 mon",
                       prec='float',compression=9)
  nctasmax_pv<-ncvar_def('pv_tasmax','',list(zdim,xdim,ydim),-9999,
                          longname="p-value between EOS and pre eos max air temperature 0-3 mon",
                          prec='float',compression=9)
  
  ncout<-nc_create(outfile,list(ncpr_cor,nctas_cor,nctasmax_cor,ncpr_pv,nctas_pv,nctasmax_pv))
  ncvar_put(ncout,varid=ncpr_cor,pr_rel[1:4,,])
  ncvar_put(ncout,varid=ncpr_pv,pr_rel[5:8,,])
  ncvar_put(ncout,varid=nctas_cor,tas_rel[1:4,,])
  ncvar_put(ncout,varid=nctas_pv,tas_rel[5:8,,])
  ncvar_put(ncout,varid=nctasmax_cor,tasmax_rel[1:4,,])
  ncvar_put(ncout,varid=nctasmax_pv,tasmax_rel[5:8,,])
  nc_close(ncout)
}

### get the correlation for each model
cor_cal_model<-function(mod){
  cmip_pheno<-list.files(paste("./pheno_CMIP5/",mod,"/pheno",sep=""),full.names = T,recursive = T,pattern = ".nc")
  cmip_cli<-list.files(paste("./pheno_CMIP5_climate/",mod,"/",sep=""),full.names = T,recursive = T,pattern = "climate.nc")
  ## get the xdim, ydim
  nc_temp<-nc_open(cmip_pheno[1])
  xdim<-nc_temp$dim[["lon"]]
  ydim<-nc_temp$dim[["lat"]]
  nc_close(nc_temp)
  ##### get the eos date and climate eos with 0,1,2,3 months of lags
  ###  function for calculate the correaltions with lags together.
  pr_rel<-get_eos_var_correlation("pre_pr",cmip_pheno,cmip_cli)
  tas_rel<-get_eos_var_correlation("pre_tas",cmip_pheno,cmip_cli)
  tasmax_rel<-get_eos_var_correlation("pre_tasmax",cmip_pheno,cmip_cli)
  ### export the data to output
  #xydim<-dim(pr_rel)[2:3]
  outfile<-paste("./analysis/correlation_CMIP5/climate_correlation_",mod,".nc",sep="")
  export_nc(pr_rel,tas_rel,tasmax_rel,outfile,xdim,ydim)
}

for (i in 1:8){
  cor_cal_model(models[i])
}









