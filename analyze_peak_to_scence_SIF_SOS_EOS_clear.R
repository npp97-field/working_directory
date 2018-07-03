#############
### calcuate the peak to senescence mean SIF or NDVI for each year
library(ncdf4)
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")
#setwd("/Users/yzhang/Project/SIF_phenology/")
## readin the SIF files
sif_files<-list.files("/rigel/glab/users/zy2309/DATA/clear_inst_SIF_4day_HD_BISE/", full.names = T,pattern=".nc")
sif_pheno<-list.files("./pheno_hd_fixed_threshold_clear/",full.names = T,pattern = ".nc")
## readin the SIF pos eos dates.
sos_file<-'./analysis/clear_daily_SOS_30N_fixed_stat.nc'
pos_file<-"./analysis/clear_daily_POS_30N_fixed_stat.nc"
eos_file<-"./analysis/clear_daily_EOS_30N_fixed_stat.nc"

##### get average SOS, EOS and POS
sos_f<-nc_open(sos_file)
sos<-round(ncvar_get(sos_f,varid="MEAN")*92+0.5)
nc_close(sos_f)
eos_f<-nc_open(eos_file)
eos<-round(ncvar_get(eos_f,varid="MEAN")*92+0.5)
nc_close(eos_f)
pos_f<-nc_open(pos_file)
pos<-round(ncvar_get(pos_f,varid="MEAN")*92+0.5)
nc_close(pos_f)
dim(sos)<-c(86400,1)
dim(eos)<-c(86400,1)
dim(pos)<-c(86400,1)

get_mean_peak2sene<-function(xts){
  if (is.na(xts[93]*xts[94]))
    return(NA)
  else{
    st<-max(xts[93],1)
    end<-min(xts[94],92)
    return(mean(xts[st:end],na.rm=T))
  }
}

### for each year, read in the sif_files
get_late_growingseason_SIF<-function(year){
  year_csif_f<-sif_files[substr(basename(sif_files),16,19)==year]
  csifin<-nc_open(year_csif_f)
  csif_year<-ncvar_get(csifin,varid="clear_daily_sif")[,241:360,]
  nc_close(csifin)
  dim(csif_year)<-c(86400,92)
  mean_csif_peak2sene<-apply(cbind(csif_year,pos,eos),1,get_mean_peak2sene)
  dim(mean_csif_peak2sene)<-c(720,120)
  
  ####
  mean_csif_start2peak<-apply(cbind(csif_year,sos,pos),1,get_mean_peak2sene)
  dim(mean_csif_start2peak)<-c(720,120)
  
  ####
  mean_csif_start2sene<-apply(cbind(csif_year,sos,eos),1,get_mean_peak2sene)
  dim(mean_csif_start2sene)<-c(720,120)
  
  ####
  mean_csif_annual<-apply(csif_year,1,sum,na.rm=T)/94
  dim(mean_csif_annual)<-c(720,120)
  
  ####
  mean_csif_peak<-apply(cbind(csif_year,pos-3,pos+3),1,get_mean_peak2sene)
  dim(mean_csif_peak)<-c(720,120)
  
  year_pheno<-sif_pheno[substr(basename(sif_pheno),35,38)==year]
  phenoin<-nc_open(year_pheno,write=T)
  xdim2<-phenoin$dim[["longitude"]]
  ydim2<-phenoin$dim[["latitude"]]
  
  p2s_csif<-ncvar_def("p2s_csif",'',list(xdim2,ydim2),-9999,compression=9)
  s2p_csif<-ncvar_def("s2p_csif",'',list(xdim2,ydim2),-9999,compression=9)
  s2s_csif<-ncvar_def("s2s_csif",'',list(xdim2,ydim2),-9999,compression=9)
  ann_csif<-ncvar_def("ann_csif",'',list(xdim2,ydim2),-9999,compression=9)
  p_csif<-ncvar_def("p_csif",'',list(xdim2,ydim2),-9999,compression=9)
  tryCatch({
    ncvar_add(phenoin,v = p2s_csif)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    ncvar_add(phenoin,v = s2p_csif)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    ncvar_add(phenoin,v = s2s_csif)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    ncvar_add(phenoin,v = ann_csif)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    ncvar_add(phenoin,v = p_csif)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  nc_close(phenoin)
  phenoin<-nc_open(year_pheno,write=T)
  ncvar_put(phenoin,p2s_csif,mean_csif_peak2sene)
  ncvar_put(phenoin,s2p_csif,mean_csif_start2peak)
  ncvar_put(phenoin,s2s_csif,mean_csif_start2sene)
  ncvar_put(phenoin,ann_csif,mean_csif_annual)
  ncvar_put(phenoin,p_csif,mean_csif_peak)
  nc_close(phenoin)
}

for (year in 2001:2016){
  print(year)
  get_late_growingseason_SIF(year)
}

################################
### calcualte the climate data for two growing periods
# 
### calcuate the peak to senescence mean SIF or NDVI for each year
library(ncdf4)
#setwd("/Users/yzhang/Project/SIF_phenology/")
## readin the SIF files
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology")
sif_pheno<-list.files("./pheno_hd_fixed_threshold_clear/",full.names = T,pattern = ".nc")

sos_file<-'./analysis/clear_daily_SOS_30N_fixed_stat.nc'
pos_file<-"./analysis/clear_daily_POS_30N_fixed_stat.nc"
eos_file<-"./analysis/clear_daily_EOS_30N_fixed_stat.nc"

##### get average EOS and SOS
eos_f<-nc_open(eos_file)
eos<-ncvar_get(eos_f,varid="MEAN")
xdim2<-eos_f$dim[["longitude"]]
ydim2<-eos_f$dim[["latitude"]]
nc_close(eos_f)
sos_f<-nc_open(sos_file)
sos<-ncvar_get(sos_f,varid="MEAN")
nc_close(sos_f)
pos_f<-nc_open(pos_file)
pos<-ncvar_get(pos_f,varid="MEAN")
nc_close(pos_f)

dim(eos)<-c(86400,1)
dim(sos)<-c(86400,1)
dim(pos)<-c(86400,1)

get_cru_var<-function(var_name){
  #get the monthly climate data from 2000 to 2016
  var_f<-list.files('/rigel/glab/users/zy2309/DATA/CRU_TS401/',full.names = T,pattern = paste(var_name,".dat.nc",sep=""))
  var_in1<-nc_open(var_f[1])
  var1<-ncvar_get(var_in1,var_name)[,241:360,109:120]
  nc_close(var_in1)
  dim(var1)<-c(86400,12)
  
  var_in2<-nc_open(var_f[2])
  var2<-ncvar_get(var_in2,var_name)[,241:360,]
  nc_close(var_in2)
  dim(var2)<-c(86400,120)
  
  var_in3<-nc_open(var_f[3])
  var3<-ncvar_get(var_in3,var_name)[,241:360,]
  nc_close(var_in3)
  dim(var3)<-c(86400,72)
  var<-cbind(var1,var2,var3)
  return(var)
}

#### interpolate monthly temperature to daily
spline_interpolate<-function(x){
  nyears<-length(x)/12
  if (sum(is.na(x))>100){
    return(rep(NA,nyears*365))
  }
  
  x_in<-(1:length(x)-0.5)/12
  x_out<-(1:(365*nyears)-0.5)/365
  spl_temp<-spline(x=x_in,y=x,xout=x_out)$y
  return(spl_temp)
}


tmean<-get_cru_var("tmp")
precip<-get_cru_var('pre')

temp_d<-apply(tmean,1,spline_interpolate)
temp_day<-t(temp_d)

radnc<-nc_open("/rigel/glab/users/zy2309/DATA/bess_hd_monthly_PAR.nc")
rad<-ncvar_get(radnc,"PAR")
par<-rad[,,241:360]
dim(par)<-c(204,86400)
par_d<-apply(t(par),1,spline_interpolate)
par_day<-t(par_d)

get_total<-function(xts){
  ### xts  first 13 is the variable for previous Dec to this Dec (13 in total)
  ### xts  the 14 and 15 are the sos and eos
  if (is.na(xts[14]*xts[15]))
    return(NA)
  else{
    st_abs<-xts[14]*12
    end_abs<-xts[15]*12
    st<-ceiling(st_abs)+1
    end<-ceiling(end_abs)+1
    if (st==end){
      sum_var = xts[st]*(end_abs-st_abs)
    }else if(end==st+1){
      sum_var= xts[st]*(1-st_abs%%1)+xts[end]*(end_abs%%1)
    }else{
      sum_var<-sum(xts[(st+1):(end-1)],na.rm=T)+xts[st]*(1-st_abs%%1)+xts[end]*(end_abs%%1)
    }
    return(sum_var)
  }
}

get_mean<-function(xts){
  #this is for temperature, modified into daily
  if (is.na(xts[397]*xts[398]))
    return(NA)
  else{
    st_abs<-round(xts[397]*365)+30
    end_abs<-round(xts[398]*365)+30
    ave_var<-mean(xts[st_abs:end_abs],na.rm=T)
    return(ave_var)
  }
}


calculate_climate<-function(year){
  #get precipitation for this year +one month from previous year
  precip_year<-precip[,((year-2000)*12):((year-1999)*12)]
  tmean_year<-temp_day[,((year-2000)*365-30):((year-1999)*365)]
  par_year<-par_day[,((year-2000)*365-30):((year-1999)*365)]
  
  prec_dat<-cbind(precip_year,sos,eos)
  s2e_precip<-apply(prec_dat,1,get_total)
  
  temp_dat<-cbind(tmean_year,sos,eos)
  s2e_temper<-apply(temp_dat,1,get_mean)
  
  par_dat<-cbind(par_year,sos,eos)
  s2e_parbess<-apply(par_dat,1,get_mean)
  ##################################
  prec_dat<-cbind(precip_year,sos,pos)
  s2p_precip<-apply(prec_dat,1,get_total)
  
  temp_dat<-cbind(tmean_year,sos,pos)
  s2p_temper<-apply(temp_dat,1,get_mean)
  
  prec_dat<-cbind(par_year,sos,pos)
  s2p_parbess<-apply(par_dat,1,get_mean)
  ## calculate the second half temp and precipitation
  temp_dat<-cbind(tmean_year,pos,eos)
  p2e_temper<-apply(temp_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,pos,eos)
  p2e_precip<-apply(prec_dat,1,get_total)
  
  par_dat<-cbind(par_year,pos,eos)
  p2e_parbess<-apply(par_dat,1,get_mean)
  ##################################
  ##----------------------------------------------
  ### calculate the average 1 month pre-season temperature
  temp_dat<-cbind(tmean_year,sos-1/12,sos)
  presos1mon_temp<-apply(temp_dat,1,get_mean)
  
  temp_dat<-cbind(tmean_year,eos-1/12,eos)
  preeos1mon_temp<-apply(temp_dat,1,get_mean)
  
  temp_dat<-cbind(tmean_year,pos-1/24,pos+1/24)
  peak_temp<-apply(temp_dat,1,get_mean)
  
  par_dat<-cbind(par_year,sos-1/12,sos)
  presos1mon_par<-apply(par_dat,1,get_mean)
  
  par_dat<-cbind(par_year,eos-1/12,eos)
  preeos1mon_par<-apply(par_dat,1,get_mean)
  
  par_dat<-cbind(par_year,pos-1/24,pos+1/24)
  peak_par<-apply(par_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,sos-1/12,sos)
  presos1mon_prec<-apply(prec_dat,1,get_total)
  
  prec_dat<-cbind(precip_year,eos-1/12,eos)
  preeos1mon_prec<-apply(prec_dat,1,get_total)
  
  prec_dat<-cbind(precip_year,pos-1/24,pos+1/24)
  peak_prec<-apply(prec_dat,1,get_total)
  ##************************************************
  ### write data to pheno
  #year_pheno<-sif_pheno[substr(basename(sif_pheno),33,36)==year]
  #phenoin<-nc_open(year_pheno)
  
  
  s2e_prec<-ncvar_def("s2e_prec",'',list(xdim2,ydim2),-9999,compression=9)
  p2e_prec<-ncvar_def("p2e_prec",'',list(xdim2,ydim2),-9999,compression=9)
  s2p_prec<-ncvar_def("s2p_prec",'',list(xdim2,ydim2),-9999,compression=9)
  s2e_temp<-ncvar_def("s2e_temp",'',list(xdim2,ydim2),-9999,compression=9)
  s2p_temp<-ncvar_def("s2p_temp",'',list(xdim2,ydim2),-9999,compression=9)
  p2e_temp<-ncvar_def("p2e_temp",'',list(xdim2,ydim2),-9999,compression=9)
  s2e_par<-ncvar_def("s2e_par",'',list(xdim2,ydim2),-9999,compression=9)
  s2p_par<-ncvar_def("s2p_par",'',list(xdim2,ydim2),-9999,compression=9)
  p2e_par<-ncvar_def("p2e_par",'',list(xdim2,ydim2),-9999,compression=9)
  p_temp<-ncvar_def("p_temp","",list(xdim2,ydim2),-9999,compression=9)
  pre_start_temp<-ncvar_def("pre_start_temp",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end_temp<-ncvar_def("pre_end_temp",'',list(xdim2,ydim2),-9999,compression=9)
  p_par<-ncvar_def("p_par","",list(xdim2,ydim2),-9999,compression=9)
  pre_start_par<-ncvar_def("pre_start_par",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end_par<-ncvar_def("pre_end_par",'',list(xdim2,ydim2),-9999,compression=9)
  p_prec<-ncvar_def("p_prec","",list(xdim2,ydim2),-9999,compression=9)
  pre_start_prec<-ncvar_def("pre_start_prec",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end_prec<-ncvar_def("pre_end_prec",'',list(xdim2,ydim2),-9999,compression=9)
  #tryCatch({
  #   ncvar_add(phenoin,v = grow_precip)
  #   ncvar_add(phenoin,v = grow_temp)
  # },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # tryCatch({
  #   ncvar_add(phenoin,v = start2peak_precip)
  #   ncvar_add(phenoin,v = start2peak_temp)
  # },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # 
  ##----------------------------------------------
  
  
  # tryCatch({
  #   ncvar_add(phenoin,v = pre_start_temp)
  #   ncvar_add(phenoin,v = pre_end_temp)
  #   ncvar_add(phenoin,v = peak2end_temp)
  #   ncvar_add(phenoin,v = peak2end_prec)
  # },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  ##************************************************
  #nc_close(phenoin)
  year_climate<-paste("./pheno_hd_fixed_threshold_climate/clear/",year,"_climate_csif_clear_daily.nc",sep="")
  if (file.exists(year_climate)){
    file.remove(year_climate)
  }
  phenoin<-nc_create(year_climate,list(s2e_prec,s2e_temp,s2p_prec,s2p_temp,p2e_temp,p2e_prec,
                                       pre_start_temp,pre_end_temp,p_temp,s2e_par,s2p_par,p2e_par,
                                       pre_start_par,pre_end_par,p_par,pre_start_prec,pre_end_prec,p_prec))
  ncvar_put(phenoin,s2e_prec,s2e_precip)
  ncvar_put(phenoin,p2e_prec,p2e_precip)
  ncvar_put(phenoin,s2p_prec,s2p_precip)
  ncvar_put(phenoin,s2e_temp,s2e_temper)
  ncvar_put(phenoin,s2p_temp,s2p_temper)
  ncvar_put(phenoin,p2e_temp,p2e_temper)
  ncvar_put(phenoin,s2e_par,s2e_parbess)
  ncvar_put(phenoin,s2p_par,s2p_parbess)
  ncvar_put(phenoin,p2e_par,p2e_parbess)
  ##----------------------------------------------
  ncvar_put(phenoin,pre_start_temp,presos1mon_temp)
  ncvar_put(phenoin,pre_end_temp,preeos1mon_temp)
  ncvar_put(phenoin,p_temp,peak_temp)
  ncvar_put(phenoin,pre_start_par,presos1mon_par)
  ncvar_put(phenoin,pre_end_par,preeos1mon_par)
  ncvar_put(phenoin,p_par,peak_par)
  ncvar_put(phenoin,pre_start_prec,presos1mon_prec)
  ncvar_put(phenoin,pre_end_prec,preeos1mon_prec)
  ncvar_put(phenoin,p_prec,peak_prec)
  ##************************************************
  nc_close(phenoin)
}

for (year in 2001:2016){
  print(year)
  calculate_climate(year)
}


#####################################
#######   analysis the partial corrlation between cSIF growing season with SOS when Tmean and precip fixed.
#####################################

## get the average of SOS, CSIF, Tmean and Precip
library(ppcor)
library(ncdf4)
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")
sif_pheno<-list.files("./pheno_hd_fixed_threshold_clear/",full.names = T,pattern = ".nc")
sif_climate<-list.files("./pheno_hd_fixed_threshold_climate/clear/",full.names = T,pattern = ".nc")

p2s_csif<-array(NA,dim=c(86400,16))
s2p_csif<-array(NA,dim=c(86400,16))
s2s_csif<-array(NA,dim=c(86400,16))
ann_csif<-array(NA,dim=c(86400,16))
sos<-array(NA,dim=c(86400,16))
eos<-array(NA,dim=c(86400,16))
s2e_temp<-array(NA,dim=c(86400,16))
s2e_prec<-array(NA,dim=c(86400,16))
s2e_par<-array(NA,dim=c(86400,16))
s2p_temp<-array(NA,dim=c(86400,16))
s2p_prec<-array(NA,dim=c(86400,16))
s2p_par<-array(NA,dim=c(86400,16))
p2e_temp<-array(NA,dim=c(86400,16))
p2e_prec<-array(NA,dim=c(86400,16))
p2e_par<-array(NA,dim=c(86400,16))
pre_start_temp<-array(NA,dim=c(86400,16))
pre_end_temp<-array(NA,dim=c(86400,16))
p_temp<-array(NA,dim=c(86400,16))
pre_start_par<-array(NA,dim=c(86400,16))
pre_end_par<-array(NA,dim=c(86400,16))
p_par<-array(NA,dim=c(86400,16))
pre_start_prec<-array(NA,dim=c(86400,16))
pre_end_prec<-array(NA,dim=c(86400,16))
p_prec<-array(NA,dim=c(86400,16))


for (year in 2001:2016){
  year_pheno<-sif_pheno[substr(basename(sif_pheno),35,38)==year]
  phenoin<-nc_open(year_pheno)
  p2s_csif[,year-2000]<-ncvar_get(phenoin,'p2s_csif')
  s2p_csif[,year-2000]<-ncvar_get(phenoin,'s2p_csif')
  s2s_csif[,year-2000]<-ncvar_get(phenoin,'s2s_csif')
  ann_csif[,year-2000]<-ncvar_get(phenoin,'ann_csif')
  sos[,year-2000]<-ncvar_get(phenoin,'SOS')*365+2.5
  eos[,year-2000]<-ncvar_get(phenoin,'EOS')*365+2.5
  xdim<-phenoin$dim[["longitude"]]
  ydim<-phenoin$dim[["latitude"]]
  nc_close(phenoin)
}


for (year in 2001:2016){
  year_climate<-sif_climate[substr(basename(sif_climate),1,4)==year]
  climatein<-nc_open(year_climate)
  
  s2e_temp[,year-2000]<-ncvar_get(climatein,'s2e_temp')
  s2e_prec[,year-2000]<-ncvar_get(climatein,'s2e_prec')
  s2e_par[,year-2000]<-ncvar_get(climatein,'s2e_par')
  s2p_temp[,year-2000]<-ncvar_get(climatein,'s2p_temp')
  s2p_prec[,year-2000]<-ncvar_get(climatein,'s2p_prec')
  s2p_par[,year-2000]<-ncvar_get(climatein,'s2p_par')
  p2e_temp[,year-2000]<-ncvar_get(climatein,'p2e_temp')
  p2e_prec[,year-2000]<-ncvar_get(climatein,'p2e_prec')
  p2e_par[,year-2000]<-ncvar_get(climatein,'p2e_par')
  p_temp[,year-2000]<-ncvar_get(climatein,'p_temp')
  pre_start_temp[,year-2000]<-ncvar_get(climatein,'pre_start_temp')
  pre_end_temp[,year-2000]<-ncvar_get(climatein,'pre_end_temp')
  p_par[,year-2000]<-ncvar_get(climatein,'p_par')
  pre_start_par[,year-2000]<-ncvar_get(climatein,'pre_start_par')
  pre_end_par[,year-2000]<-ncvar_get(climatein,'pre_end_par')
  
  p_prec[,year-2000]<-ncvar_get(climatein,'p_prec')
  pre_start_prec[,year-2000]<-ncvar_get(climatein,'pre_start_prec')
  pre_end_prec[,year-2000]<-ncvar_get(climatein,'pre_end_prec')
  nc_close(climatein)
}

calculate_partial_correlation<-function(dat){
  if (sum(is.na(dat))>=1)
    return(c(NA,NA))
  n<-length(dat)
  dim(dat)<-c(16,n/16)
  pcor_re<-pcor.test(dat[,1],dat[,2],dat[,3:(n/16)],method="pearson")
  return(c(pcor_re$estimate,pcor_re$p.value))
}

calculate_correlation<-function(dat){
  if (sum(is.na(dat))>=1)
    return(c(NA,NA))
  dim(dat)<-c(16,2)
  cor_re<-cor.test(dat[,1],dat[,2],method="pearson")
  return(c(cor_re$estimate,cor_re$p.value))
}

##### calculate correlation
cal_cor<-function(data,fileout){
  data[data< -990]<-NA
  cor_data<-apply(data,1,calculate_correlation)
  coef<-cor_data[1,]
  coef[is.nan(coef)]<- -999.9
  dim(coef)<-c(720,120)
  pv<-cor_data[2,]
  pv[is.nan(pv)]<- -999.9
  dim(pv)<-c(720,120)
  
  cor_coef<-ncvar_def("cor_coef",'',list(xdim,ydim),-999.9,prec="double",compression=9)
  cor_pv<-ncvar_def("cor_pv",'',list(xdim,ydim),-999.9,prec="double",compression=9)
  ncout<-nc_create(fileout,list(cor_coef,cor_pv))
  ncvar_put(ncout,cor_coef,coef)
  ncvar_put(ncout,cor_pv,pv)
  nc_close(ncout)
}

##### calculate partial correlation
cal_pcor<-function(data,fileout){
  data[data< -990]<-NA
  pcor_data<-apply(data,1,calculate_partial_correlation)
  pcoef<-pcor_data[1,]
  pcoef[is.nan(pcoef)]<- -999.9
  dim(pcoef)<-c(720,120)
  ppv<-pcor_data[2,]
  ppv[is.nan(ppv)]<- -999.9
  dim(ppv)<-c(720,120)
  
  pcor_coef<-ncvar_def("pcor_coef",'',list(xdim,ydim),-999.9,prec="double",compression=9)
  pcor_pv<-ncvar_def("pcor_pv",'',list(xdim,ydim),-999.9,prec="double",compression=9)
  ncout<-nc_create(fileout,list(pcor_coef,pcor_pv))
  ncvar_put(ncout,pcor_coef,pcoef)
  ncvar_put(ncout,pcor_pv,ppv)
  nc_close(ncout)
}


###### calculate the SOS effect on SIF EOS, etc
if (T){
  ##########calculate the pcor between sos and csif_P2S
  
  eos_csif<-cbind(sos,p2s_csif,s2e_temp,s2e_prec,s2e_par)
  fileout<-"./analysis/correlation_clear/pcor_sos_csif_p2s.nc"
  cal_pcor(eos_csif,fileout)
  ###########################
  eos_csif<-cbind(sos,p2s_csif)
  fileout<-"./analysis/correlation_clear/cor_sos_csif_p2s.nc"
  cal_cor(eos_csif,fileout)
  
  ##########calculate the pcor between sos and csif_S2P
  
  sos_csif<-cbind(sos,s2p_csif,s2p_temp,s2p_prec,s2p_par)
  nc_out_f2<-"./analysis/correlation_clear/pcor_sos_csif_s2p.nc"
  cal_pcor(sos_csif,nc_out_f2)
  
  #######################
  sos_csif<-cbind(sos,s2p_csif)
  nc_out_f2<-"./analysis/correlation_clear/cor_sos_csif_s2p.nc"
  cal_cor(sos_csif,nc_out_f2)
  
  ##########calculate the pcor between sos and csif_S2S
  
  sos_csif<-cbind(sos,s2s_csif,s2e_temp,s2e_prec)
  nc_out_f2<-"./analysis/correlation_clear/pcor_sos_csif_s2s.nc"
  cal_pcor(sos_csif,nc_out_f2)
  
  
  #######################
  sos_csif<-cbind(sos,s2s_csif)
  nc_out_f2<-"./analysis/correlation_clear/cor_sos_csif_s2s.nc"
  cal_cor(sos_csif,nc_out_f2)
  
  
  ##########calculate the pcor between sos and eos
  sos_eos<-cbind(sos,eos,s2e_temp,s2e_prec,s2e_par)
  nc_out_f3<-"./analysis/correlation_clear/pcor_sos_eos.nc"
  cal_pcor(sos_eos,nc_out_f3)
  
  ####################################
  sos_eos<-cbind(sos,eos)
  nc_out_f3<-"./analysis/correlation_clear/cor_sos_eos.nc"
  cal_cor(sos_eos,nc_out_f3)
  
  
}



###### calculate the controlling factors of SOS EOS
if (T){
  ##########calculate the cor between sos and eos and pre season temp
  sos_eos<-cbind(sos,pre_start_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_sos_pre_temp.nc"
  cal_cor(sos_eos,nc_out_f3)
  
  ####################################
  sos_eos<-cbind(eos,pre_end_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_eos_pre_temp.nc"
  cal_cor(sos_eos,nc_out_f3)
  
  ############ sos and late season temperature
  sos_eos<-cbind(sos,pre_end_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_sos_prefall_temp.nc"
  cal_cor(sos_eos,nc_out_f3)
  
  ####################################
  sos_eos<-cbind(sos,p2e_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_sos_p2e_temp.nc"
  cal_cor(sos_eos,nc_out_f3)
  
  ####################################
  sos_eos<-cbind(sos,s2p_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_sos_s2p_temp.nc"
  cal_cor(sos_eos,nc_out_f3)
  
}


###### temperature feedbacks
if (T){
  #################################### spring temperature and fall temperature
  pre_start_pre_end<-cbind(pre_start_temp,pre_end_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_pre_start_pre_end_temp.nc"
  cal_cor(pre_start_pre_end,nc_out_f3)
  
  s2p_p2e<-cbind(s2p_temp,p2e_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_s2p_p2e_temp.nc"
  cal_cor(s2p_p2e,nc_out_f3)
  
  
  ############################### growing season temp precip and growing season CSIF
  s2s_sif_temp<-cbind(s2s_csif,s2e_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_s2s_sif_s2s_temp.nc"
  cal_cor(s2s_sif_temp,nc_out_f3)
  
  
  s2s_sif_prec<-cbind(s2s_csif,s2e_prec)
  nc_out_f3<-"./analysis/correlation_clear/cor_s2s_sif_s2s_precip.nc"
  cal_cor(s2s_sif_prec,nc_out_f3)
  
  ############################### pre EOS and peak2end temp
  s2s_sif_prec<-cbind(pre_end_temp,p2e_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_preEOS_p2e_temp.nc"
  cal_cor(s2s_sif_prec,nc_out_f3)
  
  ############################### pre SOS and start2peak temp
  s2s_sif_prec<-cbind(pre_start_temp,s2p_temp)
  nc_out_f3<-"./analysis/correlation_clear/cor_preSOS_s2p_temp.nc"
  cal_cor(s2s_sif_prec,nc_out_f3)
  
  ############################### s2p sif and p2e sif
  s2s_sif_prec<-cbind(s2p_csif,p2s_csif)
  nc_out_f3<-"./analysis/correlation_clear/cor_s2p_p2s_csif.nc"
  cal_cor(s2s_sif_prec,nc_out_f3)
  
}

###### climate inter-correlation
if (T){
  var_dat<-cbind(pre_start_temp,pre_start_prec)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_pre_start_prec_pre_temp.nc"
  cal_cor(var_dat,fileout)
  
  var_dat<-cbind(pre_start_temp,pre_start_par)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_pre_start_par_pre_temp.nc"
  cal_cor(var_dat,fileout)
  
  var_dat<-cbind(pre_start_par,pre_start_prec)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_pre_start_par_pre_prec.nc"
  cal_cor(var_dat,fileout)
  
  var_dat<-cbind(pre_end_temp,pre_end_prec)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_pre_end_prec_pre_temp.nc"
  cal_cor(var_dat,fileout)
  
  var_dat<-cbind(pre_end_temp,pre_end_par)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_pre_end_par_pre_temp.nc"
  cal_cor(var_dat,fileout)
  
  var_dat<-cbind(pre_end_par,pre_end_prec)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_pre_end_par_pre_prec.nc"
  cal_cor(var_dat,fileout)
  
  var_dat<-cbind(s2e_temp,s2e_prec)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_s2e_prec_temp.nc"
  cal_cor(var_dat,fileout)
  
  var_dat<-cbind(s2e_temp,s2e_par)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_s2e_par_temp.nc"
  cal_cor(var_dat,fileout)
  
  var_dat<-cbind(s2e_par,s2e_prec)
  fileout<-"./analysis/correlation_clear/climate_intercor/cor_s2e_par_prec.nc"
  cal_cor(var_dat,fileout)
  
}
###### radiation feedback
sos_eos<-cbind(sos,p2e_par)
nc_out_f3<-"./analysis/correlation_clear/cor_sos_p2e_par.nc"
cal_cor(sos_eos,nc_out_f3)

sos_eos<-cbind(eos,p2e_par)
nc_out_f3<-"./analysis/correlation_clear/cor_eos_p2e_par.nc"
cal_cor(sos_eos,nc_out_f3)

sos_eos<-cbind(p2s_csif,p2e_par)
nc_out_f3<-"./analysis/correlation_clear/cor_p2e_sif_par.nc"
cal_cor(sos_eos,nc_out_f3)


sos_eos<-cbind(eos,p2e_par,p2e_temp,p2e_prec)
nc_out_f3<-"./analysis/correlation_clear/pcor_eos_p2e_par.nc"
cal_pcor(sos_eos,nc_out_f3)

sos_eos<-cbind(p2s_csif,p2e_par,p2e_temp,p2e_prec)
nc_out_f3<-"./analysis/correlation_clear/pcor_p2e_sif_par.nc"
cal_pcor(sos_eos,nc_out_f3)





