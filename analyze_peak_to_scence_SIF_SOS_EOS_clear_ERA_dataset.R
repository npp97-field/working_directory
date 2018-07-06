#############
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


#### interpolate 8-day temperature to daily
spline_interpolate<-function(x){
  nyears<-length(x)/46
  if (sum(is.na(x))>100){
    return(rep(NA,365*nyears))
  }
  
  x_in<-(1:length(x)-0.5)/46
  x_out<-(1:(365*nyears)-0.5)/365
  spl_temp<-spline(x=x_in,y=x,xout=x_out)$y
  return(spl_temp)
}

get_era_data<-function(varname){
  #var should be PAR PREC or TEMP
  dir<-paste("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/data/daily_data_ERA/",varname,sep="")
  files<-list.files(dir,full.names = T)
  dat<-array(NA,dim=c(86400,365*length(files)))
  for (i in 1:length(files)){
    nc_in<-nc_open(files[i])
    var<-ncvar_get(nc_in,varname)
    nc_close(nc_in)
    vardat<-var[,241:360,]
    dim(vardat)<-c(86400,365)
    dat[,(365*i-364):(365*i)]<-vardat
  }
  return(dat)
}


### get climate data
precip_day<-get_era_data("PREC")*2*1000 ###convert meter/halfday to mm/day
temp_day<-get_era_data("TEMP")-273.15   ###convert K to deg C
par_day<-get_era_data("PAR")/43200      ###convert J/m2/halfday to W/m2      

get_mean<-function(xts){
  #this is for temperature, modified into daily
  if (is.na(xts[457]*xts[458]))
    return(NA)
  else{
    st_abs<-ceiling(xts[457]*365)+90
    end_abs<-ceiling(xts[458]*365)+90
    ave_var<-mean(xts[st_abs:end_abs],na.rm=T)
    return(ave_var)
  }
}


calculate_climate<-function(year){
  #get precipitation for this year +90 days (3 month) from previous year
  precip_year<-precip_day[,((year-2000)*365-90):((year-1999)*365)]
  temp_year<-temp_day[,((year-2000)*365-90):((year-1999)*365)]
  par_year<-par_day[,((year-2000)*365-90):((year-1999)*365)]
  
  ### growing season
  prec_dat<-cbind(precip_year,sos,eos)
  s2e_d_precip<-apply(prec_dat,1,get_mean)
  
  temp_dat<-cbind(temp_year,sos,eos)
  s2e_d_temp<-apply(temp_dat,1,get_mean)
  
  par_dat<-cbind(par_year,sos,eos)
  s2e_d_par<-apply(par_dat,1,get_mean)
  
  #### first half_growing season
  prec_dat<-cbind(precip_year,sos,pos)
  s2p_d_precip<-apply(prec_dat,1,get_mean)
  
  temp_dat<-cbind(temp_year,sos,pos)
  s2p_d_temp<-apply(temp_dat,1,get_mean)
  
  par_dat<-cbind(par_year,sos,pos)
  s2p_d_par<-apply(par_dat,1,get_mean)
  
  ## calculate the second half temp and precipitation
  temp_dat<-cbind(temp_year,pos,eos)
  p2e_d_temp<-apply(temp_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,pos,eos)
  p2e_d_precip<-apply(prec_dat,1,get_mean)
  
  par_dat<-cbind(par_year,pos,eos)
  p2e_d_par<-apply(par_dat,1,get_mean)
  ##################################
  ##----------------------------------------------
  ### calculate the average 1 month pre-season temperature
  temp_dat<-cbind(temp_year,sos-1/12,sos)
  presos1mon_temp<-apply(temp_dat,1,get_mean)
  
  temp_dat<-cbind(temp_year,eos-1/12,eos)
  preeos1mon_temp<-apply(temp_dat,1,get_mean)
  
  temp_dat<-cbind(temp_year,pos-1/24,pos+1/24)
  peak_temp<-apply(temp_dat,1,get_mean)
  
  par_dat<-cbind(par_year,sos-1/12,sos)
  presos1mon_par<-apply(par_dat,1,get_mean)
  
  par_dat<-cbind(par_year,eos-1/12,eos)
  preeos1mon_par<-apply(par_dat,1,get_mean)
  
  par_dat<-cbind(par_year,pos-1/24,pos+1/24)
  peak_par<-apply(par_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,sos-1/12,sos)
  presos1mon_prec<-apply(prec_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,eos-1/12,eos)
  preeos1mon_prec<-apply(prec_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,pos-1/24,pos+1/24)
  peak_prec<-apply(prec_dat,1,get_mean)
  
  ### calculate the average 2 month pre-season temperature
  temp_dat<-cbind(temp_year,sos-2/12,sos)
  presos2mon_temp<-apply(temp_dat,1,get_mean)
  
  temp_dat<-cbind(temp_year,eos-2/12,eos)
  preeos2mon_temp<-apply(temp_dat,1,get_mean)
  
  
  par_dat<-cbind(par_year,sos-2/12,sos)
  presos2mon_par<-apply(par_dat,1,get_mean)
  
  par_dat<-cbind(par_year,eos-2/12,eos)
  preeos2mon_par<-apply(par_dat,1,get_mean)
  
  
  prec_dat<-cbind(precip_year,sos-2/12,sos)
  presos2mon_prec<-apply(prec_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,eos-2/12,eos)
  preeos2mon_prec<-apply(prec_dat,1,get_mean)
  
  
  ### calculate the average 3 month pre-season temperature
  temp_dat<-cbind(temp_year,sos-3/12,sos)
  presos3mon_temp<-apply(temp_dat,1,get_mean)
  
  temp_dat<-cbind(temp_year,eos-3/12,eos)
  preeos3mon_temp<-apply(temp_dat,1,get_mean)
  
  
  par_dat<-cbind(par_year,sos-3/12,sos)
  presos3mon_par<-apply(par_dat,1,get_mean)
  
  par_dat<-cbind(par_year,eos-3/12,eos)
  preeos3mon_par<-apply(par_dat,1,get_mean)
  
  
  prec_dat<-cbind(precip_year,sos-3/12,sos)
  presos3mon_prec<-apply(prec_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,eos-3/12,eos)
  preeos3mon_prec<-apply(prec_dat,1,get_mean)
  
  
  ### calculate the average half month pre-season temperature
  temp_dat<-cbind(temp_year,sos-1/24,sos)
  presos0mon_temp<-apply(temp_dat,1,get_mean)
  
  temp_dat<-cbind(temp_year,eos-1/24,eos)
  preeos0mon_temp<-apply(temp_dat,1,get_mean)
  
  
  par_dat<-cbind(par_year,sos-1/24,sos)
  presos0mon_par<-apply(par_dat,1,get_mean)
  
  par_dat<-cbind(par_year,eos-1/24,eos)
  preeos0mon_par<-apply(par_dat,1,get_mean)
  
  
  prec_dat<-cbind(precip_year,sos-1/24,sos)
  presos0mon_prec<-apply(prec_dat,1,get_mean)
  
  prec_dat<-cbind(precip_year,eos-1/24,eos)
  preeos0mon_prec<-apply(prec_dat,1,get_mean)
  
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

  pre_start1_temp<-ncvar_def("pre_start1_temp",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end1_temp<-ncvar_def("pre_end1_temp",'',list(xdim2,ydim2),-9999,compression=9)

  p_par<-ncvar_def("p_par","",list(xdim2,ydim2),-9999,compression=9)
  pre_start1_par<-ncvar_def("pre_start1_par",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end1_par<-ncvar_def("pre_end1_par",'',list(xdim2,ydim2),-9999,compression=9)
  p_prec<-ncvar_def("p_prec","",list(xdim2,ydim2),-9999,compression=9)
  pre_start1_prec<-ncvar_def("pre_start1_prec",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end1_prec<-ncvar_def("pre_end1_prec",'',list(xdim2,ydim2),-9999,compression=9)
  
  pre_start2_temp<-ncvar_def("pre_start2_temp",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end2_temp<-ncvar_def("pre_end2_temp",'',list(xdim2,ydim2),-9999,compression=9)

  
  pre_start2_par<-ncvar_def("pre_start2_par",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end2_par<-ncvar_def("pre_end2_par",'',list(xdim2,ydim2),-9999,compression=9)
  
  pre_start2_prec<-ncvar_def("pre_start2_prec",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end2_prec<-ncvar_def("pre_end2_prec",'',list(xdim2,ydim2),-9999,compression=9)
  
  pre_start3_temp<-ncvar_def("pre_start3_temp",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end3_temp<-ncvar_def("pre_end3_temp",'',list(xdim2,ydim2),-9999,compression=9)

  
  pre_start3_par<-ncvar_def("pre_start3_par",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end3_par<-ncvar_def("pre_end3_par",'',list(xdim2,ydim2),-9999,compression=9)
  
  pre_start3_prec<-ncvar_def("pre_start3_prec",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end3_prec<-ncvar_def("pre_end3_prec",'',list(xdim2,ydim2),-9999,compression=9)
  
  pre_start0_temp<-ncvar_def("pre_start0_temp",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end0_temp<-ncvar_def("pre_end0_temp",'',list(xdim2,ydim2),-9999,compression=9)

  
  pre_start0_par<-ncvar_def("pre_start0_par",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end0_par<-ncvar_def("pre_end0_par",'',list(xdim2,ydim2),-9999,compression=9)
  
  pre_start0_prec<-ncvar_def("pre_start0_prec",'',list(xdim2,ydim2),-9999,compression=9)
  pre_end0_prec<-ncvar_def("pre_end0_prec",'',list(xdim2,ydim2),-9999,compression=9)
  
  
  year_climate<-paste("./pheno_hd_fixed_threshold_climate/clear_era/",year,"_climate_csif_clear_daily.nc",sep="")
  if (file.exists(year_climate)){
    file.remove(year_climate)
  }
  phenoin<-nc_create(year_climate,list(s2e_temp,s2p_temp,p2e_temp,
                                       s2e_prec,s2p_prec,p2e_prec,
                                       s2e_par,s2p_par,p2e_par,
                                       pre_start1_temp,pre_end1_temp,p_temp,

                                       pre_start1_par,pre_end1_par,p_par,
                                       pre_start1_prec,pre_end1_prec,p_prec,
                                       pre_start2_temp,pre_end2_temp,

                                       pre_start2_par,pre_end2_par,
                                       pre_start2_prec,pre_end2_prec,
                                       pre_start3_temp,pre_end3_temp,

                                       pre_start3_par,pre_end3_par,
                                       pre_start3_prec,pre_end3_prec,
                                       pre_start0_temp,pre_end0_temp,

                                       pre_start0_par,pre_end0_par,
                                       pre_start0_prec,pre_end0_prec
  ))
  ncvar_put(phenoin,s2e_prec,s2e_d_precip)
  ncvar_put(phenoin,p2e_prec,p2e_d_precip)
  ncvar_put(phenoin,s2p_prec,s2p_d_precip)
  ncvar_put(phenoin,s2e_temp,s2e_d_temp)
  ncvar_put(phenoin,s2p_temp,s2p_d_temp)
  ncvar_put(phenoin,p2e_temp,p2e_d_temp)

  ncvar_put(phenoin,s2e_par,s2e_d_par)
  ncvar_put(phenoin,s2p_par,s2p_d_par)
  ncvar_put(phenoin,p2e_par,p2e_d_par)
  ##----------------------------------------------
  ncvar_put(phenoin,pre_start1_temp,presos1mon_temp)
  ncvar_put(phenoin,pre_end1_temp,preeos1mon_temp)
  ncvar_put(phenoin,p_temp,peak_temp)

  ncvar_put(phenoin,pre_start1_par,presos1mon_par)
  ncvar_put(phenoin,pre_end1_par,preeos1mon_par)
  ncvar_put(phenoin,p_par,peak_par)
  ncvar_put(phenoin,pre_start1_prec,presos1mon_prec)
  ncvar_put(phenoin,pre_end1_prec,preeos1mon_prec)
  ncvar_put(phenoin,p_prec,peak_prec)
  
  ##----------------------------------------------
  ncvar_put(phenoin,pre_start2_temp,presos2mon_temp)
  ncvar_put(phenoin,pre_end2_temp,preeos2mon_temp)

  ncvar_put(phenoin,pre_start2_par,presos2mon_par)
  ncvar_put(phenoin,pre_end2_par,preeos2mon_par)
  ncvar_put(phenoin,pre_start2_prec,presos2mon_prec)
  ncvar_put(phenoin,pre_end2_prec,preeos2mon_prec)
  ##----------------------------------------------
  ncvar_put(phenoin,pre_start3_temp,presos3mon_temp)
  ncvar_put(phenoin,pre_end3_temp,preeos3mon_temp)

  ncvar_put(phenoin,pre_start3_par,presos3mon_par)
  ncvar_put(phenoin,pre_end3_par,preeos3mon_par)
  ncvar_put(phenoin,pre_start3_prec,presos3mon_prec)
  ncvar_put(phenoin,pre_end3_prec,preeos3mon_prec)
  ##----------------------------------------------
  ncvar_put(phenoin,pre_start0_temp,presos0mon_temp)
  ncvar_put(phenoin,pre_end0_temp,preeos0mon_temp)

  ncvar_put(phenoin,pre_start0_par,presos0mon_par)
  ncvar_put(phenoin,pre_end0_par,preeos0mon_par)
  ncvar_put(phenoin,pre_start0_prec,presos0mon_prec)
  ncvar_put(phenoin,pre_end0_prec,preeos0mon_prec)
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
#sif_climate<-list.files("./pheno_hd_fixed_threshold_climate/clear_era/",full.names = T,pattern = ".nc")
sif_climate<-list.files("./pheno_hd_fixed_threshold_climate/clear_era/",full.names = T,pattern = ".nc")

p2e_csif<-array(NA,dim=c(86400,16))
s2p_csif<-array(NA,dim=c(86400,16))
s2e_csif<-array(NA,dim=c(86400,16))
ann_csif<-array(NA,dim=c(86400,16))
sos<-array(NA,dim=c(86400,16))
eos<-array(NA,dim=c(86400,16))

s2e_temp<-array(NA,dim=c(86400,16))
s2p_temp<-array(NA,dim=c(86400,16))
p2e_temp<-array(NA,dim=c(86400,16))

s2e_prec<-array(NA,dim=c(86400,16))
s2p_prec<-array(NA,dim=c(86400,16))
p2e_prec<-array(NA,dim=c(86400,16))

s2e_par<-array(NA,dim=c(86400,16))
s2p_par<-array(NA,dim=c(86400,16))
p2e_par<-array(NA,dim=c(86400,16))

pre_start1_temp<-array(NA,dim=c(86400,16))
pre_end1_temp<-array(NA,dim=c(86400,16))
p_temp<-array(NA,dim=c(86400,16))

pre_start1_par<-array(NA,dim=c(86400,16))
pre_end1_par<-array(NA,dim=c(86400,16))
p_par<-array(NA,dim=c(86400,16))

pre_start1_prec<-array(NA,dim=c(86400,16))
pre_end1_prec<-array(NA,dim=c(86400,16))
p_prec<-array(NA,dim=c(86400,16))

pre_start2_temp<-array(NA,dim=c(86400,16))
pre_end2_temp<-array(NA,dim=c(86400,16))

pre_start2_par<-array(NA,dim=c(86400,16))
pre_end2_par<-array(NA,dim=c(86400,16))

pre_start2_prec<-array(NA,dim=c(86400,16))
pre_end2_prec<-array(NA,dim=c(86400,16))

pre_start3_temp<-array(NA,dim=c(86400,16))
pre_end3_temp<-array(NA,dim=c(86400,16))

pre_start3_par<-array(NA,dim=c(86400,16))
pre_end3_par<-array(NA,dim=c(86400,16))

pre_start3_prec<-array(NA,dim=c(86400,16))
pre_end3_prec<-array(NA,dim=c(86400,16))

pre_start0_temp<-array(NA,dim=c(86400,16))
pre_end0_temp<-array(NA,dim=c(86400,16))

pre_start0_par<-array(NA,dim=c(86400,16))
pre_end0_par<-array(NA,dim=c(86400,16))

pre_start0_prec<-array(NA,dim=c(86400,16))
pre_end0_prec<-array(NA,dim=c(86400,16))


for (year in 2001:2016){
  year_pheno<-sif_pheno[substr(basename(sif_pheno),35,38)==year]
  phenoin<-nc_open(year_pheno)
  p2e_csif[,year-2000]<-ncvar_get(phenoin,'p2s_csif')
  s2p_csif[,year-2000]<-ncvar_get(phenoin,'s2p_csif')
  s2e_csif[,year-2000]<-ncvar_get(phenoin,'s2s_csif')
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
  p_par[,year-2000]<-ncvar_get(climatein,'p_par')
  p_prec[,year-2000]<-ncvar_get(climatein,'p_prec')

  s2e_prec[,year-2000]<-ncvar_get(climatein,'s2e_prec')
  s2p_prec[,year-2000]<-ncvar_get(climatein,'s2p_prec')
  p2e_prec[,year-2000]<-ncvar_get(climatein,'p2e_prec')
  
  s2e_par[,year-2000]<-ncvar_get(climatein,'s2e_par')
  s2p_par[,year-2000]<-ncvar_get(climatein,'s2p_par')
  p2e_par[,year-2000]<-ncvar_get(climatein,'p2e_par')
  
  pre_start1_temp[,year-2000]<-ncvar_get(climatein,'pre_start1_temp')
  pre_end1_temp[,year-2000]<-ncvar_get(climatein,'pre_end1_temp')
  p_temp[,year-2000]<-ncvar_get(climatein,'p_temp')

  
  pre_start1_par[,year-2000]<-ncvar_get(climatein,'pre_start1_par')
  pre_end1_par[,year-2000]<-ncvar_get(climatein,'pre_end1_par')
  p_par[,year-2000]<-ncvar_get(climatein,'p_par')
  
  pre_start1_prec[,year-2000]<-ncvar_get(climatein,'pre_start1_prec')
  pre_end1_prec[,year-2000]<-ncvar_get(climatein,'pre_end1_prec')
  p_prec[,year-2000]<-ncvar_get(climatein,'p_prec')
  
  pre_start2_temp[,year-2000]<-ncvar_get(climatein,'pre_start2_temp')
  pre_end2_temp[,year-2000]<-ncvar_get(climatein,'pre_end2_temp')

  
  pre_start2_par[,year-2000]<-ncvar_get(climatein,'pre_start2_par')
  pre_end2_par[,year-2000]<-ncvar_get(climatein,'pre_end2_par')
  
  pre_start2_prec[,year-2000]<-ncvar_get(climatein,'pre_start2_prec')
  pre_end2_prec[,year-2000]<-ncvar_get(climatein,'pre_end2_prec')
  
  pre_start3_temp[,year-2000]<-ncvar_get(climatein,'pre_start3_temp')
  pre_end3_temp[,year-2000]<-ncvar_get(climatein,'pre_end3_temp')

  
  pre_start3_par[,year-2000]<-ncvar_get(climatein,'pre_start3_par')
  pre_end3_par[,year-2000]<-ncvar_get(climatein,'pre_end3_par')
  
  pre_start3_prec[,year-2000]<-ncvar_get(climatein,'pre_start3_prec')
  pre_end3_prec[,year-2000]<-ncvar_get(climatein,'pre_end3_prec')
  
  pre_start0_temp[,year-2000]<-ncvar_get(climatein,'pre_start0_temp')
  pre_end0_temp[,year-2000]<-ncvar_get(climatein,'pre_end0_temp')

  
  pre_start0_par[,year-2000]<-ncvar_get(climatein,'pre_start0_par')
  pre_end0_par[,year-2000]<-ncvar_get(climatein,'pre_end0_par')
  
  pre_start0_prec[,year-2000]<-ncvar_get(climatein,'pre_start0_prec')
  pre_end0_prec[,year-2000]<-ncvar_get(climatein,'pre_end0_prec')
  
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
  eos_csif<-cbind(sos,p2e_csif,s2e_temp,s2e_prec,s2e_par)
  fileout<-"./analysis/correlation_clear_era/pcor_sos_csif_p2e.nc"
  cal_pcor(eos_csif,fileout)
  ###########################
  eos_csif<-cbind(sos,p2e_csif)
  fileout<-"./analysis/correlation_clear_era/cor_sos_csif_p2e.nc"
  cal_cor(eos_csif,fileout)
  
  ##########calculate the pcor between sos and csif_S2P
  
  sos_csif<-cbind(sos,s2p_csif,s2e_temp,s2p_prec,s2p_par)
  fileout<-"./analysis/correlation_clear_era/pcor_sos_csif_s2p.nc"
  cal_pcor(sos_csif,fileout)
  
  #######################
  sos_csif<-cbind(sos,s2p_csif)
  fileout<-"./analysis/correlation_clear_era/cor_sos_csif_s2p.nc"
  cal_cor(sos_csif,fileout)
  
  ##########calculate the pcor between sos and csif_S2S
  sos_csif<-cbind(sos,s2e_csif,s2e_temp,s2e_prec)
  fileout<-"./analysis/correlation_clear_era/pcor_sos_csif_s2e.nc"
  cal_pcor(sos_csif,fileout)
  
  
  #######################
  sos_csif<-cbind(sos,s2e_csif)
  fileout<-"./analysis/correlation_clear_era/cor_sos_csif_s2e.nc"
  cal_cor(sos_csif,fileout)
  
  
  ##########calculate the pcor between sos and eos
  sos_eos<-cbind(sos,eos,s2e_temp,s2e_prec,s2e_par)
  fileout<-"./analysis/correlation_clear_era/pcor_sos_eos.nc"
  cal_pcor(sos_eos,fileout)
  
  ####################################
  sos_eos<-cbind(sos,eos)
  fileout<-"./analysis/correlation_clear_era/cor_sos_eos.nc"
  cal_cor(sos_eos,fileout)
}



###### calculate the controlling factors of SOS EOS
if (T){
  ##########calculate the cor between sos and eos and pre season temp
  sos_eos<-cbind(sos,pre_start1_temp)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre1_temp.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start2_temp)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre2_temp.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start3_temp)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre3_temp.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start0_temp)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre0_temp.nc"
  cal_cor(sos_eos,fileout)

  ####################################     Prec
  sos_eos<-cbind(sos,pre_start1_prec)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre1_prec.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start2_prec)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre2_prec.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start3_prec)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre3_prec.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start0_prec)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre0_prec.nc"
  cal_cor(sos_eos,fileout)
  
  
  ####################################     PAR
  sos_eos<-cbind(sos,pre_start1_par)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre1_par.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start2_par)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre2_par.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start3_par)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre3_par.nc"
  cal_cor(sos_eos,fileout)
  sos_eos<-cbind(sos,pre_start0_par)
  fileout<-"./analysis/correlation_clear_era/cor_sos_pre0_par.nc"
  cal_cor(sos_eos,fileout)
  
  ############ eos and late season temperature
  eos_eos<-cbind(eos,pre_end1_temp)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre1_temp.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end2_temp)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre2_temp.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end3_temp)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre3_temp.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end0_temp)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre0_temp.nc"
  cal_cor(eos_eos,fileout)
  
  ####################################     Prec
  eos_eos<-cbind(eos,pre_end1_prec)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre1_prec.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end2_prec)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre2_prec.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end3_prec)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre3_prec.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end0_prec)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre0_prec.nc"
  cal_cor(eos_eos,fileout)
  
  
  ####################################     PAR
  eos_eos<-cbind(eos,pre_end1_par)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre1_par.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end2_par)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre2_par.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end3_par)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre3_par.nc"
  cal_cor(eos_eos,fileout)
  eos_eos<-cbind(eos,pre_end0_par)
  fileout<-"./analysis/correlation_clear_era/cor_eos_pre0_par.nc"
  cal_cor(eos_eos,fileout)
  
}

# 
# ###### temperature feedbacks
# if (T){
#   #################################### spring temperature and fall temperature
#   pre_start_pre_end<-cbind(pre_start_temp,pre_end_temp)
#   nc_out_f3<-"./analysis/correlation_clear_era/cor_pre_start_pre_end_temp.nc"
#   cal_cor(pre_start_pre_end,nc_out_f3)
#   
#   s2p_p2e<-cbind(s2p_temp,p2e_temp)
#   nc_out_f3<-"./analysis/correlation_clear_era/cor_s2p_p2e_temp.nc"
#   cal_cor(s2p_p2e,nc_out_f3)
#   
#   
#   ############################### growing season temp precip and growing season CSIF
#   s2s_sif_temp<-cbind(s2s_csif,s2e_temp)
#   nc_out_f3<-"./analysis/correlation_clear_era/cor_s2s_sif_s2s_temp.nc"
#   cal_cor(s2s_sif_temp,nc_out_f3)
#   
#   
#   s2s_sif_prec<-cbind(s2s_csif,s2e_prec)
#   nc_out_f3<-"./analysis/correlation_clear_era/cor_s2s_sif_s2s_precip.nc"
#   cal_cor(s2s_sif_prec,nc_out_f3)
#   
#   ############################### pre EOS and peak2end temp
#   s2s_sif_prec<-cbind(pre_end_temp,p2e_temp)
#   nc_out_f3<-"./analysis/correlation_clear_era/cor_preEOS_p2e_temp.nc"
#   cal_cor(s2s_sif_prec,nc_out_f3)
#   
#   ############################### pre SOS and start2peak temp
#   s2s_sif_prec<-cbind(pre_start_temp,s2p_temp)
#   nc_out_f3<-"./analysis/correlation_clear_era/cor_preSOS_s2p_temp.nc"
#   cal_cor(s2s_sif_prec,nc_out_f3)
#   
#   ############################### s2p sif and p2e sif
#   s2s_sif_prec<-cbind(s2p_csif,p2s_csif)
#   nc_out_f3<-"./analysis/correlation_clear_era/cor_s2p_p2s_csif.nc"
#   cal_cor(s2s_sif_prec,nc_out_f3)
# }
# 
# 


