#############
### calcuate the peak to senescence mean SIF or NDVI for each year
library(ncdf4)
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")
#setwd("/Users/yzhang/Project/SIF_phenology/")
## readin the SIF files
sif_files<-list.files("/rigel/glab/users/zy2309/DATA/all_daily_SIF_4day_HD/", full.names = T,pattern=".nc")
sif_pheno<-list.files("./pheno_hd_fixed_threshold/",full.names = T,pattern = ".nc")
## readin the SIF pos eos dates.
sos_file<-'./analysis/SIF_SOS_30N_fixed_stat.nc'
pos_file<-"./analysis/SIF_POS_30N_var_stat.nc"
eos_file<-"./analysis/SIF_EOS_30N_fixed_stat.nc"

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
  else
    return(mean(xts[xts[93]:xts[94]],na.rm=T))
}

### for each year, read in the sif_files
get_late_growingseason_SIF<-function(year){
  year_csif_f<-sif_files[substr(basename(sif_files),20,23)==year]
  csifin<-nc_open(year_csif_f)
  csif_year<-ncvar_get(csifin,varid="all_daily_sif")[,241:360,]
  nc_close(csifin)
  dim(csif_year)<-c(86400,92)
  mean_csif_peak2sene<-apply(cbind(csif_year,pos,eos),1,get_mean_peak2sene)
  dim(mean_csif_peak2sene)<-c(720,120)
  
  ####
  mean_csif_start2peak<-apply(cbind(csif_year,sos,pos),1,get_mean_peak2sene)
  dim(mean_csif_start2peak)<-c(720,120)
  
  year_pheno<-sif_pheno[substr(basename(sif_pheno),19,22)==year]
  phenoin<-nc_open(year_pheno,write=T)
  xdim2<-phenoin$dim[["longitude"]]
  ydim2<-phenoin$dim[["latitude"]]
  
  p2s_csif<-ncvar_def("p2s_csif",'',list(xdim2,ydim2),-9999,compression=9)
  s2p_csif<-ncvar_def("s2p_csif",'',list(xdim2,ydim2),-9999,compression=9)
  tryCatch({
    ncvar_add(phenoin,v = p2s_csif)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    ncvar_add(phenoin,v = s2p_csif)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  nc_close(phenoin)
  phenoin<-nc_open(year_pheno,write=T)
  ncvar_put(phenoin,p2s_csif,mean_csif_peak2sene)
  ncvar_put(phenoin,s2p_csif,mean_csif_start2peak)
  nc_close(phenoin)
}

for (year in 2003:2016){
  print(year)
  get_late_growingseason_SIF(year)
}

################################
### calcualte the anomaly of SOS and p2s csif
# 
setwd("/rigel/glab/users/zy2309/")
sif_pheno<-list.files("./PROJECT/SIF_phenology/pheno_hd_fixed_threshold/",full.names = T,pattern = ".nc")

sos_file<-"./PROJECT/SIF_phenology/analysis/SIF_SOS_30N_fixed_stat.nc"
eos_file<-"./PROJECT/SIF_phenology/analysis/SIF_EOS_30N_fixed_stat.nc"

##### get average EOS and SOS
eos_f<-nc_open(eos_file)
eos<-round(ncvar_get(eos_f,varid="MEAN")*12+0.5)
nc_close(eos_f)
sos_f<-nc_open(sos_file)
sos<-round(ncvar_get(sos_f,varid="MEAN")*12+0.5)
nc_close(sos_f)
dim(eos)<-c(86400,1)
dim(sos)<-c(86400,1)

get_cru_var<-function(var_name){
  var_f<-list.files('./DATA/CRU_TS401/',full.names = T,pattern = paste(var_name,".dat.nc",sep=""))
  var_in1<-nc_open(var_f[1])
  var1<-ncvar_get(var_in1,var_name)[,241:360,]
  nc_close(var_in1)
  dim(var1)<-c(86400,120)
  
  var_in2<-nc_open(var_f[2])
  var2<-ncvar_get(var_in2,var_name)[,241:360,]
  nc_close(var_in2)
  dim(var2)<-c(86400,72)
  var<-cbind(var1,var2)
  return(var)
}


tmean<-get_cru_var("tmp")
precip<-get_cru_var('pre')


get_total<-function(xts){
  if (is.na(xts[14]*xts[15]))
    return(NA)
  else
    return(sum(xts[xts[14]:(xts[15]+1)],na.rm=T))
}

get_mean<-function(xts){
  if (is.na(xts[14]*xts[15]))
    return(NA)
  else
    return(mean(xts[xts[14]:(xts[15]+1)],na.rm=T))
}


for (year in 2003:2016){
  #get precipitation for this year +one month from previous year
  precip_year<-precip[,((year-2001)*12):((year-2000)*12)]
  tmean_year<-tmean[,((year-2001)*12):((year-2000)*12)]
  
  prec_dat<-cbind(precip_year,sos,eos)
  gowning_precip<-apply(prec_dat,1,get_total)
  
  temp_dat<-cbind(tmean_year,sos,eos)
  gowning_temp<-apply(temp_dat,1,get_mean)
  
  ### write data to pheno
  year_pheno<-sif_pheno[substr(basename(sif_pheno),19,22)==year]
  phenoin<-nc_open(year_pheno,write=T)
  xdim2<-phenoin$dim[["longitude"]]
  ydim2<-phenoin$dim[["latitude"]]
  
  grow_precip<-ncvar_def("grow_precip",'',list(xdim2,ydim2),-9999,compression=9)
  grow_temp<-ncvar_def("grow_temp",'',list(xdim2,ydim2),-9999,compression=9)
  tryCatch({
    ncvar_add(phenoin,v = grow_precip)
    ncvar_add(phenoin,v = grow_temp)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  nc_close(phenoin)
  phenoin<-nc_open(year_pheno,write=T)
  ncvar_put(phenoin,grow_precip,gowning_precip)
  ncvar_put(phenoin,grow_temp,gowning_temp)
  nc_close(phenoin)
}


#####################################
#######   analysis the partial corrlation between cSIF growing season with SOS when Tmean and precip fixed.
#####################################

## get the average of SOS, CSIF, Tmean and Precip
library(ppcor)
setwd("/rigel/glab/users/zy2309/")
sif_pheno<-list.files("./PROJECT/SIF_phenology/pheno_hd_fixed_threshold/",full.names = T,pattern = ".nc")


calculate_partial_correlation<-function(dat){
  if (sum(is.na(dat))>=1)
    return(c(NA,NA))
  dim(dat)<-c(14,4)
  pcor_re<-pcor.test(dat[,1],dat[,2],dat[,c(3,4)],method="pearson")
  return(c(pcor_re$estimate,pcor_re$p.value))
}

p2s_csif<-array(NA,dim=c(86400,14))
s2p_csif<-array(NA,dim=c(86400,14))
sos<-array(NA,dim=c(86400,14))
tmean<-array(NA,dim=c(86400,14))
precip<-array(NA,dim=c(86400,14))
eos<-array(NA,dim=c(86400,14))

for (year in 2003:2016){
  year_pheno<-sif_pheno[substr(basename(sif_pheno),19,22)==year]
  phenoin<-nc_open(year_pheno,write=T)
  p2s_csif[,year-2002]<-ncvar_get(phenoin,'p2s_csif')
  s2p_csif[,year-2002]<-ncvar_get(phenoin,'p2s_csif')
  tmean[,year-2002]<-ncvar_get(phenoin,'grow_temp')
  precip[,year-2002]<-ncvar_get(phenoin,'grow_precip')
  sos[,year-2002]<-ncvar_get(phenoin,'SOS')*365+2.5
  eos[,year-2002]<-ncvar_get(phenoin,'EOS')*365+2.5
  xdim<-phenoin$dim[["longitude"]]
  ydim<-phenoin$dim[["latitude"]]
  nc_close(phenoin)
}

##########calculate the pcor between sos and csif_P2S

eos_csif<-cbind(sos,p2s_csif,tmean,precip)
eos_csif[eos_csif< -990]<-NA

pcor_eos_csif<-apply(eos_csif,1,calculate_partial_correlation)
p_coef<-pcor_eos_csif[1,]
p_coef[is.nan(p_coef)]<- -999.9
dim(p_coef)<-c(720,120)
p_pv<-pcor_eos_csif[2,]
p_pv[is.nan(p_pv)]<- -999.9
dim(p_pv)<-c(720,120)

nc_out_f1<-"./PROJECT/SIF_phenology//analysis/pcor_sos_csif_p2s.nc"

pcor_coef<-ncvar_def("pcor_coef",'',list(xdim,ydim),-999.9,prec="double",compression=9)
pcor_pv<-ncvar_def("pcor_pv",'',list(xdim,ydim),-999.9,prec="double",compression=9)
ncout<-nc_create(nc_out_f1,list(pcor_coef,pcor_pv))
ncvar_put(ncout,pcor_coef,p_coef)
ncvar_put(ncout,pcor_pv,p_pv)
nc_close(ncout)

##########calculate the pcor between sos and csif_S2P

sos_csif<-cbind(sos,s2p_csif,tmean,precip)
sos_csif[sos_csif< -990]<-NA

pcor_sos_csif<-apply(sos_csif,1,calculate_partial_correlation)
p_coef<-pcor_sos_csif[1,]
p_coef[is.nan(p_coef)]<- -999.9
dim(p_coef)<-c(720,120)
p_pv<-pcor_sos_csif[2,]
p_pv[is.nan(p_pv)]<- -999.9
dim(p_pv)<-c(720,120)

nc_out_f2<-"./PROJECT/SIF_phenology/analysis/pcor_sos_csif_s2p.nc"

pcor_coef<-ncvar_def("pcor_coef",'',list(xdim,ydim),-999.9,prec="double",compression=9)
pcor_pv<-ncvar_def("pcor_pv",'',list(xdim,ydim),-999.9,prec="double",compression=9)
ncout<-nc_create(nc_out_f1,list(pcor_coef,pcor_pv))
ncvar_put(ncout,pcor_coef,p_coef)
ncvar_put(ncout,pcor_pv,p_pv)
nc_close(ncout)

##########calculate the pcor between sos and eos

sos_eos<-cbind(sos,eos,tmean,precip)
sos_eos[sos_eos< -990]<-NA

pcor_sos_eos<-apply(sos_eos,1,calculate_partial_correlation)
p_coef<-pcor_sos_eos[1,]
p_coef[is.nan(p_coef)]<- -999.9
dim(p_coef)<-c(720,120)
p_pv<-pcor_sos_eos[2,]
p_pv[is.nan(p_pv)]<- -999.9
dim(p_pv)<-c(720,120)

nc_out_f1<-"./PROJECT/SIF_phenology/analysis/pcor_sos_eos.nc"

pcor_coef<-ncvar_def("pcor_coef",'',list(xdim,ydim),-999.9,prec="double",compression=9)
pcor_pv<-ncvar_def("pcor_pv",'',list(xdim,ydim),-999.9,prec="double",compression=9)
ncout<-nc_create(nc_out_f1,list(pcor_coef,pcor_pv))
ncvar_put(ncout,pcor_coef,p_coef)
ncvar_put(ncout,pcor_pv,p_pv)
nc_close(ncout)

