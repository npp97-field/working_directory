### analysis of correlations
## get the average of SOS, CSIF, Tmean and Precip
library(ppcor)
library(ncdf4)
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")
sif_pheno<-list.files("./pheno_hd_fixed_threshold_clear/",full.names = T,pattern = ".nc")
dataset<-"rs"
sif_climate<-list.files(paste("./pheno_hd_fixed_threshold_climate/clear_rs/"),full.names = T,pattern = ".nc")
#sif_climate_era<-list.files("./pheno_hd_fixed_threshold_climate/clear_era/",full.names = T,pattern = ".nc")

pre_end1_tday<-array(NA,dim=c(86400,14))
pre_end1_tnight<-array(NA,dim=c(86400,14))
pre_end1_par<-array(NA,dim=c(86400,14))
pre_end1_prec<-array(NA,dim=c(86400,14))
pre_end2_tday<-array(NA,dim=c(86400,14))
pre_end2_tnight<-array(NA,dim=c(86400,14))
pre_end2_par<-array(NA,dim=c(86400,14))
pre_end2_prec<-array(NA,dim=c(86400,14))
pre_end3_tday<-array(NA,dim=c(86400,14))
pre_end3_tnight<-array(NA,dim=c(86400,14))
pre_end3_par<-array(NA,dim=c(86400,14))
pre_end3_prec<-array(NA,dim=c(86400,14))
pre_end0_tday<-array(NA,dim=c(86400,14))
pre_end0_tnight<-array(NA,dim=c(86400,14))
pre_end0_par<-array(NA,dim=c(86400,14))
pre_end0_prec<-array(NA,dim=c(86400,14))

for (year in 2003:2016){
  year_climate<-sif_climate[substr(basename(sif_climate),1,4)==year]
  climatein<-nc_open(year_climate)

  pre_end1_tday[,year-2002]<-ncvar_get(climatein,'pre_end1_tday')
  pre_end1_tnight[,year-2002]<-ncvar_get(climatein,'pre_end1_tnight')
  pre_end1_par[,year-2002]<-ncvar_get(climatein,'pre_end1_par')
  pre_end1_prec[,year-2002]<-ncvar_get(climatein,'pre_end1_prec')
  pre_end2_tday[,year-2002]<-ncvar_get(climatein,'pre_end2_tday')
  pre_end2_tnight[,year-2002]<-ncvar_get(climatein,'pre_end2_tnight')
  pre_end2_par[,year-2002]<-ncvar_get(climatein,'pre_end2_par')
  pre_end2_prec[,year-2002]<-ncvar_get(climatein,'pre_end2_prec')
  pre_end3_tday[,year-2002]<-ncvar_get(climatein,'pre_end3_tday')
  pre_end3_tnight[,year-2002]<-ncvar_get(climatein,'pre_end3_tnight')
  pre_end3_par[,year-2002]<-ncvar_get(climatein,'pre_end3_par')
  pre_end3_prec[,year-2002]<-ncvar_get(climatein,'pre_end3_prec')
  pre_end0_tday[,year-2002]<-ncvar_get(climatein,'pre_end0_tday')
  pre_end0_tnight[,year-2002]<-ncvar_get(climatein,'pre_end0_tnight')
  pre_end0_par[,year-2002]<-ncvar_get(climatein,'pre_end0_par')
  pre_end0_prec[,year-2002]<-ncvar_get(climatein,'pre_end0_prec')
  xdim<-climatein$dim[["longitude"]]
  ydim<-climatein$dim[["latitude"]]
  nc_close(climatein)
}

calculate_correlation<-function(dat){
  if (sum(is.na(dat))>=1)
    return(c(NA,NA))
  n<-length(dat)
  dim(dat)<-c(n/2,2)
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

#correlation between PREC and LSTday at 0,1,2,3 lags

cli_eos<-cbind(pre_end1_tday,pre_end1_prec)
fileout<-"./analysis/correlation_clear_rs/cor_pre1_prec_tday.nc"
cal_cor(cli_eos,fileout)

cli_eos<-cbind(pre_end2_tday,pre_end2_prec)
fileout<-"./analysis/correlation_clear_rs/cor_pre2_prec_tday.nc"
cal_cor(cli_eos,fileout)

cli_eos<-cbind(pre_end3_tday,pre_end3_prec)
fileout<-"./analysis/correlation_clear_rs/cor_pre3_prec_tday.nc"
cal_cor(cli_eos,fileout)

cli_eos<-cbind(pre_end0_tday,pre_end0_prec)
fileout<-"./analysis/correlation_clear_rs/cor_pre0_prec_tday.nc"
cal_cor(cli_eos,fileout)

#########################
sif_climate<-list.files("./pheno_hd_fixed_threshold_climate/clear_era/",full.names = T,pattern = ".nc")
pre_end1_temp<-array(NA,dim=c(86400,16))
pre_end1_par<-array(NA,dim=c(86400,16))
pre_end1_prec<-array(NA,dim=c(86400,16))
pre_end2_temp<-array(NA,dim=c(86400,16))
pre_end2_par<-array(NA,dim=c(86400,16))
pre_end2_prec<-array(NA,dim=c(86400,16))
pre_end3_temp<-array(NA,dim=c(86400,16))
pre_end3_par<-array(NA,dim=c(86400,16))
pre_end3_prec<-array(NA,dim=c(86400,16))
pre_end0_temp<-array(NA,dim=c(86400,16))
pre_end0_par<-array(NA,dim=c(86400,16))
pre_end0_prec<-array(NA,dim=c(86400,16))

for (year in 2001:2016){
  year_climate<-sif_climate[substr(basename(sif_climate),1,4)==year]
  climatein<-nc_open(year_climate)
  
  pre_end1_temp[,year-2000]<-ncvar_get(climatein,'pre_end1_temp')
  pre_end1_par[,year-2000]<-ncvar_get(climatein,'pre_end1_par')
  pre_end1_prec[,year-2000]<-ncvar_get(climatein,'pre_end1_prec')
  pre_end2_temp[,year-2000]<-ncvar_get(climatein,'pre_end2_temp')
  pre_end2_par[,year-2000]<-ncvar_get(climatein,'pre_end2_par')
  pre_end2_prec[,year-2000]<-ncvar_get(climatein,'pre_end2_prec')
  pre_end3_temp[,year-2000]<-ncvar_get(climatein,'pre_end3_temp')
  pre_end3_par[,year-2000]<-ncvar_get(climatein,'pre_end3_par')
  pre_end3_prec[,year-2000]<-ncvar_get(climatein,'pre_end3_prec')
  pre_end0_temp[,year-2000]<-ncvar_get(climatein,'pre_end0_temp')
  pre_end0_par[,year-2000]<-ncvar_get(climatein,'pre_end0_par')
  pre_end0_prec[,year-2000]<-ncvar_get(climatein,'pre_end0_prec')
  xdim<-climatein$dim[["longitude"]]
  ydim<-climatein$dim[["latitude"]]
  nc_close(climatein)
}
#correlation between PREC and LSTday at 0,1,2,3 lags

cli_eos<-cbind(pre_end1_temp,pre_end1_prec)
fileout<-"./analysis/correlation_clear_era/cor_pre1_prec_temp.nc"
cal_cor(cli_eos,fileout)

cli_eos<-cbind(pre_end2_temp,pre_end2_prec)
fileout<-"./analysis/correlation_clear_era/cor_pre2_prec_temp.nc"
cal_cor(cli_eos,fileout)

cli_eos<-cbind(pre_end3_temp,pre_end3_prec)
fileout<-"./analysis/correlation_clear_era/cor_pre3_prec_temp.nc"
cal_cor(cli_eos,fileout)

cli_eos<-cbind(pre_end0_temp,pre_end0_prec)
fileout<-"./analysis/correlation_clear_era/cor_pre0_prec_temp.nc"
cal_cor(cli_eos,fileout)

