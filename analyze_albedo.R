library(ncdf4)
library(ppcor)
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology")
sif_pheno<-list.files("./pheno_hd_fixed_threshold/",full.names = T,pattern = ".nc")

sos_file<-'./analysis/all_daily_SOS_30N_fixed_stat.nc'
pos_file<-"./analysis/all_daily_POS_30N_fixed_stat.nc"
eos_file<-"./analysis/all_daily_EOS_30N_fixed_stat.nc"

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

get_mean<-function(xts){
  #this is for temperature, modified into daily
  if (is.na(xts[366]*xts[367]))
    return(NA)
  else{
    st_abs<-max(round(xts[366]*365),1)
    end_abs<-min(round(xts[367]*365),365)
    ave_var<-mean(xts[st_abs:end_abs],na.rm=T)
    return(ave_var)
  }
}


calculate_albedo<-function(year){
  #get precipitation for this year +one month from previous year
  var_f<-list.files('/rigel/glab/users/zy2309/DATA/MCD43C3_HD/',full.names = T,pattern = paste(year,".*.nc",sep=""))
  var_in1<-nc_open(var_f[1])
  bsa<-ncvar_get(var_in1,"BSA")[,241:360,]
  wsa<-ncvar_get(var_in1,"WSA")[,241:360,]
  nc_close(var_in1)
  dim(bsa)<-c(86400,365)
  dim(wsa)<-c(86400,365)
  asa<-(bsa+wsa)/2
  ###########################
  albedo_dat<-cbind(bsa,sos,eos)
  s2e_b_sa<-apply(albedo_dat,1,get_mean)
  
  albedo_dat<-cbind(bsa,sos,pos)
  s2p_b_sa<-apply(albedo_dat,1,get_mean)
  
  albedo_dat<-cbind(bsa,pos,eos)
  p2e_b_sa<-apply(albedo_dat,1,get_mean)
  ###########################
  albedo_dat<-cbind(wsa,sos,eos)
  s2e_w_sa<-apply(albedo_dat,1,get_mean)
  
  albedo_dat<-cbind(wsa,sos,pos)
  s2p_w_sa<-apply(albedo_dat,1,get_mean)
  
  albedo_dat<-cbind(wsa,pos,eos)
  p2e_w_sa<-apply(albedo_dat,1,get_mean)
  ###########################
  albedo_dat<-cbind(asa,sos,eos)
  s2e_a_sa<-apply(albedo_dat,1,get_mean)
  
  albedo_dat<-cbind(asa,sos,pos)
  s2p_a_sa<-apply(albedo_dat,1,get_mean)
  
  albedo_dat<-cbind(asa,pos,eos)
  p2e_a_sa<-apply(albedo_dat,1,get_mean)
  ###########################
  albedo_dat<-cbind(wsa,pos-15/365,pos+15/365)
  p_w_sa<-apply(albedo_dat,1,get_mean)
  
  albedo_dat<-cbind(bsa,pos-15/365,pos+15/365)
  p_b_sa<-apply(albedo_dat,1,get_mean)
  
  albedo_dat<-cbind(asa,pos-15/365,pos+15/365)
  p_a_sa<-apply(albedo_dat,1,get_mean)
  
  s2e_wsa<-ncvar_def("s2e_wsa",'',list(xdim2,ydim2),-9999,compression=9)
  s2e_bsa<-ncvar_def("s2e_bsa",'',list(xdim2,ydim2),-9999,compression=9)
  s2e_asa<-ncvar_def("s2e_asa",'',list(xdim2,ydim2),-9999,compression=9)
  s2p_wsa<-ncvar_def("s2p_wsa",'',list(xdim2,ydim2),-9999,compression=9)
  s2p_bsa<-ncvar_def("s2p_bsa",'',list(xdim2,ydim2),-9999,compression=9)
  s2p_asa<-ncvar_def("s2p_asa",'',list(xdim2,ydim2),-9999,compression=9)
  p2e_wsa<-ncvar_def("p2e_wsa",'',list(xdim2,ydim2),-9999,compression=9)
  p2e_bsa<-ncvar_def("p2e_bsa",'',list(xdim2,ydim2),-9999,compression=9)
  p2e_asa<-ncvar_def("p2e_asa",'',list(xdim2,ydim2),-9999,compression=9)
  p_wsa<-ncvar_def("p_wsa",'',list(xdim2,ydim2),-9999,compression=9)
  p_bsa<-ncvar_def("p_bsa",'',list(xdim2,ydim2),-9999,compression=9)
  p_asa<-ncvar_def("p_asa",'',list(xdim2,ydim2),-9999,compression=9)

  year_albedo<-paste("./pheno_hd_fixed_threshold_albedo/",year,"_albedo_csif_all_daily.nc",sep="")
  if (file.exists(year_albedo)){
    file.remove(year_albedo)
  }
  phenoin<-nc_create(year_albedo,list(s2e_wsa,s2e_bsa,s2e_asa,s2p_wsa,s2p_bsa,s2p_asa,p2e_wsa,p2e_bsa,p2e_asa,
                                       p_wsa,p_bsa,p_asa))
  ncvar_put(phenoin,s2e_wsa,s2e_w_sa)
  ncvar_put(phenoin,s2e_bsa,s2e_b_sa)
  ncvar_put(phenoin,s2e_asa,s2e_a_sa)
  ncvar_put(phenoin,p2e_wsa,p2e_w_sa)
  ncvar_put(phenoin,p2e_bsa,p2e_b_sa)
  ncvar_put(phenoin,p2e_asa,p2e_a_sa)
  ncvar_put(phenoin,s2p_wsa,s2p_w_sa)
  ncvar_put(phenoin,s2p_bsa,s2p_b_sa)
  ncvar_put(phenoin,s2p_asa,s2p_a_sa)
  ncvar_put(phenoin,p_wsa,p_w_sa)
  ncvar_put(phenoin,p_bsa,p_b_sa)
  ncvar_put(phenoin,p_asa,p_a_sa)
  ##************************************************
  nc_close(phenoin)
}

for (year in 2001:2016){
  print(year)
  calculate_albedo(year)
}

#############################################
########### relationship analysis
###########################################

library(ppcor)
library(ncdf4)
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")
sif_pheno<-list.files("./pheno_hd_fixed_threshold/",full.names = T,pattern = ".nc")
sif_climate<-list.files("./pheno_hd_fixed_threshold_climate/",full.names = T,pattern = ".nc")
sif_albedo<-list.files("./pheno_hd_fixed_threshold_albedo/",full.names = T,pattern = ".nc")

p2s_csif<-array(NA,dim=c(86400,16))
s2p_csif<-array(NA,dim=c(86400,16))
s2s_csif<-array(NA,dim=c(86400,16))
ann_csif<-array(NA,dim=c(86400,16))
p_csif<-array(NA,dim=c(86400,16))
sos<-array(NA,dim=c(86400,16))
eos<-array(NA,dim=c(86400,16))

p2s_asa<-array(NA,dim=c(86400,16))
s2p_asa<-array(NA,dim=c(86400,16))
s2s_asa<-array(NA,dim=c(86400,16))
p_asa<-array(NA,dim=c(86400,16))

s2e_temp<-array(NA,dim=c(86400,16))
s2e_prec<-array(NA,dim=c(86400,16))
s2p_temp<-array(NA,dim=c(86400,16))
s2p_prec<-array(NA,dim=c(86400,16))
p2e_temp<-array(NA,dim=c(86400,16))
p2e_prec<-array(NA,dim=c(86400,16))
pre_start_temp<-array(NA,dim=c(86400,16))
pre_end_temp<-array(NA,dim=c(86400,16))
p_temp<-array(NA,dim=c(86400,16))

for (year in 2001:2016){
  year_pheno<-sif_pheno[substr(basename(sif_pheno),33,36)==year]
  phenoin<-nc_open(year_pheno)
  p2s_csif[,year-2000]<-ncvar_get(phenoin,'p2s_csif')
  s2p_csif[,year-2000]<-ncvar_get(phenoin,'s2p_csif')
  s2s_csif[,year-2000]<-ncvar_get(phenoin,'s2s_csif')
  ann_csif[,year-2000]<-ncvar_get(phenoin,'ann_csif')
  p_csif[,year-2000]<-ncvar_get(phenoin,'p_csif')
  sos[,year-2000]<-ncvar_get(phenoin,'SOS')*365+2.5
  eos[,year-2000]<-ncvar_get(phenoin,'EOS')*365+2.5
  xdim<-phenoin$dim[["longitude"]]
  ydim<-phenoin$dim[["latitude"]]
  nc_close(phenoin)
}

for (year in 2001:2016){
  year_pheno<-sif_albedo[substr(basename(sif_albedo),1,4)==year]
  phenoin<-nc_open(year_pheno)
  p2s_asa[,year-2000]<-ncvar_get(phenoin,'p2e_asa')
  s2p_asa[,year-2000]<-ncvar_get(phenoin,'s2p_asa')
  s2s_asa[,year-2000]<-ncvar_get(phenoin,'s2e_asa')
  p_asa[,year-2000]<-ncvar_get(phenoin,'p_asa')
  nc_close(phenoin)
}

for (year in 2001:2016){
  year_climate<-sif_climate[substr(basename(sif_climate),1,4)==year]
  climatein<-nc_open(year_climate)
  
  s2e_temp[,year-2000]<-ncvar_get(climatein,'s2e_temp')
  s2e_prec[,year-2000]<-ncvar_get(climatein,'s2e_prec')
  s2p_temp[,year-2000]<-ncvar_get(climatein,'s2p_temp')
  s2p_prec[,year-2000]<-ncvar_get(climatein,'s2p_prec')
  p2e_temp[,year-2000]<-ncvar_get(climatein,'p2e_temp')
  p2e_prec[,year-2000]<-ncvar_get(climatein,'p2e_prec')
  pre_start_temp[,year-2000]<-ncvar_get(climatein,'pre_start_temp')
  pre_end_temp[,year-2000]<-ncvar_get(climatein,'pre_end_temp')
  p_temp[,year-2000]<-ncvar_get(climatein,'p_temp')
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


##########calculate the pcor between sos and csif_P2S

eos_csif<-cbind(s2p_csif,s2p_asa)
fileout<-"./analysis/correlation/cor_s2p_csif_albedo.nc"
cal_cor(eos_csif,fileout)

eos_csif<-cbind(p2s_csif,p2s_asa)
fileout<-"./analysis/correlation/cor_p2e_csif_albedo.nc"
cal_cor(eos_csif,fileout)

dat<-cbind(p_csif,p_asa)
fileout<-"./analysis/correlation/cor_peak_csif_albedo.nc"
cal_cor(dat,fileout)

dat<-cbind(p_csif,p_temp)
fileout<-"./analysis/correlation/cor_peak_csif_peak_temp.nc"
cal_cor(dat,fileout)

dat<-cbind(sos,p_temp)
fileout<-"./analysis/correlation/cor_sos_peak_temp.nc"
cal_cor(dat,fileout)