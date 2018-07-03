######
### use preseason temperature precipitation and PAR to predict SOS and EOS

# for each pixel, we need to fit a multiple linear regresion model, and get
# the coefficient map for temperature, precipitation and PAR. In addition, 
# we need the model explanation power.

library(ncdf4)

#### the function to get the MVR coeff and model Performance
mvr<-function(x){
  ### input vector should have 16(year)*4(variables) elements
  ## of the following order: sos(eos), temp, rad, precip
  if (sum(is.na(x))>10){
    return(rep(NA,5))
  }
  dat<-x
  dim(dat)<-c(16,4)
  dat<-as.data.frame(dat)
  names(dat)<-c("pheno",'temp','rad','precip')
  lm_pred<-lm(pheno~temp+rad+precip,data = dat)
  r_square<-summary(lm_pred)$r.square
  coeff<-coefficients(lm_pred)  # intercep, temp, rad, precip.
  return(c(coeff,r_square))
}
mvr5<-function(x){
  ### input vector should have 16(year)*4(variables) elements
  ## of the following order: sos(eos), temp, rad, precip
  if (sum(is.na(x))>10){
    return(rep(NA,6))
  }
  dat<-x
  n<-length(dat)
  dim(dat)<-c(n/5,5)
  dat<-as.data.frame(dat)
  names(dat)<-c("pheno",'temp','rad','precip','pheno2')
  lm_pred<-lm(pheno~temp+rad+precip+pheno2,data = dat)
  r_square<-summary(lm_pred)$r.square
  coeff<-coefficients(lm_pred)  # intercep, temp, rad, precip.
  return(c(coeff,r_square))
}

## export to nc files
export_nc<-function(data,fileout){
  num_coef<-dim(data)[1]
  vardim<-ncdim_def(name = 'var_num',units = "",longname = "var_names as intercept, par, temp, prec",val=1:(num_coef-1))
  pcor_coef<-ncvar_def("pcor_coef",'',list(xdim,ydim,vardim),-999.9,prec="double",compression=9)
  pcor_pv<-ncvar_def("pcor_Rsquare",'',list(xdim,ydim),-999.9,prec="double",compression=9)
  ncout<-nc_create(fileout,list(pcor_coef,pcor_pv))
  pcoef<-t(data[1:(num_coef-1),])
  dim(pcoef)<-c(720,120,(num_coef-1))
  rsq<-data[num_coef,]
  dim(rsq)<-c(720,120)
  ncvar_put(ncout,pcor_coef,pcoef)
  ncvar_put(ncout,pcor_pv,rsq)
  nc_close(ncout)
}

## export to nc files
export_cor_nc<-function(data,fileout){
  #vardim<-ncdim_def(name = 'var_num',units = "",longname = "var_names as intercept, par, temp, prec",val=1:4)
  pcor_coef<-ncvar_def("pcor_coef",'',list(xdim,ydim),-999.9,prec="double",compression=9)
  #pcor_pv<-ncvar_def("pcor_Rsquare",'',list(xdim,ydim),-999.9,prec="double",compression=9)
  ncout<-nc_create(fileout,list(pcor_coef))
  dim(data)<-c(720,120)
  ncvar_put(ncout,pcor_coef,data)
  nc_close(ncout)
}

#### maybe I will use partial least squre model for the prediction. 


####


##read in the data and prepare for the regression
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")
sif_pheno<-list.files("./pheno_hd_fixed_threshold_clear/",full.names = T,pattern = ".nc")
sif_climate<-list.files("./pheno_hd_fixed_threshold_climate/clear/",full.names = T,pattern = ".nc")

sos<-array(NA,dim=c(86400,16))
eos<-array(NA,dim=c(86400,16))

pre_start_temp<-array(NA,dim=c(86400,16))
pre_end_temp<-array(NA,dim=c(86400,16))
pre_start_par<-array(NA,dim=c(86400,16))
pre_end_par<-array(NA,dim=c(86400,16))
pre_start_prec<-array(NA,dim=c(86400,16))
pre_end_prec<-array(NA,dim=c(86400,16))

for (year in 2001:2016){
  year_pheno<-sif_pheno[substr(basename(sif_pheno),35,38)==year]
  phenoin<-nc_open(year_pheno)
  sos[,year-2000]<-ncvar_get(phenoin,'SOS')*365+2.5
  eos[,year-2000]<-ncvar_get(phenoin,'EOS')*365+2.5
  xdim<-phenoin$dim[["longitude"]]
  ydim<-phenoin$dim[["latitude"]]
  nc_close(phenoin)
}

for (year in 2001:2016){
  year_climate<-sif_climate[substr(basename(sif_climate),1,4)==year]
  climatein<-nc_open(year_climate)
  
  pre_start_temp[,year-2000]<-ncvar_get(climatein,'pre_start_temp')
  pre_end_temp[,year-2000]<-ncvar_get(climatein,'pre_end_temp')

  pre_start_par[,year-2000]<-ncvar_get(climatein,'pre_start_par')
  pre_end_par[,year-2000]<-ncvar_get(climatein,'pre_end_par')

  pre_start_prec[,year-2000]<-ncvar_get(climatein,'pre_start_prec')
  pre_end_prec[,year-2000]<-ncvar_get(climatein,'pre_end_prec')
  nc_close(climatein)
}

# ### for the sos regression,  this time, we do not detrend for all variables
# sos_dat<-cbind(sos,pre_start_temp,pre_start_par,pre_start_prec)
# reg1<-apply(sos_dat,1,mvr)
# export_nc(reg1,"./analysis/multivariate_regression/sos_with_trend.nc")
# ######
# ### for the eos regression,  this time, we do not detrend for all variables
# eos_dat<-cbind(eos,pre_end_temp,pre_end_par,pre_end_prec)
# reg2<-apply(eos_dat,1,mvr)
# export_nc(reg2,"./analysis/multivariate_regression/eos_with_trend.nc")
# ######
# 
# ### for the sos regression,  this time, we do not detrend for all variables plus pre-year EOS
# sos_dat<-cbind(sos[,2:16],pre_start_temp[,2:16],pre_start_par[,2:16],pre_start_prec[,2:16],eos[,1:15])
# reg3<-apply(sos_dat,1,mvr5)
# export_nc(reg3,"./analysis/multivariate_regression/sos_with_eos.nc")
# ######
# ### for the eos regression,  this time, we do not detrend for all variables plus current-year SOS
# eos_dat<-cbind(eos,pre_end_temp,pre_end_par,pre_end_prec,sos)
# reg4<-apply(eos_dat,1,mvr5)
# export_nc(reg4,"./analysis/multivariate_regression/eos_with_sos.nc")
# ######
# 
# diff_sos<-reg3[6,]-reg1[5,]
# export_cor_nc(diff_sos,"./analysis/multivariate_regression/sos_diff.nc")
# diff_eos<-reg4[6,]-reg2[5,]
# export_cor_nc(diff_eos,"./analysis/multivariate_regression/eos_diff.nc")

# ######
library(pracma)
dsos<-apply(sos,1,detrend)
deos<-apply(eos,1,detrend)
dpre_start_temp<-apply(pre_start_temp,1,detrend)
dpre_end_temp<-apply(pre_end_temp,1,detrend)
dpre_start_par<-apply(pre_start_par,1,detrend)
dpre_end_par<-apply(pre_end_par,1,detrend)
dpre_start_prec<-apply(pre_start_prec,1,detrend)
dpre_end_prec<-apply(pre_end_prec,1,detrend)

### for the sos regression,  this time, we do not detrend for all variables
sos_dat<-rbind(dsos,dpre_start_temp,dpre_start_par,dpre_start_prec)
reg1<-apply(sos_dat,2,mvr)
export_nc(reg1,"./analysis/multivariate_regression/dsos_with_trend.nc")
######
### for the eos regression,  this time, we do not detrend for all variables
eos_dat<-rbind(deos,dpre_end_temp,dpre_end_par,dpre_end_prec)
reg2<-apply(eos_dat,2,mvr)
export_nc(reg2,"./analysis/multivariate_regression/deos_with_trend.nc")
######

### for the sos regression,  this time, we do not detrend for all variables plus pre-year EOS
sos_dat<-rbind(dsos[2:16,],dpre_start_temp[2:16,],dpre_start_par[2:16,],dpre_start_prec[2:16,],deos[1:15,])
reg3<-apply(sos_dat,2,mvr5)
export_nc(reg3,"./analysis/multivariate_regression/dsos_with_eos.nc")
######
### for the eos regression,  this time, we do not detrend for all variables plus current-year SOS
eos_dat<-rbind(deos,dpre_end_temp,dpre_end_par,dpre_end_prec,dsos)
reg4<-apply(eos_dat,2,mvr5)
export_nc(reg4,"./analysis/multivariate_regression/deos_with_sos.nc")
######

diff_sos<-reg3[6,]-reg1[5,]
export_cor_nc(diff_sos,"./analysis/multivariate_regression/dsos_diff.nc")
diff_eos<-reg4[6,]-reg2[5,]
export_cor_nc(diff_eos,"./analysis/multivariate_regression/deos_diff.nc")








