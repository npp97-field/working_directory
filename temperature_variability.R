
library(raster)
library(ncdf4)
library(rgdal)

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
### calculate the mean annual temperature
if (F){
  mean_annual_temprature<-apply(tmean,1,mean,na.rm=T)
  lat = seq(30.25,89.75,0.5)
  long = seq(-179.75,179.75,0.5)
  
  ydim<- ncdim_def("latitude",units = "deg",lat)
  xdim<- ncdim_def('longitude',units='deg',long)
  mean_temp<-ncvar_def("mean_annual_temp",'',list(xdim,ydim),-9999,compression=9)
  temp_nc<-nc_create("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/analysis/mean_annual_temp.nc",list(mean_temp))
  ncvar_put(temp_nc,mean_temp,mean_annual_temprature)
  nc_close(temp_nc)
}
###########---------------------------------

precip<-get_cru_var('pre')

temp_d<-apply(tmean,1,spline_interpolate)
temp_day<-t(temp_d)
dim(sd_var)<-c(720,120,12)
dim(mean_var)<-c(720,120,12)


radnc<-nc_open("/rigel/glab/users/zy2309/DATA/bess_hd_monthly_PAR.nc")
rad<-ncvar_get(radnc,"PAR")
par<-rad[,,241:360]
dim(par)<-c(204,86400)
par_d<-apply(t(par),1,spline_interpolate)
par_day<-t(par_d)


#var<-tmean
var<-t(par)

sd_var<-array(NA,dim=c(86400,12))
mean_var<-array(NA,dim=c(86400,12))

for (i in 1:12){
  mon_var<-var[,i+0:16*12]
  sd_var[,i]<-apply(mon_var,1,sd,na.rm=T)
  mean_var[,i]<-apply(mon_var,1,mean,na.rm=T)
}
dim(sd_var)<-c(720,120,12)
dim(mean_var)<-c(720,120,12)

lat = seq(30.25,89.75,0.5)
long = seq(-179.75,179.75,0.5)
mon = 1:12
ydim<- ncdim_def("lat",units = "deg",lat)
xdim<- ncdim_def('long',units='deg',long)
zdim<- ncdim_def("mon",units='',mon)

# mean_temp<-ncvar_def("mean_temp",'',list(xdim,ydim,zdim),-9999,compression=9)
# sd_temp<-ncvar_def("sd_temp",'',list(xdim,ydim,zdim),-9999,compression=9)
# 
# temp_nc<-nc_create("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/analysis/monthly_temp.nc",list(mean_temp,sd_temp))
# ncvar_put(temp_nc,sd_temp,sd_var)
# ncvar_put(temp_nc,mean_temp,mean_var)
# nc_close(temp_nc)

mean_par<-ncvar_def("mean_par",'',list(xdim,ydim,zdim),-9999,compression=9)
sd_par<-ncvar_def("sd_par",'',list(xdim,ydim,zdim),-9999,compression=9)

temp_nc<-nc_create("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/analysis/monthly_par.nc",list(mean_par,sd_par))
ncvar_put(temp_nc,sd_par,sd_var)
ncvar_put(temp_nc,mean_par,mean_var)
nc_close(temp_nc)



