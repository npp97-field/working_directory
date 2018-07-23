###### calculate mean annual temperature and aridity index during 2070 to 2100
library(ncdf4)
models <-c("CanESM2","CMCC-CESM","GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR",
           "MIROC-ESM","MPI-ESM2-LR")
# first, list all files for one model
setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology/")


#### calculate delta
calc_delta<-function(ta,SVP){
  delta<-4098*SVP/(ta+237.3)^2
  return(delta)
}
calc_svp<-function(ta){
  SVP<-0.6108*10^(7.5*ta/(237.3+ta))
  return(SVP)
}

#calculate VPD
calc_vpd<-function(rh,SVP){
  VPD<-((100-rh)/100)*SVP
  return(VPD)
}

#convert wind speed
calculate_PET<-function(u10,Rn,ta,rh){
  u2=u10*(log(128)/log(661.3))
  gamma_const <- 66
  SVP<-calc_svp(ta)
  vpd<-calc_vpd(rh,SVP)
  delta<-calc_delta(ta,SVP)
  PET<-(0.408*delta*Rn+gamma_const*900/(ta+273.16)*u2*vpd)/
    (delta+gamma_const*(1+0.34*u2))
  return(PET)
}


calculate_model_pet<-function(model){
  files<-list.files("./CMIP5/rcp85/",pattern=model,full.names = T)
  ncin<-nc_open(files[1])
  hfls<-ncvar_get(ncin,'hfls')
  ncin<-nc_open(files[2])
  hfss<-ncvar_get(ncin,'hfss')
  ncin<-nc_open(files[3])
  v1<-ncin$var[[1]]
  hurs<-ncvar_get(ncin,v1)
  ncin<-nc_open(files[5])
  sfcWind<-ncvar_get(ncin,'sfcWind')
  ncin<-nc_open(files[6])
  tas<-ncvar_get(ncin,'tas')
  lon<-ncin$dim[["lon"]]
  lat<-ncin$dim[["lat"]]
  yearmonth<-ncin$dim[["yearmonth"]]
  nc_close(ncin)
  PET<-calculate_PET(sfcWind,hfls+hfss,tas,hurs)
  pet<-ncvar_def("pet",'mm/day',list(lon,lat,yearmonth),prec="float",
                 longname = "potential evaportranspiration",compression = 9)
  ncout<-nc_create(paste("./analysis/CMIP5_ANA/CMIP5_PET/rcp85_",model,'_2071_2100.nc',sep=""),pet)
  ncvar_put(ncout,varid = "pet",PET)
  nc_close(ncout)
  return(PET)
}

get_annual_values<-function(var){
  vdim<-dim(var)
  dim(var)<-c(vdim[1:2],12,30)
  annual_var<-apply(var,c(1,2,4),mean)
  return(annual_var)
}

get_multi_year<-function(annual_var){
  vdim<-dim(annual_var)
  multiyear_avg<-apply(annual_var,c(1,2),mean)
  return(multiyear_avg)
}

calculate_AI<-function(model){
  files<-list.files("./CMIP5/rcp85/",pattern=model,full.names = T)
  pet<-calculate_model_pet(model)
  ncin<-nc_open(files[4])
  pr<-ncvar_get(ncin,'pr')
  lon<-ncin$dim[["lon"]]
  lat<-ncin$dim[["lat"]]
  year<-ncdim_def(name = "year",units = '',vals = 2071:2100)
  annual_pr<-get_annual_values(pr)
  annual_pet<-get_annual_values(pet)
  avg_pr<-get_multi_year(annual_pr)
  avg_pet<-get_multi_year(annual_pet)
  annual_AI<-annual_pr/annual_pet
  average_AI<-avg_pr/avg_pet
  mean_ai<-ncvar_def("avg_AI",'',list(lon,lat),prec="float",
                     longname = "multi-year average Aridity index",compression = 9)
  annual_ai<-ncvar_def("annual_AI",'',list(lon,lat,year),prec="float",
                     longname = "yearly Aridity index",compression = 9)
  nc_out<-nc_create(filename = paste("./analysis/CMIP5_ANA/CMIP5_AI/rcp85_AI_",model,'_2071_2100.nc',sep=""),
                    list(mean_ai,annual_ai))
  ncvar_put(nc_out,varid=mean_ai,average_AI)
  ncvar_put(nc_out,varid=annual_ai,annual_AI)
  nc_close(nc_out)
}

for (i in 1:7){
  calculate_AI(models[i])
}



