##### data
###### calcualte the aridity index based on CRU-ts PET and precipitation
library(ncdf4)
library(raster)

setwd("/rigel/glab/users/zy2309/DATA/CRU_TS401/")


get_annual_average<-function(varname){
  var_files<-list.files("./",pattern=varname,full.names = T)[2:3]
  ncin<-nc_open(var_files[1])
  var_data1<-ncvar_get(ncin,varid=varname)[,241:360,]
  dim(var_data1)<-c(120*720,dim(var_data1)[3])
  nc_close(ncin)
  ncin<-nc_open(var_files[2])
  var_data2<-ncvar_get(ncin,varid=varname)[,241:360,]
  dim(var_data2)<-c(120*720,dim(var_data2)[3])
  nc_close(ncin)
  var_dat<-cbind(var_data1,var_data2)
  mean_var<-apply(var_dat,1,mean)
  return(mean_var)
}



export_nc<-function(outdata,outfile){
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
  
  ncpet<-ncvar_def('pet','mm/day',list(dimlong,dimlat),-9999,longname="potential evaporation",prec='float',compression=9)
  ncpre<-ncvar_def('pre','mm/month',list(dimlong,dimlat),-9999,longname="precipiation",prec='float',compression=9)
  ncai<-ncvar_def('ai','',list(dimlong,dimlat),-9999,longname="aridity index",prec='float',compression=9)
  
  if(file.exists(outfile)){
    file.remove(outfile)
  }
  
  ncout<-nc_create(outfile,list(ncpet,ncpre,ncai))
  
  dim(outdata)<-c(720,120,3)
  ncvar_put(ncout,varid=ncpet,outdata[,,1])
  ncvar_put(ncout,varid=ncpre,outdata[,,2])
  ncvar_put(ncout,varid=ncai,outdata[,,3])
  
  nc_close(ncout)
}

pet<-get_annual_average("pet")
pre<-get_annual_average("pre")
ai<-pet/pre/30*100
#dim(ai)<-c(120,720)
out_put<-cbind(pet,pre,ai)

outfile<-"/rigel/glab/users/zy2309/PROJECT/SIF_phenology/analysis/climate/ai_cru.nc"
export_nc(out_put,outfile)


