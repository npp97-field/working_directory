##### orgnize the CMIP5 data
library(ncdf4)
library(abind)
setwd("/rigel/glab/users/zy2309/")


get_historical<-function(f_name,varname){
  if (length(f_name)==1){
    varf<-nc_open(f_name)
    var<-ncvar_get(varf,varname)
    nmon<-dim(var)[3]
  }else{
    for (i in 1:length(f_name)){
      varf<-nc_open(f_name[i])
      vartemp<-ncvar_get(varf,varname)
      xyzdim<-dim(var)
      if (i == 1){
        var<-vartemp
      }else{
        var<-abind(var,vartemp,along=3)
      }
    }
    nmon<-dim(var)[3]
  }
  if (basename(dirname(f_name[1]))=="CMCC-CESM"){
    ydim<-ncdim_def("lat",units = "deg",vals = (24:-23)*3.75-1.875)
    xdim<-ncdim_def("lon",units = "deg",vals = (0:95)*3.75)
  }else{
    ydim<-varf$dim[["lat"]]
    xdim<-varf$dim[["lon"]]
  }
  

  if (basename(dirname(f_name[1]))=="HadGEM2-ES"){
    offset = 12
  }else{
    offset = 11
  }
  
  nc_close(varf)
  if (varname=="gpp"){
    which
    selvar<-var[,,(nmon-239-offset):(nmon-offset)]*1000*86400 ##convert to g C m-2 day-1
    unitvar<-"g C m-2 day-1"
    lname<-"gross primary production"
  }else if(varname=="pr"){
    selvar<-var[,,(nmon-239-offset):(nmon-offset)]*86400 ##convert to mm day-1
    unitvar<-"mm day-1"
    lname<-"total precipitation"
  }else if(varname=='tas'){
    selvar<-var[,,(nmon-239-offset):(nmon-offset)]+273.16 ##convert deg C
    unitvar<-"deg C"
    lname<-"surface air temperature"
  }else if(varname=='tasmax'){
    selvar<-var[,,(nmon-239-offset):(nmon-offset)]+273.16 ##convert deg C
    unitvar<-"deg C"
    lname<-"maximum surface air temperature"
  }

  time_ym<-ncdim_def(name = "yearmonth",units = "yyyymm",longname = "year and month",
                     vals = rep(1985:2004,each=12)*100+rep(1:12,20))
  var_out<-ncvar_def(name = varname, units = unitvar, missval = -9999,longname = lname,dim = list(xdim,ydim,time_ym))
  nc_out<-nc_create(paste("./PROJECT/SIF_phenology/CMIP5/historical/",
                  varname,"_",basename(dirname(f_name[1])),"_historical_198501_200412.nc",sep=""),var_out)
  ncvar_put(nc_out,varid = var_out,selvar)
  nc_close(nc_out)
}

####-------------   CanESM2
varf<-list.files("./DATA/CMIP5/CanESM2/",pattern="gpp_.*historical",full.names = T)
get_historical(varf,"gpp")
varf<-list.files("./DATA/CMIP5/CanESM2/",pattern="pr_.*historical",full.names = T)
get_historical(varf,"pr")
varf<-list.files("./DATA/CMIP5/CanESM2/",pattern="tas_.*historical",full.names = T)
get_historical(varf,"tas")
varf<-list.files("./DATA/CMIP5/CanESM2/",pattern="tasmax_.*historical",full.names = T)
get_historical(varf,"tasmax")

####-------------   GFDL-ESM2M  
varf<-list.files("./DATA/CMIP5/GFDL-ESM2M/",pattern="gpp_.*historical",full.names = T)
get_historical(varf[25:29],"gpp")
varf<-list.files("./DATA/CMIP5/GFDL-ESM2M/",pattern="pr_.*historical",full.names = T)
get_historical(varf[25:29],"pr")
varf<-list.files("./DATA/CMIP5/GFDL-ESM2M/",pattern="tas_.*historical",full.names = T)
get_historical(varf[25:29],"tas")
varf<-list.files("./DATA/CMIP5/GFDL-ESM2M/",pattern="tasmax_.*historical",full.names = T)
get_historical(varf[25:29],"tasmax")

####-------------   HadGEM2-ES 
varf<-list.files("./DATA/CMIP5/HadGEM2-ES/",pattern="gpp_.*historical",full.names = T)
get_historical(varf,"gpp")
varf<-list.files("./DATA/CMIP5/HadGEM2-ES/",pattern="pr_.*historical",full.names = T)
get_historical(varf,"pr")
varf<-list.files("./DATA/CMIP5/HadGEM2-ES/",pattern="tas_.*historical",full.names = T)
get_historical(varf,"tas")
varf<-list.files("./DATA/CMIP5/HadGEM2-ES/",pattern="tasmax_.*historical",full.names = T)
get_historical(varf,"tasmax")

####-------------   IPSL-CM5A-LR
varf<-list.files("./DATA/CMIP5/IPSL-CM5A-LR/",pattern="gpp_.*historical",full.names = T)
get_historical(varf,"gpp")
varf<-list.files("./DATA/CMIP5/IPSL-CM5A-LR/",pattern="pr_.*historical",full.names = T)
get_historical(varf,"pr")
varf<-list.files("./DATA/CMIP5/IPSL-CM5A-LR/",pattern="tas_.*historical",full.names = T)
get_historical(varf,"tas")
varf<-list.files("./DATA/CMIP5/IPSL-CM5A-LR/",pattern="tasmax_.*historical",full.names = T)
get_historical(varf,"tasmax")

####-------------   MIROC-ESM
varf<-list.files("./DATA/CMIP5/MIROC-ESM/",pattern="gpp_.*historical",full.names = T)
get_historical(varf,"gpp")
varf<-list.files("./DATA/CMIP5/MIROC-ESM/",pattern="pr_.*historical",full.names = T)
get_historical(varf,"pr")
varf<-list.files("./DATA/CMIP5/MIROC-ESM/",pattern="tas_.*historical",full.names = T)
get_historical(varf,"tas")
varf<-list.files("./DATA/CMIP5/MIROC-ESM/",pattern="tasmax_.*historical",full.names = T)
get_historical(varf,"tasmax")

####-------------   MPI-ESM2-LR
varf<-list.files("./DATA/CMIP5/MPI-ESM2-LR/",pattern="gpp_.*historical",full.names = T)
get_historical(varf,"gpp")
varf<-list.files("./DATA/CMIP5/MPI-ESM2-LR/",pattern="pr_.*historical",full.names = T)
get_historical(varf,"pr")
varf<-list.files("./DATA/CMIP5/MPI-ESM2-LR/",pattern="tas_.*historical",full.names = T)
get_historical(varf,"tas")
varf<-list.files("./DATA/CMIP5/MPI-ESM2-LR/",pattern="tasmax_.*historical",full.names = T)
get_historical(varf,"tasmax")

####-------------   NorESM1-ME
varf<-list.files("./DATA/CMIP5/NorESM1-ME/",pattern="gpp_.*historical",full.names = T)
get_historical(varf,"gpp")
varf<-list.files("./DATA/CMIP5/NorESM1-ME/",pattern="pr_.*historical",full.names = T)
get_historical(varf,"pr")
varf<-list.files("./DATA/CMIP5/NorESM1-ME/",pattern="tas_.*historical",full.names = T)
get_historical(varf,"tas")
varf<-list.files("./DATA/CMIP5/NorESM1-ME/",pattern="tasmax_.*historical",full.names = T)
get_historical(varf,"tasmax")

####-------------   CMCC-CESM
varf<-list.files("./DATA/CMIP5/CMCC-CESM/",pattern="gpp_.*historical",full.names = T)
get_historical(varf,"gpp")
varf<-list.files("./DATA/CMIP5/CMCC-CESM/",pattern="pr_.*historical",full.names = T)
get_historical(varf[28:32],"pr")
varf<-list.files("./DATA/CMIP5/CMCC-CESM/",pattern="tas_.*historical",full.names = T)
get_historical(varf[28:32],"tas")
varf<-list.files("./DATA/CMIP5/CMCC-CESM/",pattern="tasmax_.*historical",full.names = T)
get_historical(varf[28:32],"tasmax")

