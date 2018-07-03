######
library(ncdf4)
#library(rgdal)
modis_lst<-list.files("/rigel/glab/users/zy2309/DATA/MYD11C2_HD/",full.names=T)
all_data<-array(NA,dim=c(720,360,690))
all_datan<-array(NA,dim=c(720,360,690))
for (i in 1:15){
  ncf<-nc_open(modis_lst[i])
  dat<-ncvar_get(ncf,"LST_day")
  datn<-ncvar_get(ncf,"LST_night")
  all_data[,,(i*46-45):(i*46)]<-dat
  all_datan[,,(i*46-45):(i*46)]<-datn
  nc_close(ncf)
}

doy<-1:690
latmin<- -90
latmax<- 90
latd<- 0.5
lonmin<- -180
lonmax<- 180
lond<- 0.5

lat<- seq(latmin+latd/2,latmax-latd/2,latd)
long<-seq(lonmin+lond/2,lonmax-lond/2,lond)

dimlat<-ncdim_def('latitude','deg',lat)
dimlong<-ncdim_def('longitude','deg',long)
dimdoy<-ncdim_def("doy","",doy)
nclst<-ncvar_def('LST_Day','degree C',list(dimlong,dimlat,dimdoy),-9999,longname="LST_Day",prec='float',compression=9)
nclstn<-ncvar_def('LST_Night','degree C',list(dimlong,dimlat,dimdoy),-9999,longname="LST_Night",prec='float',compression=9)
ncout<-nc_create('/rigel/glab/users/zy2309/PROJECT/SIF_phenology/data/MODIS_LST.nc',list(nclst,nclstn))
ncvar_put(ncout,varid=nclst,all_data)
ncvar_put(ncout,varid=nclstn,all_datan)
nc_close(ncout)