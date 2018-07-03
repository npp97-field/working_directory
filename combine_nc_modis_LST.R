######
library(ncdf4)
library(rgdal)
modis_lst<-list.files("/rigel/glab/users/zy2309/DATA/MYD11C2_HD/",full.names=T)
all_data<-array(NA,dim=c(720,360,644))
for (i in 2:15){
  ncf<-nc_open(modis_lst[i])
  dat<-ncvar_get(ncf,"LST_day")
  all_data[,,(i*46-91):(i*46-46)]<-dat
  nc_close(ncf)
}

doy<-1:644
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
ncout<-nc_create('/rigel/glab/users/zy2309/PROJECT/SIF_phenology/data/MODIS_LST.nc',nclst)
ncvar_put(ncout,varid=nclst,all_data)
nc_close(ncout)