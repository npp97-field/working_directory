###### resample and get average temp and aridity index
library(ncdf4)
library(raster)


### read the AI data from CRU-TS 4.0
ncin<-nc_open("/Users/yzhang/Project/SIF_phenology/data/ai_cru_P_PET.nc")
current_AI<-ncvar_get(ncin,varid = "ai")
AI_cur<-raster(apply(current_AI,1,rev))
extent(AI_cur)<-c(-180,180,30,90)
lon<-ncin$dim[["longitude"]]
lat<-ncin$dim[['latitude']]

### read the global land mask
ncin<-nc_open("/Users/yzhang/Project/SIF_phenology/data/global_land_mask.nc")
lc_mask<-ncvar_get(ncin,varid = "mask")
lmask_ras<-raster(apply(lc_mask,1,rev))
extent(lmask_ras)<-c(-180,180,-90,90)
nc_close(ncin)
###function to get temperature data and AI data and resample to 0.5 degree resolution

resample_to_hd<-function(file){
  ncin<-nc_open(file)
  var1<-ncin$var[[1]]
  dat<-ncvar_get(ncin,varid = var1)
  lat<-ncin$dim[["lat"]]$vals
  lon<-ncin$dim[["lon"]]$vals
  dimext<-c(length(lon),length(lat))
  lon<-c(lon[(dimext[1]/2+1):(dimext[1])]-360,lon[1:(dimext[1]/2)])
  dat<-rbind(dat[(dimext[1]/2+1):dimext[1],],dat[1:(dimext[1]/2),])
  dlat<-lat[2]-lat[1]
  dlon<-lon[2]-lon[1]
  lonmin<-lon[1]-dlon/2
  lonmax<-lon[dimext[1]]+dlon/2
  latmax<-max(lat)+dlat/2
  latmin<-min(lat)-dlat/2
  if (length(grep("CMCC-CESM",file))==1){
    ras_dat<-raster(apply(dat,1,rev))
  }else{
    ras_dat<-raster(apply(dat,1,rev))
  }
  extent(ras_dat)<-c(lonmin,lonmax,latmin,latmax)
  model_AI<-resample(ras_dat,lmask_ras,'bilinear')
  return(model_AI)
}

model_AIs<-list()
files<-list.files("/Users/yzhang/Project/SIF_phenology/analysis/CMIP5_AI/",pattern="AI",full.names = T)[c(1,3:7)]
for (i in 1:6){
  model_AIs[[i]]<-resample_to_hd(files[i])
}

AI_stack<-stack(model_AIs)
mme_AI<-mean(AI_stack)
landAI<-mme_AI*lmask_ras
northAI<-crop(landAI,c(-180,180,30,90))
northai_val<-getValues(northAI)
dim(northai_val)<-c(720,120)
northai_val<-t(apply(northai_val,1,rev))

futureAI<-ncvar_def('rcp85ai','',dim = list(lon,lat),longname = "aridity index from multi model ensemble mean")
currentAI<-ncvar_def('currentai','',dim = list(lon,lat),longname = "aridity index from CRU-TS 4.0")
nc_out<-nc_create("/Users/yzhang/Project/SIF_phenology/analysis/future_AI.nc",list(futureAI,currentAI))
ncvar_put(nc_out,futureAI,northai_val)
ncvar_put(nc_out,currentAI,current_AI)
nc_close(nc_out)

model_tas<-list()
files<-list.files("/Users/yzhang/Project/SIF_phenology/analysis/CMIP5_AI/",pattern="tas_",full.names = T)[c(1,3:7)]
for (i in 1:6){
  model_tas[[i]]<-resample_to_hd(files[i])
}


################### mean annual temperature
ncin<-nc_open("/Users/yzhang/Project/SIF_phenology/data/mean_annual_temp.nc")
current_t<-ncvar_get(ncin,varid = "mean_annual_temp")
t_cur<-raster(apply(current_t,1,rev))
extent(AI_cur)<-c(-180,180,30,90)
lon<-ncin$dim[["longitude"]]
lat<-ncin$dim[['latitude']]



t_stack<-stack(model_tas)
mme_t<-mean(t_stack)
landt<-mme_t*lmask_ras
northt<-crop(landt,c(-180,180,30,90))
northt_val<-getValues(northt)
dim(northt_val)<-c(720,120)
northt_val<-t(apply(northt_val,1,rev))

futuret<-ncvar_def('rcp85_tas','',dim = list(lon,lat),longname = "mean annual temperature from multi model ensemble mean")
currentt<-ncvar_def('current_tas','',dim = list(lon,lat),longname = "mean annual temperature from CRU-TS 4.0")
nc_out<-nc_create("/Users/yzhang/Project/SIF_phenology/analysis/future_t.nc",list(futuret,currentt))
ncvar_put(nc_out,futuret,northai_val)
ncvar_put(nc_out,currentt,current_t)
nc_close(nc_out)

