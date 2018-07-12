###
# calcualte the monthly OCO-2 SIF for each biome type and 
# 4-day CSIF for each biome type

library(ncdf4)
setwd("/rigel/glab/users/zy2309/")
ocosif_f<-list.files("./DATA/OCO-2_v8/grid/monthly_clear_no_orbit_HD/",full.names = T)
#ocosif_f<-list.files("/Users/yzhang/Project/SIF_phenology/data/monthly_clear_no_orbit_HD/",full.names = T)
lc_nc<-nc_open("./PROJECT/SIF_phenology/analysis/North_mcd12c1_landcover1_majority.nc")
#lc_nc<-nc_open("/Users/yzhang/Project/SIF_phenology/data/North_mcd12c1_landcover1_majority.nc")
lc_data<-ncvar_get(lc_nc,varid = "lc1")
nc_close(lc_nc)
mask_nc<-nc_open("./PROJECT/SIF_phenology/analysis/North_barren_mask.nc")
mask_data<-ncvar_get(mask_nc,varid="barren")
# mask_nc<-nc_open("./PROJECT/SIF_phenology/analysis/global_land_mask.nc")
# mask_data<-ncvar_get(mask_nc,varid="mask")
nc_close(mask_nc)
lc_data1<- lc_data
lc_data1[lc_data>=1&lc_data<=5]<-1  #forest
lc_data1[lc_data>=6&lc_data<=8]<-2  # woody land
lc_data1[lc_data>=9&lc_data<=11]<-3 # grass wet sav
lc_data1[lc_data==12|lc_data==14]<-4 #cropland
lc_data1[lc_data1>4|lc_data<1]<-NA

##### read the oco-2 sif data for each month from 2014 to 2017
oco_time_series<-as.data.frame(array(NA,dim=c(12*4,11)))
names(oco_time_series)<-c("yearmonth","all","Forest","Woodland","Grassland","Cropland",
                          "all_sd","Forest_sd","Woodland_sd","Grassland_sd","Cropland_sd")
ym<-(as.numeric(substr(basename(ocosif_f),12,13))-14)*12+as.numeric(substr(basename(ocosif_f),14,15))
oco_time_series[,1]<-rep(2014:2017,each=12)*100+rep(1:12,4)
for (i in 1:length(ocosif_f)){
  oco_nc<-nc_open(ocosif_f[i])
  oco_dat<-ncvar_get(oco_nc,varid = "sif_757nm_daily")
  nc_close(oco_nc)
  north_oco<-oco_dat[,241:360]
  masked_oco<-north_oco*mask_data
  oco_time_series[ym[i],2]<-mean(masked_oco,na.rm=T)
  oco_time_series[ym[i],7]<-sd(masked_oco,na.rm=T)
  for (j in 1:4){
    biome_mask<-lc_data1
    biome_mask[biome_mask!=j]=NA
    lc_dat<-masked_oco*biome_mask
    oco_time_series[ym[i],j+2]<-mean(lc_dat,na.rm=T)
    oco_time_series[ym[i],j+7]<-sd(lc_dat,na.rm=T)
  }
}
write.csv(oco_time_series,"./PROJECT/SIF_phenology/analysis/csif_stat/OCO_2_SIF.csv",row.names = F)
##### for csif oco-2 data was used to mask all the pixels for comparisons
csif_f<-list.files("./DATA/bg_csif/clear_inst_SIF_4day_HD_BISE/",full.names = T)[14:17]
csif_time_series<-as.data.frame(array(NA,dim=c(92*4,11)))
names(csif_time_series)<-c("yearmonth","all","Forest","Woodland","Grassland","Cropland",
                          "all_sd","Forest_sd","Woodland_sd","Grassland_sd","Cropland_sd")

for (i in 1:4){
  csif_nc<-nc_open(csif_f[i])
  csif_year<-ncvar_get(csif_nc,varid='clear_daily_sif')
  nc_close(csif_nc)
  for (m in 1:92){
    # oco_f<-ocosif_f[ym==i*12-12+ceiling(m*4/30.7)]
    # if (length(oco_f)==0){
    #   next
    # }
    # oco_nc<-nc_open(oco_f)
    # oco_dat<-ncvar_get(oco_nc,varid = "sif_757nm_daily")[,241:360]
    # nc_close(oco_nc)
    # oco_dat[!is.na(oco_dat)]=1
    # 
    north_csif<-csif_year[,241:360,m]
    masked_csif<-north_csif*mask_data
    csif_time_series[i*92-92+m,2]<-mean(masked_csif,na.rm=T)
    csif_time_series[i*92-92+m,7]<-sd(masked_csif,na.rm=T)
    for (j in 1:4){
      biome_mask<-lc_data1
      biome_mask[biome_mask!=j]=NA
      lc_dat<-masked_csif*biome_mask
      csif_time_series[i*92-92+m,j+2]<-mean(lc_dat,na.rm=T)
      csif_time_series[i*92-92+m,j+7]<-sd(lc_dat,na.rm=T)
    }
  }
}
write.csv(csif_time_series,"./PROJECT/SIF_phenology/analysis/csif_stat/CSIF.csv",row.names = F)

##### EVI and NDVI
vi_f<-list.files("./DATA/MOD13C1_VI_HD/",full.names = T)[14:17]
ndvi_time_series<-as.data.frame(array(NA,dim=c(23*4,11)))
evi_time_series<-as.data.frame(array(NA,dim=c(23*4,11)))
names(ndvi_time_series)<-c("yearmonth","all","Forest","Woodland","Grassland","Cropland",
                           "all_sd","Forest_sd","Woodland_sd","Grassland_sd","Cropland_sd")
names(evi_time_series)<-c("yearmonth","all","Forest","Woodland","Grassland","Cropland",
                           "all_sd","Forest_sd","Woodland_sd","Grassland_sd","Cropland_sd")
for (i in 1:4){
  vi_nc<-nc_open(vi_f[i])
  ndvi_year<-ncvar_get(vi_nc,varid='NDVI')
  evi_year<-ncvar_get(vi_nc,varid='EVI')
  nc_close(vi_nc)
  for (m in 1:23){
    # oco_f<-ocosif_f[ym==i*12-12+ceiling(m*16/30.7)]
    # if (length(oco_f)==0){
    #   next
    # }
    # oco_nc<-nc_open(oco_f)
    # oco_dat<-ncvar_get(oco_nc,varid = "sif_757nm_daily")[,241:360]
    # nc_close(oco_nc)
    # oco_dat[!is.na(oco_dat)]=1
    
    north_ndvi<-ndvi_year[m,,241:360]
    north_ndvi[north_ndvi<0]<-NA
    masked_ndvi<-north_ndvi*mask_data
    north_evi<-evi_year[m,,241:360]
    north_evi[north_evi<0]<-NA
    masked_evi<-north_evi*mask_data
    ndvi_time_series[i*23-23+m,2]<-mean(masked_ndvi,na.rm=T)
    ndvi_time_series[i*23-23+m,7]<-sd(masked_ndvi,na.rm=T)
    evi_time_series[i*23-23+m,2]<-mean(masked_evi,na.rm=T)
    evi_time_series[i*23-23+m,7]<-sd(masked_evi,na.rm=T)
    for (j in 1:4){
      biome_mask<-lc_data1
      biome_mask[biome_mask!=j]=NA
      lc_dat<-masked_ndvi*biome_mask
      ndvi_time_series[i*23-23+m,j+2]<-mean(lc_dat,na.rm=T)
      ndvi_time_series[i*23-23+m,j+7]<-sd(lc_dat,na.rm=T)
      lc_dat<-masked_evi*biome_mask
      evi_time_series[i*23-23+m,j+2]<-mean(lc_dat,na.rm=T)
      evi_time_series[i*23-23+m,j+7]<-sd(lc_dat,na.rm=T)
    }
  }
}
write.csv(ndvi_time_series,"./PROJECT/SIF_phenology/analysis/csif_stat/ndvi.csv",row.names = F)
write.csv(evi_time_series,"./PROJECT/SIF_phenology/analysis/csif_stat/evi.csv",row.names = F)

