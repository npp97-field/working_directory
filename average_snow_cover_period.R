#### get average snow cover period.
library(ncdf4)
setwd("/rigel/glab/users/zy2309/")
snow_files<-list.files("./DATA/MOD10C1_HD/",full.names = T)[2:17]

identidy_snow_period<-function(vec){
  if (sum(is.na(vec))>180){
    return(rep(NA,2))
  }
  vec[c(1,365)]<-100
  ## the first half is 1-210 corresponding to 1:150
  first<-rev(vec[1:210])
  ## the second half is 211-365 corresponding to 151:230
  second<-vec[211:365]
  start<-211-which(first>20)[1]
  end<-which(second>5)[1]+210
  return(c(start,end))
}

annual_snow_end<-array(NA,dim=c(720,120,16))
annual_snow_start<-array(NA,dim=c(720,120,16))
for (i in 1:16){
  ncin<-nc_open(snow_files[i])
  snow_day<-ncvar_get(ncin,"snow_cover")[,241:360,]
  annual_snow<-apply(snow_day,c(1,2),identidy_snow_period)
  annual_snow_end[,,i]<-annual_snow[1,,]
  annual_snow_start[,,i]<-annual_snow[2,,]
}

average_snow_end<-apply(annual_snow_end,c(1,2),mean,na.rm=T)
average_snow_start<-apply(annual_snow_start,c(1,2),mean,na.rm=T)

save(average_snow_end,average_snow_start,annual_snow_end,annual_snow_start,
     file="./PROJECT/SIF_phenology/analysis/snow_data.RData")

