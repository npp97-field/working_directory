####### this is the site leve comparison of VI derived from multiple dataset and method 
####### MCD12Q2
setwd("/Users/yzhang/Project/SIF_phenology/")
mcd12files<-list.files("./retrieve_site/MCD12Q2/",full.names = T)

##site_list used 
site_list<-read.csv("./retrieve_site/sites_2000.csv",stringsAsFactors = F)
phenology<-read.csv("./retrieve_site/analysis/site_phenology_0.3_N30.csv",stringsAsFactors = F)
phenology_mcd<-phenology[!is.na(phenology$site),1:6]
names(phenology_mcd)[5:6]<-c("sos_mcd12","eos_mcd12")
sites_all<-unique(phenology_mcd$site)
for (i in 1:length(sites_all)){
  site_name<-sites_all[i]
  mcd_data<-read.csv(mcd12files[substr(basename(mcd12files),1,6)==site_name],stringsAsFactors = F)
  years<-unique(phenology_mcd$year[phenology_mcd$site==site_name])
  #mcd_data[is.na(mcd_data)]<--10000
  for (y in years){
    time_critcal<-as.POSIXct("2000-01-01 00:00")+
      as.difftime(as.numeric(mcd_data[as.numeric(substr(mcd_data$date,1,4))==y,4:11]*24),units="hours")
    time_as_doy<-as.numeric(format(time_critcal,"%j"))
    if (is.na(time_as_doy[2])){
      phenology_mcd$sos_mcd12[phenology_mcd$site==site_name&phenology_mcd$year==y]<-time_as_doy[6]
    }else{
      phenology_mcd$sos_mcd12[phenology_mcd$site==site_name&phenology_mcd$year==y]<-time_as_doy[2]
    }
    if(is.na(time_as_doy[4])){
      phenology_mcd$eos_mcd12[phenology_mcd$site==site_name&phenology_mcd$year==y]<-time_as_doy[8]
    }else{
      phenology_mcd$eos_mcd12[phenology_mcd$site==site_name&phenology_mcd$year==y]<-time_as_doy[4]
    }
  }
}
write.csv(phenology_mcd,"./retrieve_site/analysis/site_phenology_mcd_N30.csv",row.names = F)
