####### this is the site leve comparison of VI derived from multiple dataset and method 
####### MCD12Q2
setwd("/Users/yzhang/Project/SIF_phenology/")
mcd12files<-list.files("./retrieve_site/MCD12Q2/",full.names = T)

##site_list used 
####################################################################
########   this is for the MCD12Q2 data analysis
#################
site_list<-read.csv("./retrieve_site/sites_2000.csv",stringsAsFactors = F)
phenology<-read.csv("./retrieve_site/analysis/site_phenology_0.3_N30_all_daily.csv",stringsAsFactors = F)
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

##site_list used 
####################################################################
########   this is for the NDVI based PM method

library(zoo)
library(outliers)
##############################################
remove_outliers<-function(vits){
  numobs<-length(vits)/23
  qtiles<-c()
  for (i in 1:10){
    qtiles[i]<-quantile(vits,i*0.02,na.rm=T)
  }
  qtilediff<-qtiles[2:10]-qtiles[1:9]
  thresh_hold<-(qtiles[max(which(qtilediff>0.1))+1])
  vits[vits<thresh_hold]<-NA
  return(vits)
}



####### retrieve pehnology using the P-M method for each site
get_two_threshold<-function(vits){
  ### get the mean seasonal cycle for each pixel.
  numobs<-length(vits)/23
  #min_thresh<-quantile(vits,0.1)
  #vits[vits< min_thresh]<-NA
  if (sum(is.na(vits))>20*numobs){
    return(c(NA,NA))
  }
  ts_fill<-na.fill(vits,fill = "extend")
  dim(ts_fill)<-c(23,numobs)
  meants<-apply(ts_fill,1,median,na.rm=T)
  ### calcualte the NDVI ratio
  ndvi_ratio<-(meants[2:23]-meants[1:22])/meants[1:22]
  #get max (positive) for sos, and min (negative) for eos
  return(meants[c(which.max(ndvi_ratio),which.min(ndvi_ratio)+1)])
}

extract_phenology<-function(ts_thresh){
  if (is.na(sum(ts_thresh[24:25]))){
    return(c(NA,NA,NA))
  }
  sos_thresh<-ts_thresh[24]
  eos_thresh<-ts_thresh[25]
  #### fit the ts with a 6 order polynomial
  doy<-0.5:22.5
  if (sum(!is.na(ts_thresh[1:23]))<=6){
    return(c(NA,NA,NA))
  }
  ts_thresh[1:23]<-na.fill(ts_thresh[1:23],fill="extend")
  ts_fit<-lm(ts_thresh[1:23]~poly(doy,6,raw=T))
  b<-coef(ts_fit)
  predicted<-b[1] + poly(1:365/16,6,raw=T)%*%b[-1]
  ### get the date
  peak<-max(predicted[30:330],na.rm=T)
  
  if (is.infinite(peak)){
    return(c(NA,NA,NA))
  }
  pos<-which(predicted==peak)[1]
  increase<-c(rep(1,pos),rep(0,365-pos))
  decrease<-c(rep(0,pos-1),rep(1,365-pos+1))
  if (peak<sos_thresh){
    sos_acc<-NA
  }else{
    sos_int<-max(which(predicted<sos_thresh&increase))
    sos_acc<-sos_int+(sos_thresh-predicted[sos_int])/(predicted[sos_int+1]-predicted[sos_int])  
  }
  if (peak<eos_thresh){
    eos_acc<-NA
  }else{
    eos_int<-min(which(predicted<eos_thresh&decrease))-1
    eos_acc<-eos_int+(predicted[eos_int]-eos_thresh)/(predicted[eos_int]-predicted[eos_int+1])
  }
  return(c(sos_acc,pos,eos_acc))
}
#################
site_list<-read.csv("./retrieve_site/sites_2000.csv",stringsAsFactors = F)
phenology<-read.csv("./retrieve_site/analysis/site_phenology_0.3_N30_all_daily.csv",stringsAsFactors = F)
phenology_VI<-phenology[!is.na(phenology$site),1:8]
names(phenology_VI)[7:8]<-c("sos_VI_PM","peak_VI_PM")
phenology_VI[,7:8]<-NA
phenology_VI$eos_VI_PM<-NA

#### get all the EC_site files
com_f<-list.files('./retrieve_site/combined_all/',full.names = T)
sites_pheno<-unique(phenology_VI$site)
for (i in 1:dim(site_list)[1]){
  site_d<-read.csv(com_f[substr(basename(com_f),1,6)==sites_pheno[i]])
  site_years<-phenology_VI$year[phenology_VI$site==sites_pheno[i]]
  nobs<-dim(site_d)[1]/4
  ndvidatasite<-site_d$ndvi
  ndvidatasite[site_d$cloud==1|site_d$snow==1]<-NA
  #ndvi_dat_temp<-site_year_data$ndvi
  ndviqa_check<-ndvidatasite[(1:nobs)*4-3]
  vi_time<-site_d$TIMESTAMP[(1:nobs)*4-3]
  good_ndvi<-remove_outliers(ndviqa_check)
  ndvi_linear_fill<-na.fill(c(good_ndvi[(nobs-22):nobs],good_ndvi,good_ndvi[1:23]),fill = "extend")[24:(nobs+23)]
  #plot(ndvi_linear_fill,type='l')
  ndvi_threshold<-get_two_threshold(ndvi_linear_fill)
  for (y in 1:length(site_years)){
    ndvi_year_dat<-ndvi_linear_fill[floor(vi_time/10000)==site_years[y]]
    phenology_VI[phenology_VI$site==sites_pheno[i]&phenology_VI$year==site_years[y],7:9]<-
      extract_phenology(c(ndvi_year_dat, ndvi_threshold))
  }
}
##### the code is not good , need fix tomorrow.
write.csv(phenology_VI,"./retrieve_site/analysis/site_phenology_VI_PM.csv",row.names = F)


