library(bootstrap)
setwd("/Users/yzhang/Project/SIF_phenology/")
ncin<-nc_open("./data/ai_cru_P_PET.nc")
aridity_mask<-ncvar_get(ncin,"ai")
aridity_mask[aridity_mask<0.6]<-NA
aridity_mask[!is.na(aridity_mask)]<-1
ncin<-nc_open("./data/North_barren_mask.nc")
barren_mask<-ncvar_get(ncin,"barren")

ncin<-nc_open("./data/tree_height.nc")
tree_height_raw<-ncvar_get(ncin,"Band1")
tree_height<-tree_height_raw*barren_mask*aridity_mask
nc_close(ncin)
# plot(nc2ae(tree_height*barren_mask*aridity_mask))
load("./analysis/snow_data.RData")
# plot(nc2ae(average_snow_start*barren_mask*aridity_mask))
# plot(nc2ae(average_snow_end*barren_mask*aridity_mask))
# plot(nc2ae((average_snow_start-average_snow_end)*barren_mask*aridity_mask))
snow_free<-(average_snow_start-average_snow_end)*barren_mask*aridity_mask

#### get the sensitivity for spring and autumn and snow_free and tree height.
phen_files<-list.files(paste("./analysis/correlation_clear_era/",sep=""),
                       pattern=paste("^reg",sep=""),full.names = T)[c(4,1)]
ncin<-nc_open(phen_files[1])
sensitivity_sos<-ncvar_get(ncin,"regress_coef")[2,,]
nc_close(ncin)
ncin<-nc_open(phen_files[2])
sensitivity_eos<-ncvar_get(ncin,"regress_coef")[2,,]
nc_close(ncin)

sensitivity<-as.data.frame(cbind(as.vector(sensitivity_sos),as.vector(sensitivity_eos),
                                 as.vector(tree_height),as.vector(snow_free)))
names(sensitivity)<-c("sos","eos","tree_height","snow_free_days")
dat<-sensitivity[complete.cases(sensitivity),]

#### get the spring temperature bined by snow period
# first bin the data together and create 
### for snow free days, the bin sizes are decided by the equal numbers  

theta<-function(x){mean(x)}
perc95<-function(x){quantile(x,0.95)}
perc05<-function(x){quantile(x,0.05)}

snow_intervals<-quantile(as.vector(snow_free),0:10/10,na.rm=T)
tree_height_intervals<-quantile(as.vector(tree_height),0:10/10,na.rm=T)

snow_ind<-quantile(as.vector(snow_free),1:10/10-0.05,na.rm=T)
tree_height_ind<-quantile(as.vector(tree_height),1:10/10-0.05,na.rm=T)
### for the snow cover

ci_sos<-as.data.frame(array(NA,dim=c(10,8)))
names(ci_sos)<-c("snow_ind",'snow_mean',"snow_ci_low","snow_ci_high",
             'height_ind','height_mean','height_ci_low',"height_ci_high")
ci_sos$snow_ind<-(snow_intervals[1:10]+snow_intervals[2:11])/2
ci_sos$height_ind<-(tree_height_intervals[1:10]+tree_height_intervals[2:11])/2
ci_eos<-ci_sos

for (i in 1:10){
  snow_bin<-dat[dat$snow_free_days<snow_intervals[i+1]&dat$snow_free_days>=snow_intervals[i],]
  ci_sos$snow_ci_low[i]<-bootstrap(snow_bin$sos,1000,theta, func=perc05)$func.thetastar
  ci_sos$snow_ci_high[i]<-bootstrap(snow_bin$sos,1000,theta, func=perc95)$func.thetastar
  ci_sos$snow_mean[i]<-mean(snow_bin$sos)

  ci_eos$snow_ci_low[i]<-bootstrap(snow_bin$eos,1000,theta, func=perc05)$func.thetastar
  ci_eos$snow_ci_high[i]<-bootstrap(snow_bin$eos,1000,theta, func=perc95)$func.thetastar
  ci_eos$snow_mean[i]<-mean(snow_bin$eos)
}

for (i in 2:10){
  tree_height_bin<-dat[dat$tree_height<tree_height_intervals[i+1]&dat$tree_height>=tree_height_intervals[i],]
  ci_sos$height_ci_low[i] <-  bootstrap(tree_height_bin$sos,1000,theta, func=perc05)$func.thetastar
  ci_sos$height_ci_high[i] <-  bootstrap(tree_height_bin$sos,1000,theta, func=perc95)$func.thetastar
  ci_sos$height_mean[i]<-mean(tree_height_bin$sos)
  ## for eos
  ci_eos$height_ci_low[i] <-  bootstrap(tree_height_bin$eos,1000,theta, func=perc05)$func.thetastar
  ci_eos$height_ci_high[i] <-  bootstrap(tree_height_bin$eos,1000,theta, func=perc95)$func.thetastar
  ci_eos$height_mean[i]<-mean(tree_height_bin$eos)
}

save(ci_sos,ci_eos,file='./analysis/temperature_sensitivity_snow_treeheight.RData')


