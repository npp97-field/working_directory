# # ####################
# # ### second, get the ec-sites gpp
# # ### add NDVI and EVI into comparison
# # ###
# #
# setwd("/Users/yzhang/Project/SIF_phenology/")
# site_list<-read.csv("retrieve_site/Site_information_tier1.csv",stringsAsFactors = F)
# VI_f<-list.files("./retrieve_site/VIs/",full.names = T)
# ndvi_data<-as.data.frame(array(NA,dim=c(length(VI_f),166+1)))
# names(ndvi_data)<-c("date",site_list$SITE_ID)
# evi_data<-as.data.frame(array(NA,dim=c(length(VI_f),166+1)))
# names(evi_data)<-c("date",site_list$SITE_ID)
# snow_data<-as.data.frame(array(NA,dim=c(length(VI_f),166+1)))
# names(snow_data)<-c("date",site_list$SITE_ID)
# aerosol_data<-as.data.frame(array(NA,dim=c(length(VI_f),166+1)))
# names(aerosol_data)<-c("date",site_list$SITE_ID)
# cloud_data<-as.data.frame(array(NA,dim=c(length(VI_f),166+1)))
# names(cloud_data)<-c("date",site_list$SITE_ID)
# date_i<-substr(basename(VI_f),11,17)
# ndvi_data[,1]<-date_i
# evi_data[,1]<-date_i
# snow_data[,1]<-date_i
# aerosol_data[,1]<-date_i
# cloud_data[,1]<-date_i
# 
# for (i in 1:length(VI_f)){
#   print(VI_f[i])
#   vi_dat<-read.csv(VI_f[i],header = F)
#   ndvi_data[i,2:167]<-vi_dat[,1]
#   evi_data[i,2:167]<-vi_dat[,2]
#   snow_data[i,2:167]<-vi_dat[,3]
#   aerosol_data[i,2:167]<-vi_dat[,4]
#   cloud_data[i,2:167]<-vi_dat[,5]
# }
# write.csv(ndvi_data,'./retrieve_site/MOD13C1/ndvi_data.csv',row.names = F)
# write.csv(evi_data,'./retrieve_site/MOD13C1/evi_data.csv',row.names = F)
# write.csv(snow_data,'./retrieve_site/MOD13C1/snow_data.csv',row.names = F)
# write.csv(aerosol_data,'./retrieve_site/MOD13C1/aerosol_data.csv',row.names = F)
# write.csv(cloud_data,'./retrieve_site/MOD13C1/cloud_data.csv',row.names = F)


###-------------------------------------------------------------------
###  get DOY nc file
# ###
# setwd("/Users/yzhang/Project/SIF_phenology/")
# #### first, read the data of CSIF from nc files
# #site_list<-read.csv("./retrieve_site/sites_passed_heterogeneity_test.csv",stringsAsFactors = F)
# site_list<-read.csv("./retrieve_site/Site_information_tier1.csv",stringsAsFactors = F)
# ncfiles<-list.files("./retrieve_site/retrieved_CSIF/",full.names = T)
# site_no<-dim(site_list)[1]
# csif<-as.data.frame(array(NA,dim=c(length(ncfiles),site_no+1)))
# names(csif)<-c("DATE",site_list$SITE_ID)
# csif$DATE<-unlist(strsplit(basename(ncfiles),"[.]"))[(1:1288)*6-1]
# for (i in 1:1288){
#   ncf<-nc_open(ncfiles[i])
#   csif[i,2:(site_no+1)]<-ncvar_get(ncf,"CSIF")
#   nc_close(ncf)
# }
# write.csv(csif,"./retrieve_site/sites166_with_CSIF.csv",row.names = F)
# ###

###--------------------------------------------------------------------
###       
setwd("/Users/yzhang/Project/SIF_phenology/")
gpp_f<-list.files("./retrieve_site/EC_4day/",full.names = T)
sif<-read.csv("./retrieve_site/sites166_with_CSIF.csv",stringsAsFactors = F)
sites<-gsub('\\.','-',names(sif)[2:(dim(sif)[2])])
names(sif)<-c("DATE",sites)
site_list<-read.csv("./retrieve_site/sites_30N.csv",stringsAsFactors = F)
site_list<-site_list[site_list$LOCATION_LAT>=30,]
site_no<-dim(site_list)[1]

ndvi_data<-read.csv('./retrieve_site/MOD13C1/ndvi_data.csv',stringsAsFactors = F)
evi_data<-read.csv('./retrieve_site/MOD13C1/evi_data.csv',stringsAsFactors = F)
ndvi_data[,2:167]<-ndvi_data[,2:167]/10000
evi_data[,2:167]<-evi_data[,2:167]/10000
snow_data<-read.csv('./retrieve_site/MOD13C1/snow_data.csv',stringsAsFactors = F)
aerosol_data<-read.csv('./retrieve_site/MOD13C1/aerosol_data.csv',stringsAsFactors = F)
cloud_data<-read.csv('./retrieve_site/MOD13C1/cloud_data.csv',stringsAsFactors = F)

igbp<-unique(site_list$IGBP)
for (i in igbp){
  dir.create(paste("./retrieve_site/graph_ts/",i,sep=""))
}


for (i in 1:site_no){
  site_GPP_data<-read.csv(gpp_f[substr(basename(gpp_f),1,6)==site_list$SITE_ID[i]],stringsAsFactors = F)
  TIME_DOY<-format(strptime(site_GPP_data$TIMESTAMP,"%Y%m%d"),"%Y%j")
  gpp_dat<-cbind(TIME_DOY,site_GPP_data)
  site_SIF_dat<-subset(sif,select=c("DATE",site_list$SITE_ID[i]))
  ##### with snow 0   with aerosol 0  with cloud 1
  site_ndvi<-subset(ndvi_data,select=c("date",gsub('-','.',site_list$SITE_ID[i])))
  site_vi_dat<-cbind(site_ndvi,evi_data[,gsub('-','.',site_list$SITE_ID[i])],
                     snow_data[,gsub('-','.',site_list$SITE_ID[i])],
                     aerosol_data[,gsub('-','.',site_list$SITE_ID[i])],
                     cloud_data[,gsub('-','.',site_list$SITE_ID[i])])
  names(site_vi_dat)<-c("date",'ndvi','evi','snow','aerosol','cloud')
  site_vi_dat[!site_vi_dat$aerosol&site_vi_dat$cloud,c(2,3)]<-NA
  
  combined1<-merge(site_SIF_dat,gpp_dat,by.x="DATE",by.y="TIME_DOY")
  combined<-merge(combined1,site_vi_dat,by.x="DATE",by.y='date',all.x=T)
  ##### plot the graph here
  if (dim(combined)[1]==0|sum(!is.na(combined[,2]))==0|sum(combined$NEE_VUT_REF_QC>=0.8)==0)
    next
  combined[combined[,2]< -800,2]<-NA
  pdf(paste("./retrieve_site/graph_ts/",site_list$IGBP[i],'/',site_list$SITE_ID[i],'_',site_list$IGBP[i],
            "_GPP_SIF_graph.pdf",sep=""),width=12,height=8)
  par(fig=c(0,0.7,0.5,1),mar=c(3.5,3.5,1,2),mgp=c(3,0.5,0))
  time_ts<-strptime(combined$DATE,"%Y%j")
  bad_obs<-combined[combined$NEE_VUT_REF_QC<0.8,]
  good_obs<-combined[combined$NEE_VUT_REF_QC>=0.8,]
  maxy<-max(quantile(combined[,2],0.97,na.rm = T)*1.2,quantile(combined$GPP_NT_VUT_REF,0.97,na.rm = T)/20*1.2)
  plot(time_ts,combined[,2],ylim=c(-0.05,maxy),type="l",col="gold",xlab="",ylab="")
  legend("bottomleft",c("CSIF","GPP","badGPP"),col=c("gold",'blue',adjustcolor("red",alpha.f = 0.4)),
         lty=c(1,NA,NA),pch=c(NA,1,1),cex=0.9)
  par(new=T)
  plot(time_ts,combined$GPP_NT_VUT_REF,col=adjustcolor("blue",alpha.f = 0.4),
       cex=0.3,ylim=c(-1,maxy*20),axes=F,xlab="",ylab="")
  points(strptime(bad_obs$DATE,"%Y%j"),bad_obs$GPP_NT_VUT_REF,col=adjustcolor("red",alpha.f = 0.4),cex=0.3)
  axis(4)
  mtext(expression("CSIF"[all_daily]),side=2,line=2)
  
  par(fig=c(0.7,1,0.5,1),mar=c(3.5,3.5,1,1),new=T)
  plot(good_obs[,2],good_obs$GPP_NT_VUT_REF,xlim=c(-0.05,maxy),ylim=c(-1,maxy*20),xlab="",ylab="")
  mtext(expression("CSIF"[all_daily]),side=1,line=2)
  mtext(expression("GPP_NT"),side=2,line=2)
  reg<-lm(good_obs$GPP_NT_VUT_REF~good_obs[,2]-1)
  abline(reg,col="red",lty=2)
  text(maxy,0.1*maxy*20,pos=2,substitute(paste("R"^2~"=",a,sep=""),list(a=round(summary(reg)$r.square,2))))
  text(maxy,0,pos=2,substitute(paste("y=",a,"x",sep=""),list(a=round(reg$coefficients[1],2))))
  
  par(fig=c(0,0.7,0,0.5),mar=c(3.5,3.5,1,2),mgp=c(3,0.5,0),new=T)
  combined_vi<-combined[!is.na(combined$snow),]
  snow_ed<-combined_vi[combined_vi$snow==1,]
  snow_free<-combined_vi[combined_vi$snow==0,]
  maxy2<-quantile(combined_vi$ndvi,0.97,na.rm = T)*1.2
  time_ts2<-strptime(combined_vi$DATE,"%Y%j")
  plot(time_ts2,combined_vi$ndvi,col=adjustcolor("darkgreen",alpha.f = 0.4),pch=15,
       cex=0.3,ylim=c(0.1,maxy2),xlab="",ylab="",type='l')
  points(strptime(snow_ed$DATE,"%Y%j"),snow_ed$ndvi,col=adjustcolor("red",alpha.f = 0.4),pch=15,cex=0.3)
  lines(strptime(combined_vi$DATE,"%Y%j"),combined_vi$evi,col="darkviolet",pch=16,cex=0.3)
  points(strptime(snow_ed$DATE,"%Y%j"),snow_ed$evi,col=adjustcolor("red",alpha.f = 0.4),pch=16,cex=0.3)
  legend("bottomleft",c("NDVI","EVI","GPP"),col=c("darkgreen",'darkviolet','blue'),
         lty=c(1,1,NA),pch=c(NA,NA,1),cex=0.9)
  par(new=T)
  plot(time_ts,combined$GPP_NT_VUT_REF,col=adjustcolor("blue",alpha.f = 0.4),
       cex=0.3,ylim=c(-1,maxy*20),axes=F,xlab="",ylab="")
  points(strptime(bad_obs$DATE,"%Y%j"),bad_obs$GPP_NT_VUT_REF,col=adjustcolor("red",alpha.f = 0.4),cex=0.3)
  axis(4)
  
  par(fig=c(0.7,1,0,0.5),mar=c(3.5,3.5,1,1),new=T)
  plot(snow_free$ndvi,snow_free$GPP_NT_VUT_REF,xlim=c(-0.05,maxy2),ylim=c(-1,maxy*20),
       xlab="",ylab="",col="darkgreen")
  mtext(expression("NDVI/EVI"),side=1,line=2)
  mtext(expression("GPP_NT"),side=2,line=2)
  reg<-lm(snow_free$GPP_NT_VUT_REF~snow_free$ndvi-1)
  abline(reg,col="darkgreen",lty=2)
  points(snow_free$evi,snow_free$GPP_NT_VUT_REF,col="darkviolet")
  reg2<-lm(snow_free$GPP_NT_VUT_REF~snow_free$evi-1)
  abline(reg2,col="darkviolet",lty=2)
  text(maxy2,0.1*maxy*20,pos=2,substitute(paste("R"^2~"=",a,sep=""),list(a=round(summary(reg)$r.square,2))))
  text(maxy2,0,pos=2,substitute(paste("y=",a,"x",sep=""),list(a=round(reg$coefficients[1],2))))
  text(0,0.1*maxy*20,pos=4,substitute(paste("R"^2~"=",a,sep=""),list(a=round(summary(reg2)$r.square,2))))
  text(0,0,pos=4,substitute(paste("y=",a,"x",sep=""),list(a=round(reg2$coefficients[1],2))))
  
  dev.off()
  write.csv(combined,paste("./retrieve_site/combined/",site_list$SITE_ID[i],"_SIF_GPP_VI.csv",sep=""))
}

