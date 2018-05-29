#####
# get phenology for ec flux tower sites
library(phenopix)
library(zoo)
setwd("/Users/yzhang/Project/SIF_phenology/")

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

get_threshold<-function(var_ts,thresh=0.3){
  peak <- max(var_ts, na.rm=TRUE)
  trough<-max(min(var_ts, na.rm=TRUE),0)
  if (trough<0.1*peak){
    trough<-0
  }
  ampl <- peak - trough
  if (peak<=0){
    return(-9999)
  }
  thresh_val<- -9999
  thresh_val<-trough+thresh*ampl
  return(thresh_val)
}


extract_thresh<-function(var_ts,thresh_val){
  #var_ts<-coredata(in_var_ts)
  n <- length(var_ts)
  avg <- mean(var_ts, na.rm=TRUE)
  
  peak <- max(var_ts[3:80], na.rm=TRUE)
  trough<-max(min(var_ts, na.rm=TRUE),0)
  if (trough<0.1*peak){
    trough<-0
  }
  ampl <- peak - trough
  if (peak<=0){
    return(c(-9999,-9999,-9999))
  }
  
  increase<-c(rep(1,which(var_ts==peak)),rep(0,n-which(var_ts==peak)))
  decrease<-c(rep(0,which(var_ts==peak)-1),rep(1,n-which(var_ts==peak)+1))
  
  sos_acc<- -9999
  eos_acc<- -9999
  
  sos_int<-max(which(var_ts<thresh_val&increase))
  sos_acc<-sos_int+(thresh_val-var_ts[sos_int])/(var_ts[sos_int+1]-var_ts[sos_int])  
  eos_int<-min(which(var_ts<thresh_val&decrease))-1
  eos_acc<-eos_int+(var_ts[eos_int]-thresh_val)/(var_ts[eos_int]-var_ts[eos_int+1])
  
  return(c((sos_acc-1)/n,(eos_acc-1)/n,thresh_val))
}

threshold_pct<-0.30
retrieved_pheno<-as.data.frame(array(NA,dim=c(1000,8)))
names(retrieved_pheno)<-c("site",'year','sos_gpp','eos_gpp','sos_ndvi','eos_ndvi','sos_csif','eos_csif')

site_list<-read.csv("./retrieve_site/sites_2000.csv",stringsAsFactors = F)
com_f<-list.files('./retrieve_site/combined_all/',full.names = T)
k<-0

for (i in 1:dim(site_list)[1]){
  site_data<-read.csv(com_f[substr(basename(com_f),1,6)==site_list$SITE_ID[i]])
  #site_data$GPP_NT_VUT_REF[site_data$NEE_VUT_REF_QC<0.8]<-NA
  years<-floor(site_data$DATE/1000)
  
  uyears<-unique(years[years>2000])
  
  site_new_data<-as.data.frame(array(NA,dim=c(92*length(uyears),8)))
  names(site_new_data)<-c("year","DATE","GPP_in","GPP_sp","SIF_in","SIF_sp","NDVI_in","NDVI_sp")
  
  ndvi_dat_temp<-site_data$ndvi
  ndvi_dat_temp[site_data$cloud==1|site_data$snow==1]<-NA
  site_ndvi_dat<-remove_outliers(ndvi_dat_temp[1:(23*length(uyears))*4-3])
  num_obse<-length(site_ndvi_dat)
  ndvi_linear_fill<-na.fill(c(site_ndvi_dat[(num_obse-22):num_obse],site_ndvi_dat,site_ndvi_dat[1:23]),
                            fill = "extend")[24:(num_obse+23)]
  for (y in 1:length(uyears)){
    site_year_data_temp<-site_data[years==uyears[y],]
    time_year<-uyears[y]*1000+1:92*4-3
    site_year_date_frame<-data.frame(time_year)
    site_year_data<-merge(site_year_date_frame,site_year_data_temp,by.x ='time_year',by.y='DATE',all.x=T)
    site_year_data$NEE_VUT_REF_QC[is.na(site_year_data$NEE_VUT_REF_QC)]<-0
    
    site_new_data$year[(y*92-91):(y*92)]<-uyears[y]
    site_new_data$DATE[(y*92-91):(y*92)]<-site_year_data_temp$DATE
    ### GPP
    gpp_dat<-na.fill(site_year_data$GPP_NT_VUT_REF,fill = 'extend')
    if (sum(site_year_data$NEE_VUT_REF_QC<0.8)>30){
      next
    }
    gpp_sp_fit<-smooth.spline(gpp_dat,df=9)

    ### CSIF
    csif_dat<-site_year_data[,3]
    if (sum(is.na(csif_dat))>45){
      next
    }else{
      csif_linear_fill<-na.fill(csif_dat,fill = "extend")
    }
    csif_sp_fit<-smooth.spline(csif_linear_fill,df=9)

    ###  NDVI
    ndvi_dat<-ndvi_linear_fill[years[1:(23*length(uyears))*4-3]==uyears[y]]
    ndvi_sp_fit<-smooth.spline(ndvi_dat,df=5)
    
    
    site_new_data$GPP_in[(y*92-91):(y*92)]<-gpp_sp_fit$yin
    site_new_data$GPP_sp[(y*92-91):(y*92)]<-gpp_sp_fit$y
    site_new_data$SIF_in[(y*92-91):(y*92)]<-csif_sp_fit$yin
    site_new_data$SIF_sp[(y*92-91):(y*92)]<-csif_sp_fit$y
    site_new_data$NDVI_in[(y*23-22):(y*23)*4-3]<-ndvi_sp_fit$yin
    site_new_data$NDVI_sp[(y*23-22):(y*23)*4-3]<-ndvi_sp_fit$y
  }
  gppthresh<-get_threshold(site_new_data$GPP_sp)
  sifthresh<-get_threshold(site_new_data$SIF_sp)
  ndvithresh<-get_threshold(site_new_data$NDVI_sp)
  for (y in 1:length(uyears)){
    k=k+1
    retrieved_pheno$site[k]=site_list$SITE_ID[i]
    retrieved_pheno$year[k]=uyears[y]
    
    year_data<-site_new_data[site_new_data$year==uyears[y],]
    if (sum(is.na(year_data$GPP_in))>30){
      next
    }
    retrieved_pheno[k,c(3,4)]<-extract_thresh(year_data$GPP_sp,gppthresh)[1:2]
    retrieved_pheno[k,c(5,6)]<-extract_thresh(year_data$NDVI_sp[(1:23)*4-3],ndvithresh)[1:2]
    retrieved_pheno[k,c(7,8)]<-extract_thresh(year_data$SIF_sp,sifthresh)[1:2]
  
  
    pdf(paste('./retrieve_site/graph_sos_eos_all/',site_list$SITE_ID[i],'_',uyears[y],".pdf",sep=''),width=8,height=6)
    plot(0:91/92*12,year_data$GPP_in,xlim=c(0,12),ylim=c(0, gppthresh*4),
         col="darkgreen",xlab="", ylab="",type="p")
    lines(0:91/92*12,year_data$GPP_sp,col="green")
    abline(v=retrieved_pheno[k,c(3,4)]*12, lty=c(1,2),col="green")
    abline(h=gppthresh,col="green")
    
    par(new=T)
    plot(0:91/92*12,year_data$SIF_in,xlim=c(0,12),ylim=c(0, sifthresh*4),
         col="darkblue",xlab="", ylab="",type="p")
    lines(0:91/92*12,year_data$SIF_sp,col="blue")
    abline(v=retrieved_pheno[k,c(7,8)]*12, lty=c(1,2),col="blue")
    abline(h=sifthresh,col="blue")
    
    par(new=T)
    plot(0:22/23*12,year_data$NDVI_in[(1:23)*4-3],xlim=c(0,12),
         ylim=c(0.15, quantile(year_data$NDVI_sp,0.97,na.rm=T)*1.2),
         col="darkred",xlab="", ylab="",type="p")
    lines(0:22/23*12,year_data$NDVI_sp[(1:23)*4-3],col="red")
    abline(v=retrieved_pheno[k,c(5,6)]*12, lty=c(1,2),col="red")
    dev.off()
  }
}


par(mfrow=c(1,4))

plot(retrieved_pheno$sos_gpp-floor(retrieved_pheno$sos_gpp),
     retrieved_pheno$sos_ndvi-floor(retrieved_pheno$sos_ndvi),xlim=c(0,0.6),ylim=c(0,0.6))
abline(0,1,col='red')
cor.test(retrieved_pheno$sos_gpp-floor(retrieved_pheno$sos_gpp),
     retrieved_pheno$sos_ndvi-floor(retrieved_pheno$sos_ndvi))

plot(retrieved_pheno$eos_gpp-floor(retrieved_pheno$eos_gpp),
     retrieved_pheno$eos_ndvi-floor(retrieved_pheno$eos_ndvi),xlim=c(0.4,1),ylim=c(0.4,1))
abline(0,1,col='red')
cor.test(retrieved_pheno$eos_gpp-floor(retrieved_pheno$eos_gpp),
         retrieved_pheno$eos_ndvi-floor(retrieved_pheno$eos_ndvi))

plot(retrieved_pheno$sos_gpp-floor(retrieved_pheno$sos_gpp),
     retrieved_pheno$sos_csif-floor(retrieved_pheno$sos_csif),xlim=c(0,0.6),ylim=c(0,0.6))
abline(0,1,col='red')
cor.test(retrieved_pheno$sos_gpp-floor(retrieved_pheno$sos_gpp),
         retrieved_pheno$sos_csif-floor(retrieved_pheno$sos_csif))

plot(retrieved_pheno$eos_gpp-floor(retrieved_pheno$eos_gpp),
     retrieved_pheno$eos_csif-floor(retrieved_pheno$eos_csif),xlim=c(0.4,1),ylim=c(0.4,1))
abline(0,1,col='red')
cor.test(retrieved_pheno$eos_gpp-floor(retrieved_pheno$eos_gpp),
         retrieved_pheno$eos_csif-floor(retrieved_pheno$eos_csif))


write.csv(retrieved_pheno,"./retrieve_site/analysis/site_phenology_0.3_N30_all_daily.csv",row.names = F)

