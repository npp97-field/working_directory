#####
# get phenology for ec flux tower sites
library(phenopix)
library(zoo)
setwd("/Users/yzhang/Project/SIF_phenology/")

extract_thresh<-function(in_var_ts,thresh){
  var_ts<-coredata(in_var_ts)
  n <- length(var_ts)
  avg <- mean(var_ts, na.rm=TRUE)

  peak <- max(var_ts, na.rm=TRUE)
  trough <- max(min(var_ts, na.rm=TRUE),0)
  ampl <- peak - trough
  
  increase<-c(rep(1,which(var_ts==peak)),rep(0,n-which(var_ts==peak)))
  decrease<-c(rep(0,which(var_ts==peak)-1),rep(1,n-which(var_ts==peak)+1))
  
  thresh_val<-trough+thresh*ampl
  sos_int<-max(which(var_ts<thresh_val&increase))
  sos_acc<-sos_int+(thresh_val-var_ts[sos_int])/(var_ts[sos_int+1]-var_ts[sos_int])  
  eos_int<-min(which(var_ts<thresh_val&decrease))-1
  eos_acc<-eos_int+(var_ts[eos_int]-thresh_val)/(var_ts[eos_int]-var_ts[eos_int+1])
  
  
  return(c((sos_acc-1)/n+time(in_var_ts)[1],(eos_acc-1)/n+time(in_var_ts)[1],thresh_val))
}


threshold_pct<-0.25
retrieved_pheno<-as.data.frame(array(NA,dim=c(1000,8)))
names(retrieved_pheno)<-c("site",'year','sos_gpp','eos_gpp','sos_ndvi','eos_ndvi','sos_csif','eos_csif')

site_list<-read.csv("./retrieve_site/sites_30N.csv",stringsAsFactors = F)
com_f<-list.files('./retrieve_site/combined/',full.names = T)
k<-0
for (i in 1:dim(site_list)[1]){
  site_data<-read.csv(com_f[substr(basename(com_f),1,6)==site_list$SITE_ID[i]])
  #site_data$GPP_NT_VUT_REF[site_data$NEE_VUT_REF_QC<0.8]<-NA
  years<-floor(site_data$DATE/1000)
  uyears<-unique(years)
  for (y in uyears){
    k=k+1
    site_year_data<-site_data[years==y,]
    retrieved_pheno$site[k]=site_list$SITE_ID[i]
    retrieved_pheno$year[k]=y

    ###
    gpp_dat<-site_year_data$GPP_NT_VUT_REF
    if (sum(site_year_data$NEE_VUT_REF_QC<0.8)>30){
      next
    }
    #}else{
    #  gpp_linear_fill<-na.fill(gpp_dat,fill = "extend")
    #}
    gpp_ts<-ts(gpp_dat,start=c(y,1),frequency=92)
    gpp_sp_fit<-SplineFit(gpp_ts,uncert=F,df.factor=0.10)
    
    gpp_derivatives.phenophases <- extract_thresh(gpp_sp_fit$fit$predicted,threshold_pct)
    retrieved_pheno[k,c(3,4)]<-gpp_derivatives.phenophases[c(1,2)]
    
    ###  NDVI
    ndvi_dat<-site_year_data$ndvi[1:23*4-3]
    if (sum(is.na(ndvi_dat))>11){
      next
    }else{
      ndvi_linear_fill<-na.fill(ndvi_dat,fill = "extend")
    }
    ndvi_ts<-ts(ndvi_linear_fill,start=c(y,1),frequency=23)
    ndvi_sp_fit<-SplineFit(ndvi_ts,uncert=F,df.factor=0.2)
    ndvi_derivatives.phenophases <- extract_thresh(ndvi_sp_fit$fit$predicted,threshold_pct)
    # ndvi_derivatives.phenophases <- PhenoDeriv(
    #   x=ndvi_sp_fit$fit$predicted, fit=ndvi_sp_fit$fit
    # )
    retrieved_pheno[k,c(5,6)]<-ndvi_derivatives.phenophases[c(1,2)]
    
    
    ### CSIF
    csif_dat<-site_year_data[,3]
    if (sum(is.na(csif_dat))>45){
      next
    }else{
      csif_linear_fill<-na.fill(csif_dat,fill = "extend")
    }
    #plot(y+0:91/92,csif_dat)
    csif_ts<-ts(csif_linear_fill,start=c(y,1),frequency=92)
    csif_sp_fit<-SplineFit(csif_ts,uncert=F,df.factor=0.10)
    csif_derivatives.phenophases <- extract_thresh(csif_sp_fit$fit$predicted,threshold_pct)
    retrieved_pheno[k,c(7,8)]<-csif_derivatives.phenophases[c(1,2)]
    
    pdf(paste('./retrieve_site/graph_sos_eos/',site_list$SITE_ID[i],'_',y,".pdf",sep=''),width=8,height=6)
    plot(gpp_ts,xlim=c(y,y+1),ylim=c(0, quantile(gpp_sp_fit$fit$predicted,0.97)*1.2),
         col="darkgreen",xlab="", ylab="",type="p")
    lines(gpp_sp_fit$fit$predicted,col="green")
    abline(v=gpp_derivatives.phenophases[c(1:2)], lty=c(1,2),col="green")
    abline(h=gpp_derivatives.phenophases[3])
    par(new=T)
    plot(csif_ts,xlim=c(y,y+1),ylim=c(0, quantile(csif_sp_fit$fit$predicted,0.97)*1.2),
         col="darkblue",xlab="", ylab="",type="p")
    lines(csif_sp_fit$fit$predicted,col="blue")
    abline(v=csif_derivatives.phenophases[c(1:2)], lty=c(1,2),col="blue")
    par(new=T)
    plot(ndvi_ts,xlim=c(y,y+1),ylim=c(0.15, quantile(ndvi_sp_fit$fit$predicted,0.97)*1.2),
         col="darkred",xlab="", ylab="",type="p")
    lines(ndvi_sp_fit$fit$predicted,col="red")
    abline(v=ndvi_derivatives.phenophases[c(1:2)], lty=c(1,2),col="red")
    
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


write.csv(retrieved_pheno,"./retrieve_site/analysis/site_phenology_0.25_N30.csv",row.names = F)

