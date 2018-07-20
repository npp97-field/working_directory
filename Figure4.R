##### resample and create the threshold

library(ncdf4)
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")
setwd("/Users/yzhang/Project/SIF_phenology/")
ncin<-nc_open("./data/ai_cru_P_PET_grow.nc")
ai_data<-ncvar_get(ncin,"ai")
nc_close(ncin)
ai_dat_ras<-raster(apply(ai_data,1,rev))
extent(ai_dat_ras)<-c(-180,180,30,90)
projection(ai_dat_ras)<-longlat
ncin<-nc_open("./data/mean_annual_temp.nc")
T_data<-ncvar_get(ncin,"mean_annual_temp")
nc_close(ncin)
T_dat_ras<-raster(apply(T_data,1,rev))
extent(T_dat_ras)<-c(-180,180,30,90)
projection(T_dat_ras)<-longlat

ncin<-nc_open("./data/North_mcd12c1_landcover1_majority.nc")
lc_data<-ncvar_get(ncin,"lc1")
nc_close(ncin)
lc_data[lc_data==12|lc_data==14]<-0
lc_data[lc_data!=0]=1
lc_data_ras<-raster(apply(lc_data,1,rev))
extent(lc_data_ras)<-c(-180,180,30,90)
projection(lc_data_ras)<-longlat

############
nc2raster<-function(dat,ext=1){
  size<-dim(dat)
  format_dat<-dat
  if (ext ==2){
    format_dat<-rbind(format_dat[(size[1]/2+1):size[1],],format_dat[1:(size[1]/2),])
  }else if(ext ==3){
    format_dat<-t(apply(format_dat,1,rev))
    format_dat<-rbind(format_dat[(size[1]/2+1):size[1],],format_dat[1:(size[1]/2),])
  }
  north_dat<-format_dat[,(size[2]/3*2+1):size[2]]
  latlongdat<-raster(apply(north_dat,1,rev))
  extent(latlongdat)<-c(-180,180,30,90)
  projection(latlongdat)<-longlat
  return(latlongdat)
}

#### get the maximum correlation and its corresponding index
get_max_cor_index<-function(nc_var_dat,nc_var_sig){
  xydim<-dim(nc_var_dat)[2:3]
  dim(nc_var_dat)<-c(4,xydim[1]*xydim[2])
  dim(nc_var_sig)<-c(4,xydim[1]*xydim[2])
  abs_cor<-abs(nc_var_dat)
  maxindex<-apply(abs_cor,2,which.max)  #### the index where the maximum happens
  idx <- !(sapply(maxindex, length))
  maxindex[idx]<-NA
  maxind<-unlist(maxindex)
  
  sel<-as.matrix(cbind(maxind,1:(xydim[1]*xydim[2])))
  sign_cor<-sign(nc_var_dat[sel])
  max_cor<-apply(abs_cor,2,max)*sign_cor
  max_sig<-nc_var_sig[sel]
  dim(max_cor)<-xydim
  dim(max_sig)<-xydim
  return(list(max_cor,max_sig))
}


model_files<-list.files("./analysis/correlation_CMIP5/",full.names = T)
models <-c("CanESM2","CMCC-CESM","GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR",
           "MIROC-ESM","MPI-ESM2-LR","NorESM1-ME")
###read cmip5 data
model_ext<-c(2,3,2,2,2,2,2,2)
mod_dat<-list()
for (i in 1:8){
  ncin<-nc_open(model_files[grep(models[i],model_files)])
  size<-c(ncin$dim[["lon"]]$len,round(ncin$dim[["lat"]]$len/3))
  cor_sig_val<-array(NA,dim=c(size[1]*size[2],7))
  for (j in 1:2){
    dat<-ncvar_get(ncin,ncin$var[[j]])
    sig<-ncvar_get(ncin,ncin$var[[3+j]])
    max_cor_sig<-get_max_cor_index(dat,sig)
    ####### use the function of creating AE from matrix
    ##   2 for change 0-360
    ##   3 for flip updown
    
    cor_lonlat<-nc2raster(max_cor_sig[[1]],ext=model_ext[i])
    sig_lonlat<-nc2raster(max_cor_sig[[2]]<0.05,ext=model_ext[i])
    
    ### 1 for cor_prec 2 for pval_prec 3 for cor_tmean 4 for pval_tmean
    cor_sig_val[,j*2-1]<-getValues(cor_lonlat)
    cor_sig_val[,j*2]<-getValues(sig_lonlat)
  }
  
  ai_mod<-resample(ai_dat_ras,cor_lonlat,method='bilinear')
  T_mod<-resample(T_dat_ras,cor_lonlat,method='bilinear')
  lc_mod<-resample(lc_data_ras,cor_lonlat,method="ngb")
  cor_sig_val[,5]<-getValues(ai_mod)
  cor_sig_val[,6]<-getValues(T_mod)
  cor_sig_val[,7]<-getValues(lc_mod)
  
  mod_dat[[i+2]]<-cor_sig_val[complete.cases(cor_sig_val),]
  nc_close(ncin)
}

##### read the RS and ERA based results.
ncin<-nc_open("./analysis/correlation_clear_rs/max_correlation_eos_tday.nc")
tday<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
ncin<-nc_open("./analysis/correlation_clear_rs/max_significance_eos_tday.nc")
tday_sig<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
tday[tday_sig>0.05]<-NA

ncin<-nc_open("./analysis/correlation_clear_rs/max_correlation_eos_prec.nc")
prec<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
ncin<-nc_open("./analysis/correlation_clear_rs/max_significance_eos_prec.nc")
prec_sig<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
prec[prec_sig>0.05]<-NA
rs_dat<-cbind(as.vector(prec),as.vector(prec_sig<0.05),as.vector(tday),as.vector(tday_sig<.05),
              as.vector(ai_data),as.vector(T_data),as.vector(lc_data))
mod_dat[[1]]<-rs_dat[complete.cases(rs_dat),]
#### ERA
ncin<-nc_open("./analysis/correlation_clear_era/max_correlation_eos_temp.nc")
tday<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
ncin<-nc_open("./analysis/correlation_clear_era/max_significance_eos_temp.nc")
tday_sig<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
tday[tday_sig>0.05]<-NA

ncin<-nc_open("./analysis/correlation_clear_era/max_correlation_eos_prec.nc")
prec<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
ncin<-nc_open("./analysis/correlation_clear_era/max_significance_eos_prec.nc")
prec_sig<-ncvar_get(ncin,"max_cor")
nc_close(ncin)
prec[prec_sig>0.05]<-NA
era_dat<-cbind(as.vector(prec),as.vector(prec_sig<0.05),as.vector(tday),as.vector(tday_sig<.05),
              as.vector(ai_data),as.vector(T_data),as.vector(lc_data))
mod_dat[[2]]<-era_dat[complete.cases(era_dat),]


pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Fig/Figure_4.pdf",height=6,width=5)
### I will then analyze the correlation 
par(fig=c(0,1,0.5,1),mar=c(0,4,1.5,1),mgp=c(3,0.3,0),oma=c(6,0,0,0))
plot(NA,xlim=c(1.2,20.8),ylim=c(0,1.8),xaxs="i",yaxs="i",xlab="",ylab="",axes=F)
axis(1,tck=-0.02,at=(1:10)*2,lab=rep("",10))
axis(2,tck=-0.02,at=0:6*0.3,las=2)
box()
for (i in 1:10){
  single_mod_dat<-mod_dat[[i]]
  ### for temp limited
  prec_dat<-single_mod_dat[,c(1,2,5,7)]
  prec_lim<-prec_dat[prec_dat[,2]==1&prec_dat[,1]>0&prec_dat[,4]>0,3]
  boxplot(prec_lim,at=i*2-0.3,col=ano_col_gmt[3],border='gray45',add=T,xaxs="i",yaxs='i',outline=F,axes=F)
  ### for prec limited
  temp_dat<-single_mod_dat[,c(3,4,5,7)]
  temp_lim<-temp_dat[temp_dat[,2]==1&temp_dat[,1]>0&temp_dat[,4]>0,3]
  boxplot(temp_lim,at=i*2+0.3,col=ano_col_gmt[14],border='gray45',add=T,xaxs="i",yaxs='i',outline=F,axes=F)
}

abline(h=0.6,lty=2)
abline(v=5)
mtext(side=2,line=2,"Aridity Index")
mtext(side=2,line=3.4,'a',cex=1,font=2,padj=-9,las=2)

legend(6,1.8,c("Temperature limited","Precipitation limited"),fill=ano_col_gmt[c(14,3)],cex = 0.8,bty='n')

############
par(fig=c(0,1,0,0.605),mar=c(6,4,1.5,1),mgp=c(3,0.3,0),oma=c(0,0,0,0),new=T)
plot(NA,xlim=c(1.2,20.8),ylim=c(-25,25),xaxs="i",yaxs="i",xlab="",ylab="",axes=F)
axis(1,tck=-0.02,at=(1:10)*2,lab=rep("",10))
axis(2,tck=-0.02,las=2)
box()
#axis(1, at=c(1:10)*2,lab=c("RS","ERA-Interim",models),las=2)
for (i in 1:10){
  single_mod_dat<-mod_dat[[i]]
  ### for temp limited
  prec_dat<-single_mod_dat[,c(1,2,6,7)]
  prec_lim<-prec_dat[prec_dat[,2]==1&prec_dat[,1]>0&prec_dat[,4]>0,3]
  boxplot(prec_lim,at=i*2-0.3,col=ano_col_gmt[3],border='gray45',add=T,xaxs="i",yaxs='i',outline=F,axes=F)
  ### for prec limited
  temp_dat<-single_mod_dat[,c(3,4,6,7)]
  temp_lim<-temp_dat[temp_dat[,2]==1&temp_dat[,1]>0&temp_dat[,4]>0,3]
  boxplot(temp_lim,at=i*2+0.3,col=ano_col_gmt[14],border='gray45',add=T,xaxs="i",yaxs='i',outline=F,axes=F)
}
abline(h=2.5,lty=2)
abline(v=5)
mtext(side=2,line=2,expression(paste("Mean Annual Temperature (",degree,"C)",sep="")))
text(x=1:10*2, y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels= c('MODIS-MSWEP',"ERA-Interim",models), srt=45, adj=1, xpd=TRUE)
mtext(side=2,line=3.4,'b',cex=1,font=2,padj=-9,las=2)

dev.off()


