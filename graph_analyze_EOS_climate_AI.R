#### graph and analyze the EOS climate relationship with aridity index.
## 
library("LSD")
## use the heatscatter function
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")
setwd("/Users/yzhang/Project/SIF_phenology/")
#ncin<-nc_open("./data/ai_cru_P_PET.nc")
ncin<-nc_open("./data/ai_cru_P_PET_grow.nc")
ai_data<-ncvar_get(ncin,"ai")
nc_close(ncin)

ncin<-nc_open("./data/North_mcd12c1_landcover1_majority.nc")
lc_data<-ncvar_get(ncin,"lc1")
nc_close(ncin)

ncin<-nc_open("./analysis/correlation_clear_rs/max_correlation_eos_tday.nc")
tday<-ncvar_get(ncin,"max_cor")
nc_close(ncin)

ncin<-nc_open("./analysis/correlation_clear_rs/max_correlation_eos_prec.nc")
prec<-ncvar_get(ncin,"max_cor")
nc_close(ncin)

ncin<-nc_open("./analysis/correlation_clear_rs/max_correlation_eos_par.nc")
prec<-ncvar_get(ncin,"max_cor")
nc_close(ncin)

##### get the correlation for all points, ignore land cover tyeps
lc_type<-c("ENF","EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO","URB","CNV")
id<-c(1,4,5,7,8,9,10,12)
bind_dat<-cbind(as.vector(ai_data),as.vector(tday),as.vector(lc_data))
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/ai_climate_by_veg_grow.pdf",width=8,height=8)
par(mfrow=c(3,3),mar=c(3,3,3,3),oma=c(3,3,1,1))
for (i in 1:8){
  lc_sub<-bind_dat[bind_dat[,3]==id[i],]
  heatscatter(lc_sub[,1],lc_sub[,2],xlim=c(0,1.5),ylim=c(-1,1),main="")
  mtext(side=3,line=0,lc_type[id[i]],cex=1.3,font=2)
  abline(v=0.5)
}
mtext(side=2,line=0,"Correlation",outer=T,cex=1.6)
mtext(side=1,line=0,"Aridity Index",outer=T,cex=1.6)
dev.off()


