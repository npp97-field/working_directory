####
library(maptools)
library(raster)
library(rgdal)



lc_type<-c("ENF","DBF","MF","WSA","OSH","GRA","WET")
lc_col<-c("palegreen4","yellowgreen","lightseagreen",
          "tan4","coral3","darkolivegreen",'royalblue4')

# site_list<-read.csv("./retrieve_site/sites_used.csv",stringsAsFactors = F)
# site_list<-site_list[site_list$LOCATION_LAT>=35,]
# site_no<-dim(site_list)[1]
# write.csv(site_list,"./retrieve_site/sites_35N.csv",row.names = F)
setwd("/Users/yzhang/Project/SIF_phenology/")
site_list<-read.csv("./retrieve_site/sites_30N.csv",stringsAsFactors = F)
site_phenology<-read.csv("./retrieve_site/analysis/site_phenology_0.25_N30.csv",stringsAsFactors = F)
site_phenology_used<-site_phenology[!is.na(site_phenology$sos_gpp),]

comb_f<-list.files("./retrieve_site/combined/",full.names = T)
for (i in 1:length(site_list$SITE_ID)){
  temp<-read.csv(comb_f[substr(basename(comb_f),1,6)==site_list$SITE_ID[i]],header = T)
  names(temp)[3]<-"CSIF"
  temp$igbp<-site_list$IGBP[i]
  years_used<-unique(site_phenology_used$year[site_phenology_used$site==site_list$SITE_ID[i]])
  data_used<-temp[!is.na(match(floor(temp$DATE/1000),years_used)),]
  #temp$nor_csif<-temp$CSIF/quantile(temp$CSIF,0.95,na.rm=T)
  #temp$nor_GPP<-temp$GPP_NT_VUT_REF/quantile(temp$GPP_NT_VUT_REF,0.95,na.rm=T)
  #temp$nor_ndvi<-temp$ndvi/quantile(temp$ndvi,0.95,na.rm=T)
  if (i==1){
    combined_data<-data_used
  }
  combined_data<-rbind(combined_data,data_used)
}

write.csv(combined_data,"./retrieve_site/analysis/combined_all_sites_data.csv",row.names = F)

comb_data<-read.csv("./retrieve_site/analysis/combined_all_sites_data.csv",stringsAsFactors = F)
good_obs<-comb_data[comb_data$NEE_VUT_REF_QC>0.8,]

lc_pch<-1:length(lc_type)
site_col<-lc_col[match(site_list$IGBP,lc_type)]
site_pch<-lc_pch[match(site_list$IGBP,lc_type)]

site_loc<-SpatialPoints(site_list[c(8,7)], proj4string=longlat)
proj_site<-spTransform(site_loc,ae)

longlat <-  CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
ae<-"+proj=aeqd +lat_0=90 +lon_0=-0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


plot_lat<-function(lat){
  norths<-rep(lat,181)
  easts<- -90:90*2
  
  lati <- cbind(easts, norths)	
  lati <- rbind(lati, lati[1,])	# creates a matrix with two colums
  proj_lat <- SpatialPolygons(list(Polygons(list(Polygon(lati)), 1)))	# create a polgon from this matrix
  proj4string(proj_lat) <- longlat
  replat<-spTransform(proj_lat,ae)
  return(replat)
}

plot_long<-function(lon){
  # if (lon==-180){
  #   norths<- 10:16*5
  #   easts<- rep(lon,7)
  # }else{
  # 
  # }
  norths<- 6:16*5
  easts<- rep(lon,11)
  
  long <- cbind(easts, norths)	
  long <- rbind(long, long[1,])	# creates a matrix with two colums
  proj_lon <- SpatialPolygons(list(Polygons(list(Polygon(long)), 1)))	# create a polgon from this matrix
  proj4string(proj_lon) <- longlat
  replon<-spTransform(proj_lon,ae)
  return(replon)
}

plotlatlong<-function(){
  for (i in 3:8){
    prjlat<-plot_lat(i*10)
    lines(prjlat,lty=8,col='grey50',lwd=0.5)
  }
  
  for (i in -6:5){
    prolon<-plot_long(i*30)
    lines(prolon,lty=8,col='grey50',lwd=0.5)
  }
  
  lines(repcoa,col='black',lwd=0.2)
}



border<-plot_lat(30)
lat60<-plot_lat(60)

### coastlines
coastline<-shapefile("/Users/yzhang/Data/GIS_data/global/ne_110m_coastline.shp")
cropcoast<-crop(coastline,extent(-180,180,30,90))
repcoa<-spTransform(cropcoast,ae)

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/site_loc_0.25_30N.pdf",width=11,height=11)

#### -----------------------------------------
##  a
par(fig=c(0,0.33,0.66,1),mar=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(border)
plotlatlong()
lines(repcoa,col='black',lwd=0.2)

points(proj_site,col=site_col,pch=site_pch)
mtext(side=2,line=-1.5,"a",cex=1.8,font=2,padj=-7,las=2)
# text(0,-6500000,expression(0*degree),cex=0.5)
# text(0,6500000,expression(180*degree),cex=0.5)
# text(-6700000,0,expression(90*degree*W),cex=0.5)
# text(67000000,0,expression(90*degree*E),cex=0.5)
#box()
legend("bottomright",lc_type,col=lc_col,pch=lc_pch,cex=0.7,bg = 'white')

#### -----------------------------------------
##  b

par(fig=c(0,0.33,0.33,0.66),mar=c(3.4,3.4,1.4,1.4),new=T)
plot(NA,xlim=c(0, 1.1),ylim=c(0,20),axes=F,xlab="",ylab="")
#plot(NA,xlim=c(0, 1.3),ylim=c(0,1.3),axes=F,xlab="",ylab="")
for (i in 1:7){
  biome_dat<-good_obs[good_obs$igbp==lc_type[i],]
  points(biome_dat$CSIF,biome_dat$GPP_NT_VUT_REF,col=lc_col[i],pch=lc_pch[i],cex=0.3)
  reg<-cor.test(biome_dat$CSIF,biome_dat$GPP_NT_VUT_REF)
  text(1.15,10-1.25*i,substitute(paste(italic(r),"=",a,sep=""),
                                list(a=formatC(round(reg$estimate,2),format='f',digits = 2,flag="0"))),
       pos=2,col=lc_col[i])
}
axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02,las=2)
mtext(side=1,line=1.8,expression(paste("CSIF (mW m"^-2,"nm"^-1,"sr"^-1,")",sep="")))
mtext(side=2,line=1.8,expression(paste("GPP (gC m"^-2,"day"^-1,")",sep="")))
mtext(side=2,line=1.5,'b',cex=1.8,font=2,padj=-6,las=2)
box()

#### -----------------------------------------
##  c

par(fig=c(0,0.33,0,0.33),mar=c(3.4,3.4,1.4,1.4),new=T)
plot(NA,xlim=c(0, 1.1),ylim=c(0,20),axes=F,xlab="",ylab="")
for (i in 1:7){
  biome_dat<-good_obs[good_obs$igbp==lc_type[i],]
  points(biome_dat$ndvi,biome_dat$GPP_NT_VUT_REF,col=lc_col[i],pch=lc_pch[i],cex=0.3)
  reg<-cor.test(biome_dat$ndvi,biome_dat$GPP_NT_VUT_REF)
  text(0,19-1.25*i,substitute(paste(italic(r),"=",a,sep=""),
                                list(a=formatC(round(reg$estimate,2),format='f',digits = 2,flag="0"))),
       pos=4,col=lc_col[i])
}
axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02,las=2)
mtext(side=1,line=1.8,expression(paste("NDVI",sep="")))
mtext(side=2,line=1.8,expression(paste("GPP (gC m"^-2,"day"^-1,")",sep="")))
mtext(side=2,line=1.5,'c',cex=1.8,font=2,padj=-6,las=2)
box()

#### -----------------------------------------
##  d
par(fig=c(0.33,0.66,0.66,1),mar=c(3.4,3.4,1.4,1.4),new=T)
plot(NA,xlim=c(0, 170),ylim=c(0,170),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-site_phenology_used[!is.na(match(site_phenology_used$site,site_name_biome)),]
  points((biome_pheno$sos_gpp%%1)*365+2.5,(biome_pheno$sos_csif%%1)*365+2.5,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(site_phenology_used$sos_gpp%%1,site_phenology_used$sos_csif%%1)
reg<-lm(site_phenology_used$sos_gpp%%1~site_phenology_used$sos_csif%%1-1)
text(170,170/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                             list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                  b=corr$parameter)),
     pos=2)
text(170,170/16*1.5,substitute(paste("RMSE=",a,sep=""),
                             list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                  format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("SOS"[EC]) )
mtext(side=2,line=1.8,expression("SOS"[CSIF]))
mtext(side=2,line=1.5,'d',cex=1.8,font=2,padj=-6,las=2)
#### -----------------------------------------
##  e
par(fig=c(0.33,0.66,0.33,0.66),mar=c(3.4,3.4,1.4,1.4),new=T)
plot(NA,xlim=c(150, 365),ylim=c(150,365),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-site_phenology_used[!is.na(match(site_phenology_used$site,site_name_biome)),]
  points((biome_pheno$eos_gpp%%1)*365+2.5,(biome_pheno$eos_csif%%1)*365+2.5,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(site_phenology_used$eos_gpp%%1,site_phenology_used$eos_csif%%1)
reg<-lm(site_phenology_used$eos_gpp%%1~site_phenology_used$eos_csif%%1-1)
text(365,150+215/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(365,150+215/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                    format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("EOS"[EC]) )
mtext(side=2,line=1.8,expression("EOS"[CSIF]))
mtext(side=2,line=1.5,'e',cex=1.8,font=2,padj=-6,las=2)
#### -----------------------------------------
##  f
par(fig=c(0.33,0.66,0,0.33),mar=c(3.4,3.4,1.4,1.4),new=T)
plot(NA,xlim=c(80, 320),ylim=c(80,320),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-site_phenology_used[!is.na(match(site_phenology_used$site,site_name_biome)),]
  points((biome_pheno$eos_gpp-biome_pheno$sos_gpp)*365+2.5,
         (biome_pheno$eos_csif-biome_pheno$sos_csif)*365+2.5,col=lc_col[lc],pch=lc_pch[lc])
}
gpplos<-site_phenology_used$eos_gpp-site_phenology_used$sos_gpp
csiflos<-site_phenology_used$eos_csif-site_phenology_used$sos_csif
corr<-cor.test(gpplos,csiflos)
reg<-lm(gpplos~csiflos-1)
text(320,100+220/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(320,100+220/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                    format='f',digits = 2,flag="0")),
     pos=2)
axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("LOS"[EC]) )
mtext(side=2,line=1.8,expression("LOS"[CSIF]))
mtext(side=2,line=1.5,'f',cex=1.8,font=2,padj=-6,las=2)

#### -----------------------------------------
##  g
par(fig=c(0.66,1,0.66,1),mar=c(3.4,3.4,1.4,1.4),new=T)
plot(NA,xlim=c(0, 170),ylim=c(0,170),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-site_phenology_used[!is.na(match(site_phenology_used$site,site_name_biome)),]
  points((biome_pheno$sos_gpp%%1)*365+2.5,(biome_pheno$sos_ndvi%%1)*365+2.5,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(site_phenology_used$sos_gpp%%1,site_phenology_used$sos_ndvi%%1)
reg<-lm(site_phenology_used$sos_gpp%%1~site_phenology_used$sos_ndvi%%1-1)
text(170,170/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(170,170/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                    format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("SOS"[EC]) )
mtext(side=2,line=1.8,expression("SOS"[NDVI]))
mtext(side=2,line=1.5,'g',cex=1.8,font=2,padj=-6,las=2)
#### -----------------------------------------
##  h
par(fig=c(0.66,1,0.33,0.66),mar=c(3.4,3.4,1.4,1.4),new=T)
plot(NA,xlim=c(150, 365),ylim=c(150,365),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-site_phenology_used[!is.na(match(site_phenology_used$site,site_name_biome)),]
  points((biome_pheno$eos_gpp%%1)*365+2.5,(biome_pheno$eos_ndvi%%1)*365+2.5,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(site_phenology_used$eos_gpp%%1,site_phenology_used$eos_ndvi%%1)
reg<-lm(site_phenology_used$eos_gpp%%1~site_phenology_used$eos_ndvi%%1-1)
text(365,150+215/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(365,150+215/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),format='f',digits = 2,flag="0")),
     pos=2)


axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("EOS"[EC]) )
mtext(side=2,line=1.8,expression("EOS"[NDVI]))
mtext(side=2,line=1.5,'h',cex=1.8,font=2,padj=-6,las=2)
#### -----------------------------------------
##  i
par(fig=c(0.66,1,0,0.33),mar=c(3.4,3.4,1.4,1.4),new=T)
plot(NA,xlim=c(80, 320),ylim=c(80,320),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-site_phenology_used[!is.na(match(site_phenology_used$site,site_name_biome)),]
  points((biome_pheno$eos_gpp-biome_pheno$sos_gpp)*365+2.5,
         (biome_pheno$eos_ndvi-biome_pheno$sos_ndvi)*365+2.5,col=lc_col[lc],pch=lc_pch[lc])
}

gpplos<-site_phenology_used$eos_gpp-site_phenology_used$sos_gpp
csiflos<-site_phenology_used$eos_ndvi-site_phenology_used$sos_ndvi
corr<-cor.test(gpplos,csiflos)
reg<-lm(gpplos~csiflos-1)
text(320,100+220/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(320,100+220/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("LOS"[EC]) )
mtext(side=2,line=1.8,expression("LOS"[NDVI]))
mtext(side=2,line=1.5,'i',cex=1.8,font=2,padj=-6,las=2)
dev.off()
     
