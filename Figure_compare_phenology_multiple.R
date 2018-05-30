###### figure comparison between SOS and EOS with GPP based and other 

setwd("/Users/yzhang/Project/SIF_phenology/")
lc_type<-c("ENF","DBF","MF","WSA","OSH","GRA","WET")
lc_col<-c("palegreen4","yellowgreen","lightseagreen",
          "tan4","coral3","darkolivegreen",'royalblue4')
lc_pch<-1:length(lc_type)
site_list<-read.csv("./retrieve_site/sites_2000.csv",stringsAsFactors = F)

vi_spline<-read.csv("./retrieve_site/analysis/site_phenology_0.3_N30_all_daily.csv",stringsAsFactors = F)
vi_spline<-vi_spline[!is.na(vi_spline$sos_gpp),]
vi_PM<-read.csv("./retrieve_site/analysis/site_phenology_VI_PM.csv",stringsAsFactors = F )
vi_PM[,c(7:9)]<-vi_PM[,c(7:9)]/365
mcd12<-read.csv("./retrieve_site/analysis/site_phenology_mcd_N30.csv",stringsAsFactors = F)
mcd12[,c(5:6)]<-mcd12[,c(5:6)]/365
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/FigS_phenology_VI_MCD12.pdf",width=11,height=7.5)

#par(oma=c(0,0,0,0))
par(fig=c(0,0.33,0.5,1),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0))
####---------- a
plot(NA,xlim=c(0, 220),ylim=c(0,220),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-vi_spline[!is.na(match(vi_spline$site,site_name_biome)),]
  points((biome_pheno$sos_gpp%%1)*365+2.5,(biome_pheno$sos_ndvi%%1)*365+8.5,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(vi_spline$sos_gpp%%1,vi_spline$sos_ndvi%%1)
reg<-lm(vi_spline$sos_gpp%%1~vi_spline$sos_ndvi%%1-1)
text(220,220/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(220,220/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                    format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("SOS"[EC]) )
mtext(side=2,line=1.8,expression("SOS"[NDVI-spline]))
mtext(side=2,line=1.5,'a',cex=1.8,font=2,padj=-6,las=2)




par(fig=c(0,0.33,0,0.5),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0),new=T)
####---------- b
plot(NA,xlim=c(150, 365),ylim=c(150,365),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-vi_spline[!is.na(match(vi_spline$site,site_name_biome)),]
  points((biome_pheno$eos_gpp%%1)*365+2.5,(biome_pheno$eos_ndvi%%1)*365+8.5,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(vi_spline$eos_gpp%%1,vi_spline$eos_ndvi%%1)
reg<-lm(vi_spline$eos_gpp%%1~vi_spline$eos_ndvi%%1-1)
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
mtext(side=2,line=1.8,expression("EOS"[NDVI-PM]))
mtext(side=2,line=1.5,'b',cex=1.8,font=2,padj=-6,las=2)




par(fig=c(0.33,0.66,0.5,1),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0),new=T)
####---------- c
plot(NA,xlim=c(0, 220),ylim=c(0,220),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-vi_PM[!is.na(match(vi_PM$site,site_name_biome)),]
  points((biome_pheno$sos_gpp%%1)*365+2.5,(biome_pheno$sos_VI_PM%%1)*365+8.5,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(vi_PM$sos_gpp%%1,vi_PM$sos_VI_PM%%1)
reg<-lm(vi_PM$sos_gpp%%1~vi_PM$sos_VI_PM%%1-1)
text(220,220/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(220,220/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                    format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("SOS"[EC]) )
mtext(side=2,line=1.8,expression("SOS"[NDVI-PM]))
mtext(side=2,line=1.5,'a',cex=1.8,font=2,padj=-6,las=2)




par(fig=c(0.33,0.66,0,0.5),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0),new=T)
####---------- d
plot(NA,xlim=c(150, 365),ylim=c(150,365),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-vi_PM[!is.na(match(vi_PM$site,site_name_biome)),]
  points((biome_pheno$eos_gpp%%1)*365+2.5,(biome_pheno$eos_VI_PM%%1)*365,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(vi_PM$eos_gpp%%1,vi_PM$eos_VI_PM%%1)
reg<-lm(vi_PM$eos_gpp%%1~vi_PM$eos_VI_PM%%1-1)
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
mtext(side=2,line=1.8,expression("EOS"[MCD12Q2]))
mtext(side=2,line=1.5,'d',cex=1.8,font=2,padj=-6,las=2)





par(fig=c(0.66,1,0.5,1),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0),new=T)
####---------- e
plot(NA,xlim=c(0, 220),ylim=c(0,220),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-mcd12[!is.na(match(mcd12$site,site_name_biome)),]
  points((biome_pheno$sos_gpp%%1)*365+2.5,(biome_pheno$sos_mcd12%%1)*365+8.5,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(mcd12$sos_gpp%%1,mcd12$sos_mcd12%%1)
reg<-lm(mcd12$sos_gpp%%1~mcd12$sos_mcd12%%1-1)
text(220,220/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(220,220/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                    format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("SOS"[EC]) )
mtext(side=2,line=1.8,expression("SOS"[MCD12Q2]))
mtext(side=2,line=1.5,'e',cex=1.8,font=2,padj=-6,las=2)




par(fig=c(0.66,1,0,0.5),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0),new=T)
####---------- f
plot(NA,xlim=c(150, 365),ylim=c(150,365),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-mcd12[!is.na(match(mcd12$site,site_name_biome)),]
  points((biome_pheno$eos_gpp%%1)*365+2.5,(biome_pheno$eos_mcd12%%1)*365,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(mcd12$eos_gpp%%1,mcd12$eos_mcd12%%1)
reg<-lm(mcd12$eos_gpp%%1~mcd12$eos_mcd12%%1-1)
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
mtext(side=2,line=1.8,expression("EOS"[NDVI-spline]))
mtext(side=2,line=1.5,'f',cex=1.8,font=2,padj=-6,las=2)

dev.off()


