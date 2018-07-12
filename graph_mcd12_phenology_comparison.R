###### figure comparison between SOS and EOS with GPP based and other 

setwd("/Users/yzhang/Project/SIF_phenology/")
lc_type<-c("ENF","DBF","MF","WSA","OSH","GRA","WET")
lc_col<-c("palegreen4","yellowgreen","lightseagreen",
          "tan4","coral3","darkolivegreen",'royalblue4')
lc_pch<-1:length(lc_type)
site_list<-read.csv("./retrieve_site/sites_2000.csv",stringsAsFactors = F)

mcd12<-read.csv("./retrieve_site/analysis/site_phenology_mcd_N30.csv",stringsAsFactors = F)
mcd12[,c(5:6)]<-mcd12[,c(5:6)]/365

######get the interannual variability
site_phenology_used<-mcd12[!is.na(mcd12$sos_gpp),]
sites<-as.data.frame(table(site_phenology_used$site))
iav_sites<-sites$Var1[sites$Freq>=5]

ano_pheno<-as.data.frame(array(NA,dim=c(300,6)))
names(ano_pheno)<-c("site",'year','sos_gpp',"eos_gpp","sos_mcd12","eos_mcd12")
k=1
for(i in 1:length(iav_sites)){
  site_pheno<-site_phenology_used[site_phenology_used$site==iav_sites[i],]
  ano_pheno[k:(dim(site_pheno)[1]+k-1),1:2]<-site_pheno[,c(1,2)]
  ano_pheno$sos_gpp[k:(dim(site_pheno)[1]+k-1)]<-site_pheno$sos_gpp-mean(site_pheno$sos_gpp,na.rm=T)
  ano_pheno$eos_gpp[k:(dim(site_pheno)[1]+k-1)]<-site_pheno$eos_gpp-mean(site_pheno$eos_gpp,na.rm=T)
  ano_pheno$sos_mcd12[k:(dim(site_pheno)[1]+k-1)]<-site_pheno$sos_mcd12-mean(site_pheno$sos_mcd12,na.rm=T)
  ano_pheno$eos_mcd12[k:(dim(site_pheno)[1]+k-1)]<-site_pheno$eos_mcd12-mean(site_pheno$eos_mcd12,na.rm=T)
  k=k+dim(site_pheno)[1]
}
ano_pheno<-ano_pheno[complete.cases(ano_pheno),]



pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/FigS_phenology_VI_MCD12.pdf",width=7.5,height=7.5)

par(fig=c(0,0.5,0.5,1),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0))
####---------- a
plot(NA,xlim=c(0, 220),ylim=c(0,220),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-mcd12[!is.na(match(mcd12$site,site_name_biome)),]
  points((biome_pheno$sos_gpp%%1)*365+2.5,(biome_pheno$sos_mcd12%%1)*365,col=lc_col[lc],pch=lc_pch[lc])
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
mtext(side=1,line=1.8,expression("SOP"[EC]) )
mtext(side=2,line=1.8,expression("SOS"[MCD12Q2]))
mtext(side=2,line=1.5,'a',cex=1.8,font=2,padj=-6,las=2)


par(fig=c(0,0.5,0,0.5),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0),new=T)
####---------- b
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
mtext(side=1,line=1.8,expression("EOP"[EC]) )
mtext(side=2,line=1.8,expression("EOS"[MCD12Q2]))
mtext(side=2,line=1.5,'b',cex=1.8,font=2,padj=-6,las=2)


par(fig=c(0.5,1,0.5,1),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0),new=T)
####---------- c
plot(NA,xlim=c(-40, 40),ylim=c(-40,40),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-ano_pheno[!is.na(match(ano_pheno$site,site_name_biome)),]
  points((biome_pheno$sos_gpp)*365,(biome_pheno$sos_mcd12)*365,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(ano_pheno$sos_gpp,ano_pheno$sos_mcd12)
reg<-lm(ano_pheno$sos_gpp~ano_pheno$sos_mcd12-1)
text(40,-40+80/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                               list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                    b=corr$parameter)),
     pos=2)
text(40,-40+80/16*1.5,substitute(paste("RMSE=",a,sep=""),
                               list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                    format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("SOP"[EC]) )
mtext(side=2,line=1.8,expression("SOS"[MCD12Q2]))
mtext(side=2,line=1.5,'c',cex=1.8,font=2,padj=-6,las=2)


par(fig=c(0.5,1,0,0.5),mar=c(3.4,3.4,1.4,1.4),mgp=c(3,0.3,0),new=T)
####---------- d
plot(NA,xlim=c(-40, 40),ylim=c(-40,40),axes=F,xlab="",ylab="")
for (lc in 1:length(lc_type)){
  site_name_biome<-site_list$SITE_ID[site_list$IGBP==lc_type[lc]]
  biome_pheno<-ano_pheno[!is.na(match(ano_pheno$site,site_name_biome)),]
  points((biome_pheno$eos_gpp)*365,(biome_pheno$eos_mcd12)*365,col=lc_col[lc],pch=lc_pch[lc])
}
corr<-cor.test(ano_pheno$eos_gpp,ano_pheno$eos_mcd12)
reg<-lm(ano_pheno$eos_gpp~ano_pheno$eos_mcd12-1)
text(40,-40+80/16*2.5,substitute(paste(italic(r),"=",a,'   n=',b,sep=""),
                                   list(a=formatC(round(corr$estimate,2),format='f',digits = 2,flag="0"),
                                        b=corr$parameter)),
     pos=2)
text(40,-40+80/16*1.5,substitute(paste("RMSE=",a,sep=""),
                                   list(a=formatC(sqrt(mean((reg$residuals)^2))*365,2),
                                        format='f',digits = 2,flag="0")),
     pos=2)

axis(1,tck = -0.02)
axis(2,las=2,tck = -0.02)
box()
abline(0,1,lty=3,col="grey70")
mtext(side=1,line=1.8,expression("EOP"[EC]) )
mtext(side=2,line=1.8,expression("EOS"[MCD12Q2]))
mtext(side=2,line=1.5,'d',cex=1.8,font=2,padj=-6,las=2)



dev.off()


