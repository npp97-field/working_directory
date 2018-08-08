######################################################################
#load("/Users/yzhang/Project/SIF_phenology/analysis/ano_pheno_climate.RData")
#setwd("/Users/yzhang/habanero/PROJECT/SIF_phenology/")
setwd("/Users/yzhang/PROJECT/SIF_phenology/")
load("./analysis/ano_pheno_climate.RData")
# sos_stat_out[,,1]<-sos_stat_out[,,1]*365
# sos_stat_out[,,2]<-sos_stat_out[,,2]*365
# eos_stat_out[,,1]<-eos_stat_out[,,1]*365
# eos_stat_out[,,2]<-eos_stat_out[,,2]*365
# library(LSD)
# par(mfcol=c(1,2))
# heatscatter(as.vector(out_complete[,1:16])*365,as.vector(out_complete[,17:32]),xlim=c(-7,7))
# heatscatter(as.vector(out_complete[,33:48])*365,as.vector(out_complete[,49:64]),xlim=c(-7,7))
# library(plotKML)
# data(worldgrids_pal)
# lc_col<-worldgrids_pal$IGBP[2:15]
lc_col<-c("palegreen4","green4","green2","yellowgreen","lightseagreen",
          'darkmagenta',"coral3","tan4","orange2","darkolivegreen",'royalblue4',"darkorange3")
lc_used<-c(1,4,5,7,8,9,10,11)
lc_type<-c("ENF","EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO","URB","CNV")
s_slope<-c()
e_slope<-c()

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Fig/temp_sensitivity.pdf",width=8,height=4)
par(fig=c(0,0.5,0,1),mar=c(3.5,3.5,1,1),oma=c(0.4,0.4,0.4,0.4),mgp=c(3,0.3,0))
plot(NA,xlim=c(-3,3),ylim=c(-12,12),axes=F,xaxs='i',yaxs='i',xlab="",ylab="")
mtext(side=1,line=2,expression(paste("Preseason T anomaly (",degree,"C)",sep="")))
mtext(side=2, line=2,"SOP anomaly (days)")
box()
axis(1,tck=-0.01)
axis(2,las=2,tck=-0.01)
for (j in 1:length(lc_used)){
  points(sos_stat_out[,lc_used[j],3],sos_stat_out[,lc_used[j],1],col=lc_col[lc_used[j]],cex=0.4)
  reg<-lm(sos_stat_out[,lc_used[j],1]~sos_stat_out[,lc_used[j],3])
  s_slope[j]<-reg$coefficients[2]
  range<-c(min(sos_stat_out[,lc_used[j],3],na.rm=T),max(sos_stat_out[,lc_used[j],3],na.rm=T))
  lines(range,range*reg$coefficients[2]+reg$coefficients[1],col=lc_col[lc_used[j]])
}

legend("bottomleft",lc_type[lc_used],col=lc_col[lc_used],pch=rep(1,9),bty="n",cex=0.7)
reg_all<-lm(as.vector(sos_stat_out[,lc_used,1])~as.vector(sos_stat_out[,lc_used,3]))
corr<-cor.test(as.vector(sos_stat_out[,lc_used,1]),as.vector(sos_stat_out[,lc_used,3]))
text(3,10,pos=2,substitute(paste(italic(gamma),"=",slo%+-%std),
                           list(slo=round(mean(s_slope),2),std=round(sd(s_slope),2))))
text(3,7,pos=2,substitute(paste(italic(r),"=",slo,"   P<0.001"),
                          list(slo=round(corr$estimate,2))))

par(fig=c(0.5,1,0,1),new=T)
plot(NA,xlim=c(-3,3),ylim=c(-12,12),axes=F,xaxs='i',yaxs='i',xlab="",ylab="")
mtext(side=1,line=2,expression(paste("Preseason T anomaly (",degree,"C)",sep="")))
mtext(side=2, line=2,"EOP anomaly (days)")
box()
axis(1,tck=-0.01)
axis(2,las=2,tck=-0.01)
for (j in 1:length(lc_used)){
  points(eos_stat_out[,lc_used[j],3],eos_stat_out[,lc_used[j],1],col=lc_col[lc_used[j]],cex=0.4)
  reg<-lm(eos_stat_out[,lc_used[j],1]~eos_stat_out[,lc_used[j],3])
  e_slope[j]<-reg$coefficients[2]
  range<-c(min(eos_stat_out[,lc_used[j],3],na.rm=T),max(eos_stat_out[,lc_used[j],3],na.rm=T))
  lines(range,range*reg$coefficients[2]+reg$coefficients[1],col=lc_col[lc_used[j]])
}

reg_all<-lm(as.vector(eos_stat_out[,lc_used,1])~as.vector(eos_stat_out[,lc_used,3]))
corr<-cor.test(as.vector(eos_stat_out[,lc_used,1]),as.vector(eos_stat_out[,lc_used,3]))
text(3,10,pos=2,substitute(paste(italic(gamma),"=",slo%+-%std),
                           list(slo=round(mean(e_slope),2),std=round(sd(e_slope),2))))
text(3,7,pos=2,substitute(paste(italic(r),"=",slo,"   P<0.001"),
                          list(slo=round(corr$estimate,2))))

dev.off()

