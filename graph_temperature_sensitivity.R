source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")
setwd("/Users/yzhang/Project/SIF_phenology/")
load("./analysis/temperature_sensitivity_snow_treeheight.RData")
bin2<-colorRampPalette(ano_col_gmt)(2)


pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/temperature_sensitivity_snow_tree_height.pdf",
    height=3.5,width=8)
par(fig=c(0,0.5,0,1),mar=c(3.5,3.5,1,1),mgp=c(3,0.3,0),lend=1)
plot(NA,xlim=c(0,365),ylim=c(0,3.5),xlab="",ylab="",xaxs="i",yaxs='i',axes=F)
box()
axis(1,tck=-0.01)
axis(2,las=2,tck=-0.01)
mtext(side=1,line=2,"Snow free days")
mtext(side=2,line=2,expression(paste(gamma[T[air]]," (day/",degree,"C)",sep="")))
lines(ci_sos$snow_ind,-ci_sos$snow_mean,col=bin2[1],lwd=2.8)
polygon(c(ci_sos$snow_ind,rev(ci_sos$snow_ind)),-c(ci_sos$snow_ci_low,rev(ci_sos$snow_ci_high)),
        col=adjustcolor(bin2[1],alpha.f = 0.2),border = adjustcolor("grey",alpha.f = 0.5))
lines(ci_eos$snow_ind,ci_eos$snow_mean,col=bin2[2],lwd=2.8)
polygon(c(ci_eos$snow_ind,rev(ci_eos$snow_ind)),c(ci_eos$snow_ci_low,rev(ci_eos$snow_ci_high)),
        col=adjustcolor(bin2[2],alpha.f = 0.2),border = adjustcolor("grey",alpha.f = 0.5))

legend("topright",c(expression(gamma[T[air]]^SOP),expression(gamma[T[air]]^EOP)),
       lty=c(1,1),lwd=c(2.8,2.8),col=bin2,y.intersp = 1.3)
mtext(side=2,line=2.5,"a",cex=1,font=2,padj=-11,las=2)

par(fig=c(0.5,1,0,1),new=T)
plot(NA,xlim=c(0,40),ylim=c(0,3.5),xlab="",ylab="",xaxs="i",yaxs='i',axes=F)
box()
axis(1,tck=-0.01)
axis(2,las=2,tck=-0.01)
mtext(side=1,line=2,"Tree height (m)")
mtext(side=2,line=2,expression(paste(gamma[T[air]]," (day/",degree,"C)",sep="")))
lines(ci_sos$height_ind,-ci_sos$height_mean,col=bin2[1],lwd=2.8)
polygon(c(ci_sos$height_ind,rev(ci_sos$height_ind)),-c(ci_sos$height_ci_low,rev(ci_sos$height_ci_high)),
        col=adjustcolor(bin2[1],alpha.f = 0.2),border = adjustcolor("grey",alpha.f = 0.5))
lines(ci_eos$height_ind,ci_eos$height_mean,col=bin2[2],lwd=2.8)
polygon(c(ci_eos$height_ind,rev(ci_eos$height_ind)),c(ci_eos$height_ci_low,rev(ci_eos$height_ci_high)),
        col=adjustcolor(bin2[2],alpha.f = 0.2),border = adjustcolor("grey",alpha.f = 0.5))

mtext(side=2,line=2.5,"b",cex=1,font=2,padj=-11,las=2)
dev.off()


