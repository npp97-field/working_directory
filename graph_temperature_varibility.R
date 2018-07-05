#### graph correlation between three climate variables
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")

lst_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_color.csv",header = F)/255))
ano_color_ramp<-rev(rgb(read.csv("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/R_graph/lst_ano.csv",header = F)/255))
ano_discrete_ramp<-c(ano_color_ramp[c(1,26,51,76,101)],"grey80","grey80",ano_color_ramp[c(140,165,190,215,240)])
ano_discrete_ramp_leg<-rep(ano_discrete_ramp,each=20)

setwd("/Users/yzhang/Project/SIF_phenology/")
# #ncin<-nc_open(paste("./climate_variability/clear_daily_SOS_30N_pre_start_par_stat.nc",sep=""))
# ncin<-nc_open(paste("./data/monthly_temp.nc",sep=""))
# sd_var<-ncvar_get(nc = ncin,varid = "sd_temp")
# nc_close(ncin)
# 
# pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/month_temp_var.pdf",sep=""),width=11/3*4.3,height=11)
# mon<-c("Jan",'Feb',"Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# #############################
# for (i in 1:12){
#   row = ceiling(i/4)
#   col = i-row*4+4
#   par(fig=c((col-1)/4,col/4,1-row/3,1-row/3+1/3),mar=c(0.4,0.4,0.4,0.4),oma=c(0,0,0,4),mgp=c(3,0.3,0),new=T)
#   plot(border)
#   # image(setrange(par_var,0,10),add=T,col=rev(lst_color_ramp),
#   #       axes=F,zlim=c(0,10))
#   image(setrange(nc2ae(sd_var[,,i]),0,4),add=T,col=rev(lst_color_ramp),
#         axes=F,zlim=c(0,4))
#   text(0, 5500000,mon[i],cex=1.3)
#   plotlatlong()
# 
# }
# par(fig=c(0.5,0.97,0,1),mar=c(0,0,0,0),oma=c(0,0,0,0),mgp=c(3,0.3,0))
# plot(par_var, legend.only=TRUE, col=rev(lst_color_ramp),horizontal=F,zlim=c(0,4),#zlim=c(0,10),#for std
#      legend.width=1.3, legend.shrink=0.75,
#      axis.args=list(#at=seq(-1, 1, 1/6),
#        #c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'),
#        mgp=c(3,0.2,0),tck=0.3,
#        cex.axis=1))
# 
# dev.off()
# 
# 



#ncin<-nc_open(paste("./climate_variability/clear_daily_SOS_30N_pre_start_par_stat.nc",sep=""))
ncin<-nc_open(paste("./analysis/all_daily_SOS_30N_fixed_stat.nc",sep=""))
sd_SOS<-ncvar_get(nc = ncin,varid = "SD")
nc_close(ncin)

ncin<-nc_open(paste("./analysis/all_daily_EOS_30N_fixed_stat.nc",sep=""))
sd_EOS<-ncvar_get(nc = ncin,varid = "SD")
nc_close(ncin)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/SOS_EOS_var.pdf",sep=""),width=11/3*2.1,height=11/3)
#############################


par(fig=c(0,0.5,0,1),mar=c(0.4,0.4,0.4,0.4),oma=c(0,0,0,4),mgp=c(3,0.3,0))
plot(border)
image(setrange(nc2ae(sd_SOS)*365,0,15),add=T,col=rev(lst_color_ramp),
      axes=F,zlim=c(0,15))
text(0, 5500000,"SD SOS",cex=1.3)
plotlatlong()


par(fig=c(0.5,1,0,1),mar=c(0.4,0.4,0.4,0.4),oma=c(0,0,0,4),mgp=c(3,0.3,0),new=T)
plot(border)
image(setrange(nc2ae(sd_EOS)*365,0,15),add=T,col=rev(lst_color_ramp),
      axes=F,zlim=c(0,15))
text(0, 5500000,"SD EOS",cex=1.3)
plotlatlong()


par(fig=c(0.5,0.97,0,1),mar=c(0,0,0,0),oma=c(0,0,0,0),mgp=c(3,0.3,0))
plot(par_var, legend.only=TRUE, col=rev(lst_color_ramp),horizontal=F,zlim=c(0,15),#zlim=c(0,10),#for std
     legend.width=1.3, legend.shrink=0.75,
     axis.args=list(#at=seq(-1, 1, 1/6),
       #c('-1','','-0.66','','-0.33','','0','','0.33','','0.66','','1'),
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=1))

dev.off()



ncin<-nc_open(paste("./analysis/all_daily_SOS_30N_fixed_stat.nc",sep=""))
sd_SOS<-ncvar_get(nc = ncin,varid = "SD")
nc_close(ncin)

ncin<-nc_open(paste("./analysis/all_daily_EOS_30N_fixed_stat.nc",sep=""))
sd_EOS<-ncvar_get(nc = ncin,varid = "SD")
nc_close(ncin)
