####graph the sensitivity 
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")

# cropland<-shapefile("/Users/yzhang/Project/SIF_phenology/data/crop_north.shp")
# projection(cropland)<-longlat
# repcrop<-spTransform(cropland,ae)
setwd("/Users/yzhang/Project/SIF_phenology/")
ncin<-nc_open("./data/North_mcd12c1_landcover1_majority.nc")
lc_data<-ncvar_get(ncin,"lc1")
nc_close(ncin)
lc_data1<-lc_data
lc_data1[lc_data==12|lc_data==14]<-NA #cropland
lc_data1[lc_data<1]<-NA
lc_data1[!is.na(lc_data1)]<-1


discrete_col<-rep(drought_discrete[c(1,3,5,7)],each=20)

plot_nc_var<-function(raster_ae,v_range,title="",id="",color_ramp,position,lay_out){
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),mgp=c(3,0.3,0),new=T)
  plot(border)
  image(setrange(raster_ae,v_range[1],v_range[2]),add=T,col=rev(color_ramp),
        axes=F,zlim=v_range)
  #plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
  text(0, 5500000,title,cex=1.3)
  text(10300000,0, labels = expression(paste("day/",degree,"C",sep="")), xpd = NA, srt = -90,cex=1.1)   
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.3,font=2,padj=-9.5,las=2)
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),new=T)
  plot(raster_ae, legend.only=TRUE, col=rev(color_ramp),horizontal=F,zlim=v_range,
       legend.width=1.5, legend.shrink=0.75,
       axis.args=list(at=-2:2*4,
                      mgp=c(3,0.2,0),tck=0.3,
                      cex.axis=1))
}

notion<-c("a","b","c","d","e","f","g","h")
bin64<-colorRampPalette(ano_col_gmt)(64)
vars = c("_t_",'_p_')
exp_var<-list(expression("T"[day]~"~"~"EOP"),expression("T"[day]~"~"~"LOP"),
              expression("T"[mean]~"~"~"EOP"),expression("T"[mean]~"~"~"LOP")) ###,expression("T"[night]) ,expression("PAR")
phen="lgs"
dataset<-c('rs','era')
n_row = length(vars)
n_col=2
lay_out = c(n_col,n_row)
vrange<-list(c(-8,8),c(-0.3,0.3))
bin64[c(32,33)]<-"white"
pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/figures/",
          "temperature_sensitivity_LGS_EOS.pdf",sep=''),width=8.5,height=14/4*n_row)
par(oma=c(0,0,1,1))
for (d in 1:2){
  phen_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",dataset[d],"/",sep=""),
                         pattern=paste("^reg",sep=""),full.names = T)
  for (v in 1:length(vars)){
    ncin<-nc_open(phen_files[v])
    ncdat<-ncvar_get(ncin,"regress_coef")[2,,]
    nc_close(ncin)
    
    nc_var_ae<-nc2ae(ncdat)
    #nc_sig_ae<-convertraster2points(ncsig<0.05)  #*lc_data1
    
    plot_nc_var(nc_var_ae,vrange[[1]],exp_var[[v+d*2-2]],notion[v*2-2+d],rev(bin64),
                    c(d,n_row+1-v),c(2,n_row))
  }
}
dev.off()



plot_nc_var<-function(raster_ae,v_range,title="",id="",color_ramp,position,lay_out){
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),mgp=c(3,0.3,0),new=T)
  plot(border)
  image(setrange(raster_ae,v_range[1],v_range[2]),add=T,col=rev(color_ramp),
        axes=F,zlim=v_range)
  text(0, 5500000,title,cex=1.3)
  text(10300000,0, labels = expression(paste("day/mm",sep="")), xpd = NA, srt = -90,cex=1.1)   
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.3,font=2,padj=-9.5,las=2)
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),new=T)
  plot(raster_ae, legend.only=TRUE, col=rev(color_ramp),horizontal=F,zlim=v_range,
       legend.width=1.5, legend.shrink=0.75,
       axis.args=list(at=-1:8*0.1,
                      mgp=c(3,0.2,0),tck=0.3,
                      cex.axis=1))
}

############
exp_var<-list(expression("Prec."~"~"~"EOP"),expression("Prec."~"~"~"LOP"),
              expression("Prec."~"~"~"EOP"),expression("Prec."~"~"~"LOP")) ###,expression("T"[night]) ,expression("PAR")
phen="lgs"
dataset<-c('rs','era')
vrange<-list(c(-8,8),c(-0.1,0.4))
new_col<-colorRampPalette(ano_col_gmt)(32)
pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/figures/",
          "precipitation_sensitivity_LGS_EOS.pdf",sep=''),width=8.5,height=14/4*n_row)
par(oma=c(0,0,1,1))
for (d in 1:2){
  phen_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",dataset[d],"/",sep=""),
                         pattern=paste("^reg",sep=""),full.names = T)
  for (v in 1:length(vars)){
    ncin<-nc_open(phen_files[v])
    ncdat<-ncvar_get(ncin,"regress_coef")[3,,]
    nc_close(ncin)
    nc_var_ae<-nc2ae(ncdat)
    #nc_sig_ae<-convertraster2points(ncsig<0.05)  #*lc_data1
    plot_nc_var(nc_var_ae,vrange[[2]],exp_var[[v+d*2-2]],notion[v*2-2+d],rev(new_col[13:32]),
                c(d,n_row+1-v),c(2,n_row))
  }
}
dev.off()


