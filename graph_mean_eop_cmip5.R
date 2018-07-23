#### plot the mean EOP for models
setwd("/Users/yzhang/Project/SIF_phenology/analysis/pheno_CMIP5/")
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")
##### convert the coreatlion to ae project.
nc2ae<-function(dat,ext=1){
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
  aedat<-projectRaster(latlongdat,crs = ae,res = c(50000,50000),method = "bilinear")
  return(aedat)
}

plot_nc_var<-function(raster_ae,v_range,title="",id="",color_ramp,position,lay_out){
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
  plot(border)
  image(setrange(raster_ae,v_range[1],v_range[2]),add=T,col=rev(color_ramp),
        axes=F,zlim=v_range)
  #plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
  text(0, 5500000,title,cex=1.3)
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.3,font=2,padj=-6.5,las=2)
}




model_ext<-c(2,3,2,2,2,2,2,2)
notion<-c("a","b","c","d","e","f","g","h")
model_files<-list.files("./",pattern="pheno_mean",full.names = T,recursive = T)
models<-c("CanESM2","CMCC-CESM","GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR",
          "MIROC-ESM","MPI-ESM2-LR","NorESM1-ME")
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Fig/Fig_S_CMIP5_mean_EOP2.pdf",height=6.5,width=12)
par(oma=c(0,0,2,3))
for (i in c(1:8)){
  ncin<-nc_open(model_files[grep(models[i],model_files)])
  dat<-ncvar_get(ncin,ncin$var[[3]])*365
  ####### use the function of creating AE from matrix
  ##   2 for change 0-360
  ##   3 for flip updown
  eop_ae<-nc2ae(dat,ext=model_ext[i])
  plot_nc_var(eop_ae,v_range = c(200,320),title = models[i],color_ramp = rev(lst_color_ramp),
              id = notion[i],position = c(i-floor((i-1)/4)*4,2-floor((i-1)/4)),lay_out = c(4,2))
  #mtext(,side=3,line=0.5,cex=1.2)
  nc_close(ncin)
}


par(fig=c(0.7,1,0.3,0.7),oma=c(0.4,0.4,0.4,2),
    mar=c(1,0.4,0.4,2.0),new=T)
plot(eop_ae, legend.only=TRUE, col=rev(lst_color_ramp),horizontal=F,zlim=c(200,320),
     legend.width=1.5, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=1))
#text(100000,0, labels = expression(italic(r)), xpd = NA, srt = -90,cex=1.1)
dev.off()



