##### graph the model results.
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


plot_nc_var_sig<-function(raster_ae,sig_ae,v_range,title="",id="",color_ramp,position,lay_out){
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),mgp=c(3,0.3,0),new=T)
  plot(border)
  image(setrange(raster_ae,v_range[1],v_range[2]),add=T,col=rev(color_ramp),
        axes=F,zlim=v_range)
  #plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
  points(sig_ae,pch=15,cex=0.2)
  text(0, 5500000,title,cex=1.3)
  text(10300000,0, labels = expression(italic(r)), xpd = NA, srt = -90,cex=1.1)   
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.3,font=2,padj=-9.5,las=2)
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),new=T)
  plot(raster_ae, legend.only=TRUE, col=rev(color_ramp),horizontal=F,zlim=v_range,
       legend.width=1.5, legend.shrink=0.75,
       axis.args=list(
         mgp=c(3,0.2,0),tck=0.3,
         cex.axis=1))
}

##### convert the coreatlion to ae project.
nc2ae<-function(dat,ext=1){
  size<-dim(dat)
  north_dat<-dat[,(size[2]/3*2+1):size[2]]
  if (ext !=1){
    north_dat<-rbind(north_dat[(size[1]/2+1):size[1],],north_dat[1:(size[1]/2),])
  }
  latlongdat<-raster(apply(north_dat,1,rev))
  extent(latlongdat)<-c(-180,180,30,90)
  projection(latlongdat)<-longlat
  aedat<-projectRaster(latlongdat,crs = ae,res = c(50000,50000),method = "bilinear")
  return(aedat)
}



#### convert the significance to ae project to point.
### pval 1 for significant, pval==0 for non-significant
convertraster2points<-function(pval, ext=1){
  size<-dim(pval)
  north_pval<-pval[,(size[2]/3*2+1):size[2]]
  if (ext!=1){
    north_pval<-rbind(north_pval[(size[1]/2+1):size[1],],north_pval[1:(size[1]/2),])
  }
  ras_pval<-raster(apply(north_pval,1,rev))
  extent(ras_pval)<-c(-180,180,30,90)
  projection(ras_pval)<-longlat
  #coarse_pval<-aggregate(ras_pval, fact=4, fun=mean,na.rm=T)
  sig_2deg<-ras_pval> 0.5
  sig_2deg[sig_2deg==0]<-NA
  sig_point<-coordinates(sig_2deg)[!is.na(values(sig_2deg)),]
  ae_sig_point<-project(sig_point,proj = ae)
  return(ae_sig_point)
}

plot_nc_var_sig<-function(raster_ae,sig_ae,v_range,title="",id="",color_ramp,position,lay_out){
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),mgp=c(3,0.3,0),new=T)
  plot(border)
  image(setrange(raster_ae,v_range[1],v_range[2]),add=T,col=rev(color_ramp),
        axes=F,zlim=v_range)
  #plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
  points(sig_ae,pch=15,cex=0.2)
  text(0, 5500000,title,cex=1.3)
  text(10300000,0, labels = expression(italic(r)), xpd = NA, srt = -90,cex=1.1)   
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.3,font=2,padj=-9.5,las=2)
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),new=T)
  plot(raster_ae, legend.only=TRUE, col=rev(color_ramp),horizontal=F,zlim=v_range,
       legend.width=1.5, legend.shrink=0.75,
       axis.args=list(
         mgp=c(3,0.2,0),tck=0.3,
         cex.axis=1))
}


models<-c("CanESM2","CMCC-CESM","GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR",
          "MIROC-ESM","MPI-ESM2-LR","NorESM1-ME")
varnames <- c("pr","tas","tasmax")
exp_var<-list(expression("Precip."),expression("T"[mean]),expression("T"[max]))
notion<-c("a","b","c","d","e","f","g","h","i")
model_files<-list.files("./analysis/correlation_CMIP5/",full.names = T)
for (i in 1:8){
  ncin<-nc_open(model_files[grep(models[i],model_files)])
  for (j in 1:3){
    dat<-ncvar_get(ncin,ncin$var[[j]])
    sig<-ncvar_get(ncin,ncin$var[[1+3]])
    max_cor_sig<-get_max_cor_index(dat,sig)
    ####### use the function of creating AE from matrix
    cor_ae<-nc2ae(max_cor_sig[[1]],ext=2)
    sig_ae<-convertraster2points(max_cor_sig[[2]]<0.05,ext=2)
    plot_nc_var_sig(cor_ae,sig_ae,v_range = c(-1,1),title = exp_var[[j]],color_ramp = rev(ano_col_gmt),
                    id = notion[i*3-8+j],position = c(i,j),lay_out = c(3,3))
  }

}







