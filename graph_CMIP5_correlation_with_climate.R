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
      mar=c(1,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
  plot(border)
  image(setrange(raster_ae,v_range[1],v_range[2]),add=T,col=rev(color_ramp),
        axes=F,zlim=v_range)
  #plot(repcrop[c(1:288,290:291),],density=50,border=NA,add=T)
  points(sig_ae,pch=15,cex=0.2)
  text(0, 5500000,title,cex=1.3)
  
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.3,font=2,padj=-6.5,las=2)
}

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



#### convert the significance to ae project to point.
### pval 1 for significant, pval==0 for non-significant
convertraster2points<-function(pval, ext=1){
  size<-dim(pval)
  format_pval<-pval
  if (ext==2){
    format_pval<-rbind(format_pval[(size[1]/2+1):size[1],],format_pval[1:(size[1]/2),])
  }else if(ext==3){
    format_pval<-t(apply(format_pval,1,rev))
    format_pval<-rbind(format_pval[(size[1]/2+1):size[1],],format_pval[1:(size[1]/2),])
  }
  north_pval<-format_pval[,(size[2]/3*2+1):size[2]]

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


models<-c("CanESM2","CMCC-CESM","GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR",
          "MIROC-ESM","MPI-ESM2-LR","NorESM1-ME")
varnames <- c("pr","tas","tasmax")
exp_var<-list(expression("T"[max]),expression("T"[mean]),expression("Precip."))
notion<-c("a1","a2","a3","b1","b2","b3",
          "c1","c2","c3","d1","d2","d3",
          "e1","e2","e3","f1","f2","f3",
          "g1","g2","g3","h1","h2","h3")
model_files<-list.files("./analysis/correlation_CMIP5/",full.names = T)

model_ext<-c(2,3,2,2,2,2,2,2)
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Fig/Fig_S_CMIP5_correlation.pdf",height=8.5,width=24)
par(oma=c(0,0,2,3))
for (i in c(1:8)){
  ncin<-nc_open(model_files[grep(models[i],model_files)])
  for (j in 1:3){
    if (i==8&j==1){
      position<-c(i,4-j)
      lay_out<-c(8,3)
      par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
                (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
          mar=c(1,0.4,0.4,0.4),mgp=c(3,0.3,0),new=T)
      plot.new()
      mtext(models[i],side=3,line=0.5,cex=1.2)
      next
    }
    dat<-ncvar_get(ncin,ncin$var[[4-j]])
    sig<-ncvar_get(ncin,ncin$var[[7-j]])
    max_cor_sig<-get_max_cor_index(dat,sig)
    ####### use the function of creating AE from matrix
    ##   2 for change 0-360
    ##   3 for flip updown
    cor_ae<-nc2ae(max_cor_sig[[1]],ext=model_ext[i])
    sig_ae<-convertraster2points(max_cor_sig[[2]]<0.05,ext=model_ext[i])
    plot_nc_var_sig(cor_ae,sig_ae,v_range = c(-1,1),title = exp_var[[j]],color_ramp = rev(ano_col_gmt),
                    id = notion[i*3-3+j],position = c(i,4-j),lay_out = c(8,3))
    if (j==1){
      mtext(models[i],side=3,line=0.5,cex=1.2)
    }
  }
  nc_close(ncin)
}


par(fig=c(0.7,1,0.3,0.7),oma=c(0.4,0.4,0.4,2),
    mar=c(1,0.4,0.4,2.0),new=T)
plot(cor_ae, legend.only=TRUE, col=ano_col_gmt,horizontal=F,zlim=c(-1,1),
     legend.width=1.5, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=1))
#text(100000,0, labels = expression(italic(r)), xpd = NA, srt = -90,cex=1.1)
dev.off()






