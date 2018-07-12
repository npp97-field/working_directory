###plot the correlation between Tmean and Prec, Tday and Prec

#### ----------------------------------------------
# the maximum correlation lag is retrieved and the correspond correlation is graphed
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
  text(10300000,0, labels = expression("lag (months)"), xpd = NA, srt = -90,cex=1.1)   
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.3,font=2,padj=-9.5,las=2)
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(1,0.4,0.4,2.0),new=T)
  plot(raster_ae, legend.only=TRUE, col=rev(color_ramp),horizontal=F,zlim=v_range,
       legend.width=1.5, legend.shrink=0.75,
       axis.args=list(at=0:3,label=c("0.5","1","2","3"),
                      mgp=c(3,0.2,0),tck=0.3,
                      cex.axis=1))
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
notion<-c("a","b","c","d","e","f","g","h")
convertraster2points<-function(pval){
  ras_pval<-raster(apply(pval,1,rev))
  extent(ras_pval)<-c(-180,180,30,90)
  projection(ras_pval)<-longlat
  coarse_pval<-aggregate(ras_pval, fact=4, fun=mean,na.rm=T)
  sig_2deg<-coarse_pval> 0.5
  sig_2deg[sig_2deg==0]<-NA
  sig_point<-coordinates(sig_2deg)[!is.na(values(sig_2deg)),]
  ae_sig_point<-project(sig_point,proj = ae)
  return(ae_sig_point)
}

### export nc
exportnc<-function(dat,fileout){
  longi<-seq(-179.75,179.75,0.5)
  lati<-seq(30.25,89.75,0.5)
  xdim <-ncdim_def("longitude",units = "",vals = longi)
  ydim<-ncdim_def("latitude",units = "",vals = lati)
  dim(dat)<-c(720,120)
  dat_nc<-ncvar_def("max_cor",'',list(xdim,ydim),-999.9,prec="double",compression=9)
  ncout<-nc_create(fileout,list(dat_nc))
  ncvar_put(ncout,dat_nc,dat)
  nc_close(ncout) 
}


vars = c("tday",'temp')
dataset = c("rs","era")
exp_var<-list(expression(paste(italic(r)["T"[day],"Prec"],sep="")),
              expression(paste(italic(r)["T"[mean],"Prec"],sep=""))) 
v_tag<-c("0.5-month","1-month","2-month","3-month")
n_col = 2
n_row = 4
lay_out = c(n_col,n_row)

pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/figures/climate_inter_correlation.pdf",sep=''),
    width=8.5,height=14)
par(oma=c(0,0,1,1))
for (v in 1:length(vars)){
  climate_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",dataset[v],"/",sep=""),
                            pattern=paste("^cor_pre",sep=""),full.names = T)
  nc_var_dat<-array(NA,dim=c(720*120,4))
  nc_var_sig<-array(NA,dim=c(720*120,4))
  for (i in 1:length(climate_files)){
    ncin<-nc_open(climate_files[i])
    ncdat<-ncvar_get(ncin,"cor_coef")
    ncsig<-ncvar_get(ncin,"cor_pv")
    nc_close(ncin)
    nc_var_ae<-nc2ae(ncdat*lc_data1)
    nc_sig_ae<-convertraster2points(ncsig<0.05*lc_data1)
    
    plot_nc_var_sig(nc_var_ae,nc_sig_ae,c(-1,1),v_tag[i],notion[i*2-2+v],rev(ano_col_gmt),
                    c(v,n_row+1-i),c(2,n_row))
  }
}
dev.off()
