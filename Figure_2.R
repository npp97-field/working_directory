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


vars = c("tday",'prec')
exp_var<-list(expression("T"[day]),expression("Precip.")) ###,expression("T"[night]) ,expression("PAR")
phen="eos"
dataset="rs"
n_col = 1
n_row = length(vars)  
lay_out = c(n_col,n_row)
phen_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",dataset,"/",sep=""),
                       pattern=paste("^cor_",phen,"_pre",sep=""),full.names = T)
pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/figures/",phen,"_",
          dataset,"_max_lag.pdf",sep=''),width=8.5,height=14/4*n_row)
par(oma=c(0,0,1,1))
for (v in 1:length(vars)){
  var_file<-phen_files[grep(vars[v],basename(phen_files))]
  nc_var_dat<-array(NA,dim=c(720*120,4))
  nc_var_sig<-array(NA,dim=c(720*120,4))
  for (i in 1:length(var_file)){
    ncin<-nc_open(var_file[i])
    ncdat<-ncvar_get(ncin,"cor_coef")
    ncsig<-ncvar_get(ncin,"cor_pv")
    nc_close(ncin)
    nc_var_dat[,i]<-ncdat
    nc_var_sig[,i]<-ncsig
  }
  abs_cor<-abs(nc_var_dat)
  maxindex<-apply(abs_cor,1,which.max)  #### the index where the maximum happens
  idx <- !(sapply(maxindex, length))
  maxindex[idx]<-NA
  maxind<-unlist(maxindex)
  #
  sel<-as.matrix(cbind(1:86400,maxind))
  sign_cor<-sign(nc_var_dat[sel])
  max_cor<-apply(abs_cor,1,max)*sign_cor
  max_sig<-nc_var_sig[sel]
  dim(maxind)<-c(720,120)
  dim(max_cor)<-c(720,120)
  dim(max_sig)<-c(720,120)
  
  nc_var_ae<-nc2ae(max_cor*lc_data1)
  nc_max_ae<-nc2ae(maxind*lc_data1)-1
  nc_sig_ae<-convertraster2points(max_sig<0.05*lc_data1)
  
  plot_nc_var_sig(nc_var_ae,nc_sig_ae,c(-1,1),exp_var[[v]],notion[v*2-1],rev(ano_col_gmt),
                  c(1,n_row+1-v),c(2,n_row))
  plot_nc_var(nc_max_ae,c(-0.5,3.5),exp_var[[v]],notion[v*2],rev(discrete_col),
              c(2,n_row+1-v),c(2,n_row))
  data_all<-cbind(as.vector(maxind*lc_data1),as.vector(sign_cor*lc_data1),as.vector(max_sig*lc_data1))
  data_sign<-data_all[complete.cases(data_all),]
  #neg_sig<-data_sign[data_sign[,2]<0|data_sign[,3]<0.05,1]
  pos_sig<-data_sign[data_sign[,2]>0|data_sign[,3]<0.05,1]
  ct<-table(pos_sig)
  lbls<-c("","","","")
  par(fig = c(0.55,0.72,1-v/n_row,1-v/n_row+0.17),new=T)
  pie(ct,labels = lbls,col=drought_discrete[c(1,3,5,7)],clockwise=T)
}

dev.off()
