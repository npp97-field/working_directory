#### analyze the controlling factor of SOS and EOS
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")


discrete_col<-rep(drought_discrete[c(1,3,5,7)],each=20)

plot_nc_var<-function(raster_ae,v_range,title="",id="",color_ramp,position,lay_out){
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(0.4,0.4,0.4,3.4),mgp=c(3,0.3,0),new=T)
  plot(border)
  image(setrange(raster_ae,v_range[1],v_range[2]),add=T,col=rev(color_ramp),
        axes=F,zlim=v_range)
  text(0, 5500000,title,cex=1.3)
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.8,font=2,padj=-7,las=2)
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(0.4,0.4,0.4,0.8),new=T)
  plot(raster_ae, legend.only=TRUE, col=rev(color_ramp),horizontal=F,zlim=v_range,
       legend.width=1.5, legend.shrink=0.75,
       axis.args=list(
         mgp=c(3,0.2,0),tck=0.3,
         cex.axis=1))
}

plot_nc_var_sig<-function(raster_ae,sig_ae,v_range,title="",id="",color_ramp,position,lay_out){
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(0.4,0.4,0.4,3.4),mgp=c(3,0.3,0),new=T)
  plot(border)
  image(setrange(raster_ae,v_range[1],v_range[2]),add=T,col=rev(color_ramp),
        axes=F,zlim=v_range)
  points(sig_ae,pch=15,cex=0.2)
  text(0, 5500000,title,cex=1.3)
  plotlatlong()
  mtext(side=2,line=-1.5,id,cex=1.8,font=2,padj=-7,las=2)
  par(fig=c((position[1]-1)/lay_out[1],(position[1]/lay_out[1]),
            (position[2]-1)/lay_out[2],(position[2]/lay_out[2])),
      mar=c(0.4,0.4,0.4,0.8),new=T)
  plot(raster_ae, legend.only=TRUE, col=rev(color_ramp),horizontal=F,zlim=v_range,
       legend.width=1.5, legend.shrink=0.75,
       axis.args=list(
         mgp=c(3,0.2,0),tck=0.3,
         cex.axis=1))
}
## plot_nc_var input (raster_ae, v_range (two value), title="", id="",
##                    color_ramp, position (col, row),  lay_out (col,row))
#compare the maximum correlation
graph_lags<-function(dataset){
  if (dataset=="rs"){
    vars = c("tday","tnight",'prec','par')
  }else{
    vars = c("temp","prec","par")
  }
  n_col = 1
  n_row = length(vars)  
  lay_out = c(n_col,n_row)
  sos_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",dataset,"/",sep=""),
                        pattern="cor_sos_pre",full.names = T)

  for (v in 1:length(vars)){
    var_file<-sos_files[grep(vars[v],basename(sos_files))]
    pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/figures/SOS_",
              dataset,'_',vars[v],".pdf",sep=''),width=4.5,height=14)
    for (i in 1:length(var_file)){
      ncin<-nc_open(var_file[i])
      ncdat<-ncvar_get(ncin,"cor_coef")
      nc_var_ae<-nc2ae(ncdat)
      plot_nc_var(nc_var_ae,c(-1,1),paste("SOS_",dataset,'_',vars[v],"_lag",i-1,sep=''),notion[i],rev(ano_col_gmt),
                  c(1,5-i),c(1,4))
    }
    dev.off()
  }
  
  eos_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",dataset,"/",sep=""),
                        pattern="cor_eos_pre",full.names = T)
  
  for (v in 1:length(vars)){
    var_file<-eos_files[grep(vars[v],basename(eos_files))]
    pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/figures/EOS_",
              dataset,'_',vars[v],".pdf",sep=''),width=4.5,height=14)
    for (i in 1:length(var_file)){
      ncin<-nc_open(var_file[i])
      ncdat<-ncvar_get(ncin,"cor_coef")
      nc_var_ae<-nc2ae(ncdat)
      plot_nc_var(nc_var_ae,c(-1,1),paste("EOS_",dataset,'_',vars[v],"_lag",i-1,sep=''),notion[i],rev(ano_col_gmt),
                  c(1,5-i),c(1,4))
    }
    dev.off()
  }
}

graph_lags('rs')
graph_lags('era')

#### ----------------------------------------------
# the maximum correlation lag is retrieved and the correspond correlation is graphed
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


graph_maximum_lags<-function(dataset,phen="sos"){
  if (dataset=="rs"){
    vars = c("tday","tnight",'prec','par')
    exp_var<-list(expression("T"[day]),expression("T"[night]),expression("Precip."),expression("PAR"))
  }else{
    vars = c("temp","prec","par")
    exp_var<-list(expression("T"[mean]),expression("Precip."),expression("PAR"))
  }
  n_col = 1
  n_row = length(vars)  
  lay_out = c(n_col,n_row)
  phen_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",dataset,"/",sep=""),
                        pattern=paste("^cor_",phen,"_pre",sep=""),full.names = T)
  pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/figures/",phen,"_",
            dataset,"_max_lag.pdf",sep=''),width=8.5,height=14/4*n_row)
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
    
    nc_var_ae<-nc2ae(max_cor)
    nc_max_ae<-nc2ae(maxind)-1
    nc_sig_ae<-convertraster2points(max_sig<0.05)
    
    ##export maximum correaltion
    exportnc(max_cor,paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",
                              dataset,"/max_correlation_",phen,"_",vars[v],".nc",sep=""))

    plot_nc_var_sig(nc_var_ae,nc_sig_ae,c(-1,1),exp_var[[v]],notion[v*2-1],rev(ano_col_gmt),
                    c(1,n_row+1-v),c(2,n_row))
    plot_nc_var(nc_max_ae,c(-0.5,3.5),exp_var[[v]],notion[v*2],rev(discrete_col),
                    c(2,n_row+1-v),c(2,n_row))
  }
  dev.off()
}

graph_maximum_lags(dataset = "rs")
graph_maximum_lags(dataset = "era") 

graph_maximum_lags(dataset = "rs","eos")
graph_maximum_lags(dataset = "era","eos") 








