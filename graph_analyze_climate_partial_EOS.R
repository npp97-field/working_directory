####graph analyze climate and EOS partial correlation

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
notion<-c("a","b","c","d","e","f","g","h")
graph_lags<-function(dataset){
  if (dataset=="rs"){
    vars = c("tday",'prec','par')
    exp_var<-list(expression("T"[day]),expression("Precip."),expression("PAR"))
  }else{
    vars = c("temp","prec","par")
    exp_var<-list(expression("T"[mean]),expression("Precip."),expression("PAR"))
  }
  n_col = 1
  n_row = length(vars)  
  lay_out = c(n_col,n_row)
  
  eos_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_",dataset,"/",sep=""),
                        pattern="pcor_eos_pre",full.names = T)
  pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/figures/pcor_EOS_",
            dataset,".pdf",sep=''),width=4.5,height=14/4*n_row)
  for (v in 1:length(vars)){
    var_file<-eos_files[grep(vars[v],basename(eos_files))]
    ncin<-nc_open(var_file)
    ncdat<-ncvar_get(ncin,"pcor_coef")
    nc_var_ae<-nc2ae(ncdat)
    ncsig<-ncvar_get(ncin,"pcor_pv")
    nc_sig_ae<-convertraster2points(ncsig<0.05)
    nc_close(ncin)
    plot_nc_var_sig(nc_var_ae,nc_sig_ae,c(-1,1),exp_var[[v]],notion[v],rev(ano_col_gmt),
                    c(1,n_row+1-v),c(1,n_row))
  }
  dev.off()
}
graph_lags("rs")
graph_lags("era")
