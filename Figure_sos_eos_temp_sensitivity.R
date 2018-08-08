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
ncin<-nc_open("./data/ai_cru_P_PET.nc")
aridity_mask<-ncvar_get(ncin,"ai")
aridity_mask[aridity_mask<0.6]<-NA
aridity_mask[!is.na(aridity_mask)]<-1
aridity_shade<-ncvar_get(ncin,"ai")
aridity_shade[aridity_shade>0.6|is.na(lc_data1)]<-NA
aridity_shade[!is.na(aridity_shade)]<-1
nc_close(ncin)

discrete_col<-rep(drought_discrete[c(1,3,5,7)],each=20)

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
  mtext(side=2,line=-1.5,id,cex=1.3,font=2,padj=-9.5,las=2)

}

notion<-c("a","b","c","d","e","f","g","h")
bin64<-colorRampPalette(ano_col_gmt)(64)
vars = c("_t_")
exp_var<-list(expression(gamma[T[air]]^SOP),expression(gamma[T[air]]^EOP)) ###,expression("T"[night]) ,expression("PAR")
dataset<-c('rs','era')
n_row = length(vars)
n_col=2
lay_out = c(n_col,n_row)
vrange<-list(c(-3,3),c(-3,3))
bin64[c(32,33)]<-"white"
pdf(paste("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/Fig/",
          "temperature_sensitivity_SOS_EOS.pdf",sep=''),width=8.5,height=4.2)
par(oma=c(0.4,0.4,0.4,4))
for (d in 1:2){
  phen_files<-list.files(paste("/Users/yzhang/Project/SIF_phenology/analysis/correlation_clear_era/",sep=""),
                         pattern=paste("^reg",sep=""),full.names = T)[c(4,1)]
  ncin<-nc_open(phen_files[d])
  ncdat<-ncvar_get(ncin,"regress_coef")[2,,]
  nc_close(ncin)
  
  nc_var_ae<-nc2ae(ncdat*lc_data1*aridity_mask)
  #nc_sig_ae<-convertraster2points(ncsig<0.05)  #*lc_data1

  aridity_ae<-nc2ae(aridity_shade)
  plot_nc_var(nc_var_ae,vrange[[1]],exp_var[[d]],notion[d],rev(bin64),
              c(d,n_row),c(2,n_row))
  image(aridity_ae,col=adjustcolor("grey50",alpha.f = 0.1),add=T)
}

text(9200000,0, labels = expression(paste("day/",degree,"C",sep="")), xpd = NA, srt = -90,cex=1.1)   
par(fig=c(0.5,1,0,1),oma=c(0.4,0.4,0.4,1.4),
    mar=c(1,0.4,0.4,1),new=T)
plot(nc_var_ae, legend.only=TRUE, col=bin64,horizontal=F,zlim=vrange[[1]],
     legend.width=1.5, legend.shrink=0.75,
     axis.args=list(
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=1))
dev.off()

