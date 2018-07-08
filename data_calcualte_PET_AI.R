##### data
###### calcualte the aridity index based on CRU-ts PET and precipitation
library(ncdf4)
library(raster)

setwd("/rigel/glab/users/zy2309/DATA/CRU_TS401/")


get_annual_average<-function(varname){
  var_files<-list.files("./",pattern=varname,full.names = T)[2:3]
  ncin<-nc_open(var_files[1])
  var_data1<-ncvar_get(ncin,varid=varname)[,241:360,]
  dim(var_data1)<-c(120*720,dim(var_data1)[3])
  nc_close(ncin)
  ncin<-nc_open(var_files[2])
  var_data2<-ncvar_get(ncin,varid=varname)[,241:360,]
  dim(var_data2)<-c(120*720,dim(var_data2)[3])
  nc_close(ncin)
  var_dat<-cbind(var_data1,var_data2)
  mean_var<-apply(var_dat,1,mean)
  return(mean_var)
}



export_nc<-function(outdata,outfile){
  latmin<- 30
  latmax<- 90
  latd<- 0.5
  lonmin<- -180
  lonmax<- 180
  lond<- 0.5
  
  lat<- seq(latmin+latd/2,latmax-latd/2,latd)
  long<-seq(lonmin+lond/2,lonmax-lond/2,lond)
  
  dimlat<-ncdim_def('latitude','deg',lat)
  dimlong<-ncdim_def('longitude','deg',long)
  
  ncpet<-ncvar_def('pet','mm/day',list(dimlong,dimlat),-9999,longname="potential evaporation",prec='float',compression=9)
  ncpre<-ncvar_def('pre','mm/month',list(dimlong,dimlat),-9999,longname="precipiation",prec='float',compression=9)
  ncai<-ncvar_def('ai','',list(dimlong,dimlat),-9999,longname="aridity index",prec='float',compression=9)
  
  if(file.exists(outfile)){
    file.remove(outfile)
  }
  
  ncout<-nc_create(outfile,list(ncpet,ncpre,ncai))
  
  dim(outdata)<-c(720,120,3)
  ncvar_put(ncout,varid=ncpet,outdata[,,1])
  ncvar_put(ncout,varid=ncpre,outdata[,,2])
  ncvar_put(ncout,varid=ncai,outdata[,,3])
  
  nc_close(ncout)
}

pet<-get_annual_average("pet")
pre<-get_annual_average("pre")
ai<-pre/(pet*30)
#dim(ai)<-c(120,720)
out_put<-cbind(pet,pre,ai)

outfile<-"/rigel/glab/users/zy2309/PROJECT/SIF_phenology/analysis/climate/ai_cru_P_PET.nc"
export_nc(out_put,outfile)
######## calcualte the average growing season AI
library(ncdf4)
library(raster)


setwd("/rigel/glab/users/zy2309/PROJECT/SIF_phenology")
sos_file<-'./analysis/clear_daily_SOS_30N_fixed_stat.nc'
pos_file<-"./analysis/clear_daily_POS_30N_fixed_stat.nc"
eos_file<-"./analysis/clear_daily_EOS_30N_fixed_stat.nc"

##### get average EOS and SOS
eos_f<-nc_open(eos_file)
eos<-ncvar_get(eos_f,varid="MEAN")
xdim2<-eos_f$dim[["longitude"]]
ydim2<-eos_f$dim[["latitude"]]
nc_close(eos_f)
sos_f<-nc_open(sos_file)
sos<-ncvar_get(sos_f,varid="MEAN")
nc_close(sos_f)

dim(eos)<-c(86400,1)
dim(sos)<-c(86400,1)


setwd("/rigel/glab/users/zy2309/DATA/CRU_TS401/")
get_total<-function(xts){
  ### xts  the 13 and 14 are the sos and eos
  if (is.na(xts[13]*xts[14]))
    return(NA)
  else{
    st_abs<-xts[13]*12
    end_abs<-xts[14]*12
    st<-ceiling(st_abs)
    end<-ceiling(end_abs)
    if (st==end){
      sum_var = xts[st]*(end_abs-st_abs)
    }else if(end==st+1){
      sum_var= xts[st]*(1-st_abs%%1)+xts[end]*(end_abs%%1)
    }else{
      sum_var<-sum(xts[(st+1):(end-1)],na.rm=T)+xts[st]*(1-st_abs%%1)+xts[end]*(end_abs%%1)
    }
    return(sum_var)
  }
}

get_growingseason_average<-function(varname){
  var_files<-list.files("./",pattern=varname,full.names = T)[2:3]
  ncin<-nc_open(var_files[1])
  var_data1<-ncvar_get(ncin,varid=varname)[,241:360,]
  dim(var_data1)<-c(120*720,dim(var_data1)[3])
  nc_close(ncin)
  ncin<-nc_open(var_files[2])
  var_data2<-ncvar_get(ncin,varid=varname)[,241:360,]
  dim(var_data2)<-c(120*720,dim(var_data2)[3])
  nc_close(ncin)
  var_dat<-cbind(var_data1,var_data2)

  # growing season average
  growingseason_ave<-array(NA,dim=c(86400,16))
  for (i in 1:16){
    year_data<-var_dat[,(i*12-11):(i*12)]
    combined<-cbind(year_data,sos,eos)
    growingseason_ave[,i]<-apply(combined,1,get_total)
  }
  mean_grow<-apply(growingseason_ave,1,mean,na.rm=T)
  return(mean_grow)
}

##### need to calcualte the PET and P during the growing season
pet_grow<-get_growingseason_average("pet")
p_grow<-get_growingseason_average("pre")

ai<-p_grow/(pet_grow*30)
#dim(ai)<-c(120,720)
out_put<-cbind(pet_grow,p_grow,ai)

outfile<-"/rigel/glab/users/zy2309/PROJECT/SIF_phenology/analysis/climate/ai_cru_P_PET_grow.nc"
export_nc(out_put,outfile)

######graph ai data
source("/Users/yzhang/Documents/GitHub/Phenology_SIF/graph_ae_tools.R")

plot_nc_var<-function(raster_ae,v_range,title="",id="",color_ramp,position=c(1,1),lay_out=c(1,1)){
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

setwd("/Users/yzhang/Project/SIF_phenology/")
ncin<-nc_open("./data/ai_cru_P_PET.nc")
ai_data<-ncvar_get(ncin,"ai")
nc_close(ncin)
col_ram<-colorRampPalette(drought_discrete)(15)
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/AI_continuous.pdf",width=5,height=3.9)
ae_ai<-nc2ae(ai_data)
plot_nc_var(ae_ai,c(0,1.5),title="Aridity Index",color_ramp = jet.color)
dev.off()
ncin<-nc_open("./data/ai_cru_P_PET_grow.nc")
ai_data<-ncvar_get(ncin,"ai")
nc_close(ncin)
col_ram<-colorRampPalette(drought_discrete)(15)
pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/AI_grow_continuous.pdf",width=5,height=3.9)
ae_ai<-nc2ae(ai_data)
plot_nc_var(ae_ai,c(0,1.5),title="Aridity Index",color_ramp = jet.color)
dev.off()


pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/2018_SIF_phenology/AI_classification.pdf",width=5,height=3.9)
ae_ai<-nc2ae(ai_data)
color_r = rep(jet.color[c(1,25,50,75,100)],each=20)
ae_ai[ae_ai>=0.65]=5
ae_ai[ae_ai>=0.5&ae_ai<0.65]=4
ae_ai[ae_ai>0.2&ae_ai<=0.5]=3
ae_ai[ae_ai>0.05&ae_ai<=0.2]=2
ae_ai[ae_ai<=0.05]=1
plot_nc_var(ae_ai,c(0.5,5.5),title="Aridity Index",color_ramp = color_r)

dev.off()
