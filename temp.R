csif<-read.csv("/Users/yzhang/Project/SIF_phenology/analysis/csif_stat/CSIF.csv")
osif<-read.csv("/Users/yzhang/Project/SIF_phenology/analysis/csif_stat/OCO_2_SIF.csv")
ndvi<-read.csv('/Users/yzhang/Project/SIF_phenology/analysis/csif_stat/ndvi.csv')
evi<-read.csv('/Users/yzhang/Project/SIF_phenology/analysis/csif_stat/evi.csv')
csif_yd<-rep(2014:2017,each=92)*1000+rep(1:92*4,4)-2
csif_yd[c(92,184,276,368)]<-csif_yd[c(92,184,276,368)]-1
csif_date<-strptime(csif_yd,format = "%Y%j")

osif_yd<-osif$yearmonth*100+15
osif_date<-strptime(osif_yd,format= "%Y%m%d")
 
vi_yd<-rep(2014:2017,each=23)*1000+rep(1:23*16,4)-8
vi_date<-strptime(vi_yd,format = "%Y%j")

plot(NA,ylim=c(0,0.3),xlim=c(as.POSIXct(csif_date[1]),as.POSIXct(csif_date[368])))
lines(osif_date,osif$all)
lines(csif_date,csif$all,col='red')
lines(vi_date,ndvi$all/20000)
