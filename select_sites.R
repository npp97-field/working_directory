setwd("/Users/yzhang/Project/SIF_phenology")

site_list<-read.csv("./retrieve_site/flux_sites_with_lc_and_EVI_stat.csv")

mean_val<-abs(site_list$central_mean_VI-site_list$all_mean_VI)
mean_ind<-abs(site_list$central_mean_VI-site_list$all_mean_VI)/abs(site_list$all_mean_VI)
sd_val<-abs(site_list$central_mean_SD-site_list$all_mean_SD)
sd_ind<-abs(site_list$central_mean_SD-site_list$all_mean_SD)/abs(site_list$all_mean_SD)
max_val<-abs(site_list$central_max_VI-site_list$all_mean_max)
max_ind<-abs(site_list$central_max_VI-site_list$all_mean_max)/abs(site_list$all_mean_max)
min_val<-abs(site_list$central_min_VI-site_list$all_mean_min)
min_ind<-abs(site_list$central_min_VI-site_list$all_mean_min)/abs(site_list$all_mean_min)

sel_sites<-site_list[((mean_ind<0.2)|(mean_val<500))&((sd_ind<0.2)|(sd_val<500))
                     &((max_ind<0.2)|(max_val<500))
                     &site_list$central_max_VI>2000&site_list$all_mean_max>2000,]
write.csv(sel_sites,"./retrieve_site/sites_passed_heterogeneity_test.csv")
