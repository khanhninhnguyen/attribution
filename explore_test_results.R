# this program is to explore the test result 
source(paste0(path_code_att,"plot_sic_diff_series.R"))

Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))
all.var = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_var.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
all.arima = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_autocorrelation.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)

# catch the problematic cases  --------------------------------------------

## case of G-G' = 0 --------------------------------------------------------
# data.GG = Total.res[which(abs(Total.res$`tG-G'`)>1.96 & abs(Total.res$`tG-G'`)<3.1 & Total.res$distance <25),]
data.GG = Total.res[which(abs(Total.res$`tG-G'`)<1.96 ),]
all.case.p = data.GG$station
for (i in c(481:length(all.case.p))) {
  plot_six(all.case.p[i])
}
### variance plot 
dat.p = data.frame(mean.msd = (all.var$`mean.G-G'`)/(all.var$`mean.G-E`),
                   sig = as.factor(Total.res$`G-G'`),
                   distance = Total.res$distance)
ggplot(data = dat.p, aes(x = sig, y = distance))+theme_bw()+
  geom_boxplot()
### length plot  
dat.p = data.frame(n = Total.res$n1+ Total.res$n2, sig = as.factor(Total.res$`G-G'`))
ggplot(data = dat.p, aes(x = sig, y = n))+theme_bw()+
  geom_boxplot()

## case of E-E' != 0 --------------------------------------------------------
data.EE = Total.res[which(abs(Total.res$`tE-E'`)>3),]
all.case.p = data.EE$station
for (i in c(1:length(all.case.p))) {
  plot_six(all.case.p[i])
}

dat.p = data.frame(mean.msd = (all.var$`mean.E-E'`)/(all.var$`mean.G-E`), sig = as.factor(Total.res$`E-E'`))
ggplot(data = dat.p, aes(x = sig, y = mean.msd))+theme_bw()+
  geom_boxplot()

dat.p = data.frame(n = Total.res$n1+ Total.res$n2, sig = as.factor(Total.res$`E-E'`))
ggplot(data = dat.p, aes(x = sig, y = n))+theme_bw()+
  geom_boxplot()

## case of G'-E' != 0 --------------------------------------------------------
data.G1E1 = Total.res[which(Total.res$`G'-E'`!=0),]
all.case.p = data.G1E1$station
for (i in c(1:length(all.case.p))) {
  plot_six(all.case.p[i])
}

dat.p = data.frame(mean.msd = (all.var$`mean.G'-E'`)/(all.var$`mean.G-E`), sig = as.factor(Total.res$`G'-E'`))
ggplot(data = dat.p, aes(x = sig, y = mean.msd))+theme_bw()+
  geom_boxplot()



# case of distance 50-100 -------------------------------------------------
lim.list = lapply(seq(0,175, 25), function(x) list(lower = x, upper = x+25))

for (i in c(1:length(lim.list))) {
  limit.dist = lim.list[[i]]
  data.GG = Total.res[which(Total.res$distance > limit.dist$lower & Total.res$distance < limit.dist$upper),]
  data.GG = data.GG[,c(paste0(c("jump", "t", ""), list.name.test[2]), names(Total.res)[19:24])]
  data.GG = cbind( data.GG, all.var[which(Total.res$distance >limit.dist$lower & Total.res$distance <limit.dist$upper), paste0(c("mean.", "range"), list.name.test[2])])
  data.GG = cbind( data.GG, all.arima[which(Total.res$distance >limit.dist$lower & Total.res$distance <limit.dist$upper), paste0(c("phi.", "theta"), list.name.test[2])])
  data.GG$arima = data.GG$`phi.G-G'` + data.GG$`thetaG-G'`
  data.GG[,-c(9, 12, 13)] = abs(data.GG[,-c(9, 12, 13)])
  data.GG$n = data.GG$n1+data.GG$n2
  dat.p = data.GG[,c("G-G'", "distance",  "ver.distance", "mean.G-G'", "rangeG-G'", "phi.G-G'","thetaG-G'", "arima", "n" )]
  colnames(dat.p)[c(1,3:8)] = c("sig", "ver.dist", "mean.var", "range.var", "phi","theta", "phi+theta")
  print(table(dat.p$sig))
  # dat.p[,-1] <- as.data.frame(scale(dat.p[,-1]))
  dat.p[,-1] <- as.data.frame(lapply(dat.p[,-1], min_max_norm))
  dat.p = reshape2::melt(dat.p, id = "sig")
  dat.p$sig = as.factor(dat.p$sig)
  p = ggplot(data = dat.p, aes(x = variable, y = value, col = sig))+theme_bw()+
    geom_boxplot(lwd = 0.3, outlier.size = 0.3, width=0.5)+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank())
  ggsave(paste0(path_results,"attribution0/understand",limit.dist$lower, "GG.jpg" ), plot = p, width = 8, height = 4, units = "cm", dpi = 300)
}

list.case.p = data.GG[which(data.GG$`G-G'`==0),"station"]
for (i in c(1:length(list.case.p))) {
  plot_six_details(name.case = list.case.p[i])
}


#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#apply Min-Max normalization to first four columns in iris dataset


# cases of configuration 35 

data.35 = Total.res[which(Total.res$config==35),]
all.case.p = data.35$station
for (i in c(1:length(all.case.p))) {
  plot_six(all.case.p[i])
  # plot_six_details(all.case.p[i])
}

a = read.table(file = paste0(path_results, "validation/34-1BM_BJ1062.txt"), header = TRUE, stringsAsFactors = FALSE)
# a = a[which(a$valid==0),]
list.detection =  as.Date(a$detected,format="%Y-%m-%d")
month.year = as.Date(list.detection, d,format="%Y-%m")
data.p = data.frame(date = list.detection, month.year = month.year, 
                    month = month(as.POSIXlt(list.detection, format = "%Y-%m-%d")),
                    year = year(as.POSIXlt(list.detection, format = "%Y-%m-%d")))

p = data.p %>%
  ggplot(aes(x = month)) + theme_bw()+
  geom_histogram(bins = 12,  binwidth = 0.5) +
  scale_x_continuous(breaks = seq(1,12,1), limits = c(0,13))+
  theme(axis.text.x = element_text(size = 5), 
        axis.text.y = element_text(size = 5),
        legend.text=element_text(size=4),
        axis.title = element_text(size = 5), 
        legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 5), 
        plot.subtitle = element_text(size = 5),
        legend.title=element_blank(), 
        legend.position = "none")+
  facet_wrap(~ year, ncol = 4)
  
ggsave(paste0(path_results,"attribution0/dist.detection.jpg" ), plot = p, width = 20, height = 15, units = "cm", dpi = 1200)



a = get(load(file = paste0(path_results, "attribution/check_homo/six_diff_series_rm_crenel_restricted_closed_brp_10year_NGL.RData")))
plot(a$`kour.2018-01-05.koug`$era.era)
plot_six("vill.2001-03-26.madr")



# contradiction  --------------------------------------------------------
contradicted = Total.res[which(Total.res$config==0),]
contradicted.5 = contradicted[,c(14:18)]
a = contradicted.5 %>% group_by_all() %>% summarise(COUNT = n())
##  with distance  --------------------------------------------------------

dat.p = Total.res[,c(19:21)]
dat.p$contra = ifelse(dat.p$config!=0, 0, 1)
contra.ratio = sapply(seq(0, 190, 10), function(x){
  b = dat.p[which(dat.p$distance>0 & dat.p$distance<(x+10)),]
  return( length(which(b$contra==1))/ nrow(b))
})

plot(seq(10, 200, 10), contra.ratio, xlab ="Distance (km)", ylab = "Rate of contradiction")
dat.p$contra = as.factor(dat.p$contra)
dat.p$n = Total.res$n1+Total.res$n2
dat.p = dat.p[which(dat.p$config==0),]
ggplot(dat.p, aes(x = distance, y = n, col = contra))+
  theme_bw()+
  geom_point()

## case of contradiction --------------------------------------------------------
data.EE = Total.res[which(Total.res$config==0 & Total.res$`G-G'`==0),]
all.case.p = data.EE$station
for (i in c(1:length(all.case.p))) {
  plot_six(all.case.p[i])
}
b = inner_join(a[17,], Total.res, by = c( "G-G'" , "G-E'",  "E-E'" , "G'-E'", "G'-E" ))
e = d[which(Total.res$station %in% b$station == TRUE)]


