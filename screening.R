source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))
# This function is used for data screening ---------------------------------
# first normalized data
window.thres = 2
data.cr = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
name.var = list.test[2]
dist.mean <- data.frame(matrix(NA, nrow = 0, ncol = (length(seq(-30,30,0.05))-1)))
data.all <- list()

for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  data.i = data.cr[[i]]
  data.i <- tidyr::complete(data.i, date = seq(min(data.i$date), max(data.i$date), by = "day"))
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  before =  data.i[which(data.i$date <= breakpoint),]
  after =  data.i[which(data.i$date > breakpoint),]
  if(nrow(na.omit(before)) > 30){
    ind.sta = which(before$date == max(breakpoint %m+% years(-1), min(before$date)))
    if(length(ind.sta)>0){
      bef.norm.all <- list()
      for (k in c(1:6)) {
        bef.norm = one.step.norm(before, name.var = list.test[k], estimator = "Sca", length.wind = 30) 
        names(bef.norm) <- NULL
        bef.norm.all[[list.test[k]]] <- bef.norm
      }
      bef.norm.all[["date"]] = before$date
      bef.norm.all <- as.data.frame(bef.norm.all)
      bef.norm.all <- bef.norm.all[c(ind.sta:length(bef.norm)),]
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$bef <- bef.norm.all
    }
  }
  if(nrow(na.omit(after)) > 30){
    ind.end = which(after$date == min(breakpoint %m+% years(1), max(after$date)))
    if(length(ind.end) > 0){
      aft.norm.all <- list()
      for (k in c(1:6)) {
        aft.norm = one.step.norm(after, name.var = list.test[k], estimator = "Sca", length.wind = 30) 
        names(aft.norm) <- NULL
        aft.norm.all[[list.test[k]]] <- aft.norm
      }
      aft.norm.all[["date"]] = after$date
      aft.norm.all <- as.data.frame(aft.norm.all)
      aft.norm.all = aft.norm.all[c(1:ind.end),]
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$aft <- aft.norm.all
    }
  }
}

save(data.all, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"normalized2.RData"))

# Study the distribution of all data sets ----------------------------

dat = get(load( file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"normalized.RData")))
# concatenate all data
all.dat <- c()
zero <- data.frame(matrix(NA, ncol = 2, nrow = length(dat)))
for (r in 1:length(dat)) {
  if(length(dat[[r]]$bef$gps.gps)>30){
    ind1 = length(dat[[r]]$bef$gps.gps)-30
    all.dat = c(all.dat, dat[[r]]$bef$gps.gps[30:ind1])
  } 
  if(length(dat[[r]]$aft$gps.gps)>30){
    ind2 = length(dat[[r]]$aft$gps.gps)-30
    all.dat = c(all.dat, dat[[r]]$aft$gps.gps[30:ind2])
  }
  if (length(dat[[r]]$bef$gps.gp) >0){
    zero[r,1] <- length(which(dat[[r]]$bef$gps.gps<=0 & dat[[r]]$bef$gps.gps>-0.02))
  }
  if (length(dat[[r]]$aft$gps.gp) >0){
    zero[r,2] <- length(which(dat[[r]]$aft$gps.gps<=0 & dat[[r]]$aft$gps.gps>-0.02))
  }
}
x = na.omit(all.dat)

a = hist(x,
         main = "Histogram of normalized data",
         xlab = "",
         breaks = 1000,
         xlim=c(-5, 5),
         prob = TRUE)

y <- x[-which(x>-0.008 & x<0.008)] 
fit = Mclust(y, G=2, model="V")
param = fit$parameters
# # add the normal density 
# x2 <- a$breaks
# fun1 <- dnorm(x2, mean = 0, sd = 1)
# lines(x2, fun1, col = 2, lwd = 2)
# QQ PLOT 
y <- y[-which(abs(y)>8)] 
jpeg(paste0(path_results,"attribution/" , "histogram of all data.jpeg"),
     width = 3000, height = 2000,res = 300)
  data.p = data.frame(x = y)
  colors <- c("N(0,1)" = "#CD3700", "Fit1" = "#FF3E96", "Fit2" = "#6495ED", "Fit12" = "#000000")
  p <- ggplot(data.p, aes(x=x)) +   theme_bw()+ 
    geom_histogram(aes(y =..density..), 
                   breaks = seq(-5, 5, by = 0.1),
                   fill="#69b3a2",
                   color="#e9ecef", 
                   alpha=0.9)+
    # xlim(c(-7.5,7.5))+
    labs(y = "Density", x = " GPS-GPS' ",
         subtitle = "Normalized with median + ScaleTau / Fit by Mclust ")+
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size=15,face="bold"),
          plot.subtitle = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
  p <- p +
    geom_function(fun = function(x) {
      param$pro[1]*dnorm(x = x,mean = param$mean[1], sd = sqrt(param$variance$scale[1])) +
        param$pro[2]*dnorm(x = x,mean = param$mean[2] , sd = sqrt(param$variance$scale[2]))
      }, aes(col = "Fit12"))+
  geom_function(fun = function(x) dnorm(x = x,mean = 0, sd =1), aes(col = "N(0,1)"))+
  geom_function(fun = function(x) param$pro[1]*dnorm(x = x,mean = param$mean[1], sd = sqrt(param$variance$scale[1])), aes(col = "Fit2"))+
  geom_function(fun = function(x) param$pro[2]*dnorm(x = x,mean = param$mean[2], sd = sqrt(param$variance$scale[2])), aes(col = "Fit1"))+
  scale_color_manual(values = colors)
  
  print(p)
dev.off()

n= length(y)
a =  rnorm(n,mean = param$mean[2] , sd = sqrt(param$variance$scale[2]))
ind.out = rbinom(n, 1, 0.2299717)
a[which(ind.out!=0)] <- rnorm(length(which(ind.out!=0)),mean = param$mean[1], sd = sqrt(param$variance$scale[1])) 

library(tidyverse)

set.seed(10)
dat <- data.frame(Observed = y, 
                  mixture = a, 
                  normal = rnorm(n, 0 ,1))

plot_data <- map_dfr(names(dat)[-1], ~as_tibble(qqplot(dat[[.x]], dat$Observed, plot.it = FALSE)) %>% 
                       mutate(id = .x))
jpeg(paste0(path_results,"attribution/" , "QQplots.jpeg"),
     width = 3000, height = 1500,res = 300)
ggplot(plot_data, aes(x, y, color = id)) + theme_bw()+
  geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~id)+
  labs(y = "", x = " Observed ")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()

# CDF PLOT  

df <- data.frame(x = c(y, a, rnorm(n, 0 ,1)), val=factor(rep(c("Observed", "Mixture", "Normal"), c(n,n,n))))
jpeg(paste0(path_results,"attribution/" , "PPplots.jpeg"),
     width = 3000, height = 1500,res = 300)
ggplot(df, aes(x, colour = val, linetype = val)) + theme_bw()+
  stat_ecdf()+
  labs(y = "CDF", x = " ")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()






# compute the skewness of data --------------------------------------------
dat = get(load( file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"normalized.RData")))


