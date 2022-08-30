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

dat = get(load( file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"normalized2.RData")))
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

# a = hist(x, 
#          main = "Histogram of normalized data", 
#          xlab = "", 
#          breaks = 200,
#          xlim=c(-5, 5),
#          prob = TRUE)
# # add the normal density 
# x2 <- a$breaks
# fun1 <- dnorm(x2, mean = 0, sd = 1)
# lines(x2, fun1, col = 2, lwd = 2)
jpeg(paste0(path_results,"attribution/" , "histogram of all data.jpeg"),
     width = 3000, height = 2000,res = 300)
  data.p = data.frame(x = x)
  p <- ggplot(data.p, aes(x=x)) +   theme_bw()+ 
    geom_histogram(aes(y =..density..),
                   breaks = seq(-5, 5, by = 0.1),
                   fill="#69b3a2",
                   color="#e9ecef", 
                   alpha=0.9)+
    # xlim(c(-5,5))+
    labs(y = "Density", x = " GPS-GPS' ",
         subtitle = "Normalized with trimmed mean (10%) + ScaleTau ")
  p <- p + stat_function(fun = dnorm, args = list(mean = -0.009158155, sd = sqrt(3.1588313)), color = "orange")+
    stat_function(fun = dnorm, args = list(mean =0.002219531 , sd = sqrt(0.7382873)), color = "orange")
  print(p)
dev.off()




