source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))

##### Study in the raw data #####
# This function is used for data screening
# first normalized data ---------------------------------
window.thres = 10
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
        bef.norm = one.step.norm(before, name.var = list.test[k], estimator = "Sca", length.wind = 30, loes = 0) 
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
        aft.norm = one.step.norm(after, name.var = list.test[k], estimator = "Sca", length.wind = 30, loes = 0) 
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
save(data.all, file = paste0(path_results,"attribution/data.all_1year_", nearby_ver,"paired.normalized.RData"))

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
fit = Mclust(y, G=3, model="V")
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
 
res = rep(NA, length(dat))
for (i in c(1:length(dat))) {
  y = dat[[i]]$bef$gps.gps
  n = length(na.omit(y))
  if(n>100){
    numer = sum((y - mean(y, na.rm = TRUE))**3, na.rm = TRUE)
    denom = (n-1)*(robustbase::scaleTau2(na.omit(y)))**3
    sk = numer/denom
    res[i] <- sk
  }
}
boxplot(res)







#### Study on the pair data ####
# read the paired data (after - before)

data.p = get(load(file = paste0(path_results,"attribution/data.all_1year_", nearby_ver,"paired.normalized.RData")))
list1 <- list()
list2 <- list()
for (i in c(1:length(data.p))) {
  data.i = data.p[[i]]
  data.i <- tidyr::complete(data.i, date = seq(min(data.i$date), max(data.i$date), by = "day"))
  a <- screen.O(Y = data.i, name.var = "gps.era", method = 'sigma', iter = 1, estimator = "Sca", fix.thres = 0, loes = 0)
  a1 <- screen.O(Y = data.i, name.var = "gps.era", method = 'sigma', iter = 1, estimator = "Sca", fix.thres = 0, loes = 1)
  
  list1[[i]] <- a$point.rm
  list2[[i]] <- a1$point.rm
  
}

d1 <- sapply(c(1:length(data.p)), function(x) length(list1[[x]]))
d2 <- sapply(c(1:length(data.p)), function(x) length(list2[[x]]))

# concatenate all data
all.dat <- c()
zero <- data.frame(matrix(NA, ncol = 2, nrow = length(dat)))
for (r in 1:length(dat)) {
  if(length(dat[[r]]$bef$gps.era1)>30){
    all.dat = c(all.dat, dat[[r]]$bef$gps.era1)
  } 
  if (length(dat[[r]]$bef$gps.era1) >0){
    zero[r,1] <- length(which(dat[[r]]$bef$gpa.era1<=0 & dat[[r]]$bef$gpa.era1>-0.02))
  }
}
x = na.omit(all.dat)

a = hist(x,
         main = "Histogram of normalized data",
         xlab = "",
         breaks = 1000,
         xlim=c(-5, 5),
         prob = TRUE)

y <- x[-which(x==0)] 
fit = Mclust(y, G=2, model="V")
param = fit$parameters
text1 = paste0("pro",  ": ", paste(round(param$pro, digits = 2), collapse = " & "), 
               ",   mean", ": ", paste(round(param$mean, digits = 2), collapse = " & "),
              ",   variance", ": ", paste(round(param$variance$scale, digits = 2), collapse = " & "))

# # add the normal density 
# x2 <- a$breaks
# fun1 <- dnorm(x2, mean = 0, sd = 1)
# lines(x2, fun1, col = 2, lwd = 2)
# QQ PLOT 
y <- y[-which(abs(y)>8)] 
jpeg(paste0(path_results,"attribution/" , "histogram of all paired normalized data.jpeg"),
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
       title = "Normalized with median + ScaleTau / Fit by Mclust ", 
       subtitle = text1)+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"),
        plot.title = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
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


# FINAL SCREENING - RUN ONLY THIS SECTION  --------------------------------------------------------

# read data : can be paired or not 
window.thres = 10

data.in = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
data.all <- list()
list.outlier <- list()
sd.all <- list()
for (i in c(1:length(data.in))) {
  # read data
  case.name = names(data.in)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  data.i = as.data.frame(data.in[[i]])
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  begin = breakpoint - (window.thres*365-1)
  fin = breakpoint + window.thres*365
  data.i = tidyr::complete(data.i, date = seq(begin, fin, by = "day"))
  
  # screen only series has data 
  before1 = na.omit(data.i[which(data.i$date <= breakpoint),])
  after1 = na.omit(data.i[which(data.i$date > breakpoint),])
  condi = 0
  if(nrow(before1)>100){
    before = data.i[which(data.i$date <= breakpoint),]
    bef.all = list()
    bef.outlier = list()
    bef.sd = list()
    
    for (k in c(1:6)) {
      bef.scr <- screen.O(Y = before, name.var = list.test[k], method = 'sigma', global.mu = 0, iter = 1, estimator = "Sca", fix.thres = 0, loes = 0, loes.method = 0)
      bef.all[[list.test[k]]] <- bef.scr$data
      bef.outlier[[list.test[k]]] <- bef.scr$point.rm
      bef.sd[[list.test[k]]] <- bef.scr$sd.est[[length(bef.scr$sd.est)]]
    }
    bef.all[["date"]] = before$date
    bef.all = as.data.frame(bef.all)
    rownames(bef.all) <- NULL
    condi = condi+1
  }
  if(nrow(after1) > 100){
    after = data.i[which(data.i$date > breakpoint),]
    aft.all = list()
    aft.outlier = list()
    aft.sd = list()
    
    for (k in c(1:6)) {
      aft.scr <- screen.O(Y = after, name.var = list.test[k], method = 'sigma',  global.mu = 0, iter = 1, estimator = "Sca", fix.thres = 0, loes = 0, loes.method = 0)
      aft.all[[list.test[k]]] = aft.scr$data
      aft.outlier[[list.test[k]]] = aft.scr$point.rm
      aft.sd[[list.test[k]]] = aft.scr$sd.est[[length(aft.scr$sd.est)]]
    }
    aft.all[["date"]] = after$date
    aft.all = as.data.frame(aft.all)
    rownames(aft.all) <- NULL
    condi = condi+1
  }
  if(condi == 2){
    data.all[[case.name]] <- rbind(bef.all, aft.all)
    list.outlier[[case.name]] <- list(bef = bef.outlier, aft = aft.outlier)
    sd.all[[case.name]] <- list(bef = bef.sd, aft = aft.sd)
  }
}
save(data.all, file = paste0(path_results,"attribution/data.all_", window.thres,"years_", nearby_ver,"screened.RData"))
save(list.outlier, file = paste0(path_results,"attribution/list.outlier_",  window.thres,"years_", nearby_ver,"screened.RData"))
save(sd.all, file = paste0(path_results,"attribution/sd.all_",  window.thres,"years_", nearby_ver,"screened.RData"))

  
# SCREEN A REAL TIME SERIES -----------------------------------------------
case.name = c("albh.2015-12-29.bcns")
station.ref = substr(case.name ,start = 1, stop = 4)
station.near = substr(case.name ,start = 17, stop = 20)
breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")

Y = data.in[[case.name]]
y = Y[which(Y$date > (breakpoint-364) & Y$date <= breakpoint),]
var.name = "gps.era"

met1 <- screen.O(Y = y, name.var = var.name, method = 'sigma', iter = 1, estimator = "Sca", fix.thres = 0, 
                 loes = 0, loes.method ="", global.mu = 0)

# visualize 
list.method = list( met1)

for (j in c(1:length(list.method))) {
  method = list.method[[j]]
  mu = as.data.frame(method$mu.est, col.names = paste("mu", c(1:length(method$mu.est))))
  sd1 = as.data.frame(method$sd.est, col.names = paste("sd", c(1:length(method$mu.est))))
  data.i = cbind(date = y$date, mu)
  for(x in c(1:length(method$normalized))) {
    data.i[paste("lower",x, sep = "")] = mu[,x] - 3* sd1[,x]
    data.i[paste("upper",x, sep = "")] = mu[,x] + 3* sd1[,x]
    series.rm <- rep(NA, length(mu[,x]))
    if(x < length(method$normalized)){
      series.rm[method$point.rm[[x]]] <- y[method$point.rm[[x]], var.name]
    }
    data.i[paste("rm",x, sep = "")] = series.rm 
  }
  data.i$y = unlist(y[var.name])
  
  Y.i <- data.i
  a = reshape2::melt(Y.i, id = "date")
  gr = c(c(1:length(method$normalized)), rep(c(1:length(method$normalized)), each = 3), 0)
  gr[c(1:length(method$normalized))*3+length(method$normalized)] <- 10
  a$col <- as.factor(rep(gr, each = 365 ))
  b = a[which(a$col == 10),]
  d = a[which(a$col != 10),]
  b$col <- as.factor(rep(c(1:length(method$normalized)), each = 365 ))
  p <- ggplot(data = d, aes(x = date, y = value, group = variable))+
    geom_line(aes(colour= col))+ theme_bw()
  p <- p + geom_point( data = b, aes(y = value, colour= col), size=2)
  p <- p + theme(axis.text = element_text(size = 15),
                 axis.title = element_text(size=15,face="bold"),
                 plot.subtitle = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
  jpeg(paste0(path_results,"attribution/" , case.name, ".screened.", j, ".jpeg"),
       width = 3000, height = 1500,res = 300)
 
  print(p)
  dev.off()
  
}

# first iteration
# iter = 1
# method = met1
# norm1 = method$normalized[[iter]]
# mu1 = method$mu.est[[iter]]
# sd1 = method$sd.est[[iter]]
# point.rm = method$point.rm[[iter]]
# 
# a = unlist(y[var.name])
# plot(a, type = "l", ylim = c(-6, 4))
# lines(mu1)
# lines(mu1 + 3*sd1)
# lines(mu1 - 3*sd1)











