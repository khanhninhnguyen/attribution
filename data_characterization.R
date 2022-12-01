#### This function is used for the data characterization. This is step after pairing and screening data ####
source(paste0(path_code_att, "newUsed_functions.R"))

# Heteroskedastic ---------------------------------------------------------

# from OLS resdiual 
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365

all.cases.name = names(dat)
all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
unique.ind = match(unique(all.cases.ind), all.cases.ind )
gps.era.dat = dat[unique.ind]
data.test = gps.era.dat 
n = length(data.test)
sd.all <- rep(NA, n)
# offset.all <- rep(NA, n)
trend.all = data.frame(matrix(NA, ncol = 2, nrow = n))
for (k in c(1:n)) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")
  
  # Contruction of the dataset 
  Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    # mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>% 
    mutate(complete.time=1:(2*one.year*win.thres)) %>% 
    mutate(Xt=complete.time-one.year*win.thres/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  Data.bef = Data.mod[c(1:(one.year*win.thres)),]
  Data.aft = Data.mod[-c(1:(one.year*win.thres)),]
  trend.bef = lm(signal ~ ., data =  Data.bef)
  trend.aft = lm(signal ~ ., data =  Data.aft)
  trend.all[k,] <- c(trend.bef$coefficients[2], trend.aft$coefficients[2])
  # res.hac.1step <- Test_OLS_vcovhac_1step(Data.mod)
  # offset.all[k] = res.hac.1step$fit.hac$Estimate[2]
}

# from the moving window

sd.all= get(load( file = paste0(path_results,"attribution/sd.all_",  win.thres,"years_", nearby_ver,"screened.RData")))
# filter only gps-era unique
sd.gpsera= sd.all[unique.ind]
mean.bef = sapply(c(1:length(sd.gpsera)), function(x) mean(unlist(sd.gpsera[[x]]$bef[name.series])^2, na.rm = TRUE))
delta.bef = sapply(c(1:length(sd.gpsera)), function(x){
  a = unlist(sd.gpsera[[x]]$bef[name.series])^2
  (max(a, na.rm = TRUE) -min(a, na.rm = TRUE))/2
}) 
mean.aft = sapply(c(1:length(sd.gpsera)), function(x) mean(unlist(sd.gpsera[[x]]$aft[name.series])^2, na.rm = TRUE))
delta.aft = sapply(c(1:length(sd.gpsera)), function(x) {
  a = unlist(sd.gpsera[[x]]$aft[name.series])^2
  (max(a, na.rm = TRUE) -min(a, na.rm = TRUE))/2
})
all.var = data.frame(var.m.bef = mean.bef, var.m.aft = mean.aft, 
                     var.d.bef = delta.bef, var.d.aft = delta.aft, 
                     r.var.d.bef = delta.bef/mean.bef, r.var.d.aft = delta.aft/mean.aft)
summary(all.var)
# fit the variance 

# autocorrelation  --------------------------------------------------------


# relative trend ----------------------------------------------------------

library(tidyverse)   
library(attempt)
library(nlme)

source(paste0(path_code_att, "newUsed_functions.R"))
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365

all.cases.name = names(dat)
all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
unique.ind = match(unique(all.cases.ind), all.cases.ind )
gps.era.dat = dat[unique.ind]

data.test = gps.era.dat
list.ind = c(1:length(data.test))
tot.res <- data.frame(matrix(NA, ncol = 11, nrow = length(list.ind)))
n = length(data.test)
len = data.frame(matrix(NA, ncol = 2, nrow = n))
Res <- list()
trend.all = data.frame(matrix(NA, ncol = 2, nrow = n))
p.all = data.frame(matrix(NA, ncol = 2, nrow = n))

for (k in list.ind) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")
  
  Bef = Y.with.NA[c(1:(one.year*win.thres)),]
  Aft = Y.with.NA[c((one.year*win.thres+1) :(one.year*win.thres*2)),]
  # Contruction of the dataset before breaks
  Data.bef <- Bef %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(complete.time=1:(one.year*win.thres)) %>% 
    mutate(Xt=complete.time-one.year*win.thres/2) %>%
    dplyr::select(-date)
  # Contruction of the dataset after breaks
  Data.aft <- Aft %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(complete.time=1:(one.year*win.thres)) %>% 
    mutate(Xt=complete.time-one.year*win.thres/2) %>%
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.bef <- Data.bef %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
    eval(parse(text=paste0("Data.aft <- Data.aft %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.bef <- Data.bef %>% dplyr::select(-complete.time)
  Data.aft <- Data.aft %>% dplyr::select(-complete.time)
  
  # test HAC sith significant vars 
  bef.test = Test_OLS_vcovhac1(Data.bef)
  aft.test = Test_OLS_vcovhac1(Data.aft)
  # remove insignificant
  print(k)
  trend.all[k,] <- c(bef.test$fit.ols$coefficients[2], aft.test$fit.ols$coefficients[2])
  p.all[k,] <- c(bef.test$fit.hac$`Pr(>|t|)`[2], aft.test$fit.hac$`Pr(>|t|)`[2])
  len[k,] = c( nrow(na.omit(Data.bef)), nrow(na.omit(Data.aft)))
}
colnames(tot.res) = names(Data.mod)
trend.abs = tot.res$Xt
# compute the relative trend: t(gps-era)/mean(gps)
four.series = get(load(file = paste0(path_results,"attribution/four_main_series_",win.thres,"year_", nearby_ver,".RData")))
all.names = names(four.series)
mean.gps = rep(NA, length(list.ind))
for (i in list.ind) {
  mean.gps[i] = mean(four.series[[names(data.test)[i]]]$GPS, na.rm = TRUE)
}

trend.rel = trend.all/mean.gps

# sort only the case have length >1000 

bef.t = trend.all$X1[which(len$X1>1000)]
aft.t = trend.all$X2[which(len$X2>1000)]

a = c(bef.t, aft.t)
b = c(len$X1[which(len$X1>1000)], len$X2[which(len$X2>1000)])
hist(a, main = "Histogram of the trend of the homogenous segment when n >1000", breaks=100)
plot(b,a, ylab = "Absolute trend", xlab = "Length")
p.all

bef.p = p.all$X1[which(len$X1>1000)]
aft.p = p.all$X2[which(len$X2>1000)]
a = c(bef.p, aft.p)
b = c(len$X1[which(len$X1>1000)], len$X2[which(len$X2>1000)])
hist(a, main = "Histogram of the p value of trend when n >1000", xaxt='n',breaks = seq(0, 1, 0.05), xlim = c(0,1))
axis(side = 1, at=seq(0,1, 0.05), labels=seq(0,1, 0.05))

which(side = 1, len$X1>1000 & p.all$X1 <0.05)

bef.p = trend.all$X1[which(len$X1>1000 & p.all$X1 < 0.05)]
aft.p = trend.all$X2[which(len$X2>1000 &  p.all$X2 < 0.05)]
a = c(bef.p, aft.p)


