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
offset.all <- rep(NA, n)
for (k in c(1:n)) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")
  
  # Contruction of the dataset 
  Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>% 
    mutate(complete.time=1:(2*one.year*win.thres)) %>% 
    mutate(Xt=complete.time-one.year*win.thres/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  res.hac.1step <- Test_OLS_vcovhac_1step(Data.mod)
  offset.all[k] = res.hac.1step$fit.hac$Estimate[2]
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
all.var = data.frame(mean.bef, mean.aft, delta.bef/mean.bef, delta.aft/mean.aft)
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
Res <- list()
for (k in list.ind) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")
  
  # Contruction of the dataset 
  Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>% 
    mutate(complete.time=1:(2*one.year*win.thres)) %>% 
    mutate(Xt=complete.time-one.year*win.thres/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod = Data.mod %>% dplyr::select(-complete.time)
  res.hac.1step = Test_OLS_vcovhac_1step(Data.mod)
  tot.res[k,] = res.hac.1step$fit.ols$coefficients
  print(k)
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

a = trend.abs/mean.gps

a = lm(signal~Jump, data = Data.mod)








