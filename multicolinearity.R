# this prog is used to investigate the multicolinearity problem 
# inspect by the autocorrelation between regressors

# inspect the variance inflation index 
# this prog used to apply the ancova/fgls on the real data 
library(tidyverse)   
library(attempt)
library(nlme)

source(paste0(path_code_att, "newUsed_functions.R"))
dat  = get(load(file = paste0(path_results,"attribution/data.all_1year_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365

all.cases.name = names(dat)
all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
unique.ind = match(unique(all.cases.ind), all.cases.ind )
gps.era.dat = dat[unique.ind]

data.test = gps.era.dat
list.ind = c(1:length(data.test))
tot.res <- data.frame(matrix(NA, ncol = 3, nrow = length(list.ind)))

for (k in list.ind) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")

  # Contruction of the dataset 
  Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(Jump=c(rep(0,one.year),rep(1,one.year))) %>% 
    mutate(complete.time=1:(2*one.year)) %>% 
    mutate(Xt=complete.time-one.year/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  res.hac.1step <- Test_OLS_vcovhac_1step(Data.mod)
  res.hac <- Test_OLS_vcovhac(Data.mod)
  tot.res[k,2] = car::vif(res.hac.1step$fit.ols)[1]
  ind.jump = which(rownames(res.hac$fit.hac) == "Jump")
  if(length(rownames(res.hac$fit.hac)) >2){
    if( length(ind.jump) !=0){
      tot.res[k,3] <- car::vif(res.hac$fit.ols)[ind.jump]
    }else{tot.res[k,3] =  NA }
  }else{tot.res[k,3] =  -1 }
  
  # res.hac <- Test_OLS_vcovhac(Data.mod)
  ########################
  # compute the correlation between regressors: jump and trend 
  tot.res[k,1] = res.hac.1step$vcov.para[3,2]/(sqrt(res.hac.1step$vcov.para[2,2]*res.hac.1step$vcov.para[3,3]))
  # tot.res[[name.dataset]] <- res.hac
  # tot.res[[name.dataset]] <- list(hac = round(res.fgls$res.gls, digits = 5),
  #                                 predicted = res.fgls$predicted)
  #                                 fgls = round(res.fgls$res.gls, digits = 5))
  print(k)
}

colnames(tot.res) <- c("corr", "VIF.full", "VIF.selected")
hist(tot.res[,1] , breaks =50, main = "Histogram of the covariance", xlab = "")
hist(tot.res[,2] , breaks =200, main = "Histogram of the VIF with all variables", xlab = "", xlim = c(0,50))
hist(tot.res[,3] , breaks =50, main = "Histogram of the VIF after variable selection", xlab = "", xlim = c(0,50))
save(tot.res, file = paste0(path_results, "attribution/multicolinear.2years.RData"))


