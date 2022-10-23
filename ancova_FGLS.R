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


tot.res <- list()
for (i in c(1:length(gps.era.dat))) {
  name.dataset = names(gps.era.dat)[i]
  Y.with.NA = gps.era.dat[[i]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")
  # if(which(Y.with.NA$date==date.detected.break)!=one.year) next
  
  # Contruction of the dataset 
  Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(Jump=c(rep("Left",one.year),rep("Right",one.year))) %>% 
    mutate(complete.time=1:(2*one.year)) %>% 
    mutate(Xt=complete.time-one.year/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  ########################
  # Test on the OLS estimate with the HAC covariance and selection of the 
  # significant parameters
  res.hac <- Test_OLS_vcovhac(Data.mod)
  ####################################
  # Test on the FGLS estimate with a ARMA(1,1) model
  res.fgls <- Test_FGLS(Data.mod)

  tot.res[[name.dataset]] <- list(hac = round(res.hac$fit.hac, digits = 5), 
                                  fgls = round(res.fgls$res.gls, digits = 5))
}

save(tot.res, file = paste0(path_results, "GPS.ERA.mean_test.R"))







