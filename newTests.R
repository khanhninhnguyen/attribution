rm(list=ls())
library(tidyverse)   
library(padr)
library(attempt)
library(nlme)


path.main <- "/Users/lebarbier/Desktop/Boulot/Theses&Stages/Theses/NinhNguyen/Test_CP_corr_hetero/"
setwd(path.main)
source("Used_functions.R")

name.dataset <-"data.all_1years_NGL.auck.2005-11-07.whng.RData"
name.series <- "gps.era"

            ####################################
            #Importation of the data, add the NA 
load(name.dataset)
Y.with.NA <- pad(Y)


            ######
            # Construction of the dataset for the regression on the series name.series

# We test if the break is in the middle of the data. 
# If it is not true, we stop
date.detected.break=str_extract_all(name.dataset, pattern = str_c("([0-9]{4})[- \\.]",  "([0-9]{2})[- \\.]","([0-9]{2})"))[[1]]
one.year=365
stop_if(which(Y.with.NA$date==date.detected.break), ~.x!=one.year)

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
Data.mod1 <- Data.mod1 %>% dplyr::select(-complete.time)
  
                  ########################
                  # Test on the OLS estimate with the HAC covariance and selection of the 
                  # significant parameters
      
res.hac <- Test_OLS_vcovhac(Data.mod)
res.hac$fit.hac
                    
# plot of the solution
date=Y.with.NA$date
plot.predict(res.hac,Data.mod,one.year,date)


                  ####################################
                  # Test on the FGLS estimate with a ARMA(1,1) model

res.fgls <- Test_FGLS(Data.mod)
res.fgls$res.gls

# plot of the solution
plot.predict(res.fgls,Data.mod,one.year,date)
