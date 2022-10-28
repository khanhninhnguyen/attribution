rm(list=ls())
library(tidyverse)   
# library(padr)
library(attempt)
library(nlme)
library(car)


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
  mutate(Jump=c(rep(0,one.year),rep(1,one.year))) %>% 
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
res.hac$fit.hac
fit.hac

# correlation between the beta_hat
cor.beta.hac=cov2cor(res.hac$vcov.para)
corrplot(cor.beta.hac)

                  ########################
                  # Test on the OLS estimate with the HAC covariance and selection of the 
                  # significant parameters

# But the estimators are biased. The right solution is the Beta_GLS with the estimated covariance matrix of Y

#X <- model.matrix(res.hac$fit.ols)
#mat.n <- X%*%t(X)
#mat.p <- t(X)%*%X
#Omega=(solve(mat.n))%*%X%*%mat.p%*%res.hac$vcov.para%*%mat.p%*%t(X)%*%solve(mat.n)

# Indice of colinearity
vif(res.hac$fit.ols)
#The rule of thumb is that if VIF_i>10 then multicollinearity is high[5] (a cutoff of 5 is also commonly used[6]) Wikipedia
#If the variance inflation factor of a predictor variable were 5.27 (√5.27 = 2.3), this means that the standard error for the coefficient of that predictor variable is 2.3 times larger than if that predictor variable had 0 correlation with the other predictor variables.

                    
# plot of the solution
date=Y.with.NA$date
predicted=res.hac$fit.ols$fitted.values
plot.predict(predicted,Data.mod,one.year,date)

# 
#                   ####################################
#                   # Test on the FGLS estimate with a ARMA(1,1) model
# 
# res.fgls <- Test_FGLS(Data.mod)
# res.fgls$res.gls
# 
# # plot of the solution
# plot.predict(res.fgls,Data.mod,one.year,date)

a = dat1$`alic.2011-06-08.20na`[c(3286:4015),]

win.thres = 1
# Contruction of the dataset 
Data.mod <- a %>% dplyr::select(name.series,date) %>%
  rename(signal=name.series) %>% 
  mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>% 
  mutate(complete.time=1:(2*one.year*win.thres)) %>% 
  mutate(Xt=complete.time-one.year*win.thres/2) %>% 
  dplyr::select(-date)
for (i in 1:4){
  eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
}
Data.mod <- Data.mod %>% dplyr::select(-complete.time)
res.hac <- Test_OLS_vcovhac(Data.mod)
res.hac1 <- Test_OLS_vcovhac_1step(Data.mod)


