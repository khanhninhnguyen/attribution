# this function is ancova test for the change in mean 

ancova.test <- function(Y){
  one.year=365
  time.series <- seq(1-one.year/2,one.year/2,1) 
  rg.break = 365
  signal.without.NA = Y$gps.era[which(is.na(Y$gps.era) == FALSE)]
  
  Data.mod <- Y %>% dplyr::select(gps.era,date) %>%
    mutate(signal=gps.era) %>% dplyr::select(-gps.era) %>% 
    slice((rg.break-one.year+1):(rg.break+one.year)) %>% 
    mutate(Group=c(rep("Left",one.year),rep("Right",one.year))) %>% 
    mutate(Xt=rep(time.series,2)) %>% dplyr::select(-date)
  
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=rep(cos(i*time.series*(2*pi)/one.year),2),sin",i,"=rep(sin(i*time.series*(2*pi)/one.year),2))")))
  }
  
  # Test on the OLS estimate with the HAC covariance and selection of the 
  # significant parameters
  ####################################
  
  # OLS estimates
  fit.ols <- lm(signal~ Group+Xt+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=Data.mod)
  X <- model.matrix(fit.ols)
  Y <- Data.mod$signal
  list.coeff.names <- c("Intercept","jump",names(fit.ols$coefficients)[3:dim(X)[2]])
  X.names <- list.coeff.names
 
  # HAC
  vcov.para=sandwich::NeweyWest(fit.ols, lag = 1)
  
  # Test with vcov (HAC)
  fit.signal.hac=lmtest::coeftest(fit.ols,df=Inf,vcov.=vcov.para)[, ] %>% as.data.frame()
  
  # Selection parameter by parameter
  fit.signal.hac.without.intercept=fit.signal.hac[2:dim(X)[2],]
  
  return(fit.signal.hac.without.intercept$`Pr(>|z|)`)
  pval.test.without.intercept <-fit.signal.hac.without.intercept$`Pr(>|z|)`
  max.pval.test <- max(pval.test.without.intercept)
  threshold <- 0.01
  while (max.pval.test>threshold){
    # which is the most unsignificant
    rg.max.pval.test <- which(max.pval.test==pval.test.without.intercept)
    names.max <- row.names(fit.signal.hac.without.intercept)[rg.max.pval.test]
    names.max

    X <- X[,!(X.names %in% names.max)]
    X.names <- X.names[!(X.names %in% names.max)]

    fit.select <- c()
    fit.select <- lm(signal.without.NA~-1+X)
    X <- model.matrix(fit.select)
    colnames(X) <- X.names
  #
  #
    vcov.para.select=sandwich::NeweyWest(fit.select, lag = 1)
    fit.signal.hac=lmtest::coeftest(fit.select,df=Inf,vcov.=vcov.para.select)[, ] %>% as.data.frame()
    row.names(fit.signal.hac) <- X.names
    fit.signal.hac.without.intercept=fit.signal.hac[2:dim(X)[2],]
    pval.test.without.intercept <-fit.signal.hac.without.intercept$`Pr(>|z|)`
    max.pval.test <- max(pval.test.without.intercept)
  }

  # return(fit.signal.hac)

  # plot of the solution
  # signal.predicted <- rep(NA,length(Data.mod$signal))
  # signal.predicted[!((1:length(Data.mod$signal)) %in% na.Y)] <- fit.select$fitted.values
  #
  # plot(Y.with.NA$date,Data.mod$signal,type="n",xlab="dates",ylab="gps-era")
  #
  # Left <- 1:one.year
  # lines(Y.with.NA$date[Left],Data.mod$signal[Left],col="red", type="l",xlab="time",ylab="signal")
  # lines(Y.with.NA$date[Left],signal.predicted[Left],col="red")
  #
  # Right <- (one.year+1):(2*one.year)
  # lines(Y.with.NA$date[Right],Data.mod$signal[Right],col="blue")
  # lines(Y.with.NA$date[Right],signal.predicted[Right],col="blue")

}
re = rep(NA, length(dat))
for (i in c(1:length(dat))) {
  a = ancova.test(dat[[i]])
  re[i] <- a[1]
}





rm(list=ls())
setwd("/Users/lebarbier/Desktop/Boulot/Theses&Stages/Theses/NinhNguyen/Test_CP_corr_hetero/")
library(tidyverse)   
library(padr)

#######
# Data importation
load("data.all_1years_NGL.auck.2005-11-07.whng.RData")
#load("data.all_1years_NGL.gope.2009-05-08.zdib.RData")
head(Y)
dim(Y)
Y.without.NA <- Y
signal.without.NA <-Y.without.NA$gps.era   #?

# Adding the missing data corresponding to missing date
Y.with.NA <- pad(Y)
na.Y <- which(is.na(Y.with.NA$gps.era))
n.Y <- length(Y)

######
# Model construction
date.detected.break="2005-11-07"
#date.detected.break="2009-05-08"

rg.break <- which(Y.with.NA$date==date.detected.break)
one.year=365
time.series <- seq(1-one.year/2,one.year/2,1) 

Data.mod <- Y.with.NA %>% dplyr::select(gps.era,date) %>%
  mutate(signal=gps.era) %>% dplyr::select(-gps.era) %>% 
  slice((rg.break-one.year+1):(rg.break+one.year)) %>% 
  mutate(Group=c(rep("Left",one.year),rep("Right",one.year))) %>% 
  mutate(Xt=rep(time.series,2)) %>% dplyr::select(-date)

for (i in 1:4){
  eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=rep(cos(i*time.series*(2*pi)/one.year),2),sin",i,"=rep(sin(i*time.series*(2*pi)/one.year),2))")))
}
#head(Data.mod)

#######
# Test on the OLS estimate with the HAC covariance and selection of the 
# significant parameters
####################################

# OLS estimates
fit.ols <- lm(signal~ Group+Xt+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=Data.mod)
X <- model.matrix(fit.ols)
Y <- Data.mod$signal
list.coeff.names <- c("Intercept","jump",names(fit.ols$coefficients)[3:dim(X)[2]])
X.names <- list.coeff.names

# HAC
vcov.para=sandwich::NeweyWest(fit.ols, lag = 1)

# Test with vcov (HAC)
fit.signal.hac=lmtest::coeftest(fit.ols,df=Inf,vcov.=vcov.para)[, ] %>% as.data.frame()

# Selection parameter by parameter
fit.signal.hac.without.intercept=fit.signal.hac[2:dim(X)[2],]

pval.test.without.intercept <-fit.signal.hac.without.intercept$`Pr(>|z|)`
max.pval.test <- max(pval.test.without.intercept)
threshold <- 0.01
while (max.pval.test>threshold){
  # which is the most unsignificant
  rg.max.pval.test <- which(max.pval.test==pval.test.without.intercept)
  names.max <- row.names(fit.signal.hac.without.intercept)[rg.max.pval.test]
  names.max
  
  X <- X[,!(X.names %in% names.max)]
  X.names <- X.names[!(X.names %in% names.max)]
  
  fit.select <- c()
  fit.select <- lm(signal.without.NA~-1+X)
  X <- model.matrix(fit.select)
  colnames(X) <- X.names
  
  
  vcov.para.select=sandwich::NeweyWest(fit.select, lag = 1)
  fit.signal.hac=lmtest::coeftest(fit.select,df=Inf,vcov.=vcov.para.select)[, ] %>% as.data.frame()
  row.names(fit.signal.hac) <- X.names
  fit.signal.hac.without.intercept=fit.signal.hac[2:dim(X)[2],]
  pval.test.without.intercept <-fit.signal.hac.without.intercept$`Pr(>|z|)`
  max.pval.test <- max(pval.test.without.intercept)
}

fit.signal.hac


# plot of the solution
signal.predicted <- rep(NA,length(Data.mod$signal))
signal.predicted[!((1:length(Data.mod$signal)) %in% na.Y)] <- fit.select$fitted.values

plot(Y.with.NA$date,Data.mod$signal,type="n",xlab="dates",ylab="gps-era")

Left <- 1:one.year
lines(Y.with.NA$date[Left],Data.mod$signal[Left],col="red", type="l",xlab="time",ylab="signal")
lines(Y.with.NA$date[Left],signal.predicted[Left],col="red")

Right <- (one.year+1):(2*one.year)
lines(Y.with.NA$date[Right],Data.mod$signal[Right],col="blue")
lines(Y.with.NA$date[Right],signal.predicted[Right],col="blue")


# 
# 
# # OLS estimates
# fit.ols <- lm(signal~ -1+Group+Xt+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=Data.mod)
# X <- model.matrix(fit.ols)
# Y <- Data.mod$signal
# list.coeff.names <- names(fit.ols$coefficients)
# 
# #summary(fit.signal)
# 
# # HAC
# vcov.para=sandwich::NeweyWest(fit.ols, lag = 1)
# 
# # Test
# 
# 
# #fit.signalb <- lm(signal~Group+Xt+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=Data.mod)
# #fit.signal.hac=lmtest::coeftest(fit.signalb,df=Inf,vcov.=vcov.para)[, ]
# coeff.esti <- fit.ols$coefficients
# n.coeff <- length(coeff.esti)
# list.coeff.names.test <- c("jump",names(coeff.esti)[3:n.coeff])
# 
# # Construction of the matrix C for each test
# C.matrix <- c()
# C.diag <- diag(1,n.coeff-2,n.coeff-2)
# C.matrix <- cbind(matrix(0,nrow=n.coeff-2,ncol=2),C.diag)  
# C.matrix <- rbind(c(-1,1,rep(0,n.coeff-2)),C.matrix)
# names.col.C.matrix <- list.coeff.names
# names.row.C.matrix <- list.coeff.names.test
# 
# 
# test.hand <- function(coeff.esti,C,vcov.para){
#   stat.test <- (t(C)%*%coeff.esti)/sqrt(t(C)%*%vcov.para%*% C)
#   pval.test <- 2*pnorm(abs(stat.test),lower.tail=FALSE)
#   return(list(stat.test=stat.test,pval.test=pval.test))
# }
# 
# 
# res.test=purrr::map(1:length(list.coeff.names.test),~test.hand(coeff.esti,C.matrix[.x,],vcov.para)) %>% bind_rows() %>% mutate(coeff.esti=C.matrix%*%coeff.esti)%>% as.data.frame()
# row.names(res.test)=list.coeff.names.test
# res.test <- res.test %>% relocate(coeff.esti,stat.test,pval.test)
# 
# max.pval.test <- max(res.test$pval.test)
# threshold <- 0.01
# while (max.pval.test>threshold){
#   # which is the most unsignificant
#   rg.max.pval.test <- which(max.pval.test==res.test$pval.test)
#   names.max <- row.names(res.test)[rg.max.pval.test]
#   
#   names.max
#   # New matrix X, coeff.esti, C.matrix ...
#   if (names.max=="jump"){
#     names.max.both=c("GroupLeft","GroupRight")
#     X <- X[,!(list.coeff.names %in% names.max.both)]
#     
#     fit.select <- c()
#     fit.select <- lm(signal.without.NA~-1+X)
#     vcov.para.select=sandwich::NeweyWest(fit.select, lag = 1)
#     
#     coeff.esti.select <- fit.select$coefficients
#     names(coeff.esti.select) <- colnames(X)
#     
#     list.coeff.names <- list.coeff.names[!(list.coeff.names %in% names.max.both)]
#     list.coeff.names.test <- list.coeff.names.test[!(list.coeff.names.test %in% names.max)]
#     
#     C.matrix <- C.matrix[,!(names.col.C.matrix %in% names.max.both)]
#     C.matrix <- C.matrix[!(names.row.C.matrix %in% names.max),]
#     names.col.C.matrix <- names.col.C.matrix[!(names.col.C.matrix %in% names.max.both)]
#     names.row.C.matrix <- names.row.C.matrix[!(names.row.C.matrix %in% names.max)]
#   } else {
#     
#   X <- X[,!(list.coeff.names %in% names.max)]
#   
#   fit.select <- c()
#   fit.select <- lm(signal.without.NA~-1+X)
#   vcov.para.select=sandwich::NeweyWest(fit.select, lag = 1)
#   
#   coeff.esti.select <- fit.select$coefficients
#   names(coeff.esti.select) <- colnames(X)
#     
#   list.coeff.names <- list.coeff.names[!(list.coeff.names %in% names.max)]
#   list.coeff.names.test <- list.coeff.names.test[!(list.coeff.names.test %in% names.max)]
#   
#   C.matrix <- C.matrix[,!(names.col.C.matrix %in% names.max)]
#   C.matrix <- C.matrix[!(names.row.C.matrix %in% names.max),]
#   names.col.C.matrix <- names.col.C.matrix[!(names.col.C.matrix %in% names.max)]
#   names.row.C.matrix <- names.row.C.matrix[!(names.row.C.matrix %in% names.max)]
#   }
#   
#   res.test=purrr::map(1:length(list.coeff.names.test),~test.hand(coeff.esti.select,C.matrix[.x,],vcov.para.select)) %>% 
#     bind_rows() %>% 
#     mutate(coeff.esti.select=C.matrix%*%coeff.esti.select) %>% 
#     as.data.frame()
#   
#   row.names(res.test)=list.coeff.names.test
#   res.test <- res.test %>% relocate(coeff.esti.select,stat.test,pval.test)
#   max.pval.test <- max(res.test$pval.test)
#   
# }
# 
# res.test


