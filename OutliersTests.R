rm(list=ls())
setwd("/Users/lebarbier/Desktop/Boulot/Theses&Stages/Theses/NinhNguyen/MixtureModelForOutliers/")
source("UsedFunctions.R")
library(tidyverse)   
          ########################@
          #### Simulation

# Parameters
n             = 500   #length of the series
prob.outliers = 0.01
size.outliers = 5 
P.true             = 2  # P can be 1, 2 or 3

# Simulated series
SimData = SimulatedSeries(n,P.true,prob.outliers,size.outliers)
Y=SimData$Y
cluster.true=SimData$cluster.true
# Plot
plot(Y,col=cluster.true,pch=16,cex=0.8)
hist(Y,breaks = 20)

Z=(Y-mean(Y))/sd(Y)
hist(Z,breaks = 20,freq=FALSE)
zseq=seq(min(Z),max(Z),0.01)
fnorm=dnorm(zseq,0,1)
lines(zseq,fnorm,col=p)

Zb=(Y-median(Y))/mad(Y)
hist(Zb,breaks = 20,freq=FALSE)
zseq=seq(min(Zb),max(Zb),0.01)
fnorm=dnorm(zseq,0,1)
lines(zseq,fnorm,col=p)


              ######
              #### Real data
load("auck.2005-11-07.whng.RData")
dates=Y$date
Y=Y$gps.gps
n=length(Y)

rg.na=which(is.na(Y))
Y=Y[-rg.na]
dates=dates[-rg.na]
n=length(Y)
# Plot
plot(dates[1:200],Y[1:200],pch=16,cex=0.8)
hist(Y,breaks = 30)

x=Y
out <- LocScaleB(x = x,  k = 3, method='MAD')
out$pars
out$bounds
out$outliers

plot(dates,Y,pch=16,cex=0.8)
points(dates[out$outliers],Y[out$outliers],col="red")

            ########
###############@
##### TESTS on real data

load("data.all_1years_NGL.gope.2009-05-08.zdib.paired.RData")
#data.all_1years_NGL.auck.2005-11-07.whng.paired.RData
head(Y)
y=Y$gps.era
dates=Y$date
n=length(y)
times=1:n

plot(dates,y,pch=16,cex=0.8)


residus=y-median(y)
plot(dates,residus^2)
lissage=loess(residus^2~times,degree=2) # revoir pour une meilleure estimation
lines(dates,lissage$fitted,col="red")

z=residus/sqrt(lissage$fitted)
plot(dates,z)
hist(z,breaks = 20)

prob.z=dnorm(z,0,1)
f=2
OS=log(prob.z)^(2*f)
plot(dates,OS)
ScOS=(OS-min(OS))/(max(OS)-min(OS))*10
plot(dates,ScOS)
thresh=7
abline(h=thresh,col="red")
rg=which(ScOS>thres)

plot(dates,z,pch=16,cex=0.8)
points(dates[rg],z[rg],col="red")

plot(dates,y,pch=16,cex=0.8)
points(dates[rg],y[rg],col="red")

#mean(y,na.rm=TRUE)
#y[c(50,60,80)]=NA
reg=lm(z~1)
### OLS
print(coef(summary(reg)))
### Test with variance correction
library(sandwich)
V_HAC <- sandwich::NeweyWest(reg, lag = 1)
print(lmtest::coeftest(reg, sandwich::NeweyWest(reg, lag = 1))[, ])

            ########
            ### Tests 
library(outliers)
subset <- Y[1:30]
test <- dixon.test(subset)
test

library(EnvStats)
k=prob.outliers*n
test <- rosnerTest(Y,k=k)
test
num.outliers <- test$all.stats$Obs.Num[test$all.stats$Outlier==TRUE]
plot(Y,pch=16,cex=0.8)
points(num.outliers,Y[num.outliers],col="red")

install.packages("~/Downloads/OutlierDetection_0.1.1.tar.gz", repos = NULL, type = "source",dependencies = TRUE)
library(OutlierDetection)
test <- dens(Y, k = 0.05 * nrow(Y))


test <- OutlierDetection(x, k = 0.05*nrow(Y))
test

library(ldbod)
scores <- ldbod(Y, k = c(10,20), nsub = 0.10*length(Y), method = c('lof','lpdf'))
plot(Y)
top5outliers <- Y[order(scores$lof[,2],decreasing=TRUE)[1:5]]  ###?
points(top5outliers,col=2)   ##?

library(fdaoutlier)
msplot(dts = as.matrix(Y))  #????


library(outliertree)
### example dataset with interesting outliers
data(hypothyroid)
### fit the model and get a print of outliers
model <- outlier.tree(hypothyroid,
                      outliers_print=10,
                      save_outliers=TRUE,
                      nthreads=1)


car::outlierTest()


library(extremevalues)
test <- getOutliers(Y, method="I")
testI <- getOutliersI(Y, rho=c(1,1), FLim=c(0.1,0.9), distribution="normal")
testII <- getOutliersII(Y, alpha=c(0.05, 0.05), FLim=c(0.1, 0.9),
              distribution="normal", returnResiduals=TRUE)
outlierPlot(Y,test,mode="qq")
outlierPlot(Y,testII,mode="residual")


# #######
# Simulation without correlation
n=300
delta=5
s2=c(rep(0.1,100),rep(0.3,100),rep(0.1,100))
times=1:n

y=delta+rnorm(n,0,sqrt(s2)) 
data_mat <- data.frame(y,times)
plot(y)

residus=y-mean(y)
plot(residus^2)
lissage=loess(residus^2~times,degree=2)    # revoir pour une meilleure estimation
lines(lissage$fitted,col="red")

res=lm(y~1,data=data_mat,weights = lissage$fitted)
summary(res)

#####
# Avec correlation
delta=1
rho=0.7
y1=arima.sim(n=100,list(ar=rho),sd=0.1)
y2=arima.sim(n=100,list(ar=rho),sd=0.4)
y3=arima.sim(n=100,list(ar=rho),sd=0.1)
y=c(y1,y2,y3)+delta
plot(y)


reg=lm(y~1)
### OLS
print(coef(summary(reg)))
### Test with variance correction
library(sandwich)
V_HAC <- sandwich::NeweyWest(reg, lag = 1)
print(lmtest::coeftest(reg, sandwich::NeweyWest(reg, lag = 1))[, ])


# FGLS with constant variance
library(orcutt)
reg_CORC <- orcutt::cochrane.orcutt(reg)
print(coef(summary(reg_CORC)))
print(reg_CORC$rho)



### # Estimation of the heteroscedastic variance



#data.test=data.frame(y=y,times=times)
#
#rho.hat=0.7

result_acf=acf(y-mean(y)) # If the noise is heteroscdestic??
rho.hat=result_acf$acf[,,1][2]
rho.hat

# 1. Whitening
n.y=length(y)
y.prime=y[2:n.y]-rho.hat*y[1:(n.y-1)]
n=length(y.prime)
times=1:n
data.test=data.frame(y.prime=y.prime,times=times)

residus=y.prime-mean(y.prime)

# 2. estimation of sigma^2_t as a function of t using a smoothing
#library(gam)
lissage=loess(residus^2~times)    # revoir pour une meilleure estimation
plot(residus^2)
lines(lissage$fitted,col="red")
lissage2=lowess(residus^2~times)
lines(lissage2,col="blue")
esti.s2=lissage$fitted

#?lowess

# nulity test for delta

delta.prime.hat=(sum(y.prime/esti.s2))/(sum(1/esti.s2))

#aa=c(rep(0.1,49),rep(0.5,50),rep(0.1,50))
#delta.prime.hat=(sum(y.prime/aa))/(sum(1/aa))

delta.prime.hat
#delta*(1-rho.hat)

T.stat=delta.prime.hat/sqrt(sum(1/esti.s2))
T.stat
p.val.test=2*pnorm(-abs(T.stat))
p.val.test




library(mvoutlier)
Z=data.frame(x=1:length(y),y=y)
a=dd.plot(Z)
outliers.y=which(a$outliers==TRUE)
plot(y)
x=1:length(y)
points(x[outliers.y],y[outliers.y],col="red")

library(psych)
d2=outlier(Z,plot=TRUE)
plot(Z$x,Z$y,bg=c("yellow","blue")[(d2 > 15)+1],pch=21)

library(Routliers)
res=outliers_mad(y)
plot_outliers_mad(res,x=x)

# 
r=c()
n=length(y)
for (i in 1:n){
  mu.i=mean(y[-i])
  r[i]=y[i]-mu.i
}
r=r/sd(r)  #RES residus externes standadisés
qt(0.975,n)
plot(abs(r))
abline(h=qt(0.975,n))
rg=which(abs(r)>qt(0.975,n))

plot(y)
points(x[rg],y[rg],col="red")
