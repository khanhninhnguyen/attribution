rm(list=ls())
setwd("/Users/lebarbier/Desktop/Boulot/Theses&Stages/Theses/NinhNguyen/MixtureModelForOutliers/")
source("UsedFunctions.R")
library(tidyverse) 



########################@
#### Simulation

# Parameters
n             = 500   #length of the series
prob.outliers = 0.1
size.outliers = 5 
P.true        = 3  # P can be 1, 2 or 3

# Simulated series
SimData = SimulatedSeries(n,P.true,prob.outliers,size.outliers)
Y=SimData$Y
cluster.true=SimData$cluster.true
# Plot
plot(Y,col=cluster.true,pch=16,cex=0.8)
hist(Y,breaks = 20)

######
#### Real data
load(paste0(path_results,"attribution/auck.2005-11-07.whng.RData"))
Y = Y[which(Y$date > as.Date("2005-11-07") & Y$date < as.Date("2006-11-07")),]
dates=Y$date
Y=Y$gps.gps
n=length(Y)

rg.na=which(is.na(Y))
Y=Y[-rg.na]
dates=dates[-rg.na]
n=length(Y)
# Plot
plot(dates,Y,pch=16,cex=0.8)

hist(Y,breaks = 30,freq=FALSE)
yseq=seq(min(Y),max(Y),0.01)
fnorm=dnorm(yseq,0,1)
lines(yseq,fnorm,col="red")
ft=dt(yseq,n-1)
lines(yseq,ft,col="blue")

#########
# Gaussian mixture model with same variance

#### 1. N(0,1) and the others variances=1

option.init="CAH"

Pmax=10  # Pmax>1. 
lvincP_imp=c()
PHI_imp=list()
TAU_imp=list()
empty_imp=c()
dv_imp=c()
lvincP_imp[1]=-0.5*n*log(2*pi)- 0.5*sum(Y^2) # log-lik for P=1: the mean is imposed to 0

Pseq=2:Pmax
for (P in Pseq){
  Out.EM.init_imp=EM.init_imp(Y,P,option.init)
  Out.EM_imp =EM.algo_imp(Y,Out.EM.init_imp$phi,P,Out.EM.init_imp$Id.cluster1)
  PHI_imp[[P]]= Out.EM_imp$phi
  TAU_imp[[P]]=Out.EM_imp$tau
  empty_imp[P] = Out.EM_imp$empty
  dv_imp[P] = Out.EM_imp$dv
  lvincP_imp[P]=Out.EM_imp$lvinc
}

pen.bic_imp=(2*((1:Pmax)-1))*log(n)
BIC.c_imp=2*lvincP_imp-pen.bic_imp
plot(2*lvincP_imp)
lines(BIC.c_imp,col="red")


#mu.est=mean(Y)
#s2.est=var(Y)
#a=-0.5*n*log((2*pi))- 0.5*sum(Y^2)
#b=-0.5*n*log((2*pi*s2.est))- 0.5*sum((Y-mu.est)^2)/s2.est
### Choice of P
P.best_imp=which(BIC.c_imp==max(BIC.c_imp))
phi_imp=PHI_imp[[P.best_imp]]
tau_imp=TAU_imp[[P.best_imp]]
P.best_imp

# #         #### EM algorithm for a fixed P
# P=2
# Out.EM.init_imp=EM.init_imp(Y,P,option.init)
# Out.EM_imp =EM.algo_imp(Y,Out.EM.init_imp$phi,P,Out.EM.init_imp$Id.cluster1)
# PHI_imp[[P]]= Out.EM_imp$phi
# TAU_imp[[P]]=Out.EM_imp$tau
# empty_imp[P] = Out.EM_imp$empty
# dv_imp[P] = Out.EM_imp$dv
# lvincP_imp[P]=Out.EM_imp$lvinc

#### Classification

cluster_imp = apply(tau_imp,1,which.max)
#table(cluster.true,cluster)
table(cluster_imp)

#### Plot

plot(Y,col=cluster_imp,pch=16,cex=0.8)

yseq=seq(min(Y),max(Y),0.001)
f.per.cluster=matrix(0,ncol=length(yseq),nrow=P)
for (p in 1:P){
  f.per.cluster[p,]=dnorm(yseq,phi_imp$mu[p],sqrt(phi_imp$s2[p]))
}
hist(Y,breaks = 50,freq=FALSE)
for (p in 1:P){
  lines(yseq,f.per.cluster[p,],col=p)
}


#### 2. gaussians with same variance=1
option.init="CAH"

Pmax=10  # Pmax>1. 
lvincP=c()
PHI=list()
TAU=list()
empty=c()
dv=c()
# for one group
mu=mean(Y)
lvincP[1]=-0.5*n*log(2*pi)- 0.5*sum((Y-mu)^2) # log-lik for P=1: the mean is imposed to 0

#s2=var(Y)*(n-1)/n
#lvincP[1]=-0.5*n*log(2*pi*s2)- 0.5*sum((Y-mu)^2)/s2 # log-lik for P=1: the mean is imposed to 0


Pseq=2:Pmax
for (P in Pseq){
  Out.EM.init=EM.init(Y,P,option.init)
  Out.EM =EM.algo(Y,Out.EM.init$phi,P)
  PHI[[P]]= Out.EM$phi
  TAU[[P]]=Out.EM$tau
  empty[P] = Out.EM$empty
  dv[P] = Out.EM$dv
  lvincP[P]=Out.EM$lvinc
}

pen.bic=(2*(1:Pmax)-1)*log(n)
BIC.c=2*lvincP-pen.bic
plot(2*lvincP)
lines(BIC.c,col="red")

### Choice of P
P.best=which(BIC.c==max(BIC.c))
phi=PHI[[P.best]]
tau=TAU[[P.best]]


# #         #### EM algorithm for a fixed P
# P=3
# Out.EM.init=EM.init(Y,P,option.init="CAH")
# Out.EM =EM.algo(Y,Out.EM.init$phi,P,Out.EM.init$Id.cluster1)
# phi = Out.EM$phi
# tau=Out.EM$tau
# lvinc=Out.EM$lvinc
# empty = Out.EM$empty
# dv = Out.EM$dv


phi
#### Classification

cluster = apply(tau,1,which.max)
#table(cluster.true,cluster)
table(cluster)

#### Plot

plot(Y,col=cluster,pch=16,cex=0.8)

yseq=seq(min(Y),max(Y),0.001)
f.per.cluster=matrix(0,ncol=length(yseq),nrow=P)
for (p in 1:P){
  f.per.cluster[p,]=dnorm(yseq,phi$mu[p],sqrt(phi$s2[p]))
}
hist(Y,breaks = 20,freq=FALSE)
for (p in 1:P){
  lines(yseq,f.per.cluster[p,],col=p)
}





#### 3. Gaussians with same variances but estimated

library(mclust)

Pmax=10
res.bic <- mclustBIC(Y, G=1:Pmax,modelName="E")   # ---> donne 1 group
plot(res.bic)

vraiss=c()
for (i in 1:Pmax){
  dist.Y=dist(Y)
  Clust.cah<- hclust(dist.Y^2, method = "ward.D")  
  classif.outliers <- Mclust(Y,G=i,modelName="E",initialization=Clust.cah)
  #classif.outliers <- Mclust(Y,G=i,modelName="E")
  vraiss[i]=classif.outliers$loglik
}
classif.outliers$parameters
table(classif.outliers$classification)
#aa=apply(classif.outliers$z,1,max) 
plot(Y,col=classif.outliers$classification,pch=16,cex=0.8)



#########
# comparison between the classifications

table(cluster,cluster_imp)



#### 4. Gaussians with unequal variances with P=2
P=2
dist.Y=dist(Y)
Clust.cah<- hclust(dist.Y^2, method = "ward.D")  
classif.outliers <- Mclust(Y,G=P,modelName="V",initialization=Clust.cah)
classif.outliers$parameters

plot(Y,col=classif.outliers$classification,pch=16,cex=0.8)


#### 5. P=2, N(mu,s2) et N(mu, a*s2)

# initialisation with a mixture model with unequal variance
P=2
dist.Y=dist(Y)
Clust.cah<- hclust(dist.Y^2, method = "ward.D")  
classif.outliers <- Mclust(Y,G=P,modelName="V",initialization=Clust.cah)

mu=classif.outliers$parameters$mean
s2=classif.outliers$parameters$variance$sigmasq
prop=classif.outliers$parameters$pro

b    = order(s2)
s2   = sort(s2)
mu   = mu[b]
prop = prop[b]

phi=list()
phi=list(prop=prop,mu=mu,s2=s2)
# 
Out.EM_SameMean_PropVariance =EM.algo_SameMean_PropVariance(Y,phi,P)
phi_SameMean_PropVariance= Out.EM_SameMean_PropVariance$phi
tau_SameMean_PropVariance=Out.EM_SameMean_PropVariance$tau
empty_SameMean_PropVariance = Out.EM_SameMean_PropVariance$empty
dv_SameMean_PropVariance= Out.EM_SameMean_PropVariance$dv
lvinc_SameMean_PropVariance=Out.EM_SameMean_PropVariance$lvinc

phi_SameMean_PropVariance
#### Classification

cluster_SameMean_PropVariance = apply(tau_SameMean_PropVariance,1,which.max)
#table(cluster.true,cluster)
table(cluster_SameMean_PropVariance)

#### Plot

plot(Y,col=cluster_SameMean_PropVariance+1,pch=16,cex=0.8)

yseq=seq(min(Y),max(Y),0.001)
f.per.cluster=matrix(0,ncol=length(yseq),nrow=P)
for (p in 1:P){
  f.per.cluster[p,]=dnorm(yseq,phi_SameMean_PropVariance$mu[p],sqrt(phi_SameMean_PropVariance$s2[p]))
}
hist(Y,breaks = 50,freq=FALSE)
for (p in 1:P){
  lines(yseq,f.per.cluster[p,],col=p)
}

table(classif.outliers$classification,cluster_SameMean_PropVariance)
table(cluster.true,cluster_SameMean_PropVariance)
table(cluster.true,classif.outliers$classification)


hist(Y[cluster_SameMean_PropVariance==1],breaks = 20,freq=FALSE)
fnorm=dnorm(yseq,0,1)
lines(yseq,fnorm,col=p)

#######
# Simualtion gaussian and uniform
n=200
Y <- rnorm(n,0,1)
pos.outliers <- sample(1:n,10,FALSE)
Y0 <- runif(10, min = -8, max = 8)
Y[pos.outliers] <- X0
color=rep(1,n)
color[pos.outliers]=2

plot(Y,col=color,pch=16,cex=0.8)
hist(Y,breaks = 20,freq=FALSE)
# # Y transform
# Y.min=min(Y)
# Y.max=max(Y)
# Y.new=sqrt(Y-Y.min)/sqrt(Y.max-Y.min)
# 
# Z=(Y-mean(Y))/sd(Y)
# # Plot
# plot(Y,col=cluster.true,pch=16,cex=0.8)
# hist(Y,breaks = 20)
# plot(Z,col=cluster.true,pch=16,cex=0.8)
# hist(Z,breaks = 20)
# 
# ####? paper of Cousineau

#####
# fitting

# Take first the true values and try the ratio

Y.min=min(Y)
Y.max=max(Y)
# an oextreme outlier
a.point= Y[pos.outliers[1]]
plot(Y,col=color,pch=16,cex=0.8)
points(pos.outliers[1],a.point,col="green")


ratio.point=dnorm(a.point,0,1)/dunif(a.point, min = Y.min, max = Y.max)
ratio.point

# a moderate outlier 
a.point= Y[38]
plot(Y,col=color,pch=16,cex=0.8)
points(38,a.point,col="green")


ratio.point=dnorm(a.point,0,1)/dunif(a.point, min = Y.min, max = Y.max)
ratio.point
# a normal point 

a.point= Y[1]
plot(Y,col=color,pch=16,cex=0.8)
points(1,a.point,col="green")


ratio.point=dnorm(a.point,0,1)/dunif(a.point, min = Y.min, max = Y.max)
ratio.point

# remove points potential "good observation"

d1 = abs(diff(Y))
d2 = sapply(c(1:length(d1)), function(x) mean(d1[x:(x+1)], na.rm = T))

Q3 <- quantile(d2, .75, na.rm = TRUE)
d3 = which(d2>Q3)


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







