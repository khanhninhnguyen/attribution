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
P.true             = 2  # P can be 1, 2 or 3

# Simulated series
SimData = SimulatedSeries(n,P.true,prob.outliers,size.outliers)
Y=SimData$Y
cluster.true=SimData$cluster.true
# Plot
plot(Y,col=cluster.true,pch=16,cex=0.8)
hist(Y,breaks = 20)

              ######
              #### Real data
load(paste0(path_results,"auck.2005-11-07.whng.RData"))
dates=Y$date
Y=Y$gps.gps
n=length(Y)

rg.na=which(is.na(Y))
Y=Y[-rg.na]
dates=dates[-rg.na]
n=length(Y)
# Plot
plot(dates,Y,pch=16,cex=0.8)
hist(Y,breaks = 20)


              #### EM with imposed distribution for one group and variances=1 for P=1:Pmax
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
hist(Y,breaks = 20,freq=FALSE)
for (p in 1:P){
  lines(yseq,f.per.cluster[p,],col=p)
}


            #### EM without imposed distribution for P=1:Pmax
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





###############################
#### Test de mclust

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
