########## This function is used to compare the two method of outlier screening

rm(list=ls())
setwd("/home/knguyen/Documents/PhD/Code/from_Emilie/")
source("UsedFunctions.R")
source("/home/knguyen/Documents/PhD/Code/attribution/support_screening.R")
source("/home/knguyen/Documents/PhD/Code/attribution/sliding_variance.R")


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

a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad")
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

P.best_imp=which(BIC.c_imp==max(BIC.c_imp))
phi_imp=PHI_imp[[P.best_imp]]
tau_imp=TAU_imp[[P.best_imp]]
P.best_imp

P=2
Out.EM.init_imp=EM.init_imp(Y,P,option.init)
Out.EM_imp =EM.algo_imp(Y,Out.EM.init_imp$phi,P,Out.EM.init_imp$Id.cluster1)
PHI_imp[[P]]= Out.EM_imp$phi
TAU_imp[[P]]=Out.EM_imp$tau
empty_imp[P] = Out.EM_imp$empty
dv_imp[P] = Out.EM_imp$dv
lvincP_imp[P]=Out.EM_imp$lvinc

cluster_imp = apply(tau_imp,1,which.max)

plot(Y,col=cluster_imp,pch=16,cex=0.8)


yseq=seq(min(Y),max(Y),0.001)
f.per.cluster=matrix(0,ncol=length(yseq),nrow=P)
l=457/500
for (p in 1:P){
  if (p!=1){
    l = 43/500
  }
  f.per.cluster[p,]=l*dnorm(yseq,phi_imp$mu[p],sqrt(phi_imp$s2[p]))
}
hist(Y,breaks = 20,freq=FALSE)
lines(yseq,f.per.cluster[1,],col=1)

for (p in 2:P){
  lines(yseq,f.per.cluster[p,],col=p)
}

a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = 2.87)
emi = sort(which(cluster_imp==2))
a1 <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = 3)
obi = sort(a1$point.rm)
obi

# conclusion: emilie and Olivier threshold give the same results and same with the cluster themself 

