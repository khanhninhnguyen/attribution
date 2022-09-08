########## This function is used to compare the two method of outlier screening
source(paste0(path_code_att,"UsedFunctions.R"))

source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))

library(tidyverse) 
library(mclust)
########################@
#### Simulation

# Parameters
n             = 1000   #length of the series
prob.outliers = 0.02
size.outliers = 3
P.true             = 2 # P can be 1, 2 or 3

# Simulated series ------------------
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

# conclusion: emilie and Olivier threshold give the same results and same with the cluster themself ------------------------

nb.sim = 1000
res <- data.frame(matrix(NA, nrow = nb.sim, ncol = 3))
P=3
for (i in c(1:nb.sim)) {
  # sim series with P.true groups
  set.seed(i)
  SimData = SimulatedSeries(n,P=4,prob.outliers,size.outliers)
  Y=SimData$Y
  cluster.true= which(SimData$cluster.true != 1)
  
  # classification with fixed Ptrue group 
  Out.EM.init_imp=EM.init_imp(Y,P=P,option.init="CAH")
  Out.EM_imp =EM.algo_imp(Y,Out.EM.init_imp$phi,P=P,Out.EM.init_imp$Id.cluster1)
  tau_imp = Out.EM_imp$tau 
  cluster_imp = apply(tau_imp,1,which.max)
  main.g = as.numeric(names(sort(table(cluster_imp),decreasing=TRUE)[1]))
  emi = which(cluster_imp!=main.g)
  
  # algorithm from olivier
  a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = 3)
  obi = sort(a$point.rm)
  clust_obi = rep(1, n)
  clust_obi[obi] <- 2
  mcl = Mclust(Y, G=2, model="V")
  mcl.g = as.numeric(names(sort(table(mcl$classification),decreasing=TRUE)[1]))
  nin = which(mcl$classification!=mcl.g)
  
  res[i,1] <- length(which(emi %in% cluster.true == TRUE))/length(cluster.true)
  res[i,2] <- length(which(obi %in% cluster.true == TRUE))/length(cluster.true)
  res[i,3] <- length(which(nin %in% cluster.true == TRUE))/length(cluster.true)
  
}

summary(res)

hist(rbeta(100000,100,1)*10)

plot(Y, col = cluster_imp)
plot(Y, col = mcl$classification)

clust_obi = rep(1, n)
clust_obi[obi] <- 2
plot(Y, col = clust_obi)

plot(Y, type = "l")
points(Y, col=clust_obi)
abline(h=3)
abline(h=-3)
# To see the impact of skewness on the outlier dectection 
# skewed data 
library(univOutl)
library(sn)# Azzalini skewnormal package
set.seed(7*11*13) 
test  <-  rsn(n = 1000, xi = 0, omega = 1, alpha=1)
hist(test)

nb.sim = 100
res <- data.frame(matrix(NA, nrow = nb.sim, ncol = 6))
P=3
for (i in c(1:nb.sim)) {
  # sim series with P.true groups
  set.seed(i)
  SimData = SimulatedSeries(n,P=3,prob.outliers,size.outliers)
  Y=SimData$Y
  cluster.true= which(SimData$cluster.true != 1)
  
  # classification with fixed Ptrue group 
  Out.EM.init_imp=EM.init_imp(Y,P=P,option.init="CAH")
  Out.EM_imp =EM.algo_imp(Y,Out.EM.init_imp$phi,P=P,Out.EM.init_imp$Id.cluster1)
  tau_imp = Out.EM_imp$tau 
  cluster_imp = apply(tau_imp,1,which.max)
  main.g = as.numeric(names(sort(table(cluster_imp),decreasing=TRUE)[1]))
  emi = which(cluster_imp!=main.g)
  
  # algorithm from olivier
  a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = 3)
  obi = sort(a$point.rm)
  clust_obi = rep(1, n)
  clust_obi[obi] <- 2
  # automatic clustering
  # mcl = Mclust(Y, G=3, model="V")
  # mcl.g = as.numeric(names(sort(table(mcl$classification),decreasing=TRUE)[1]))
  # nin = which(mcl$classification!=mcl.g)
  # adjust quantile 
  # adj.qua = boxB(Y,method="adjbox")
  # clust_qua = rep(1, n)
  # clust_qua[adj.qua$] <- 2
  
  res[i,1] <- length(which(emi %in% cluster.true == TRUE))/length(cluster.true)
  res[i,2] <- length(which(obi %in% cluster.true == TRUE))/length(cluster.true)
  res[i,3] <- length(which(nin %in% cluster.true == TRUE))/length(cluster.true)
  
  res[i,4] <- length(which(emi %in% cluster.true == TRUE))/length(emi)
  res[i,5] <- length(which(obi %in% cluster.true == TRUE))/length(obi)
  res[i,6] <- length(which(nin %in% cluster.true == TRUE))/length(nin)
}
summary(res)


hist(Y, breaks = 100, main =  "Histogram of simulated data")

r <- rep(NA, 1000)
for(k in c(1:1000)){
  set.seed(k)
  test  <-  rsn(n = 1000, xi = 0, omega = 1, alpha=1)
  r[k] <-  robustbase::scaleTau2(test)
}

