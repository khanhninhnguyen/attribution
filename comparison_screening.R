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
choose_model <- function(x){
  if(x == 1){ P = 5 } # fixed value of outlier
  else if(x == 2){ P = 3 } # normal outlier 
  else if(x == 3){ P = 4 } # normal overlayed outlier 
  else if(x == 4){ P = 7 } # normal overlayed outlier + arma data
  else if(x == 5){ P = 8 } # normal overlayed outlier + replaced at a threshold
  return(P)
}
model.out = 3

# Simulated series ------------------
SimData = SimulatedSeries(n,P.true,prob.outliers,size.outliers)
Y=SimData$Y
cluster.true=SimData$cluster.true
# Plot
plot(Y,col=cluster.true,pch=16,cex=0.8)
hist(Y,breaks = 20)


# conclusion: emilie and Olivier threshold give the same results and same with the cluster themself ------------------------

nb.sim = 100
res <- data.frame(matrix(NA, nrow = nb.sim, ncol = 6))
size.outliers = 3
prob.outliers = 0.1
n0 = 100
for (i in c(1:nb.sim)) {
  # sim series with P.true groups
  set.seed(i)
  SimData = SimulatedSeries(n = n0, P = P.ini, prob.outliers = prob.outliers, size.outliers = size.outliers, rho = 0, theta = 0)
  Y=SimData$Y
  cluster.true= which(SimData$cluster.true != 1)
  
  # classification with fixed Ptrue group 
  gmm.imp = GMM_imp(P = 3, Y = Y)
  emi = gmm.imp
  
  # algorithm from olivier
  a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = 3)
  obi = sort(a$point.rm)
  clust_obi = rep(1, n)
  clust_obi[obi] <- 2
  
  # 2 groups unequal variances
  gmm.unequal = GMM_sameMean_uequalvar(P = 2, Y = Y)
  nin = gmm.unequal
  # save results
  res[i,1] <- length(which(emi %in% cluster.true == TRUE))/length(cluster.true)
  res[i,2] <- length(which(obi %in% cluster.true == TRUE))/length(cluster.true)
  res[i,3] <- length(which(nin %in% cluster.true == TRUE))/length(cluster.true)
  
  res[i,4] <- length(emi)
  res[i,5] <- length(obi)
  res[i,6] <- length(nin)
  
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
  SimData = SimulatedSeries(n,P=6,prob.outliers,size.outliers)
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
  adj.qua = boxB(Y,method="adjbox")
  clust_qua = rep(1, n)
  clust_qua[adj.qua$outliers] <- 2
  
  res[i,1] <- length(which(emi %in% cluster.true == TRUE))/length(cluster.true)
  res[i,2] <- length(which(obi %in% cluster.true == TRUE))/length(cluster.true)
  res[i,3] <- length(which( clust_qua %in% cluster.true == TRUE))/length(cluster.true)
  
  res[i,4] <- length(which(emi %in% cluster.true == TRUE))/length(emi)
  res[i,5] <- length(which(obi %in% cluster.true == TRUE))/length(obi)
  res[i,6] <- length(which( clust_qua %in% cluster.true == TRUE))/length( clust_qua)
}
summary(res)


hist(Y, breaks = 100, main =  "Histogram of simulated data")

r <- rep(NA, 1000)
for(k in c(1:1000)){
  set.seed(k)
  test  <-  rsn(n = 1000, xi = 0, omega = 1, alpha=0.6)
  # r[k] <-  robustbase::scaleTau2(test)
  r[k] <-  IQR(test)/1.349
}
summary(r)


load(file=paste0(path_results, "skew.RData"))


# # Loop to check all situation  ------------------------------------------
nb.sim = 100
res <- data.frame(matrix(NA, nrow = nb.sim, ncol = 6))
size.outliers = 3
prob.outliers = 0.1
n0 = 100
outlier.mod = c(2,3)
for (j in c(1:length(outlier.mod))) {
  P.ini = choose_model(outlier.mod[j])
  
  for (i in c(1:nb.sim)) {
    # sim series with P.true groups
    set.seed(i)
    SimData = SimulatedSeries(n = n0, P = P.ini, prob.outliers = prob.outliers, size.outliers = size.outliers, rho = 0, theta = 0)
    Y=SimData$Y
    cluster.true= which(SimData$cluster.true != 1)
    
    # classification with fixed Ptrue group 
    gmm.imp = GMM_imp(P = 3, Y = Y)
    emi = gmm.imp$cluster
    
    # algorithm from olivier
    a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = 3)
    obi = sort(a$point.rm)
    clust_obi = rep(1, n)
    clust_obi[obi] <- 2
    
    # 2 groups unequal variances
    gmm.unequal = GMM_sameMean_uequalvar(P = 2, Y = Y)
    nin = gmm.unequal$cluster
    
    
    # save results
    res[i,1] <- length(which(emi %in% cluster.true == TRUE))/length(cluster.true)
    res[i,2] <- length(which(obi %in% cluster.true == TRUE))/length(cluster.true)
    res[i,3] <- length(which(nin %in% cluster.true == TRUE))/length(cluster.true)
    
    res[i,4] <- length(emi)
    res[i,5] <- length(obi)
    res[i,6] <- length(nin)
    
  }
  
}












