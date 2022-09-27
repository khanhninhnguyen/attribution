########## This function is used to compare the two method of outlier screening
source(paste0(path_code_att,"UsedFunctions.R"))

source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))

library(tidyverse) 
library(mclust)
library(extremevalues)
library(EnvStats)

########################@
#### Simulation

# Parameters
# n             = 1000   #length of the series
# prob.outliers = 0.02
# size.outliers = 3
choose_model <- function(x){
  if(x == 1){ P = 5 } # fixed value of outlier
  else if(x == 2){ P = 3 } # normal outlier 
  else if(x == 3){ P = 4 } # normal overlayed outlier 
  else if(x == 4){ P = 7 } # normal overlayed outlier + arma data
  else if(x == 5){ P = 8 } # normal overlayed outlier + replaced at a threshold
  else if(x == 6){ P = 1 } # no outlier
  return(P)
}
# model.out = 3

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
n0 = 300
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
nb.sim = 1000
size.outliers = 5
prob.outliers = 0.1
n0 = 300
outlier.mod = c(1)
res.tot <- list()
res.all1 <- list()

for (j in c(1:length(outlier.mod))) {
  P.ini = choose_model(outlier.mod[j])
  res <- data.frame(matrix(NA, nrow = nb.sim, ncol = 20))
  res.all <- list()
  for (i in c(1:nb.sim)) {
    # sim series with P.true groups
    print(i)
    set.seed(i)
    SimData = SimulatedSeries(n = n0, P = P.ini, prob.outliers = prob.outliers, size.outliers = size.outliers, rho = 0, theta = 0)
    Y=SimData$Y
    cluster.true= which(SimData$cluster.true != 1)
    
    # classification with fixed Ptrue group 
    gmm.imp.0.1 = GMM_imp(P = 3, Y = Y, thres = 0.1)
    gmm.imp.0.2 = GMM_imp(P = 3, Y = Y, thres = 0.2)
    gmm.imp.0.3 = GMM_imp(P = 3, Y = Y, thres = 0.5)
    
    # algorithm from olivier
    a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = 3)
    obi = sort(a$point.rm)
    clust_obi = rep(1,length(Y))
    clust_obi[obi] <- 2
    
    # 2 groups unequal variances
    # gmm.unequal.0.1 = GMM_sameMean_uequalvar(P = 2, Y = Y, thres = 0.1)
    # gmm.unequal.0.2 = GMM_sameMean_uequalvar(P = 2, Y = Y, thres = 0.2)
    gmm.unequal.0.3 = GMM_sameMean_uequalvar(P = 2, Y = Y, thres = 0.5)
    
    # OS Tests
    OS = testOS(Y, thres = 5)
    
    # Tests
    testI <- getOutliersI(Y, rho=c(1,1), FLim=c(0.1,0.9), distribution="normal")
    testI.out = c(testI$iLeft, testI$iRight)
    testI.res <- rep(1, length(Y))
    if (length(testI.out) >0){
      testI.res[testI.out] <- 2
    }
    testII <- getOutliersII(Y, alpha=c(0.01, 0.01), FLim=c(0.1, 0.9),
                            distribution="normal", returnResiduals=TRUE)
    testII.out = c(testII$iLeft, testII$iRight)
    testII.res <- rep(1, length(Y))
    if (length(testII.out) >0){
      testII.res[testII.out] <- 2
    }
    # save results - TPR
    # TPR = list(gmm.imp.0.1$outliers, gmm.imp.0.1$outliers, gmm.imp.0.1$outliers,
    #            gmm.unequal.0.1$outliers, gmm.unequal.0.2$outliers, gmm.unequal.0.3$outliers,
    #            obi, OS$outliers, testI.out, testII.out)
    # for (k in c(1:length(TPR))) {
    #   res[i,k] <- length(which(TPR[[k]] %in% cluster.true == TRUE))/length(cluster.true)
    #   res[i,(k+length(TPR))] <- length(TPR[[k]])
    #   
    # }

    res.all[[i]] <-  list( cluster.true = SimData$cluster.true, 
                           # gmm.imp.0.1 = gmm.imp.0.1$cluster,
                           # gmm.imp.0.2 = gmm.imp.0.2$cluster, 
                           gmm.imp.0.3 = gmm.imp.0.3$cluster, 
                           # gmm.uneq.0.1 = gmm.unequal.0.1$cluster, 
                           # gmm.uneq.0.2 = gmm.unequal.0.2$cluster, 
                           gmm.uneq.0.3 = gmm.unequal.0.3$cluster, 
                           three.sigma = clust_obi, 
                           testOS = OS$cluster, 
                           MethodI = testI.res, 
                           MethodII = testII.res )
  }
  res.all1[[j]] <-  res.all
  res.tot[[j]] <- res
}

save(res.tot, file = paste0(path_results,"attribution/TPR", prob.outliers,"P8.RData"))
save(res.all1 , file = paste0(path_results,"attribution/comparison_screening_methods", prob.outliers,"mod", outlier.mod, size.outliers, "1.RData"))

res.load <- get(load( file = paste0(path_results,"attribution/comparison_screening_methods", prob.outliers,"mod", outlier.mod, size.outliers, "1.RData")))

res.total = res.all1[[1]]
res.total.arrange = list()
# translate name group 
for (i in c(1:length(res.total))) {
  res.i = as.data.frame(res.total[[i]])
  res.arrange <- data.frame(true.cluster = res.i$cluster.true)
  for (j in c(1:(dim(res.i)[2] -1))) {
    col.j <- rep(NA, nrow(res.i))
    g0 = as.numeric(names(sort(table(res.i[,(j+1)]),decreasing=TRUE)))
    for (k in c(1:length(g0))) {
      col.j[which(res.i[,(j+1)] == g0[k])] <- k
    }
    res.arrange[names(res.i)[j+1]] <- col.j
  }
  res.total.arrange[[i]] <- res.arrange
}

sta <- lapply(c(1:nb.sim), function (x){
  r = as.data.frame(res.total.arrange[[x]])
  val.true = r$true.cluster
  n.pos = length(which(val.true > 1))
  n.neg = length(which(val.true == 1))
  a <- sapply(c(1:(dim(res.i)[2] -1)), function (y) {
    P.true = which(val.true > 1)
    N.true = which(val.true == 1)
    P.ind = which(r[,(y+1)] > 1)
    N.ind = which(r[,(y+1)] == 1)
    TPR = length(which(P.ind %in% P.true == TRUE))/n.pos
    TNR = length(which(N.ind %in% N.true == TRUE))/n.neg
    ACC = (TPR*n.pos+TNR*n.neg)/(n.pos+n.neg) 
    c(TPR, TNR, ACC)
  })
  b = data.frame(t(a))
  colnames(b) <- c("TPR", "TNR", "ACC")
  b
})

# plot TPR, TNR and ACC
name.methods = names(res.total[[1]])[-1]
TPR.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$TPR)))
colnames(TPR.data) <- name.methods
dat <- reshape2::melt(TPR.data[,-1])
jpeg(paste0(path_results,"attribution/variances/bias.model.TPR",outlier.mod, size.outliers, prob.outliers, "601.jpeg" ),
     width = 3500, height = 1500,res = 300)
p <- ggplot(dat, aes(x=variable, y=value)) + 
  geom_boxplot()+theme_bw()  +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=10,face="bold"))+xlab("methods") + ylab("TPR")
p
dev.off()

TNR.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) 1 - sta[[x]]$TNR)))
colnames(TNR.data) <- name.methods
dat <- reshape2::melt(TNR.data[,-1])
jpeg(paste0(path_results,"attribution/variances/bias.model.TNR",outlier.mod, size.outliers, prob.outliers, "601.jpeg" ),
     width = 3500, height = 1500,res = 300)
p <- ggplot(dat, aes(x=variable, y=value)) + 
  geom_boxplot()+theme_bw()  +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=10,face="bold"))+xlab("methods") + ylab("FPR")
p
dev.off()

ACC.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$ACC)))
colnames(ACC.data) <- name.methods
dat <- reshape2::melt(ACC.data[,-1])
jpeg(paste0(path_results,"attribution/variances/bias.model.ACC",outlier.mod, size.outliers, prob.outliers, "601.jpeg" ),
     width = 3500, height = 1500,res = 300)
p <- ggplot(dat, aes(x=variable, y=value)) + 
  geom_boxplot()+theme_bw()  +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=10,face="bold"))+xlab("methods") + ylab("ACC")
p
dev.off()



# testOS ------------------------------------------------------------------
thres.t <- 2
tot <- data.frame(matrix(NA, nrow = length(thres.t), ncol = 3))
for (j in c(1:length(thres.t))) {
  res.all <- list()
  statis = data.frame(matrix(NA, ncol = 3, nrow = nb.sim))
  for (i in c(1:nb.sim)) {
    # sim series with P.true groups
    print(i)
    set.seed(i)
    SimData = SimulatedSeries(n = n0, P = choose_model(6), prob.outliers = prob.outliers, size.outliers = size.outliers, rho = 0, theta = 0)
    Y=SimData$Y
    cluster.true= which(SimData$cluster.true != 1)
    
    # OS = testOS(Y, thres = thres.t[j])
    a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = thres.t[j])
    obi = sort(a$point.rm)
    clust_obi = rep(1,length(Y))
    clust_obi[obi] <- 2
    rst.outlier = obi
    rst.cluster = clust_obi
    
    TP = length(which(rst.outlier %in% cluster.true == TRUE))/length(cluster.true)
    N.true = which(SimData$cluster.true == 1)
    N.pre = which(rst.cluster == 1)
    TN = length(which(N.pre %in% N.true))/length(N.true)
    ACC = (TP*length(cluster.true)+TN*(300 - length(cluster.true)))/300
    
    statis[i, ] <- c( TP, 1-TN, ACC)
    res.all[[i]] <- list( cluster.true, rst.cluster)
  }
  tot[j,] <- colMeans(statis)
}

colnames(tot) <- c("TPR", "FPR", "ACC")
tot$thres = thres.t
dat <- reshape2::melt(tot, id = "thres")

jpeg(paste0(path_results,"attribution/variances/bias.model.ACC",outlier.mod, size.outliers, prob.outliers, "60OS.jpeg" ),
     width = 3500, height = 1500,res = 300)
p <- ggplot(dat, aes(x= thres, y=value, col = variable)) + 
  geom_point()+theme_bw()  +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=10,face="bold"))+xlab("methods") + ylab("ACC")
p
dev.off()

boxplot(statis)

# individual case 

plot(SimData$Y, cex = 0.5)
points(cluster.true, SimData$Y[cluster.true], cex = 0.5, col = "red",pch = 19)
points(gmm.imp.0.1$outliers, SimData$Y[gmm.imp.0.1$outliers], cex = 0.6, col = "blue", pch = 0)
points(gmm.imp.0.2$outliers, SimData$Y[gmm.imp.0.2$outliers], cex = 0.6, col = "blue", pch = 2)
points(gmm.imp.0.3$outliers, SimData$Y[gmm.imp.0.3$outliers], cex = 0.6, col = "blue", pch = 3)






