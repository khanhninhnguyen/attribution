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
  else if(x == 7){ P = 9 } # normal outlier + different simulated 
  
  return(P)
}

# # Loop to check all methods ------------------------------------------
nb.sim = 1000
prob.outliers = 0.1
n0 = 300
outlier.mod = c(7)
list.gmm.imp = list()
list.gmm.une = list()
test.length = c()
list.loop = c(5)
P.ini = choose_model(outlier.mod)
data.sim = list()
# loop

for (j in c(1:length(list.loop))) {
  res <- data.frame(matrix(NA, nrow = nb.sim, ncol = 20))
  res.all <- list()
  size.outlier = list.loop[j]
  for (i in c(1:nb.sim)) {
    # sim series with P.true groups
    print(i)
    set.seed(i)
    SimData = SimulatedSeries(n = n0, P = P.ini, prob.outliers = prob.outliers, size.outliers = size.outlier, rho = 0, theta = 0)
    Y=SimData$Y
    cluster.true= which(SimData$cluster.true != 1)
    test.length = c(test.length, length(cluster.true))
    
    # classification with fixed Ptrue group 
    gmm.imp = GMM_imp(P = 3, Y = Y, thres = 0.5)
    
    # algorithm from olivier
    a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'def', iter = 0, estimator = "mad", fix.thres = 3, loes = 0)
    # a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'sigma', iter = 0, estimator = "Sca", fix.thres = 0, loes = 0)
    obi = sort(a$point.rm)
    clust_obi = rep(1,length(Y))
    clust_obi[obi] <- 2
    
    # 2 groups unequal variances
    gmm.unequal = GMM_sameMean_uequalvar(P = 2, Y = Y, thres = 0.5)
    # 
    # OS Tests
    OS = testOS(Y, thres = 9)
    # 
    # # Tests
    testI <- getOutliersI(Y, rho=c(1,1), FLim=c(0.1,0.9), distribution="normal")
    testI.out = c(testI$iLeft, testI$iRight)
    testI.res <- rep(1, length(Y))
    if (length(testI.out) >0){
      testI.res[testI.out] <- 2
    }

    res.all[[i]] <-  list( cluster.true = SimData$cluster.true,
                           gmm.imp = gmm.imp$cluster,
                           gmm.uneq = gmm.unequal$cluster,
                           three.sigma = clust_obi,
                           testOS = OS$cluster,
                           MethodI = testI.res )
    list.gmm.imp[[i]] <- gmm.imp
    list.gmm.une[[i]] <- gmm.unequal
    data.sim[[i]] <- Y
  }
  save(res.all, file = paste0(path_results,"attribution/comparison_screening_methods", prob.outliers,"mod", outlier.mod, size.outlier, "1.RData"))
  save(list.gmm.imp, file = paste0(path_results,"attribution/gmm.imp", prob.outliers,"mod", outlier.mod, size.outlier, "1.RData"))
  save(list.gmm.une, file = paste0(path_results,"attribution/gmm.unequal", prob.outliers,"mod", outlier.mod, size.outlier, "1.RData"))
  save(data.sim, file = paste0(path_results,"attribution/data.sim", prob.outliers,"mod", outlier.mod, size.outlier, "1.RData"))
}

print(table(test.length))
# adding results of GMM regarding different threshold  --------------------
size.outlier = 5
res.all <- get(load( file = paste0(path_results,"attribution/comparison_screening_methods", prob.outliers,"mod", outlier.mod, size.outlier,".RData")))
res.imp <- get(load( file = paste0(path_results,"attribution/gmm.imp", prob.outliers,"mod", outlier.mod, size.outlier,".RData")))
res.une <- get(load( file = paste0(path_results,"attribution/gmm.unequal", prob.outliers,"mod", outlier.mod, size.outlier,".RData")))
sim.dat <- get(load( file = paste0(path_results,"attribution/data.sim", prob.outliers,"mod", outlier.mod, size.outlier,".RData")))

thres.outlier = function(tau, thres){
  cluster_imp0 = apply(tau, 1, which.max)
  main.g = as.numeric(names(sort(table(cluster_imp0),decreasing=TRUE)[1]))
  cluster_imp = sapply(tau[,main.g], function(x) ifelse(x>thres, 1, 2) ) # change here to modify how to choose 
  outliers = which(cluster_imp >1)
  # return(cluster_imp)
  return(outliers)
}
thres.list = seq(0.1,0.4,0.1)
for (i in c(1:length(res.all))) {
  tau_imp = as.data.frame(res.imp[[i]]$tau)
  tau_une = as.data.frame(res.une[[i]]$tau)
  for (k in c(1:length(thres.list))) {
    thres = thres.list[k]
    res.all[[i]][[paste0("gmm.imp.",thres)]] <- thres.outlier(tau = tau_imp, thres = thres)
    res.all[[i]][[paste0("gmm.uneq.",thres)]] <- thres.outlier(tau = tau_une, thres = thres)
  }
}


# run again with estimated variance in three sigma rule -------------------

# res.load <- get(load( file = paste0(path_results,"attribution/comparison_screening_methods", prob.outliers,"mod", outlier.mod = 5, ".RData")))
# 
# 
# nb.sim = 1000
# size.outliers = 5
# prob.outliers = 0.1
# n0 = 300
# outlier.mod = c(1, 2, 5, 6)
# res.tot <- list()
# res.all.O <- list()
# 
# for (j in c(1:length(outlier.mod))) {
#   P.ini = choose_model(outlier.mod[j])
#   res <- data.frame(matrix(NA, nrow = nb.sim, ncol = 20))
#   res.all <- list()
#   for (i in c(1:nb.sim)) {
#     # sim series with P.true groups
#     print(i)
#     set.seed(i)
#     SimData = SimulatedSeries(n = n0, P = P.ini, prob.outliers = prob.outliers, size.outliers = size.outliers, rho = 0, theta = 0)
#     Y=SimData$Y
#     cluster.true= which(SimData$cluster.true != 1)
#     
#     # algorithm from olivier
#     a <- screen.O(Y = data.frame(y = Y), name.var = "y", method = 'sigma', iter = 0, estimator = "Sca", fix.thres = 0, loes = 0)
#     obi = sort(a$point.rm)
#     clust_obi = rep(1,length(Y))
#     clust_obi[obi] <- 2
#     
#     res.all[[i]] <- clust_obi
#   }
#   res.all.O[[j]] <- res.all
# }

# analyze all result ------------------------------------------------------

# for (k in c(1:length(res.total))) {
#   res.total[[k]]$three.sigma <- res.all.O[[5]][[k]]
# }

# res.total.arrange = list()
# translate name group some time the main group is 2 but this mistake is corrected in the source code 
# for (i in c(1:length(res.total))) {
#   res.i = as.data.frame(res.total[[i]])
#   res.arrange <- data.frame(true.cluster = res.i$cluster.true)
#   for (j in c(1:(dim(res.i)[2] -1))) {
#     col.j <- rep(NA, nrow(res.i))
#     g0 = as.numeric(names(sort(table(res.i[,(j+1)]),decreasing=TRUE)))
#     for (k in c(1:length(g0))) {
#       col.j[which(res.i[,(j+1)] == g0[k])] <- k
#     }
#     res.arrange[names(res.i)[j+1]] <- col.j
#   }
#   res.total.arrange[[i]] <- res.arrange
# }


# compute statistics  -----------------------------------------------------

sta <- lapply(c(1:nb.sim), function (x){
  r = as.data.frame(res.all[[x]])
  val.true = r$cluster.true
  n.pos = length(which(val.true > 1))
  n.neg = length(which(val.true == 1))
  sim.seri = sim.dat[[x]]
  a <- sapply(c(1:(dim(r)[2] -1)), function (y) {
    P.true = which(val.true > 1)
    N.true = which(val.true == 1)
    P.ind = which(r[,(y+1)] > 1)
    N.ind = which(r[,(y+1)] == 1)
    FN.pos = which(r[,(y+1)] == 1 & val.true != 1)
    FP.pos = which(r[,(y+1)] != 1 & val.true == 1)
    TP = length(which(P.ind %in% P.true == TRUE))
    TN = length(which(N.ind %in% N.true == TRUE))
    FN = n.pos - TP
    FP = n.neg - TN
    TPR = TP/n.pos
    FPR = FP/n.neg
    TNR = TN/n.neg
    PPV = TP/(TP+FP)
    if(is.na(PPV) == TRUE){ PPV = 0}
    ACC = (TPR*n.pos+TNR*n.neg)/(n.pos+n.neg)
    F1 = (2*TP)/(2*TP + FN + FP)
    TS = TP/( TP+ FN+FP)
    if(FN >0){
      FNs = sum(sim.seri[FN.pos]^2)/FN
    }else{FNs = 0}
    if(FP >0){
      FPs = sum(sim.seri[FP.pos]^2)/FP
    }else{FPs = 0}
    c(TP, TN, FP, FN, TPR, FPR, FNs, FPs, PPV, ACC, F1, TS)
  })
  b = data.frame(t(a))
  colnames(b) <- c("TP", "TN", "FP", "FN", "TPR", "FPR", "FNs", "FPs", "PPV", "ACC", "F1", "TS")
  b
})

TP.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$TP)))
FP.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$FP)))
TN.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$TN)))
FN.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$FN)))
TPR.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$TPR)))
FPR.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$FPR)))
PPV.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$PPV)))
ACC.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$ACC)))
F1.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$F1)))
TS.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$TS)))
FPs.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$FPs)))
FNs.data = as.data.frame( t(sapply(c(1:nb.sim), function(x) sta[[x]]$FNs)))

all.data <- list(TP.data, FP.data, TN.data, FN.data, TPR.data, FPR.data, FNs.data, FPs.data, PPV.data, ACC.data, F1.data, TS.data)
mean.val = data.frame(matrix(NA, ncol = 13, nrow = 12))
sd.val = data.frame(matrix(NA, ncol = 13, nrow = 12))
for (h in c(1:length(all.data))) {
  dat = all.data[[h]]
  mean.val[h,] <- round(apply(dat,2,mean), digits = 3)
  sd.val[h,] <- round(apply(dat,2,sd), digits = 3)
}

colnam = names(res.all[[1]])[-1]
colnam[c(1,2)] <- paste(colnam[c(1,2)], ".0.5", sep = "") 
colnames(mean.val) = colnam
colnames(sd.val) = colnam
rownam = c("TP", "FP", "TN", "FN", "TPR", "FPR", "FNs", "FPs", "PPV", "ACC", "F1", "TS")
rownames(mean.val) = rownam
rownames(sd.val) = rownam

mean.data = mean.val[,order(colnames(mean.val))]
sd.data = sd.val[,order(colnames(sd.val))]


write.table(mean.data , 
            file = paste0(path_results,'attribution/mean', outlier.mod, size.outlier, ".txt"), 
            sep="\t",
            col.names = TRUE,
            row.names = TRUE,
            quote=FALSE)

write.table(sd.data, 
            file = paste0(path_results,'attribution/sd',outlier.mod, size.outlier, ".txt"), 
            sep="\t",
            col.names = TRUE,
            row.names = TRUE,
            quote=FALSE)

# plot TPR, TNR and ACC -------------
name.methods = names(res.all[[1]])[-1]
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

# individual case ----------------------------------

hist(SimData$Y, breaks = 100)

# results of gmm.imp
jpeg(paste0(path_results,"attribution/gmm.imp",outlier.mod, size.outlier, prob.outliers, ".jpeg" ),
     width = 2000, height = 1000,res = 300)
par(mar = c(2, 2, 2, 2))

plot(SimData$Y, cex = 0.8, ann=FALSE)
points(cluster.true, SimData$Y[cluster.true], cex = 0.8, col = "red", pch = 19)
points(gmm.imp$gmm.imp.0.1, SimData$Y[gmm.imp$gmm.imp.0.1], cex = 0.9, col = "blue", pch = 0)
points(gmm.imp$gmm.imp.0.2, SimData$Y[gmm.imp$gmm.imp.0.2], cex = 1, col = "purple", pch = 2)
points(gmm.imp$gmm.imp.0.3, SimData$Y[gmm.imp$gmm.imp.0.3], cex = 1, col = "green", pch = 3)
points(gmm.imp$gmm.imp.0.4, SimData$Y[gmm.imp$gmm.imp.0.4], cex = 1.2, col = "orange", pch = 4)
points(gmm.imp$outliers, SimData$Y[gmm.imp$outliers], cex = 1.5, col = "pink", pch = 5)

legend("bottomright", c("true", "0.1", "0.2", "0.3","0.4", "0.5"), 
       col=c("red","blue","purple", "green","orange","pink"),
       pch=c(19,0, 2, 3, 4, 5),
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n"
)
dev.off()

a = reshape2::melt(tau_imp)
colnames(a) <- c("index", "group", "value")
a$group = as.factor(a$group)
ggplot(a, aes(x=index, y=value, col = group))+geom_point()+ theme_bw()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"),
        legend.text=element_text(size=15)) + ylab("Tau")



# results of Olivier and test OS
jpeg(paste0(path_results,"attribution/OB",outlier.mod, size.outlier, prob.outliers, ".jpeg" ),
     width = 2000, height = 1000,res = 300)
par(mar = c(2, 2, 2, 2))

plot(SimData$Y, cex = 0.8, ann=FALSE)
points(cluster.true, SimData$Y[cluster.true], cex = 0.8, col = "red", pch = 19)
points(a$point.rm, SimData$Y[a$point.rm], cex = 0.9, col = "blue", pch = 0)
# points(OS$outliers, SimData$Y[OS$outliers], cex = 0.9, col = "green", pch = 2)
legend("bottomright", c("true", "3.sigma.m"), 
       col=c("red","blue"),
       pch=c(19,0),
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n"
)
dev.off()

# results of GMM.unequal 
jpeg(paste0(path_results,"attribution/gmm.une",outlier.mod, size.outlier, prob.outliers, ".jpeg" ),
     width = 2000, height = 1000,res = 300)
par(mar = c(2, 2, 2, 2))

plot(SimData$Y, cex = 0.8, ann=FALSE)
points(cluster.true, SimData$Y[cluster.true], cex = 0.8, col = "red", pch = 19)
points(gmm.unequal$gmm.imp.0.1, SimData$Y[gmm.unequal$gmm.imp.0.1], cex = 0.9, col = "blue", pch = 0)
points(gmm.unequal$gmm.imp.0.2, SimData$Y[gmm.unequal$gmm.imp.0.2], cex = 1, col = "purple", pch = 2)
points(gmm.unequal$gmm.imp.0.3, SimData$Y[gmm.unequal$gmm.imp.0.3], cex = 1, col = "green", pch = 3)
points(gmm.unequal$gmm.imp.0.4, SimData$Y[gmm.unequal$gmm.imp.0.4], cex = 1.2, col = "orange", pch = 4)
points(gmm.unequal$outliers, SimData$Y[gmm.unequal$outliers], cex = 1.5, col = "pink", pch = 5)

legend("bottomright", c("true", "0.1", "0.2", "0.3","0.4", "0.5"), 
       col=c("red","blue","purple", "green","orange","pink"),
       pch=c(19,0, 2, 3, 4, 5),
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n"
)
dev.off()
a = reshape2::melt(tau_SameMean_PropVariance)
colnames(a) <- c("index", "group", "value")
a$group = as.factor(a$group)
ggplot(a, aes(x=index, y=value, col = group))+geom_point()+ theme_bw()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"),
        legend.text=element_text(size=15)) + ylab("Tau")


# results of Loo's methods 

plot(SimData$Y, cex = 0.8)
points(cluster.true, SimData$Y[cluster.true], cex = 0.8, col = "red", pch = 19)
points(testI.out, SimData$Y[testI.out], cex = 0.9, col = "blue", pch = 0)
points(testII.out, SimData$Y[testII.out], cex = 0.9, col = "green", pch = 2)
legend("bottomright", c("true", "MethodI", "MethodII"), 
       col=c("red","blue","green"),
       pch=c(19,0, 2),
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n"
)

outlierPlot(SimData$Y,testI,mode="qq")
outlierPlot(Y,testII,mode="residual")
