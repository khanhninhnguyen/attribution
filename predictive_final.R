# Final predictive rule prog

# output: last predtive rule and prediction results 
# input: test result, name of series used for training: 
# G-E, G-G', G-E',E-E', G'-E', G'-E, remove them by name of remove variable 
# prob for each configuration, significance.level, name.version


# load package and paths --------------------------------------------------

if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")
library(ggrepel)
library(mclust)
library(forecast)
library(lmtest)
library(ggvenn)
library(fields)
library(MASS)
library(reshape2)
library(factoextra)
library(mvtnorm)
library(caret)
library(randomForest)
library(rpart)
library(rpart.plot)
library(class)
library(rsample)
library(tidyverse)

source(paste0(path_code_att,"support_predictive_rule.R"))
path_restest <- paste0(path_results,"attribution/predictive_rule/")
file_path_Results=paste0(path_results,'attribution/predictive_rule/')
# name.results <- paste0(path_restest,"stats_test_real_data_corrected_dist_fdr.txt")
name.version ="FGLS_on_real_data_t.txt"
name.results <- paste0(path_restest, name.version) # name of test result file
NbSim = 37865
significance.level = 0.01
B = 20
offset=0
GE=0
number.pop = 3
# make the truth table and probability ------------------------------------

G=c(rep(1,9), rep(0,9),rep(-1,9),rep(0,9),rep(1,9),rep(-1,9))
E=c(rep(0,9),rep(-1,9),rep(0,9),rep(1,9),rep(-1,9),rep(1,9))
a=c(rep(0,3),rep(1,3),rep(-1,3))
Gp=rep(a,6)
Ep=rep(c(0,1,-1),18)
Y=data.frame(G=G,E=E,Gp=Gp,Ep=Ep)

Z=data.frame(GE=Y$G-Y$E,GGp=Y$G-Y$Gp,GEp=Y$G-Y$Ep,EEp=Y$E-Y$Ep,GpEp=Y$Gp-Y$Ep,GpE=Y$Gp-Y$E)
knitr::kable(head(Z))
List.names.tot <- colnames(Z)

Z.trunc <- Z 
Z.trunc[Z.trunc==-2]=-1
Z.trunc[Z.trunc==2]=1
knitr::kable(head(Z.trunc))

# prob <- c( 0.1944,0.0108,0.0108,0.0108,6e-04,6e-04,0.0108,6e-04,6e-04,
#           0.0108,0.0108, 0.1944,6e-04,6e-04,0.0108,6e-04,6e-04,0.0108,
#            0.1944,0.0108,0.0108,0.0108,6e-04,6e-04,0.0108,6e-04,6e-04,
#           0.0108, 0.1944,0.0108,6e-04,0.0108,6e-04,6e-04,0.0108,6e-04,
#           9e-04,9e-04,0.0162,5e-05,5e-05,9e-04,5e-05,5e-05,9e-04,9e-04,
#          0.0162,9e-04,5e-05,9e-04,5e-05,5e-05,9e-04,5e-05)
prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
          0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
          0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
          0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)
keep.config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)
Y <- Y[keep.config,]
Z <- Z[keep.config,]
Z.trunc <- Z.trunc[keep.config,]
num.conf <- 1:length(keep.config)
row.names(Z.trunc) <- num.conf

prob <- prob[keep.config]
prob <- prob/sum(prob)

# remove the duplicated due to the selection of variable ------------------
# keep the one with break in GPS
#+ 
remove.var = "G-E"

Y1 = Y
rownames(Y1) = NULL
#+ 
rm.ind = which(list.name.test == remove.var)
Z.trunc.rm = Z.trunc[,-rm.ind]
List.names.final = list.name.test[-rm.ind]
sum(duplicated2(Z.trunc[,-rm.ind]))/2
rg.duplicate <- which(duplicated2(Z.trunc.rm)) 
keep.ind = duplicated3(Z.trunc.rm, Y1)

remove.config = rg.duplicate[which(rg.duplicate %in% keep.ind == FALSE)]
Z.trunc.final <- Z.trunc[-remove.config,]

prob.final <- prob[-remove.config] 
prob.final <- prob.final/sum(prob.final)

Z.trunc.final.code <- Z.trunc.final
Z.trunc.final.code[Z.trunc.final.code==1]=3
Z.trunc.final.code[Z.trunc.final.code==-1]=2
Z.trunc.final.code[Z.trunc.final.code==0]=1
# save the list of configurations
head(Z.trunc.final.code)
config.list <- 1:38
config.list.final <- config.list[-remove.config] 
saveRDS(Z.trunc.final.code , file = paste0(file_path_Results,"List_config.rds"))

# 
Data.Res.Test <- read.table(name.results,header = TRUE, stringsAsFactors = FALSE)
# colnames(Data.Res.Test)[4:9] <- paste0("t", list.name.test)
colnames(Data.Res.Test)[4:9] <- paste0("t", list.name) 

# Separate into learning and test dataset --------------------------------------------------------

NumDetected.per.station <- Data.Res.Test %>% group_by(main) %>% 
  summarise(ldetected.per.station=length(unique(brp))) %>% dplyr::select(ldetected.per.station)

NumNearby.per.station.detection <- Data.Res.Test %>% group_by(main,brp) %>% 
  count() %>% as.data.frame() %>% dplyr::select(n)

List.names.final = List.names.tot[-rm.ind]
for (i in  1:length(List.names.final)){
  eval(parse(text=paste0("Data.",List.names.final[i],"=Keep.Data(List.names.final[i])")))  
}

# equal prob data  --------------------------------------------------------
R = 100
Nbconfig <- nrow(Z.trunc.final)

error.test.4.methods <- matrix(NA,nrow=B,ncol=4)
NbSim <- R*Nbconfig #in order to have at least 5 samples for each configurations. 
# B <- 5
NbSimLearn <- NbSim*0.8
NbSimTest <- NbSim*0.2
set.seed(1)

for (b in 13:B){
  
  ######
  # Vfold 
  existing.pop.Learn <- 0
  existing.pop.Test <- 0
  
  while ((existing.pop.Learn!=15) & (existing.pop.Test !=15)){
    trainIndex<-createDataPartition(1:nrow(Data.GGp) , p=0.8, list = FALSE)
    
    Data.GGp.Learn <- Data.GGp %>% slice(trainIndex)
    Data.GGp.Test <- Data.GGp %>% slice(-trainIndex)
    
    Data.EEp.Learn <- Data.EEp %>% slice(trainIndex)
    Data.EEp.Test <- Data.EEp %>% slice(-trainIndex)
    
    Data.GEp.Learn <- Data.GEp %>% slice(trainIndex)
    Data.GEp.Test <- Data.GEp %>% slice(-trainIndex)
    
    Data.GpEp.Learn <- Data.GpEp %>% slice(trainIndex)
    Data.GpEp.Test <- Data.GpEp %>% slice(-trainIndex)
    
    Data.GpE.Learn <- Data.GpE %>% slice(trainIndex)
    Data.GpE.Test <- Data.GpE %>% slice(-trainIndex)
    
    
    existing.pop.Learn=sum(c(length(unique(Data.GGp.Learn$pop)),length(unique(Data.GEp.Learn$pop)),length(unique(Data.EEp.Learn$pop)),length(unique(Data.GpEp.Learn$pop)),length(unique(Data.GpE.Learn$pop)))) 
    
    existing.pop.Test=sum(c(length(unique(Data.GGp.Test$pop)),length(unique(Data.GEp.Test$pop)),length(unique(Data.EEp.Test$pop)),length(unique(Data.GpEp.Test$pop)),length(unique(Data.GpE.Test$pop))))  
  }
  
  #####
  # Construction of the pop of t.values for each series Learn and Test
  for (i in 1:length(List.names.final)){
    eval(parse(text=paste0("Pop.",List.names.final[i],".Learn=Pop.create(Data.",List.names.final[i],".Learn)")))
    eval(parse(text=paste0("Pop.",List.names.final[i],".Test=Pop.create(Data.",List.names.final[i],".Test)")))
  }
  
  
  #####
  # Bootstrap: ccontruction de DataLearn et DataTest respecting the proportion of configuration
  
  DataLearn <- c()
  DataTest <- c()
  DataLearn <- Boot("Learn",Z.trunc.final.code,NbSimLearn)
  DataTest <- Boot("Test",Z.trunc.final.code,NbSimTest)
  saveRDS(DataLearn, file = paste0(file_path_Results,"DataLearn_",b,significance.level, offset, GE, number.pop,".rds"))
  saveRDS(DataTest, file = paste0(file_path_Results,"DataTest_",b,significance.level, offset, GE, number.pop,".rds"))
  
  Res.pred <- PredRule_ET_ErrorTest(DataLearn,DataTest,b,Nbconfig)
  saveRDS(Res.pred, file = paste0(file_path_Results,"Res.pred_",b,significance.level, offset, GE, number.pop,".rds"))
  error.test.4.methods[b,] <- Res.pred$err.tot %>% unlist()
}

type.dataset = "Learn"
Z.trunc.code = Z.trunc.final.code
NbSim = NbSimLearn

B=20
error.test.4.methods <- matrix(NA,nrow=B,ncol=4)

for (b in 1:B){
  Res.pred <- readRDS(paste0(file_path_Results,"Res.pred_",b,significance.level, offset, GE, number.pop,".rds"))
  error.test.4.methods[b,] <- Res.pred$err.tot %>% unlist()
}
colnames(error.test.4.methods) <- c("lda","cart","knn","rf")

boxplot(error.test.4.methods)

print(colMeans(error.test.4.methods))
print(apply(error.test.4.methods,2,sd))
print(apply(error.test.4.methods,2,which.min))

FinalPred <- readRDS(paste0(file_path_Results,"modrf_b",b = 4,significance.level, offset, GE, number.pop,".rds"))

# diff prob data  --------------------------------------------------------
Nbconfig <- nrow(Z.trunc.final)

set.seed(1)
error.test <- rep(NA, B)
for (b in 3:B){
  
  ######
  # Vfold 
  existing.pop.Learn <- 0
  existing.pop.Test <- 0
  
  while ((existing.pop.Learn!=15) & (existing.pop.Test !=15)){
    trainIndex<-createDataPartition(1:nrow(Data.GGp) , p=0.8, list = FALSE)
    
    Data.GGp.Learn <- Data.GGp %>% slice(trainIndex)
    Data.GGp.Test <- Data.GGp %>% slice(-trainIndex)
    
    Data.EEp.Learn <- Data.EEp %>% slice(trainIndex)
    Data.EEp.Test <- Data.EEp %>% slice(-trainIndex)
    
    Data.GEp.Learn <- Data.GEp %>% slice(trainIndex)
    Data.GEp.Test <- Data.GEp %>% slice(-trainIndex)
    
    Data.GpEp.Learn <- Data.GpEp %>% slice(trainIndex)
    Data.GpEp.Test <- Data.GpEp %>% slice(-trainIndex)
    
    Data.GpE.Learn <- Data.GpE %>% slice(trainIndex)
    Data.GpE.Test <- Data.GpE %>% slice(-trainIndex)
    
    
    existing.pop.Learn=sum(c(length(unique(Data.GGp.Learn$pop)),length(unique(Data.GEp.Learn$pop)),length(unique(Data.EEp.Learn$pop)),length(unique(Data.GpEp.Learn$pop)),length(unique(Data.GpE.Learn$pop)))) 
    
    existing.pop.Test=sum(c(length(unique(Data.GGp.Test$pop)),length(unique(Data.GEp.Test$pop)),length(unique(Data.EEp.Test$pop)),length(unique(Data.GpEp.Test$pop)),length(unique(Data.GpE.Test$pop))))  
  }
  
  #####
  # Construction of the pop of t.values for each series Learn and Test
  for (i in 1:length(List.names.final)){
    eval(parse(text=paste0("Pop.",List.names.final[i],".Learn=Pop.create(Data.",List.names.final[i],".Learn)")))
    eval(parse(text=paste0("Pop.",List.names.final[i],".Test=Pop.create(Data.",List.names.final[i],".Test)")))
  }
  
  
  #####
  # Bootstrap: ccontruction de DataLearn et DataTest respecting the proportion of configuration
  
  DataLearn <- c()
  DataTest <- c()
  DataLearn <- Boot1("Learn",Z.trunc.final.code,NbSim)
  DataTest <- Boot1("Test",Z.trunc.final.code,NbSim)
  saveRDS(DataLearn, file = paste0(file_path_Results,"DataLearn_",b,significance.level, offset, GE, number.pop,".rds"))
  saveRDS(DataTest, file = paste0(file_path_Results,"DataTest_",b,significance.level, offset, GE, number.pop,".rds"))
  
  Res.pred1 <- PredRule_RF(DataLearn,DataTest,b,Nbconfig)
  saveRDS(Res.pred1, file = paste0(file_path_Results,"Res.pred_",b,significance.level, offset, GE, number.pop,".rds"))
  error.test[b] <- Res.pred1$err.tot 
}
# read he best predictive rule
tot.err <- rep(NA, B) 
for (i in c(1:B)) {
  r <- readRDS(paste0(file_path_Results,"Res.pred_",b = i,significance.level, offset, GE, number.pop,".rds"))
  tot.err[i] = r$err.tot
}
FinalPred <- readRDS(paste0(file_path_Results,"modrf_b",b = which.min(tot.err ),significance.level, offset, GE, number.pop,".rds"))

# apply the best rule to the real data  -----------------------------------

# check cases in the table
Thresh <- significance.level
p.values.i=c()
truth.vec.i=c()
Z.truth.i <- c()
for (i in 1:nrow(Data.Res.Test)){
  a <- c()
  p.values.i <- 2*pnorm(-abs(as.numeric(Data.Res.Test[i,paste0("t", List.names.final)])))
  a<- ifelse(p.values.i<Thresh,1,0)*sign(as.numeric(Data.Res.Test[i,paste0("t", List.names.final)]))
  truth.vec.i <- rbind(truth.vec.i,a)
  Z.truth.i <- c(Z.truth.i,config.list.final[which(duplicated2(rbind(Z.trunc.final[,-rm.ind],a)))[1]])
}

pred.truth <- as.data.frame(cbind(truth.vec.i,Z.truth.i))
colnames(pred.truth) <- c("code.GGp", "code.GEp" , "code.EEp", "code.GpEp","code.GpE","Z.truth")

Data.Res.Test <- cbind(Data.Res.Test,pred.truth)
RealData.x <- Data.Res.Test[,colnames(Data.Res.Test) %in% c("tGGp","tGEp","tEEp", "tGpEp","tGpE")]
colnames(RealData.x) <- List.names.final

RealData.predy <- predict(FinalPred,newdata=RealData.x) 
FinalTable <- cbind(Data.Res.Test,config.list.final[RealData.predy])
colnames(FinalTable)[which(colnames(FinalTable)=="config.list.final[RealData.predy]")] <- "pred.y"

List.main <- unique(FinalTable$main)
FinalTable$w <- 1/FinalTable$distance

Post.Prob.List <- list()
Nb.main.break=1
for (i in List.main){
  Data.tmp1 <- FinalTable %>% dplyr::filter(main==i)
  List.break.i <-unique(Data.tmp1$brp) 
  for (j in List.break.i){
    Data.tmp2 <- Data.tmp1 %>% dplyr::filter(brp==j)
    A <- c(i,j)
    Post.Prob <- c()
    Config.Pred.Post <- c()
    if (nrow(Data.tmp2) >1){
      unique.pred <-unique(Data.tmp2$pred.y) 
      for (l in 1:length(unique.pred)){
        Post.Prob[l] <- sum((Data.tmp2$pred.y==unique.pred[l])*Data.tmp2$w)/sum(Data.tmp2$w)
      }
      names(Post.Prob) <-  unique.pred
      Config.Pred.Post <- unique.pred[which.max(Post.Prob)]
    } else {
      Post.Prob[1] <- 1
      names(Post.Prob) <-  Data.tmp2$pred.y
      Config.Pred.Post <- Data.tmp2$pred.y
    }
    
    Post.Prob.List[[Nb.main.break]] <- list(MainBreak=A,PostProb=Post.Prob,Config.Pred.Post=Config.Pred.Post)
    
    
    Nb.main.break=Nb.main.break+1
  }
}

Post.Prob.List
save(FinalTable, file = paste0(path_restest, "Final.Table", significance.level, offset, GE, number.pop, ".RData"))
save(Post.Prob.List, file = paste0(path_restest, "Post.Prob.List", significance.level, offset, GE, number.pop, ".RData"))




# RUN OTHER CLASSIFIER ----------------------------------------------------

Nbconfig <- nrow(Z.trunc.final)

set.seed(1)
for (b in 1:B){
  Res.pred1 <-c()
  DataLearn <- readRDS(file = paste0(file_path_Results,"DataLearn_",b,significance.level, offset, GE, number.pop,".rds"))
  DataTest <- readRDS(file = paste0(file_path_Results,"DataTest_",b,significance.level, offset, GE, number.pop,".rds"))
  
  # Res.pred1 <- PredRule_LDA(DataLearn,DataTest,b,Nbconfig)
  # Res.pred1 <- PredRule_CART(DataLearn,DataTest,b,Nbconfig)
  Res.pred1 <- PredRule_knn(DataLearn,DataTest,b,Nbconfig)
  saveRDS(Res.pred1, file = paste0(file_path_Results,"Res.pred_",b,significance.level, offset, GE, number.pop,"knn.rds"))
}
tot.err <- rep(NA, B) 
for (i in c(1:B)) {
  r <- readRDS(paste0(file_path_Results,"Res.pred_",b = i,significance.level, offset, GE, number.pop,"CART.rds"))
  tot.err[i] = r$err.tot
}
FinalPred <- readRDS(paste0(file_path_Results,"Res.pred_",b = which.min(tot.err ),significance.level, offset, GE, number.pop,"LDA.rds"))
