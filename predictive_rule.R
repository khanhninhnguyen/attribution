predictive_rule <- function(path_results, significance.level, offset, GE, number.pop, R){
  if (!require(devtools)) install.packages("devtools")
  devtools::install_github("yanlinlin82/ggvenn")
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
  
  path_restest <- paste0(path_results,"attribution/predictive_rule/")
  file_path_Results=paste0(path_results,'attribution/predictive_rule/')
  
  ####---------
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
  # knitr::kable(head(Z.trunc))
  
  prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
            0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
            0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
            0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
            0.00225,0.00225,0.0405,0.010125,0.010125,0.00225,0.010125,0.010125,0.00225,0.00225,
            0.0405,0.00225,0.010125,0.00225,0.010125,0.010125,0.00225,0.010125)
  
  keep.config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)
  Y <- Y[keep.config,]
  Z <- Z[keep.config,]
  Z.trunc <- Z.trunc[keep.config,]
  num.conf <- 1:length(keep.config)
  row.names(Z.trunc) <- num.conf
  
  prob <- prob[keep.config]
  prob <- prob/sum(prob)
  
  sum(duplicated(Z.trunc))
  
  duplicated2 <- function(x){
    auxi <- do.call("paste",c(x,sep="_"))
    y <- table(auxi)
    y <- y>1
    res <- y[match(auxi,names(y))]
    dimnames(res) <- NULL
    return(res)
  }
  
  sum(duplicated2(Z.trunc[,-1]))/2
  rg.duplicate <-which(duplicated2(Z.trunc[,-1])) 
  cbind(Z.trunc[rg.duplicate,],prob[rg.duplicate])
  
  # original configurations
  cbind(Y[c(7,12,19,28),],c(7,12,19,28))
  # Removing of the configurations 19 and 28
  Z.trunc.final <- Z.trunc[-c(12,28),]
  dim(Z.trunc.final)
  dim(Z.trunc.final[Z.trunc.final$GE==1,])
  dim(Z.trunc.final[Z.trunc.final$GE==-1,])
  
  # New table and probabilities
  Z.trunc.final <- Z.trunc.final[,-1]
  prob.final <- prob[-c(12,28)]
  prob.final <- prob.final/sum(prob.final)
  
  # The coded Z.trunc in terms of population
  Z.trunc.final.code <- Z.trunc.final 
  Z.trunc.final.code[Z.trunc.final==1]=3
  Z.trunc.final.code[Z.trunc.final==-1]=2
  Z.trunc.final.code[Z.trunc.final==0]=1
  
  head(Z.trunc.final.code)
  config.list <- 1:38
  config.list.final <- config.list[-c(12,28)]
  saveRDS(Z.trunc.final.code , file = paste0(file_path_Results,"List_config.rds"))
  
  name.results <- paste0(path_restest,"FGLS_on_real_data_t.txt")
  Data.Res.Test <- read.table(name.results,header = TRUE, stringsAsFactors = FALSE)
  colnames(Data.Res.Test)[4:9] <- paste0("t", List.names.tot)
  
  NumDetected.per.station <- Data.Res.Test %>% group_by(main) %>% 
    summarise(ldetected.per.station=length(unique(brp))) %>% dplyr::select(ldetected.per.station)
  
  NumNearby.per.station.detection <- Data.Res.Test %>% group_by(main,brp) %>% 
    count() %>% as.data.frame() %>% dplyr::select(n)
  
  if(GE == 0){
    List.names.final <- List.names.tot[!(List.names.tot %in% "GE")]
  }
  
  Keep.Data <- function(name.series){
    Thresh <- significance.level
    names.t.series <- paste0("t",name.series)
    colnames.info <- c("main","brp","nearby")
    Data.Res <- Data.Res.Test[colnames(Data.Res.Test) %in% c(colnames.info,names.t.series)] 
    Data.Res$t <- Data.Res[,which(colnames(Data.Res) %in% names.t.series)]
    Data.Res <- Data.Res[,-which(colnames(Data.Res) %in% names.t.series)]
    Data.Res <- Data.Res %>% dplyr::mutate(p=2*pnorm(-abs(t))) %>%
      mutate(signif=ifelse(p<Thresh,1,0),
             t.sign=sign(t),
             pop=(signif==0)*1+((signif==1) & (t.sign==-1))*2+ ((signif==1) & (t.sign==1))*3,
             pop.code=(signif==0)*0+((signif==1) & (t.sign==-1))*-1+ ((signif==1) & (t.sign==1))*1)
    return(Data.Res)
  }
  
  for (i in  1:length(List.names.final)){
    eval(parse(text=paste0("Data.",List.names.final[i],"=Keep.Data(List.names.final[i])")))  
  }
  
  Nbconfig <- nrow(Z.trunc.final)
  NbSim <- R*Nbconfig #in order to have at least 5 samples for each configurations. 
  B <- 20
  NbSimLearn <- NbSim*0.8
  NbSimTest <- NbSim*0.2
  
  Pop.create <- function(dataset){
    Pop1 <- dataset %>% dplyr::filter(pop==1) %>% dplyr::select(t)  %>% unlist()
    Pop2 <- dataset %>% dplyr::filter(pop==2) %>% dplyr::select(t)  %>% unlist()
    Pop3 <- dataset %>% dplyr::filter(pop==3) %>% dplyr::select(t)  %>% unlist()
    Pop=list(Pop1,Pop2,Pop3)
    return(Pop)
  }
  
  
  Boot <- function(type.dataset,Z.trunc.code,NbSim, seed){
    config=c()
    Data.res <- c()
    Nbconfig <- nrow(Z.trunc.code)
    
    for (c in 1:Nbconfig){
      config.code <- Z.trunc.code[c,]
      N <- eval(parse(text=paste0("NbSim",type.dataset)))
      if (type.dataset=="Learn"){
        Nb.config.c <- R*0.8
      } else {Nb.config.c <- R*0.2}
      config <- c(config,rep(c,Nb.config.c))
      
      data.c <- purrr::map(List.names.final,~{
        Pop <- c()
        code.names.series <- config.code[which(colnames(Z.trunc.code) %in% .x)]
        eval(parse(text=paste0("Pop=Pop.",.x,".",type.dataset,"[[",code.names.series,"]]")))
        set.seed(seed[c])
        res.t=sample(Pop,Nb.config.c , replace = TRUE)
        return(res.t)
      }) %>% bind_cols() %>% as.data.frame()
      colnames(data.c)=List.names.final
      Data.res = rbind(Data.res,data.c)
    }
    
    Data.res <-Data.res %>% mutate(config=config)
    Data.res$config <- as.factor(Data.res$config)
    return(Data.res)
    
  }
  
  PredRule_ET_ErrorTest <- function(DataLearn,DataTest,b,Nbconfig){
    
    config.tot <- 1:Nbconfig
    ##############
    # LDA
    ##############
    modlda<-lda(config~ ., data=DataLearn, CV=FALSE)
    saveRDS(modlda, file = paste0(file_path_Results,'modlda_b',b,significance.level, offset, GE, number.pop,'.rds'))
    
    pred.lda <- predict(modlda,newdata=DataTest,na.action = na.pass)$class 
    err.lda <- mean(pred.lda!=DataTest$config)
    which.misclassif.cluster.lda <- which(pred.lda!=DataTest$config)
    misclassif.cluster.lda <- rbind(DataTest$config[which.misclassif.cluster.lda],pred.lda[which.misclassif.cluster.lda])
    never.pred.lda <- config.tot[!(config.tot %in% pred.lda)]
    
    #####
    #
    cvControl=trainControl(method="cv",number=10)
    #
    ##############
    # Cart 
    ##############
    # Tree.max <- rpart(config ~ ., data = DataLearn,method="class",control=rpart.control(minsplit=2,cp=1e-06))
    # alpha_opt <- Tree.max$cptable[which.min(Tree.max$cptable[,"xerror"]),"CP"]
    # modcart <- prune(Tree.max,cp=alpha_opt)
    # #saveRDS(modcart, file = paste0(file_path_Results,'modcart_b',b,'.rds'))
    # 
    # pred.cart <- predict(modcart,newdata=DataTest,na.action = na.pass,type="class") 
    # err.cart <- mean(pred.cart!=DataTest$config)
    # which.misclassif.cluster.cart <- which(pred.cart!=DataTest$config)
    # misclassif.cluster.cart <- rbind(DataTest$config[which.misclassif.cluster.cart],pred.cart[which.misclassif.cluster.cart])
    # never.pred.cart <- config.tot[!(config.tot %in% pred.cart)]
    # err.cart
    
    modcart = caret::train(config~ ., data = DataLearn, method = "rpart", tuneLength = 10,trControl = cvControl,metric = "Accuracy")
    saveRDS(modcart, file = paste0(file_path_Results,'modcart_b',b,significance.level, offset, GE, number.pop,'.rds'))
    
    pred.cart=predict(modcart, newdata = DataTest)
    err.cart <- mean(pred.cart!=DataTest$config)
    which.misclassif.cluster.cart <- which(pred.cart!=DataTest$config)
    misclassif.cluster.cart <- rbind(DataTest$config[which.misclassif.cluster.cart],pred.cart[which.misclassif.cluster.cart])
    never.pred.cart <- config.tot[!(config.tot %in% pred.cart)]
    err.cart
    
    
    ##############
    # knn
    ##############
    # grid <- c(1,2,3,5,10)
    # error.knn.i <- c()
    # error.knn.i <- purrr::map(grid,~{
    #   pred.knn <- knn(train=DataLearn[,1:5],test=DataTest[,1:5],cl=DataLearn$config,k=.x)
    #   error.knn <- mean(pred.knn!=DataTest$config)
    #   return(error.knn)
    # }) %>% unlist()
    # 
    # k.opt=grid[which.min(error.knn.i)]
    # pred.knn <- knn(train=DataLearn[,1:5],test=DataTest[,1:5],cl=DataLearn$config,k=k.opt)
    modknn = caret::train(config~ ., data = DataLearn, method = "knn", tuneLength = 10,trControl = cvControl,metric = "Accuracy")
    saveRDS(modknn, file = paste0(file_path_Results,'modknn_b',b,significance.level, offset, GE, number.pop,'.rds'))
    
    pred.knn=predict(modknn, newdata = DataTest)
    err.knn <- mean(pred.knn!=DataTest$config)
    which.misclassif.cluster.knn <- which(pred.knn!=DataTest$config)
    misclassif.cluster.knn <- rbind(DataTest$config[which.misclassif.cluster.knn],pred.knn[which.misclassif.cluster.knn])
    never.pred.knn <- config.tot[!(config.tot %in% pred.knn)]
    err.knn
    
    
    
    
    ##############
    # RF
    ##############
    
    # Testé avec le paramètre par défaut pour m=sqrt(p)
    # m=seq(1,4,1)
    # error.rf.i <- c()
    # error.rf.i <- purrr::map(m,~{
    #   foret <- randomForest(config~.,data=DataLearn,mtry=.x,ntree=500)
    #   pred.rf <- predict(foret,newdata=DataTest,na.action = na.pass)
    #   error.rf <- mean(pred.rf!=DataTest$config)
    #   return(error.rf)
    # }) %>% unlist()
    # m.opt=m[which.min(error.rf.i)]
    # modrf <- randomForest(config~.,data=DataLearn,mtry=m.opt)
    # #saveRDS(modrf, file = paste0(file_path_Results,'modrf_b',b,'.rds'))
    # pred.rf <- predict(modrf,newdata=DataTest)
    # err.rf <- mean(pred.rf!=DataTest$config)
    # err.rf
    
    modrf = caret::train(config~ ., data = DataLearn, method = "rf", tuneLength = 4,trControl = cvControl,metric = "Accuracy", importance=TRUE)
    saveRDS(modrf, file = paste0(file_path_Results,'modrf_b',b,significance.level, offset, GE, number.pop,'.rds'))
    pred.rf=predict(modrf, newdata = DataTest)
    err.rf <- mean(pred.rf!=DataTest$config)
    
    which.misclassif.cluster.rf <- which(pred.rf!=DataTest$config)
    misclassif.cluster.rf <- rbind(DataTest$config[which.misclassif.cluster.rf],pred.rf[which.misclassif.cluster.rf])
    never.pred.rf <- config.tot[!(config.tot %in% pred.rf)]
    
    # return
    err.tot <- list(lda=err.lda,cart=err.cart,knn=err.knn,rf=err.rf)
    
    which.misclassif.tot <- list(lda=which.misclassif.cluster.lda,cart=which.misclassif.cluster.cart,knn=which.misclassif.cluster.knn,rf=which.misclassif.cluster.rf)
    
    misclassif.tot <- list(lda=misclassif.cluster.lda,cart=misclassif.cluster.cart,knn=misclassif.cluster.knn,rf=misclassif.cluster.rf)
    
    never.pred.tot <- list(lda=never.pred.lda,cart=never.pred.cart,knn=never.pred.knn,rf=never.pred.rf)
    
    return(list(err.tot=err.tot,which.misclassif.tot=which.misclassif.tot,misclassif.tot=misclassif.tot,never.pred.tot=never.pred.tot,modlda=modlda,modcart=modcart,modknn=modknn,modrf=modrf))
    #err.tot=c(err.lda,err.cart,err.knn,err.rf)
    #return(err.tot)
    
  }
  
  error.test.4.methods <- matrix(NA,nrow=B,ncol=4)
  
  for (b in 1:B){
    
    set.seed(b)
    seed.learn = sample(c(1:10000), 36, replace = FALSE)
    set.seed(b+20)
    seed.test = sample(c(1:10000), 36, replace = FALSE)
    ######
    # Vfold 
    existing.pop.Learn <- 0
    existing.pop.Test <- 0
    
    while ((existing.pop.Learn!=15) & (existing.pop.Test !=15)){
      trainIndex<-createDataPartition(1:nrow(Data.GGp) , p=0.8, list = FALSE)
      
      Data.GGp.Learn <- Data.GGp %>% slice(trainIndex)
      Data.GGp.Test <- Data.GGp %>% slice(-trainIndex)
      
      Data.GEp.Learn <- Data.GEp %>% slice(trainIndex)
      Data.GEp.Test <- Data.GEp %>% slice(-trainIndex)
      
      Data.EEp.Learn <- Data.EEp %>% slice(trainIndex)
      Data.EEp.Test <- Data.EEp %>% slice(-trainIndex)
      
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
    DataLearn <- Boot("Learn",Z.trunc.final.code,NbSimLearn,seed.learn)
    DataTest <- Boot("Test",Z.trunc.final.code,NbSimTest, seed.test)
    saveRDS(DataLearn, file = paste0(file_path_Results,"DataLearn_",b,significance.level, offset, GE, number.pop,".rds"))
    saveRDS(DataTest, file = paste0(file_path_Results,"DataTest_",b,significance.level, offset, GE, number.pop,".rds"))
    
    Res.pred <- PredRule_ET_ErrorTest(DataLearn,DataTest,b,Nbconfig)
    saveRDS(Res.pred, file = paste0(file_path_Results,"Res.pred_",b,significance.level, offset, GE, number.pop,".rds"))
    error.test.4.methods[b,] <- Res.pred$err.tot %>% unlist()
  }
  colnames(error.test.4.methods) <- c("lda","cart","knn","rf")
  
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
  
  b = apply(error.test.4.methods,2,which.min)[4]
  FinalPred <- readRDS(paste0(file_path_Results,"modrf_b",b,significance.level, offset, GE, number.pop,".rds"))
  Thresh <- significance.level
  p.values.i=c()
  truth.vec.i=c()
  Z.truth.i <- c()
  for (i in 1:nrow(Data.Res.Test)){
    a <- c()
    p.values.i <- 2*pnorm(-abs(as.numeric(Data.Res.Test[i,5:9])))
    a<- ifelse(p.values.i<Thresh,1,0)*sign(as.numeric(Data.Res.Test[i,5:9]))
    truth.vec.i <- rbind(truth.vec.i,a)
    Z.truth.i <- c(Z.truth.i,config.list.final[which(duplicated2(rbind(Z.trunc.final,a)))[1]])
  }
  
  pred.truth <- as.data.frame(cbind(truth.vec.i,Z.truth.i))
  colnames(pred.truth) <- c("code.GGp",  "code.GEp" , "code.EEp" , "code.GpEp","code.GpE","Z.truth")
  
  Data.Res.Test <- cbind(Data.Res.Test,pred.truth)
  RealData.x <- Data.Res.Test[,colnames(Data.Res.Test) %in% c("tGGp","tGEp","tEEp","tGpEp","tGpE")]
  colnames(RealData.x) <- List.names.final
  
  RealData.predy <- predict(FinalPred,newdata=RealData.x) 
  FinalTable <- cbind(Data.Res.Test,config.list.final[RealData.predy])
  colnames(FinalTable)[which(colnames(FinalTable)=="config.list.final[RealData.predy]")] <- "pred.y"
  
  Res.List <- list()
  List.main <- unique(FinalTable$main)
  
  Nb.main.break <- 1
  
  for (i in List.main){
    Data.tmp1 <- FinalTable %>% dplyr::filter(main==i)
    List.break.i <-unique(Data.tmp1$brp) 
    for (j in List.break.i){
      Data.tmp2 <- Data.tmp1 %>% dplyr::filter(brp==j)
      A <- c(i,j)
      B <- Data.tmp2[,colnames(Data.tmp2) %in% c("tGE","tGGp" ,"tGEp","tEEp","tGpEp","tGpE")]
      C <- Data.tmp2[,colnames(Data.tmp2) %in% c("code.GGp","code.GpE","code.EEp","code.GpEp","code.GpE")]
      D <- Data.tmp2[,colnames(Data.tmp2) %in% c( "Z.truth" ,"pred.y", "distances")]
      Res.List[[Nb.main.break]] <- list(MainBreak=A,t.value=B,Code.res=C,Pred=D)
      Nb.main.break=Nb.main.break+1
    }
  }
  
  Res.List
  
  
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
  
  FinalTable$w <- FinalTable$n1+FinalTable$n2
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
  
  save(Post.Prob.List, file = paste0(path_restest, "Post.Prob.Listn", significance.level, offset, GE, number.pop, ".RData"))
  
}


# b=2
# FinalPred <- readRDS(paste0(file_path_Results,"modrf_b",b,".rds"))
# varImpPlot(FinalPred$finalModel)
# varImp(FinalPred)
# plot(varImp(FinalPred), top =5)
# importance <- randomForest::importance(FinalPred$finalModel,type=1)

# Application on the real data --------------
# Data.Res.Test <- read.table(name.results,sep="\t",header=TRUE,na.strings  ="NA")


# 
# 
# FinalPred <- readRDS(paste0(path_results,"attribution/","modrf_b",2,".rds"))
# importance <- randomForest::importance(FinalPred$finalModel,type=1)
# plot(importance)
# 
# 
# save(out,file = paste0(path_result,"validation/",nb_test,"-",criterion,"metacompa",screen.value,".RData"))
# 
# 

































