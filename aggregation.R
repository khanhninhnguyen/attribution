# test the aggregation

# Run aggregation under changes -------------------------------------------

five = get(load(file = "/home/knguyen/Documents/PhD/paper/attribution/result/attribution/Final.Table.RData"))
version = "low_dist"

FinalTable = five[which(five$distance>50),]
List.main <- unique(FinalTable$main)

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

save(Post.Prob.List, file = paste0(path_results, "attribution/Post.Prob.List", version, ".RData"))

# make the group plot  ----------------------------------------------------

plot_group <- function(dat.p){
  G = c(1,15)
  E = c(8,22)
  G.G = c(5,18)
  G.E = c(29, 34)
  G.E. = c(2,3,16,17)
  E.E = c(9,10,23,24)
  G.E.E = c(30,31,35,36)
  G1.E.E = c(12,13,14,25,27,28)
  three = c(4,6,7, 11:14,19,20,21, 25:28, 30:33, 35,36)
  
  chosen.w = "last.pre"
  dat.p$g = "others"
  dat.p$g[which(dat.p[,chosen.w] %in% G == TRUE)] = "G"
  dat.p$g[which(dat.p[,chosen.w] %in% E == TRUE)] = "E"
  dat.p$g[which(dat.p[,chosen.w] %in% G.G == TRUE)] = "G&G'"
  dat.p$g[which(dat.p[,chosen.w] %in% G.E == TRUE)] = "G&E"
  dat.p$g[which(dat.p[,chosen.w] %in% G.E. == TRUE)] = "G&E'"
  dat.p$g[which(dat.p[,chosen.w] %in% E.E == TRUE)] = "E&E'"
  dat.p$g[which(dat.p[,chosen.w] %in% G.E.E == TRUE)] = "G&E&E'"
  dat.p$g[which(dat.p[,chosen.w] %in% G1.E.E == TRUE)] = "G'&E&E'"
  
  dat.p$c =1 
  data.plot = aggregate(c~last.pre+g, dat.p, sum)
  colnames(data.plot)[1] = "config"
  data.plot$config = as.factor(data.plot$config)
  
  p <- ggplot(data.plot, aes(fill=config, y=c, x=g)) + 
    geom_bar(position="stack", stat="identity")+theme_bw()+
    geom_text(aes(label = config),
              colour = "black",  size=1.8,
              position = position_stack(vjust = 0.5)) +
    labs(y ="Count", x = "Group",  fill='Configuration') + 
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=4),
          axis.title = element_text(size = 6), legend.key.size = unit(0.2, "cm"), 
          plot.tag = element_text(size = 6),legend.title=element_text(size=5),legend.box.spacing = unit(0, "pt"),
          plot.margin = rep(unit(0,"null"),4))
  return(p)
}

prob.t = get(load(file = paste0(path_results, "attribution/Post.Prob.List", version, ".RData")))
n = length(prob.t)
last.name = sapply(c(1:n), function(x) prob.t[[x]]$MainBreak[1])
last.brp = sapply(c(1:n), function(x) prob.t[[x]]$MainBreak[2])
last.pre = sapply(c(1:n), function(x) prob.t[[x]]$Config.Pred.Post)

dat.p1 = data.frame( name = last.name, brp = last.brp, last.pre)
p = plot_group(dat.p = dat.p1)

ggsave(paste0(path_results,"attribution/config_dist", version = "prob", ".jpg" ), plot = p, width = 8.8, height = 6, units = "cm", dpi = 600)

# agrregate with probability ----------------------------------------------
FinalTable = five

prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
          0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
          0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
          0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)
keep.config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)
prob <- prob[keep.config]

List.main <- unique(FinalTable$main)
FinalTable$p = sapply(c(1:nrow(FinalTable)), function(x) prob[FinalTable$pred.y[x]])
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
        Post.Prob[l] <- sum((Data.tmp2$pred.y==unique.pred[l])*Data.tmp2$p)/sum(Data.tmp2$p)
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

save(Post.Prob.List, file = paste0(path_results, "attribution/Post.Prob.List", version, "prob.RData"))






