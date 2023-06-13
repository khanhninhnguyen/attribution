# test the aggregation

# Run aggregation under changes -------------------------------------------

five = get(load(file = "/home/knguyen/Documents/PhD/paper/attribution/result/attribution/Final.Table.RData"))
version = "low_dist"

FinalTable = five[which(five$distance>50),]

aggre1 <- function(FinalTable){
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
  
  return(Post.Prob.List)
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
  
  cols <- c("aquamarine","aquamarine1", "aquamarine2", "aquamarine3","aquamarine4",
            "cyan", "cyan1", "cyan2", "cyan3", "cyan4",
            "brown", "brown1", "brown2", "brown3", "brown4",
            "coral", "coral1", "coral2", "coral3", "coral4",
            "gold", "gold1", "gold2", "gold3", "gold4",
            paste0("gray", c(1:13)))
  
  dat.p$c =1 
  data.plot = aggregate(c~last.pre+g, dat.p, sum)
  colnames(data.plot)[1] = "config"
  data.plot$config = as.factor(data.plot$config)
  max.group = aggregate(c~g, data.plot, sum)
  max.y = ceiling(max(max.group$c/10))*10
  
  p <- ggplot(data.plot, aes(fill=config, y=c, x=g)) + 
    geom_bar(position="stack", stat="identity")+theme_bw()+
    geom_text(aes(label = config),
              colour = "black",  size=1.8,
              position = position_stack(vjust = 0.5)) +
    scale_y_continuous(breaks = seq(0,  max.y , 10), limit = c(0, max.y))+
    labs(y ="Count", x = "Group",  fill='Configuration') + 
    scale_fill_manual(breaks = c(1:38),
                         values= cols)
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

ggsave(paste0(path_results,"attribution/config_dist", version = "unequal.dist", ".jpg" ), plot = p, width = 8.8, height = 6, units = "cm", dpi = 600)

# agrregate with probability ----------------------------------------------
FinalTable = five

aggre2 <- function(FinalTable){
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
  
  return(Post.Prob.List)
}

save(Post.Prob.List, file = paste0(path_results, "attribution/Post.Prob.List", version, "prob.RData"))

# comparison 2 methods ----------------------------------------------------

l1 = aggre1(FinalTable)
l2 = aggre2(FinalTable)

n = length(l1)
dat1 = data.frame(name = sapply(c(1:n), function(x) l1[[x]]$MainBreak[1]), 
                  brp = sapply(c(1:n), function(x) l1[[x]]$MainBreak[2]),
                  last.pre = sapply(c(1:n), function(x) l1[[x]]$Config.Pred.Post))

dat2 = data.frame(name = sapply(c(1:n), function(x) l2[[x]]$MainBreak[1]), 
                  brp = sapply(c(1:n), function(x) l2[[x]]$MainBreak[2]),
                  last.pre = sapply(c(1:n), function(x) l2[[x]]$Config.Pred.Post))

both = full_join(dat1, dat2, by = c("name", "brp")) 
colnames(both)[c(3:4)] <- c("equal", "unequal")

same = both[which(both$equal== both$unequal),]
dif = both[which(both$equal!= both$unequal),]

# why they give the same small config? (1 nearby?)
small.sim = both[which(both$equal== both$unequal),] %>% 
  filter(!unequal %in% c(1, 10, 15, 23))
colnames(small.sim)[1] = "main"
full.small = inner_join(small.sim, FinalTable, by = c("main", "brp"))

data.p = full.small[,c("main", "brp", "nearby", "Z.truth", "pred.y")] %>%
  mutate(case = paste0(data.p$main, data.p$brp)) %>%
  mutate(contradict = as.factor(sapply(c(1:nrow(full.small)), function(x) ifelse(is.na(Z.truth[x]), 1, 0))))

ggplot(data = data.p, aes(x = case, fill = contradict)) + theme_bw() +
  geom_bar(width = 0.7)+
  theme(axis.text.x = element_text(angle = 90))

# what make differences?
colnames(dif)[1] = "main"

statis = inner_join(dif, FinalTable, by = c("main", "brp"))[,-c(12:13, 15:19, 22)]
statist = inner_join(statis, one[,c(1:3,18:19)], by = c("main", "brp", "nearby"))
colnames(statist)[13:16] <- c("contra.5", "pred.5", "contra.1", "pred.1") 

statistic = statist[c(1,2,5, 6:16, 3, 4)]

write.table(format(statistic, digits = 2), file = paste0(path_results, "attribution/print.cases.analyze.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# count.config = data.frame(table(e$pred.y))
# count.config.test = data.frame(table(e$Z.truth))
# count.configs = full_join(count.config, count.config.test, by = "Var1") %>% 
#   replace(is.na(.), 0) %>% 
#   mutate(truth = Freq.y) %>% 
#   mutate(pred = Freq.x - Freq.y) %>% 
#   dplyr::select(-c(Freq.x,Freq.y)) %>% 
#   mutate_at(vars(Var1), list(factor)) %>% 
#   rename("config" = "Var1") %>% 
#   reshape2::melt(id="config") %>% 
#   arrange(factor(variable, levels= c("pred", "truth"))) %>% 
#   group_by(config) %>%
#   mutate(lab_ypos = cumsum(value) - 0.5 * value)  %>% 
#   mutate_all(~na_if(., 0))
# ggplot(data = count.configs, aes(x = config, y = value)) + theme_bw()+
#   geom_col(aes(fill = variable), width = 0.7)+
#   geom_text(aes(y = lab_ypos, label = value, group =variable), color = "black") +
#   scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+
#   scale_fill_manual(values = c("#0073C2FF", "#EFC000FF"))


# PLOT FOR PAPER ----------------------------------------------------------
FinalTable = five  
l1 = aggre1(FinalTable)

n = length(l1)
dat1 = data.frame(name = sapply(c(1:n), function(x) l1[[x]]$MainBreak[1]), 
                  brp = sapply(c(1:n), function(x) l1[[x]]$MainBreak[2]),
                  last.pre = sapply(c(1:n), function(x) l1[[x]]$Config.Pred.Post))

p = plot_group(dat.p = dat1)

ggsave(paste0(path_results,"attribution/config_dist", version = "unequal.dist", ".jpg" ), plot = p, width = 8.8, height = 6, units = "cm", dpi = 600)







