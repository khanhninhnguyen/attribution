# study the contradiction
#  create the truncated table  ------
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

keep.config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)
Y <- Y[keep.config,]
Z <- Z[keep.config,]
Z.trunc <- Z.trunc[keep.config,]
save(Z.trunc, file = paste0(path_results, "attribution0/truncated.table.RData"))

# check the contradiction function ---------------
check_contradict <- function(y, table.selected){
  names(y) = NULL
  colnames(table.selected) = NULL
  res = sapply(c(1:nrow(table.selected)), function(x) identical(unlist(table.selected[x,]), y))
  ind.o = which(res==TRUE)
  out = ifelse(length(ind.o)>0, ind.o, 0)
  return(out)
}

# convert to coded table 
# read test result 
Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))
# result for all G-E

convert_coded <- function(x, significance.level, length.x){
  sapply(c(1:length(x)), function(i) ifelse( (2*pnorm(-abs(as.numeric(x[i])))) < significance.level[i], 1*sign(x[i]), 0)) 
}

#check contradiction for different level of significance 
trunc.table = get(load(file = paste0(path_results, "attribution0/truncated.table.RData")))

sig.list = seq(0.00, 0.1,0.002)[-1]
nb.contradicted = rep(NA, length(sig.list))
tot.res = data.frame(matrix(NA, ncol = length(sig.list), nrow = 494))
for (j in 1:length(sig.list)) {
  Total.coded = data.frame(matrix(NA, ncol = 6, nrow = nrow(Total.res)))
  for (i in c(1:nrow(Total.res))) {
    if(Total.res$distance[i] <25){
      sig.level = rep(0.002, 6)
    }else{
      sig.level =  c(0.002, 0.05, 0.05, 0.05, 0.05, 0.002, 0.05)
    }
    case.i = Total.res[i, c(paste0("t", list.name.test[1:6]))]
    Total.coded[i,] = convert_coded(case.i, significance.level =sig.level)
  }
  contra = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded[x,]), trunc.table))
  tot.res[,j] = contra
  nb.contradicted[j] = length(which(contra == 0))
  a = data.frame(table(contra))
  a$per = paste0(a$contra, "-", round((a$Freq*100/494), digits = 1), "%")
  png(file = paste0(path_results,"attribution/pie3",".png"))
  pie(a$Freq,a$per)
  dev.off()
}


plot(sig.list, nb.contradicted, xlab = "Significance level", ylab = "No. contradicted cases")

# plot the evolution of configurations 
list.config.test = unique(tot.res$X50)[-c(7, 22:24)]
count.all = data.frame(t(sapply(c(1:length(list.config.test)), function(x) sapply(c(1:50), function(y) length(which(tot.res[,y] == list.config.test[x]))))))
colnames(count.all)= sig.list*100
count.all$config = paste0("config.", list.config.test)

dat.p = reshape2::melt(count.all, id = "config")

ggplot(dat.p, aes( x = variable, y = value, col = config))+ theme_bw()+ geom_point(size=0.7)+
  xlab("Significance level (%)") + scale_x_discrete(breaks = seq(0, 10, 1))+
  ylab("Count")

# investigate cases, where results are changing to another configuration -----------
ind.inv = which(tot.res$X1 ==1 & tot.res$X25 != 1 & tot.res$X25 != 0)

jump.inv = Total.res[ind.inv,c(paste0("jump", list.name.test))]

dat.p = reshape2::melt(jump.inv)
ggplot(dat.p, aes(x = variable,y = value))+ theme_bw()+ geom_boxplot(size=0.7)+
  xlab("series") + 
  ylab("jump")

t.inv = Total.res[ind.inv,c(paste0("t", list.name.test))]

dat.p = reshape2::melt(t.inv)
ggplot(dat.p, aes(x = variable,y = value))+ theme_bw()+ 
  geom_boxplot(size=0.7)+ geom_hline(yintercept = 1.96)
  xlab("series") + 
  ylab("jump")








  

# barplot of distribution of configurations at different significant level --------

trunc.table = get(load(file = paste0(path_results, "attribution0/truncated.table.RData")))

sig.list = c(0.05, 0.01, 0.002, 4e-05)
nb.contradicted = rep(NA, length(sig.list))
tot.res = data.frame(matrix(NA, ncol = length(sig.list), nrow = 494))
for (j in 1:length(sig.list)) {
  Total.coded = data.frame(matrix(NA, ncol = 6, nrow = nrow(Total.res)))
  for (i in c(1:nrow(Total.res))) {
    sig.level = rep(sig.list[j], 6)
    # if(Total.res$distance[i] <25){
    #   sig.level = rep(sig.list, 6)
    # }else{
    #   sig.level =  c(0.002, 0.05, 0.05, 0.05, 0.05, 0.002, 0.05)
    # }
    case.i = Total.res[i, c(paste0("t", list.name.test[1:6]))]
    Total.coded[i,] = convert_coded(case.i, significance.level =sig.level)
  }
  contra = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded[x,]), trunc.table))
  tot.res[,j] = contra
  nb.contradicted[j] = length(which(contra == 0))
  a = data.frame(table(contra))
  b = data.frame(contra = factor(0:38), Freq =0)
  data.p = left_join(b, a, by ="contra")
  data.p$Freq.y[which(is.na(data.p$Freq.y)==TRUE)] = 0
  p = ggplot(data.frame(data.p),aes(x = contra, Freq.y*100/494))+
    geom_bar(stat="identity", width = 0.7)+ 
    theme_bw()+ 
    ylab("Percentage (%)")+
    xlab("Configuration")+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5))
  
  ggsave(paste0(path_results,"attribution0/config.table", sig.list[j], ".jpg" ), plot = p, width = 12, height = 4, units = "cm", dpi = 300)
  
}

