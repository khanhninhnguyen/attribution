# to investigate result of predictive rule
significance.level = 0.05
offset = 0
GE = 0
number.pop = 1
final.t = get(load(file = paste0(path_results, "attribution/predictive_rule/Final.Table", significance.level, offset, GE, number.pop, ".RData")))
prob.t = get(load(file = paste0(path_results, "attribution/predictive_rule/Post.Prob.List", significance.level, offset, GE, number.pop, ".RData")))
# prob.tn = get(load(file = paste0(path_results, "attribution/predictive_rule/Post.Prob.Listn", significance.level, offset, GE, number.pop, ".RData")))
tot = get(load(file = paste0(path_results,"attribution0/stats_test_real_data.RData")))
# old result fromm paper
final.t = get(load(file = "/home/knguyen/Documents/PhD/paper/attribution/result/attribution/Final.Table.RData"))
prob.t = get(load(file = "/home/knguyen/Documents/PhD/paper/attribution/result/attribution/Post.Prob.List.RData"))

# learning  ---------------------------------------------------------------

# plot t-values of each configuration 
config1 = final.t[which(final.t$pred.y==1),]
dat.p = config1[,c(5:9)]
res1 = reshape2::melt(dat.p)
res1$value = abs(res1$value)
ggplot(res1, aes(x = variable, y = value))+ geom_boxplot()+geom_hline(yintercept = 1.96)

# All 0
d = rep(0,5)

# a =final.t[which(apply(final.t[,c(13:17)], 1, function(x) all(x==d))),c(5:9,19)]

a = final.t[which(final.t$pred.y %in% c(1,31)),c(5:9,19)]
res1 = reshape2::melt(a, id = 'pred.y')
# res1=res1[which(res1$pred.y %in% c(1,4)),]
res1$pred.y = as.factor(res1$pred.y)
# res1$value = abs(res1$value)
ggplot(res1, aes(x = variable, y = value, fill = pred.y))+ geom_boxplot()+geom_hline(yintercept = 1.96)

# Accuracy of the prediction 
b=2
FinalPred <- readRDS(paste0(path_main,"paper/attribution/result/attribution/", "modrf_b",b,".rds"))
dataLearn <- readRDS(file = paste0(path_main,"paper/attribution/result/attribution/", "DataLearn_",b,".rds"))
dataTest <- readRDS(file = paste0(path_main,"paper/attribution/result/attribution/", "DataTest_",b,".rds"))

pred.rf = predict(FinalPred$finalModel, newdata = dataTest)
err.rf <- mean(pred.rf!=dataTest$config)
pred.rf = predict(FinalPred$finalModel, newdata = dataLearn)
err.rf <- mean(pred.rf!=dataLearn$config)

for (i in c(1:36)) {
  a = dataLearn[which(dataLearn$config ==i),]
  res1 = reshape2::melt(a, id = "config")
  res1$config = as.factor(res1$config)
  p <- ggplot(res1, aes(x = variable, y = value, fill = config))+ theme_bw()+ geom_boxplot()+
    geom_hline(yintercept = c(-1.96, 1.96))+
    theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
          axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
          legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
  ggsave(paste0(path_results,"attribution/config", i, ".jpg" ), plot = p, width = 10, height = 8, units = "cm", dpi = 600)
  
}


# PLOT SIMILARITY BETWEEN CONFIGURATIONS ----------------------------------

list.config = readRDS(paste0(path_results,"attribution/List_config.rds"))
name.config = as.numeric(rownames(list.config))
list.config$GE = 2
list.config$GE[which(name.config %in% c(1:14,29:33))] = 3
similar = list()
similar.c = list()

for (i in c(1:36)) {
  list.i = list.config[i,c(1:5)]
  list.wti = c(1:36)[-i]
  similarity = sapply(list.wti, function(x) sum(list.config[x,c(1:5)] == list.i))
  same.g = name.config[which(list.config$GE == list.config$GE[i])]
  l1 = name.config[list.wti[which(similarity>2)]]
  ind1 = which(l1 %in% same.g)
  similar[[i]] = l1[ind1]
  l2 = similarity[which(similarity>2)]
  similar.c[[i]] = l2[ind1]
}

dat.sim = data.frame(config = unlist(sapply(c(1:36), function(x) rep(name.config[x], length(similar[[x]])))),
                     config.sim = unlist(similar),
                     count = unlist(similar.c))
# dat.sim = dat.sim[which(dat.sim$config<16),]
dat.sim$config.sim = as.factor(dat.sim$config.sim)
p = ggplot(dat.sim, aes(fill=config.sim, y=count, x=config, col = config.sim)) + 
  geom_bar(position="stack", stat="identity")+theme_bw()+
  geom_text(aes(label = config.sim),
            colour = "black",  size=1.8,
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y ="Count") + 
  scale_x_continuous(breaks = seq(1,36,1))+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
ggsave(paste0(path_results,"attribution/similarity1.jpg" ), plot = p, width = 14, height = 8, units = "cm", dpi = 600)


# PLOT T VALUES OF CLASSES IN 1 GROUP -------------------------------------
b = 3
FinalPred <- readRDS(paste0(path_results,"attribution/predictive_rule/details/", "modrf_b",b,significance.level, offset, GE, number.pop,".rds"))
dataLearn <- readRDS(file = paste0(path_results,"attribution/predictive_rule/details/", "DataLearn_",b,significance.level, offset, GE, number.pop,".rds"))
dataTest <- readRDS(file = paste0(path_results,"attribution/predictive_rule/details/", "DataTest_",b,significance.level, offset, GE, number.pop,".rds"))
b=2
FinalPred <- readRDS(paste0(path_main,"paper/attribution/result/attribution/", "modrf_b",b,".rds"))
dataLearn <- readRDS(file = paste0(path_main,"paper/attribution/result/attribution/", "DataLearn_",b,".rds"))
dataTest <- readRDS(file = paste0(path_main,"paper/attribution/result/attribution/", "DataTest_",b,".rds"))

grp = 3
list.classes = c(1,8,15,22)

dataLearn$config = rep(name.config, each = 80)
a = dataLearn[which(dataLearn$config %in% list.classes),]
res1 = reshape2::melt(a, id = "config")
res1$config = as.factor(res1$config)
p <- ggplot(res1, aes(x = variable, y = value, fill = config))+ theme_bw()+ geom_boxplot(lwd = 0.1, outlier.size = 0.1)+
  geom_hline(yintercept = c(-1.96, 1.96), lwd = 0.1)+ labs( y = "t_value")+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
ggsave(paste0(path_results,"attribution/group",grp, "train1.jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)
# test set 
grp = 2
list.classes = c(35)
b=2
dataTest$config = rep(name.config, each = 20)
a = dataTest[which(dataTest$config %in% list.classes),]
res1 = reshape2::melt(a, id = "config")
res1$config = as.factor(res1$config)
p <- ggplot(res1, aes(x = variable, y = value, fill = config))+ theme_bw()+ geom_boxplot(lwd = 0.1, outlier.size = 0.1)+
  geom_hline(yintercept = c(-1.96, 1.96), lwd = 0.1)+ labs( y = "t_value")+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
ggsave(paste0(path_results,"attribution/group",grp, "test.jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)


# additional importance for this group
dat.imp = varImp(FinalPred$finalModel)
dat.imp$var = rownames(dat.imp)
class.imp = dat.imp[,c(which(name.config %in% list.classes), "var")]
colnames(class.imp) = c(as.character(list.classes), "var")
dat.p.imp = reshape2::melt(class.imp, id = "var")

p <- ggplot(dat.p.imp, aes(x = variable, y = value, fill = var))+ theme_bw()+ 
  geom_bar(position="dodge", stat="identity")+ labs(x = "class", y = "Importance")+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
ggsave(paste0(path_results,"attribution/group_var_Imp",grp, ".jpg" ), plot = p, width = 10, height = 8, units = "cm", dpi = 600)


# PLOT DISTRIBUTION IN THE REAL DATA --------------------------------------
list.classes = c(1, 8, 15, 22 )
a = final.t[which(final.t$pred.y %in% list.classes),c(5:9,19)]
res1 = reshape2::melt(a, id = 'pred.y')
# res1=res1[which(res1$pred.y %in% c(1,4)),]
res1$pred.y = as.factor(res1$pred.y)
# res1$value = abs(res1$value)
ggplot(res1, aes(x = variable, y = value, fill = pred.y))+ geom_boxplot()+geom_hline(yintercept = 1.96)

p <- ggplot(res1, aes(x = variable, y = value, fill = pred.y))+ theme_bw()+ 
  geom_boxplot(lwd = 0.1, outlier.size = 0.1)+
  geom_hline(yintercept = c(-1.96, 1.96), lwd = 0.1)+ labs( y = "t_value")+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
ggsave(paste0(path_results,"attribution/real_group",do.call(paste, c(as.list(list.classes), sep = "")), ".0.05.jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)



# PLOT DISTRIBUTION OF CONFIGURATION --------------------------------------

count.config = data.frame(table(final.t$pred.y))
count.config.test = data.frame(table(final.t$Z.truth))
count.configs = full_join(count.config, count.config.test, by = "Var1") %>% 
  replace(is.na(.), 0) %>% 
  mutate(truth = Freq.y) %>% 
  mutate(pred = Freq.x - Freq.y) %>% 
  dplyr::select(-c(Freq.x,Freq.y)) %>% 
  mutate_at(vars(Var1), list(factor)) %>% 
  rename("config" = "Var1") %>% 
  reshape2::melt(id="config") %>% 
  arrange(factor(variable, levels= c("pred", "truth"))) %>% 
  group_by(config) %>%
  mutate(lab_ypos = cumsum(value) - 0.5 * value)  %>% 
  mutate_all(~na_if(., 0))
 
  

dat.p = data.frame(count.config)

ggplot(final.t, aes( y = pred.y))+ theme_bw()+
  geom_bar(stat= "count")+
  scale_y_continuous(breaks = seq(1,38,1))+
  scale_x_continuous(breaks = seq(0,140,20))+
  geom_text(aes( label = scales::percent(..prop..),
                 x = ..prop..), stat= "count", hjust = -5.5) 

ggplot(data = count.configs, aes(x = config, y = value)) + theme_bw()+
  geom_col(aes(fill = variable), width = 0.7)+
  geom_text(aes(y = lab_ypos, label = value, group =variable), color = "black") +
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF"))
# PLOT SUSPICIOUS CASES ---------------------------------------------------
source(paste0(path_code_att,"plot_sic_diff_series.R"))

# sort out cases in the table and low probability
list.on.table = final.t[which(is.na(final.t$Z.truth)==TRUE),]
list.high.prob = c(3,5,7,29,31,33)
list.low.prob = list.on.table[which(list.on.table$pred.y %in% list.high.prob == TRUE),]

all.case.p = paste0(list.low.prob$main,".",as.character(list.low.prob$brp), ".", list.low.prob$nearby)
for (i in c(1:length(all.case.p))) {
  plot_six(all.case.p[i])
}

list.low.prob$case = all.case.p
all.case.p = "auck.2016-07-22.wark"


# link to ditance ---------------------------------------------------------
cor.dis = get(load(file = paste0(path_results, "attribution0/stats_test_real_data_corrected_dist.RData")))
variance = read.table(file = paste0(path_results, "attribution0/FGLS_on_real_data_var.txt"), header = TRUE, stringsAsFactors = FALSE)
data.plot = data.frame(n = cor.dis$n1+cor.dis$n2, 
                       code = sapply(c(1:nrow(final.t)), function(x) do.call("paste",c(final.t[x,13:17],sep="_"))),
                       var.rat = variance$mean.G.G./variance$mean.G.E,
                       hdist = cor.dis$distance, 
                       vdist = cor.dis$ver.distance, 
                       pred = final.t$pred.y,
                       truth = sapply(c(1:nrow(final.t)), function(x) ifelse(is.na(final.t$Z.truth[x]), 2, 1)),
                       ID = c(1:nrow(cor.dis)),
                       stringsAsFactors = FALSE)
data.plot$pred = sapply(c(1:nrow(data.plot)), function(x) ifelse(data.plot$pred[x] %in% c(22,8), data.plot$pred[x], 40))
data.plot$pred = as.factor(data.plot$pred )

data.plot$code = sapply(c(1:nrow(data.plot)), function(x) ifelse(data.plot$code[x] %in% unique(data.plot$code)[15], data.plot$code[x], 40))
data.plot$code = as.factor(data.plot$code)

data.plot$code[which(is.na(final.t$Z.truth) == FALSE)] <- 0
a = data.frame(table(data.plot$code), stringsAsFactors = FALSE)
list.con = a$Var1[which(a$Freq>8)]
data.plot$code = sapply(c(1:nrow(data.plot)), function(x) ifelse(data.plot$code[x] %in% list.con, data.plot$code[x], 40))

data.plot$truth = data.plot$truth+15
ggplot(data = data.plot, aes(x = hdist, y = var.rat, col = as.factor(truth))) + theme_bw()+
  geom_point()
  # geom_point(shape=data.plot$truth)
  
# +geom_text(aes(label=ID), nudge_y =0.5)

do.call("paste",c(final.t[x,13:17],sep="_"))

a= final.t[which(is.na(final.t$Z.truth)==TRUE),]
data.plot$j = abs(Total.res1$`jumpG-E`)
data.plot = data.plot[which(is.na(final.t$Z.truth)==TRUE),]

data.plot$g = sapply(c(1:nrow(data.plot)), function(x){
  t = abs(final.t$tGGp[x]) 
  if(t<1.96 & t>0){y =0}
  else if(t<2.58 & t>1.96){y =1}
  # else if(t<3.3 & t>2.58){y =2}
  # else if(t<3.9 & t>3.3){y =3}
  else{y =4}
})

ggplot(data = data.plot, aes(x = hdist, y = var.rat, shape = as.factor(g2), col = as.factor(g))) + theme_bw()+
  geom_point()

data.plot$g1 = sapply(c(1:nrow(data.plot)), function(x){
  t = data.plot$g[x]
  if(t<500 & t>0){y =0}
  else if(t<900 & t>500){y =1}
  # else if(t<3.3 & t>2.58){y =2}
  # else if(t<3.9 & t>3.3){y =3}
  else{y =4}
  return(y)
})
data.plot$g2 = sapply(c(1:nrow(data.plot)), function(x) ifelse(sum(abs(final.t[x, 13:17])) >1, 1,0) )

dat = data.plot
dat$j = abs(Total.res1$`jumpG-E`)

tree <- rpart(truth ~ n+hdist+var.rat+j, data = dat, method = "class")
rpart.plot(tree)

dat$g1 =  sapply(c(1:nrow(data.plot)), function(x) ifelse(dat$n[x] >652, 1, 0))
dat$g2 = sapply(c(1:nrow(data.plot)), function(x) ifelse(dat$j[x]>0.38,1,0) )
dat$g2 = as.factor(dat$g2+10)
ggplot(data = dat, aes(x = hdist, y = var.rat, shape = g2, col = as.factor(truth))) + theme_bw()+
  geom_point()

cvControl=trainControl(method="cv",number=10)
caret::train(truth ~ n+hdist+j, data = dat, method = "rpart", tuneLength = 10,trControl = cvControl)

# test --------------------------------------------------------------------
all = five
g1a = all[which(five$code.GGp==1 & five$code.GEp ==1 & five$code.EEp == 0 & five$code.GpE ==0 & five$code.GpEp==1),]
g1b = all[which(five$code.GGp==-1 & five$code.GEp ==-1 & five$code.EEp == 0 & five$code.GpE ==0 & five$code.GpEp==-1),]
g2 = all[which(five$code.GGp==0 & five$code.GEp ==0 & five$code.EEp == 0 & five$code.GpE ==0 & five$code.GpEp==0),]
g3a = all[which(five$code.GGp==-1 & five$code.GEp ==0 & five$code.EEp == 0 & five$code.GpE ==0 & five$code.GpEp==0),]
g3b = all[which(five$code.GGp==1 & five$code.GEp ==0 & five$code.EEp == 0 & five$code.GpE ==0 & five$code.GpEp==0),]
g4 = all[which(five$code.EEp!=0 & is.na(five$Z.truth)==FALSE),]
g3b = all[which(five$code.GGp==0 & five$code.GEp ==1 & five$code.EEp == 0 & five$code.GpE ==0 & five$code.GpEp==0),]

f1a = g1a[order(g1a $pred.y,decreasing=FALSE),]

f1b = g1b[order(g1b $pred.y,decreasing=FALSE),]
f2 = g2[order(g2 $pred.y,decreasing=FALSE),]
f3a = g3a[order(g3a$pred.y,decreasing=FALSE),]
f3b = g3b[order(g3b$pred.y,decreasing=FALSE),]

rownames(f1a) = NULL
rownames(f1b) = NULL
rownames(f2) = NULL

rownames(f3a) = NULL
rownames(f3b) = NULL

d = g2[,c(4,19)] %>% mutate( pred.y = as.factor(g2$pred.y))  %>%
  reshape2::melt("pred.y")
ggplot(d, aes(x = pred.y, y = abs(value)))+ theme_bw()+ 
  geom_boxplot(lwd = 0.5, outlier.size = 1)+
  labs( y = "jump in the main stations")+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))

d = Total.res1[which(five$code.GGp==0 & five$code.GEp ==0 & five$code.EEp == 0 & five$code.GpE ==0 & five$code.GpEp==0),(2:5)] %>% mutate( pred = as.factor(final.t$pred.y))  %>%
  reshape2::melt("pred")
