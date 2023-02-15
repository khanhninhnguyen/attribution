# to investigate result of predictive rule
significance.level = 0.01
offset = 0
GE = 0
number.pop =3
final.t = get(load(file = paste0(path_results, "attribution/predictive_rule/Final.Table", significance.level, offset, GE, number.pop, ".RData")))
prob.t = get(load(file = paste0(path_results, "attribution/predictive_rule/Post.Prob.List", significance.level, offset, GE, number.pop, ".RData")))
prob.tn = get(load(file = paste0(path_results, "attribution/predictive_rule/Post.Prob.Listn", significance.level, offset, GE, number.pop, ".RData")))
tot = get(load(file = paste0(path_results,"attribution0/stats_test_real_data.RData")))

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

FinalPred <- readRDS(paste0(path_main,"paper/attribution/result/attribution/", "modrf_b",b,".rds"))
dataLearn <- readRDS(file = paste0(path_main,"paper/attribution/result/attribution/", "DataLearn_",b,".rds"))
dataTest <- readRDS(file = paste0(path_main,"paper/attribution/result/attribution/", "DataTest_",b,".rds"))

grp = 1
list.classes = c(1,3,5,7,29,31,33)
b=2
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
ggsave(paste0(path_results,"attribution/group",grp, ".jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)
# test set 
grp = 1
list.classes = c(1,3,5,7,29,31,33)
b=2
dataTest$config = rep(name.config, each = 20)
a = dataLearn[which(dataTest$config %in% list.classes),]
res1 = reshape2::melt(a, id = "config")
res1$config = as.factor(res1$config)
p <- ggplot(res1, aes(x = variable, y = value, fill = config))+ theme_bw()+ geom_boxplot(lwd = 0.1, outlier.size = 0.1)+
  geom_hline(yintercept = c(-1.96, 1.96), lwd = 0.1)+ labs( y = "t_value")+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
ggsave(paste0(path_results,"attribution/group",grp, ".jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)


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
ggsave(paste0(path_results,"attribution/real_groupp",grp, ".jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)



# PLOT DISTRIBUTION OF CONFIGURATION --------------------------------------

count.config = table(final.t$pred.y)
dat.p = data.frame(count.config)

ggplot(final.t, aes( y = pred.y))+ theme_bw()+
  geom_bar(stat= "count")+
  scale_y_continuous(breaks = seq(1,38,1))+
  scale_x_continuous(breaks = seq(0,140,20))+
  geom_text(aes( label = scales::percent(..prop..),
                 x = ..prop..), stat= "count", hjust = -5.5) 


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
all.case.p = "wes2.2006-05-27.npri"
