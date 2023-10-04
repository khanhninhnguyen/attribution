# analyse the test results 
unicode_minus = function(x) sub('^-', '\U2212', format(x))

# check the contradiction function ------
G = c(1,0,-1)
E = c(1,0,-1)
G. = c(1,0,-1)
E. =c(1,0,-1)
a = expand.grid(G, E, G., E.)
colnames(a) = c("E.", "G.", "E", "G")
a = a[,c("G", "E", "G.", "E.")]
Truth = a[-which(a$G==1 & a$E==1 | a$G==0 & a$E==0 | a$G==-1 & a$E==-1),]

logic.table = data.frame(gps.era = Truth$G- Truth$E, 
                         gps.gps = Truth$G- Truth$G.,
                         gps.era1 = Truth$G- Truth$E.,
                         era.era = Truth$E- Truth$E.,
                         gps1.era1 = Truth$G.- Truth$E.,
                         gps1.era = Truth$G.- Truth$E)
truncate.table = logic.table
truncate.table.m = sapply(c(1:6), function(x) {
  a = truncate.table[,x]
  a[which(a==2)]=1
  a[which(a==-2)]=-1
  return(a)
})
truncate.table.mo = unique(truncate.table.m)
truncate.table.5 = as.data.frame(unique(truncate.table.mo[,-1]))
colnames(truncate.table.5) = colnames(logic.table)[2:6]
save(truncate.table.5, file = paste0(path_results, "attribution/truncated.table.RData"))
# check is it the same as in table overleaf
check_contradict <- function(y, table.selected){
  names(y) = NULL
  colnames(table.selected) = NULL
  res = sapply(c(1:36), function(x) identical(unlist(table.selected[x,]), y))
  ind.o = which(res==TRUE)
  out = ifelse(length(ind.o)>0, ind.o, 0)
  return(out)
}




# synthesis result --------------------------------------------------------
# significant level
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres, "years_", nearby_ver,"screened.RData")))
source(paste0(path_code_att,"FGLS.R"))

full.list = get(load( file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData")))
full.list$station = paste0(full.list$main,".",as.character(full.list$brp), ".", full.list$nearby)
full.list$nbc = sapply(c(1:nrow(full.list)), function(x) min(full.list[x,c(4:5)]))
# Reduced list  
full.list$nbc.max = sapply(c(1:nrow(full.list)), function(x) max(full.list[x,c(4:5)]))
ind.sel = which(full.list$nearby!="pama" & full.list$min.var>0.002 & full.list$nbc>200)
# ind.sel = which(full.list$nearby!="pama" & full.list$min.var>0.002 & full.list$nbc>270)
reduced.list = full.list[ind.sel,]
reduced.list = reduced.list[-8,]
rownames(reduced.list) = NULL

# reduced.list$chose[c(65,140,144,204,300,337,378,381,383,384)]=1

# check the significance of GPS-ERA------------
a = aggregate(nbc~main+brp, reduced.list, which.max)
# colnames(a)[3] = "ind"
# d = left_join(reduced.list, a, by = c("main", "brp"))
# d$GE = NA
# for (i in c(1:nrow(a))) {
#   ind = which(d$main == a$main[i] & d$brp == a$brp[i])
# 
#   d$GE[ind[a$ind[i]]] = 1
# }
# 
# reduced.list$GE = d$GE
# list.GE = reduced.list[which(is.na(reduced.list$GE)==FALSE),]
# rownames(list.GE) = NULL

a1 = list.files(path = paste0(path_results, "attribution/FGLS-GE/"))
list.GE = data.frame(station = as.character(substr(a1, 1, 20)))
list.GE = list.GE[which(list.GE$station %in% reduced.list$station == TRUE),]

t.value.GE = sapply(c(1:length(list.GE)), function(x){
  station = get(load(file = paste0(path_results,"attribution/FGLS-GE/", list.GE[x], "fgls.RData")))
  station$gps.era$t.table$`t value`[9]
} )

Total.res = data.frame(matrix(NA, ncol = 15, nrow = nrow(reduced.list)))
for (i in c(1:nrow(reduced.list))) {
  name.i = reduced.list$station[i]
  dat.i = get(load(file = paste0(path_results,"attribution/FGLS-full/", name.i, "fgls.RData")))
  jump.est = sapply(c(2:6), function(x) dat.i[[list.test[x]]]$t.table$Estimate[9])
  t.values = sapply(c(2:6), function(x) dat.i[[list.test[x]]]$t.table$`t value`[9])
  p.values = sapply(c(2:6), function(x) dat.i[[list.test[x]]]$t.table$`Pr(>|t|)`[9])
  Total.res[i,] = c(jump.est, t.values, p.values)
}

colnames(Total.res) = c( paste0("jump", list.name.test[2:6]), paste0("t", list.name.test[2:6]), paste0("p", list.name.test[2:6]))
Total.res$distance = reduced.list$distances
Total.res = cbind(Total.res, reduced.list[,c(4:5)])

# convert to coded table 
convert_coded <- function(x, thres){
  sapply(c(1:length(x)), function(i) ifelse((2*pnorm(-abs(x[i])))<thres, 1*sign(x[i]), 0)) 
}

Total.coded = data.frame(matrix(NA, ncol = 5, nrow = nrow(Total.res)))
for (i in c(1:nrow(Total.res))) {
  case.i = Total.res[i, c(paste0("t", list.name.test[2:6]))]
  Total.coded[i,] = convert_coded(case.i)
}

#check contradiction
trunc.table = get(load(file = paste0(path_results, "attribution/truncated.table.RData")))
contra = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded[x,]), trunc.table)
)
table(unlist(contra))

# plot distribution FOR PAPER ----------------------------
text1 = "Distance < 50 km"
text2 = "Distance > 50 km"
colnames(Total.coded) = list.name.test[2:6]
data1 = Total.coded[which(reduced.list$distances<50),]
# data1 = Total.coded
data.p = reshape2::melt(data1) 
data.p = rbind(data.p, data.frame(variable = rep(list.name.test[1], length(t.value.GE)), value = convert_coded(t.value.GE)))
data.p$c = 1
data.plot = aggregate(c~., data = data.p, sum)
data.plot$S= nrow(data1)
data.plot$S[which(data.plot$variable == list.name.test[1])] = length(t.value.GE)
data.plot$fre = data.plot$c/data.plot$S
data.plot$value = as.factor(data.plot$value)
data.plot$variable = factor(data.plot$variable,  levels = reoder.list.name1)
library(gtable)    
library(grid)
library(gridExtra) 
unicode_minus = function(x) sub('^-', '\U2212', format(x))
p1 <- ggplot(data.plot, aes(fill=value, y=c, x=variable)) + 
  geom_bar(position="stack", stat="identity")+theme_bw()+
  geom_text(aes(label = c),
            colour = "black",  size=1.8,
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y ="Count", tag = "(a)", subtitle = text1) + 
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))

data1 = Total.coded[which(reduced.list$distances>50),]
# data1 = Total.coded
data.p = reshape2::melt(data1) 
data.p = rbind(data.p, data.frame(variable = rep(list.name.test[1], length(t.value.GE)), value = convert_coded(t.value.GE)))
data.p$c = 1
data.plot = aggregate(c~., data = data.p, sum)
data.plot$S= nrow(data1)
data.plot$S[which(data.plot$variable == list.name.test[1])] = length(t.value.GE)
data.plot$fre = data.plot$c/data.plot$S
data.plot$value = as.factor(data.plot$value)
data.plot$variable = factor(data.plot$variable,  levels = reoder.list.name1)

p2 <- ggplot(data.plot, aes(fill=value, y=c, x=variable)) + 
  geom_bar(position="stack", stat="identity")+theme_bw()+
  geom_text(aes(label = c),
            colour = "black",  size=1.8,
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y ="Count", tag = "(b)", subtitle = text2) + 
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)

gB$widths <- gA$widths
# Arrange the two charts.

grid.newpage()
p = (grid.arrange(gA, gB,nrow = 1))

ggsave(paste0(path_results,"attribution/pop_significance_level1.jpg" ), plot = p, width = 14.4, height = 5, units = "cm", dpi = 1200)



# scale_y_discrete()+
  # scale_fill_discrete(breaks = c("1","0","-1"))

contradict = reduced.list[,c(4,5,9:11)]
contradict$contra = sapply(c(1:length(contra)), function(x) ifelse(contra[x] ==0, 1, 0))
sum(contradict$contra[which(contradict$distances<50)])
sum(contradict$contra[which(contradict$distances>50)])


a = read.table(file = paste0(path_results, "attribution/FGLS_on_real_data.txt"),header = TRUE, sep = "\t")

# arima model and variance 
arima.res = data.frame(matrix(NA, ncol = 10, nrow = nrow(reduced.list)))
var.res = data.frame(matrix(NA, ncol = 10, nrow = nrow(reduced.list)))
# for the 5 series
for (i in c(1:nrow(reduced.list))) {
  name.i = reduced.list$station[i]
  dat.i = get(load(file = paste0(path_results,"attribution/FGLS-full/", name.i, "fgls.RData")))
  arma.coefs = sapply(c(2:6), function(x) dat.i[[list.test[x]]]$coef.arma)
  
  var.inf = sapply(c(2:6), function(x){
    y = na.omit(dat.i[[list.test[x]]]$var)
    ye = round(length(y)/365)
    range. = mean(sapply(c(1:ye), function(x) (max(y[(365*x):(365*(x-1))],na.rm=TRUE) - min(y[(365*x):(365*(x-1))],na.rm=TRUE))/2 ))
    return(list(mean(y, na.rm =TRUE), range.))
  }) 
    
  arima.res[i,] = unlist(arma.coefs )
  var.res[i,] = unlist(var.inf)
}
# for G-E
arima.res.1 = data.frame(matrix(NA, ncol = 2, nrow = length(list.GE)))
var.res.1 = data.frame(matrix(NA, ncol = 2, nrow = length(list.GE)))
for (i in c(1:length(list.GE))) {
  name.i = list.GE[i]
  dat.i = get(load(file = paste0(path_results,"attribution/FGLS-GE/", name.i, "fgls.RData")))
  arma.coefs = dat.i$gps.era$coef.arma
  
  var.inf = sapply(c(1), function(x){
    y = na.omit(dat.i$gps.era$var)
    ye = round(length(y)/365)
    range. = mean(sapply(c(1:ye), function(x) (max(y[(365*x):(365*(x-1))],na.rm=TRUE) - min(y[(365*x):(365*(x-1))],na.rm=TRUE))/2 ))
    return(list(mean(y, na.rm =TRUE), range.))
  }) 
  
  arima.res.1[i,] = unlist(arma.coefs )
  var.res.1[i,] = unlist(var.inf)
}
arima.res.1.m =  data.frame(matrix(NA, ncol = 2, nrow = nrow(reduced.list)))
arima.res.1.m[which(reduced.list$station %in% list.GE == TRUE),] =  arima.res.1
all.arima = cbind(arima.res.1.m, arima.res)

var.res.1.m =  data.frame(matrix(NA, ncol = 2, nrow = nrow(reduced.list)))
var.res.1.m[which(reduced.list$station %in% list.GE == TRUE),] =  var.res.1
all.var = cbind(var.res.1.m, var.res)

colnames(all.var) = c(paste0(c("mean.","range"), rep(list.name.test,each=2)))
colnames(all.arima) = c(paste0(c("phi.","theta"), rep(list.name.test,each=2)))

write.table(format(all.var, digits=2), file = paste0(path_results, "attribution/FGLS_on_real_data_var.txt"), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(format(all.arima, digits=2), file = paste0(path_results, "attribution/FGLS_on_real_data_autocorrelation.txt"), sep = '\t', quote = FALSE, row.names = FALSE)


# update the length in the table --------------------
data.vai = data.frame(matrix(NA, ncol =8, nrow = nrow(reduced.list)))
for (i in c(1:nrow(reduced.list))) {
  name.i = reduced.list$station[i]
  dat.i = get(load(file = paste0(path_results,"attribution/FGLS-full/", name.i, "fgls.RData")))
  a = dat[[name.i]]
  b = remove_na_2sides(a, "gps.gps")
  l1 = length(na.omit(a[c(1:3650),"gps.gps"]))
  l2 = length(na.omit(a[c(3651:7300),"gps.gps"]))
  l3 = length(na.omit(a[c(2651:3650),"gps.gps"]))
  l4 = length(na.omit(a[c(3651:4649),"gps.gps"]))
  n = length(dat.i$gps.gps$fit)
  duration = length(na.omit(dat.i$gps.gps$var))
  data.vai[i,] = c(n, duration, nrow(b),l1, l2, l3, l4, name.i)
}
colnames(data.vai) = c("n","duration","duration.raw", "l1", "l2","l3","l4", "name")
data.vai[,c(1:7)] <- as.data.frame(sapply(data.vai[,c(1:7)], as.numeric)) #<- sapply is here

data.vai = cbind(data.vai, reduced.list[,c(4,5)])
data.vai$treat = 2
data.vai$treat[which(data.vai$n == (data.vai$l3+data.vai$l4))]=3
data.vai$treat[which(data.vai$n == (data.vai$l1+data.vai$l2))]=1
data.vai$l1m = NA

n12 = data.frame(matrix(NA, ncol = 2, nrow = nrow(data.vai)))
for (j in c(1:nrow(data.vai))) {
  cri = data.vai$treat[j]
  if(cri==1){
    n12[j,] = c(data.vai$l1[j], data.vai$l2[j])
  }else if(cri==3){
    n12[j,] = c(data.vai$l3[j], data.vai$l4[j])
  }else{
    l.raw = data.vai[j,c("l1","l2")]
    ind.c = which.min(data.vai[j,c("nbc1","nbc2")])
    l.raw[-ind.c] = data.vai$n[j]- l.raw[ind.c] 
    n12[j,] = l.raw
  }
}
reduced.list = cbind(reduced.list, n12)

save(reduced.list, file = paste0(path_results, "attribution/lengthlist.RData"))
# Output significant level -------------
reduced.list$t = NA
reduced.list$t[which(reduced.list$station %in% list.GE ==TRUE)] = t.value.GE
Out.res = cbind(reduced.list[,c(1:3,19)], Total.res[,c(6:10)], reduced.list[,c(17,18,11)])

colnames(Out.res)[c(4,10,11)] = c(paste0("t", list.name.test[1]), "n1", "n2")
write.table(format(Out.res, digits=2), file = paste0(path_results, "attribution/FGLS_on_real_data1.txt"), sep = '\t', quote = FALSE, row.names = FALSE)



# check the length of data old and new 
# plot the noise vs distance 

cor.dis = get(load(file = paste0(path_results, "attribution0/unlist/stats_test_real_data_corrected_dist.RData")))
variance = read.table(file = paste0(path_results, "attribution0/unlist/FGLS_on_real_data_var.txt"),
                      header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
name.test = "G'-E'"
dat.p = data.frame(mean.var = variance[,paste0("mean.", name.test)], 
                   distance = cor.dis$distance,
                   group = as.factor(sapply(c(1:nrow(variance)), function(x) ifelse(cor.dis$distance[x] >50, "far", "close"))))
p = ggplot(data = dat.p, aes(x = distance, y = mean.var))+
  theme_bw()+
  geom_point(size = 0.3)+
  labs(title = name.test)+
  ylim(c(0,12.5))+
  xlab("Distance in associated non-collocated series (km)") + ylab("Mean of MW variance")+
  theme(axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5),
      plot.title = element_text(size = 6),
      axis.title = element_text(size = 6),
      plot.tag = element_text(size = 6),
      legend.box.spacing = unit(3, "pt"),
      legend.key.size = unit(6, 'pt'),
      legend.title.align=0.5,
      plot.subtitle = element_text(size = 6)) 

ggsave(paste0(path_results,"attribution0/mean.MW.var",name.test,".jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 600)

# distribution of noise 

dat.p %>%
  ggplot(aes(x=mean.var, fill=group)) + theme_bw()+
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
  scale_fill_manual(values=c("#69b3a2", "coral")) +
  xlab(paste0("Mean of MW variance ", name.test)) +
  theme(axis.text = element_text(size = 12))
# for non-collocated 
dat.p %>%
  ggplot(aes(x=mean.var)) + theme_bw()+
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
  xlab(paste0("Mean of MW variance ", name.test)) +
  theme(axis.text = element_text(size = 12))

# variance of jump with distance 

name.test = "G-E'"
dat.p = data.frame(mean.var = Total.res[,paste0(name.test)], 
                   distance = cor.dis$distance,
                   group = as.factor(sapply(c(1:nrow(variance)), function(x) ifelse(cor.dis$distance[x] >50, "far", "close"))))

dat.p %>%
  ggplot(aes(x=mean.var, fill=group)) + theme_bw()+
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
  scale_fill_manual(values=c("#69b3a2", "coral")) +
  xlab(paste0("Mean of est variance of the jump ", name.test)) +
  theme(axis.text = element_text(size = 12))
# boxplot 

dat.plot = Total.res[,c(-1,-5)] %>% 
  mutate(group = as.factor(sapply(c(1:nrow(variance)), function(x) ifelse(cor.dis$distance[x] >50, "far", "close"))))  %>% 
  reshape2::melt(id = "group")  %>% 
  rbind(dat = data.frame(group = as.factor("close"), 
                         variable = rep(list.name.test[c(1,5)], each = 494),
                         value = c(Total.res[,1],Total.res[,5])))  

  










# DISTRIBUTION OF JUMPS, STD ERR AND T VALUES FOR 2 GROUPS ----------------

tot = get(load(file = paste0(path_results,"attribution0/stats_test_real_data.RData")))
five = get(load(file = "/home/knguyen/Documents/PhD/paper/attribution/result/attribution/Final.Table.RData"))
for (i in c(1:6)) {
  name.t = list.name.test[i]
  tot[paste0("std.err", name.t)] <- tot[paste0("jump", name.t)]/ tot[paste0("t", name.t)]
}
# remove replicated value in G-E
name.t = "G-E"
tot[which(is.na(five$tGE)==TRUE),c(paste0("std.err", name.t), paste0("jump", name.t), paste0("t", name.t))] <- NA


collo = list.name.test[c(1,5)]
noncollo = list.name.test[-c(1,5)]

param = "t"
param.in.plot = "Absolute t-value"

data.plot = tot[,paste0(param, noncollo)] %>%
  mutate(group = as.factor(sapply(c(1:494), function(x) ifelse(tot$distance[x]<50, "<50", ">50")))) %>%
  reshape2::melt(id = "group")  %>%
  rbind(dat = data.frame(group = as.factor(rep("<50", 494*2)),
                         variable = rep(paste0(param, collo), each = 494),
                         value = unlist(tot[,paste0(param, collo)]))) %>%
  `rownames<-`( NULL ) %>% 
  mutate(variable = gsub(param, "", variable)) %>% 
  mutate(variable = factor(variable, levels = reoder.list.name1))


p <- ggplot(data.plot, aes(x = variable, col = group, y = abs(value)))+ theme_bw()+ 
  geom_boxplot(lwd = 0.2, outlier.size = 0.1)+
  labs(y = paste0(param.in.plot))+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), legend.title=element_text(size=5),
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.box.spacing = unit(0, "pt"), 
        plot.margin = rep(unit(0,"null"),4)
        ) +
  guides(color = guide_legend(title = "Distance"))

ggsave(paste0(path_results,"attribution/abs_", param, "_distance.jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 1200)

# FOR PAPER ---------------------------------------------------------------

data.plot1 = tot[,paste0(param, list.name.test)] %>%
  reshape2::melt()  %>%
  `rownames<-`( NULL ) %>% 
  mutate(variable = gsub(param, "", data.plot$variable)) %>% 
  mutate(variable = factor(data.plot1$variable, levels = reoder.list.name1))

p1 <- ggplot(data.plot1, aes(x = variable, y = abs(value)))+ theme_bw()+ 
  geom_boxplot(lwd = 0.2, outlier.size = 0.1)+
  labs(y = paste0(param.in.plot))+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), legend.title=element_text(size=5),
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.box.spacing = unit(0, "pt"), 
        # plot.margin = rep(unit(0,"null"),4)
  ) 

param = "jump"
param.in.plot = "Jump"

data.plot2 = tot[,paste0(param, noncollo)] %>%
  mutate(group = as.factor(sapply(c(1:494), function(x) ifelse(tot$distance[x]<50, "<50", ">50")))) %>%
  reshape2::melt(id = "group")  %>%
  rbind(dat = data.frame(group = as.factor(rep("<50", 494*2)),
                         variable = rep(paste0(param, collo), each = 494),
                         value = unlist(tot[,paste0(param, collo)]))) %>%
  `rownames<-`( NULL ) %>% 
  mutate(variable = gsub(param, "", data.plot1$variable))  %>%
  mutate(variable = factor(data.plot2$variable, levels = reoder.list.name1))

p2 <- ggplot2(data.plot2, aes(x = variable, col = group, y = abs(value)))+ theme_bw()+ 
  geom_boxplot(lwd = 0.2, outlier.size = 0.1)+
  labs(y = paste0(param.in.plot))+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), legend.title=element_text(size=5),
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.box.spacing = unit(0, "pt"), 
        # plot.margin = rep(unit(0,"null"),4)
  ) +
  guides(color = guide_legend(title = "Distance"))


# compute statistics of the 6 series  -------------------------------------

coef = read.table( file = paste0(path_results, "attribution0/FGLS_on_real_data_autocorrelation.txt"), header = TRUE, check.names = FALSE)
five = get(load(file = "/home/knguyen/Documents/PhD/paper/attribution/result/attribution/Final.Table.RData"))

model.recons = as.data.frame(lapply(c(2:6), function(x){
  sapply(c(1:494), function(y) {
    coef.xy = coef[y,c((2*x-1): (2*x))]
    if(coef.xy[1] ==0 & coef.xy[2] ==0){ r = "White"}
    else if (coef.xy[1] !=0 & coef.xy[2] ==0){ r = "ar"}
    else if (coef.xy[1] ==0 & coef.xy[2] !=0){ r = "ma"}
    else if (coef.xy[1] !=0 & coef.xy[2] !=0){ r = "arma"}
    return(r)
  })
}))
colnames(model.recons) = list.name.test[-1]

ind = which(model.recons$`G'-E'` == "ar" & five$distance < 50)
summary(coef$`phi.G'-E'`[ind])




# For Thesis: plot the dependence of result on the distance ---------------
tot = get(load(file = paste0(path_results,"attribution0/unlist/stats_test_real_data.RData")))
five = get(load(file = "/home/knguyen/Documents/PhD/paper/attribution/result/attribution/Final.Table.RData"))
for (i in c(1:6)) {
  name.t = list.name.test[i]
  tot[paste0("std.err", name.t)] <- tot[paste0("jump", name.t)]/ tot[paste0("t", name.t)]
}
cor.dis = get(load(file = paste0(path_results, "attribution0/unlist/stats_test_real_data_corrected_dist.RData")))
variance = read.table(file = paste0(path_results, "attribution0/unlist/FGLS_on_real_data_var.txt"),
                      header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
name.test = "E-E'"
ind.GE = which(!is.na(five$tGE))
dat.p = data.frame(std.err = tot[,paste0("std.err", name.test)], 
                   t = abs(tot[,paste0("t", name.test)]),
                   mean.var = sqrt(variance[,paste0("mean.", name.test)]),
                   range.var = variance[,paste0("range", name.test)],
                   dist = tot$distance,
                   group = sapply(c(1:nrow(variance)), function(x) ifelse(cor.dis$distance[x] <250, "<50", ">50")))
# dat.p = dat.p[ind.GE,]
# r = round(cor(dat.p$std.err, dat.p$mean.var), digits = 2)

p = ggplot(data = dat.p, aes(x = dist, y = std.err))+
# p = ggplot(data = dat.p, aes(x = mean.var, y = std.err))+
  theme_bw()+
  # geom_point(size = 0.3, aes(color = group))+
  geom_point(size = 0.3)+
  labs(title = name.test)+
  ylim(c(0,0.5))+
  # xlim(c(0,3.5))+
  # annotate("text", x=3, y=0.5, label= paste0("r = ",r ), size = 2) + 
  # ylim(c(0,40))+
  xlab("Distance (km)") + ylab("Standard error of the jump")+
  # xlab("Distance between associated non-collocated series (km)") + ylab("Standard error of the jump")+
  # xlab("Mean of standard deviation") + ylab("Standard error of the jump")+
  # scale_colour_manual(values = c('<50' = 'red','>50' = 'blue'),labels = expression(d<50,d>=50))+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.title = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 6),
        legend.box.spacing = unit(3, "pt"),
        legend.key.size = unit(6, 'pt'),
        legend.title.align=0.5,
        plot.subtitle = element_text(size = 6))+
  theme(legend.position = "none")
  # guides(color = guide_legend(title = "Distance(km)"))

ggsave(paste0(path_results,"attribution/std.err", name.test, "_dist.jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 1200)
# ggsave(paste0(path_results,"attribution/std.err", name.test, "_mean.sd.jpg" ), plot = p, width = 8, height = 5, units = "cm", dpi = 1200)


