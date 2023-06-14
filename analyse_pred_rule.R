
final.t = get(load(file = paste0(path_results, "attribution/predictive_rule/details/Final.Table.0.01.RData")))
prob.t = get(load(file = paste0(path_results, "attribution/predictive_rule/details/Post.Prob.List.0.01.RData")))
prob.tn = get(load(file = paste0(path_results, "attribution/Post.Prob.Listn.RData")))
tot = get(load(file = paste0(path_results,"attribution0/stats_test_real_data.RData")))
tot = cbind(final.t[,c(1:3)], tot)
rownames(tot) = NULL
rownames(final.t) = NULL

valid = get(load(file = paste0(path_results,"validation/",nb_test.ref,"-",criterion,"metacompa",screen.value="",".RData")))
last.name = sapply(c(1:114), function(x) prob.t[[x]]$MainBreak[1])
last.brp = sapply(c(1:114), function(x) prob.t[[x]]$MainBreak[2])
last.pre = sapply(c(1:114), function(x) prob.t[[x]]$Config.Pred.Post)
valid.list = sapply(c(1:114), function(x) valid$valid[which(valid$name == last.name[x] & valid$detected == last.brp[x])])

# test with length
final.t$n = final.t$n1 + final.t$n2
last.p = sapply(vc(1:114), function(x) prob.t[[x]]$PostProb[ which(names( prob.t[[x]]$PostProb) == as.character(prob.t[[x]]$Config.Pred.Post))])


res = data.frame(main = last.name, brp = last.brp, pred = last.pre, valid = valid.list,last.p )
G = c(1:7,15:21,29:38)
res$G = 0
res$G[which(res$pred %in% G ==TRUE)] = 1
length(which(res$valid==1 & res$G==1))
sum(res$valid==1)
sum(res$G==1)

# plot statistics 
table(final.t$pred.y)
a = table(final.t$pred.y)
# GPS only 
sum(a[1:5])
sum(a[9:14])
# ERA only 
sum(a[6:8])
sum(a[15:17])

d = data.frame(a)
d$Var1 = as.numeric(names(a))
d$Freq = as.numeric(d$Freq)

last.name = sapply(c(1:114), function(x) prob.t[[x]]$MainBreak[1])
last.brp = sapply(c(1:114), function(x) prob.t[[x]]$MainBreak[2])

last.pre = sapply(c(1:114), function(x) prob.t[[x]]$Config.Pred.Post)

final.t$last = NA
for (i in c(1:length(last.name))) {
  ind = which(final.t$main == last.name[i] & final.t$brp == last.brp[i])
  final.t$last[ind[1]] = last.pre[i]
}

dat.p = data.frame( name = last.name, brp = last.brp, last.pre)
# dat.p = final.t
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
ggsave(paste0(path_results,"attribution/config_dist.test_1.jpg" ), plot = p, width = 8.8, height = 6, units = "cm", dpi = 600)

# data -------------------------------
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres = 10,"years_", nearby_ver,"screened.RData")))
shorten = final.t[,c(1:12, 18:20)]

name.i = "brus.2005-06-03.dies"
brp = as.Date(substr(name.i, 6, 15))
par(mfrow=c(1,1))    # set the plotting area into a 1*2 array
ind.case = which(all.case %in% name.i  == TRUE)

# read data
name.s = "gps.era"
datai = remove_na_2sides(datai, name.s)
brp.ind = which(datai$date == brp)
all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)
# read result FGLS 5 series 
station1 = get(load(file = paste0(path_results,"attribution/FGLS-GE/",name.i, "fgls.RData")))

datai$fit = station1[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station1[[name.s]]$t.table$Estimate[10] + station1[[name.s]]$t.table$Estimate[9]
text1 = paste0("jump = ", round(station1[[name.s]]$t.table$Estimate[9], digits = 3), 
               ", t = ", round(station1[[name.s]]$t.table$`t value`[9], digits = 3), 
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
plot(datai$date, datai[,name.s], col="gray", type = "l", ylab = name.s, xlab = "", main = text1)
abline(v = brp)
lines(datai$date,datai$fit)

name.s = "gps.gps"
# read data
datai = dat[[name.i]]
datai = remove_na_2sides(datai, name.s)
brp.ind = which(datai$date == brp)
all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)

# read result FGLS 5 series 

station = get(load(file = paste0(path_results,"attribution/FGLS-full/",name.i, "fgls.RData")))
datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 3), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 3), 
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])

plot(datai$date, datai[,name.s], col="gray", type = "l", ylab = name.s, xlab = "", main = text1)
abline(v = brp)
lines(datai$date,datai$fit)

# 
name.s = "gps.era1"
# read data
datai = dat[[name.i]]
datai = remove_na_2sides(datai, name.s)
brp.ind = which(datai$date == brp)
all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)

# read result FGLS 5 series 

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 3), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 3), 
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
plot(datai$date, datai[,name.s], col="gray", type = "l", ylab = name.s, xlab = "", main = text1)
abline(v = brp)
lines(datai$date,datai$fit)

# 
name.s = "era.era"
# read data
datai = dat[[name.i]]
datai = remove_na_2sides(datai, name.s)
brp.ind = which(datai$date == brp)
all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)

# read result FGLS 5 series 

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 3), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 3), 
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
plot(datai$date, datai[,name.s], col="gray", type = "l", ylab = name.s, xlab = "", main = text1)
abline(v = brp)
lines(datai$date,datai$fit)

# 

name.s = "gps1.era1"
# read data
datai = dat[[name.i]]
datai = remove_na_2sides(datai, name.s)
brp.ind = which(datai$date == brp)
all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)

# read result FGLS 5 series 

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 3), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 3), 
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
plot(datai$date, datai[,name.s], col="gray", type = "l", ylab = name.s, xlab = "", main = text1)
abline(v = brp)
lines(datai$date,datai$fit)

# 
name.s = "gps1.era"
# read data
datai = dat[[name.i]]
datai = remove_na_2sides(datai, name.s)
brp.ind = which(datai$date == brp)
all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)
# read result FGLS 5 series 

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 3), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 3), 
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
plot(datai$date, datai[,name.s], col="gray", type = "l", ylab = name.s, xlab = "", main = text1)
abline(v = brp)
lines(datai$date,datai$fit)

# plot(final.t$n1+final.t$n2, final.t$tEEp)
write.table(format(final.t, digits =4), file = paste0(path_results,"attribution/FinalTable.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(format(res, digits =4), file = paste0(path_results,"attribution/last_decision.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

a = left_join(final.t, res, by = c("main","brp"))


# investigate the distance 

list.s = res[which(res$pred %in% c(8,22)),]
list.s.all = inner_join(list.s, final.t, by = c("main","brp"))
list.s.all$r = d
a = aggregate(r~main+brp+pred.y, list.s.all, min)
a$case = paste0(a$main,".",a$brp)
a$pred.y = as.factor(a$pred.y)

ggplot(a, aes(x = pred.y, y = distance))+ theme_bw()+ geom_point()

a = list.s1$X1
for (i in c(2:(length(a)))){
  if(is.na(a[i])==TRUE){
    a[i] = a[i-1]
    }
}
d = list.s1$X2/a

a = aggregate(distance~main+brp+pred.y, final.t, min)

in.t = cbind( final.t[,c(1:3, 12,19)], tot[,c(4,27:28)])
a = na.omit(in.t[,c(1:2,6:7)])
colnames(a)[3:4] = c("j.g","var.g")
inves.t1 = left_join(in.t, a, by = c("main","brp"))
inves.t1$r = inves.t1$X2/inves.t1$var.g
inves.t1$pred.y = as.factor(inves.t1$pred.y)

inves.t1$rati = inves.t1$j.g/inves.t1$r
ggplot(inves.t1, aes(x = pred.y, y =rati))+geom_boxplot()

