# analyse the test results 

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
convert_coded <- function(x){
  sapply(c(1:length(x)), function(i) ifelse(abs(x[i])>1.96, 1*sign(x[i]), 0)) 
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

