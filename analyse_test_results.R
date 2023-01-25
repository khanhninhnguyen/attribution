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
win.thres = 10
full.list = get(load( file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData")))
full.list$station = paste0(full.list$main,".",as.character(full.list$brp), ".", full.list$nearby)
full.list$nbc = sapply(c(1:nrow(full.list)), function(x) min(full.list[x,c(4:5)]))
full.list$chose[c(119, 260, 278, 378, 587, 661, 728, 742, 760, 770)]=1
# check available results 
all_file = list.files(path =paste0(path_results,"attribution/FGLS-full/"))
all_station = substr(all_file, 1, 20)
ind.avai = which(full.list$station %in% all_station == TRUE)
full.list = full.list[ind.avai,]
# Reduced list  
full.list$nbc.max = sapply(c(1:nrow(full.list)), function(x) max(full.list[x,c(4:5)]))
# ind.sel = which(full.list$nearby!="pama" & full.list$min.var>0.002 & full.list$nbc>200 & full.list$nbc.max >365)
ind.sel = which(full.list$nearby!="pama" & full.list$min.var>0.002 & full.list$nbc>270)
reduced.list = full.list[ind.sel,]
rownames(reduced.list) = NULL
# reduced.list$chose[c(65,140,144,204,300,337,378,381,383,384)]=1

# check the significance of GPS-ERA
# a = aggregate(nbc~main+brp, reduced.list, which.max)
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

a = list.files(path = paste0(path_results, "attribution/FGLS-GE/"))
list.GE = data.frame(station = substr(a, 1, 20))
t.value.GE = sapply(c(1:nrow(list.GE)), function(x){
  station = get(load(file = paste0(path_results,"attribution/FGLS-GE/", list.GE$station[x], "fgls.RData")))
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

colnames(Total.res) = c( paste("jump", list.name.test[2:6]), paste("t", list.name.test[2:6]), paste("p", list.name.test[2:6]))
Total.res$distance = reduced.list$distances
Total.res = cbind(Total.res, reduced.list[,c(4:5)])

# convert to coded table 
convert_coded <- function(x){
  sapply(c(1:length(x)), function(i) ifelse(abs(x[i])>1.96, 1*sign(x[i]), 0)) 
}

Total.coded = data.frame(matrix(NA, ncol = 5, nrow = nrow(Total.res)))
for (i in c(1:nrow(Total.res))) {
  case.i = Total.res[i, c(paste("t", list.name.test[2:6]))]
  Total.coded[i,] = convert_coded(case.i)
}

#check contradiction
trunc.table = get(load(file = paste0(path_results, "attribution/truncated.table.RData")))
contra = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded[x,]), trunc.table)
)
table(unlist(contra))

# plot distribution
colnames(Total.coded) = list.name.test[2:6]
data1 = Total.coded[which(reduced.list$distances<50),]
data.p = reshape2::melt(data1) 
data.p = rbind(data.p, data.frame(variable = rep(list.name.test[1], length(t.value.GE)), value = convert_coded(t.value.GE)))
data.p$c = 1
data.plot = aggregate(c~., data = data.p, sum)
data.plot$S= nrow(data1)
data.plot$S[which(data.plot$variable == list.name.test[1])] = length(t.value.GE)
data.plot$fre = data.plot$c/data.plot$S
data.plot$value = as.factor(data.plot$value)
data.plot$variable = factor(data.plot$variable,  levels = reoder.list.name)

ggplot(data.plot, aes(x=variable, y = fre, fill = value))+theme_bw()+
  geom_col(position = position_fill()) +
  geom_text(aes(label = scales::percent(fre)),
            colour = "black",
            position = position_fill(vjust = 0.5)) +
  labs(x = NULL, y ="Series", subtitle =  "Distance <50 km")  


# scale_y_discrete()+
  # scale_fill_discrete(breaks = c("1","0","-1"))

contradict = reduced.list[,c(4,5,9:11)]
contradict$contra = sapply(c(1:length(contra)), function(x) ifelse(contra[x] ==0, 1, 0))
sum(contradict$contra[which(contradict$distances<50)])
sum(contradict$contra[which(contradict$distances>50)])

# Output -------------
reduced.list$t = NA
reduced.list$t[which(is.na(reduced.list$GE)==FALSE)] = t.value.GE
Out.res = cbind(reduced.list[,c(1:3,18)], Total.res[,c(6:10)], reduced.list[,c(4,5,11)])
Out.res$GE = reduced.list$chose
Out.res$GE[which(is.na(Out.res$GE)==FALSE)]=1

colnames(Out.res)[4] = paste("t", list.name.test[1])
format(Out.res, digits=2)
write.table(format(Out.res, digits=2), file = paste0(path_results, "attribution/FGLS_on_real_data.txt"), sep = '\t', quote = FALSE, row.names = FALSE)

