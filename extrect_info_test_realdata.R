# extract infor from the test on the real data 

convert_coded <- function(x){
  sapply(c(1:length(x)), function(i) ifelse(abs(x[i])>1.96, 1*sign(x[i]), 0)) 
}
check_contradict <- function(y, table.selected){
  names(y) = NULL
  colnames(table.selected) = NULL
  res = sapply(c(1:36), function(x) identical(unlist(table.selected[x,]), y))
  ind.o = which(res==TRUE)
  out = ifelse(length(ind.o)>0, ind.o, 0)
  return(out)
}



# list of station ---------------------------------------------------------

win.thres = 10
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



# statistics of G-E -------------------------------------------------------

a1 = list.files(path = paste0(path_results, "attribution/FGLS-GE/"))
list.GE = data.frame(station = as.character(substr(a1, 1, 20)))
ind.GE = which(list.GE$station %in% reduced.list$station == TRUE)
list.GE = list.GE[ind.GE,]

t.value.GE = sapply(c(1:length(list.GE)), function(x){
  station = get(load(file = paste0(path_results,"attribution/FGLS-GE/", list.GE[x], "fgls.RData")))
  station$gps.era$t.table$`t value`[9]
} )

jump.GE = sapply(c(1:length(list.GE)), function(x){
  station = get(load(file = paste0(path_results,"attribution/FGLS-GE/", list.GE[x], "fgls.RData")))
  station$gps.era$t.table$Estimate[9]
} )

# statistics of 5 others --------------------------------------------------

Total.res = data.frame(matrix(NA, ncol = 15, nrow = nrow(reduced.list)))
for (i in c(1:nrow(reduced.list))) {
  name.i = reduced.list$station[i]
  dat.i = get(load(file = paste0(path_results,"attribution/FGLS-full/", name.i, "fgls.RData")))
  jump.est = sapply(c(2:6), function(x) dat.i[[list.test[x]]]$t.table$Estimate[9])
  t.values = sapply(c(2:6), function(x) dat.i[[list.test[x]]]$t.table$`t value`[9])
  p.values = sapply(c(2:6), function(x) dat.i[[list.test[x]]]$t.table$`Pr(>|t|)`[9])
  Total.res[i,] = c(jump.est, t.values, p.values)
}

colnames(Total.res) = c(paste0("jump", list.name.test[2:6]), paste0("t", list.name.test[2:6]), paste0("p", list.name.test[2:6]))

# add the G_E
Total.res[,c(paste0("jump", list.name.test[1]), paste0("t", list.name.test[1]))] = NA
ind.GE1 = which(reduced.list$station %in% list.GE == TRUE)

Total.res$`jumpG-E`[ind.GE1] = jump.GE
Total.res$`tG-E`[ind.GE1] = t.value.GE

Total.res = Total.res[,c(paste0("jump", list.name.test[1:6]), paste0("t", list.name.test[1:6]))]
# add the significance code -----------------------------------------------

Total.coded = data.frame(matrix(NA, ncol = 6, nrow = nrow(Total.res)))
for (i in c(1:nrow(Total.res))) {
  case.i = Total.res[i, c(paste0("t", list.name.test[1:6]))]
  Total.coded[i,] = convert_coded(case.i)
}
colnames(Total.coded) = list.name.test
Total.res = cbind(Total.res, Total.coded)

# check contracdiction

trunc.table = get(load(file = paste0(path_results, "attribution/truncated.table.RData")))
contra = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded[x,c(2:6)]), trunc.table)
)

Total.res$config = contra
# combine all results
Total.res$distance = reduced.list$distances
# length 

last.result = read.table(file = paste0(path_results,"attribution/FGLS_on_real_data_t.txt"),
                         header = TRUE, 
                         stringsAsFactors = FALSE)

Total.res = cbind(Total.res, last.result[,c(10,11)])
Total.res$station = reduced.list$station

lengthlist = get(load(file = paste0(path_results, "attribution/lengthlist.RData")))
Total.res$n1 = lengthlist$X1
Total.res$n2 = lengthlist$X2

save(Total.res, file = paste0(path_results,"stats_test_real_data.RData"))

# 
# b = which(a$l2<1000 & a$de>100) 
# save(b, file = paste0(path_results, "attribution/add.list.RData"))

# variance 
variance = data.frame(matrix(NA, ncol = 5, nrow = nrow(reduced.list)))
for (i in c(1:nrow(reduced.list))) {
  name.i = reduced.list$station[i]
  dat.i = get(load(file = paste0(path_results,"attribution/FGLS-full/", name.i, "fgls.RData")))
  p.values = sapply(c(2:6), function(x) mean(sqrt(dat.i[[list.test[x]]]$var),na.rm =TRUE))
  variance[i,] = p.values
}





