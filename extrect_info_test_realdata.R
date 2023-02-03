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

# Add length of series 

lengthlist = get(load(file = paste0(path_results, "attribution/lengthlist.RData")))
Total.res$n1 = lengthlist$X1
Total.res$n2 = lengthlist$X2

save(Total.res, file = paste0(path_results,"stats_test_real_data.RData"))

# To send 

Total.res <- Total.res %>% dplyr::select(-c(config, station))
Total.res = cbind(reduced.list[,c(1:3)], Total.res)
write.table(Total.res, file = paste0(path_results,"stats_test_real_data.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# 
# b = which(a$l2<1000 & a$de>100) 
# save(b, file = paste0(path_results, "attribution/add.list.RData"))

# Additional information: variance and ARMA(1,1) model --------------------

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




