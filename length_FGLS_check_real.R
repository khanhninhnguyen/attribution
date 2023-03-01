# check the length of FGLS test data 
# read list of stations
win.thres = 10
# dat = get(load( file = paste0(path_results,"attribution0/data.all_", win.thres, "years_", nearby_ver,"screened.RData")))
source(paste0(path_code_att,"FGLS.R"))

full.list = get(load( file = paste0(path_results, "attribution0/list.segments.selected", win.thres,".RData")))
full.list$station = paste0(full.list$main,".",as.character(full.list$brp), ".", full.list$nearby)
full.list$nbc = sapply(c(1:nrow(full.list)), function(x) min(full.list[x,c(4:5)]))
# Reduced list  
full.list$nbc.max = sapply(c(1:nrow(full.list)), function(x) max(full.list[x,c(4:5)]))
ind.sel = which(full.list$nearby!="pama" & full.list$min.var>0.002 & full.list$nbc>200)
# ind.sel = which(full.list$nearby!="pama" & full.list$min.var>0.002 & full.list$nbc>270)
reduced.list = full.list[ind.sel,]
reduced.list = reduced.list[-8,]
rownames(reduced.list) = NULL

# read total length of 2 data sets 

data.vai = data.frame(matrix(NA, ncol = 4, nrow = nrow(reduced.list)))
for (i in c(1:nrow(reduced.list))) {
  name.i = reduced.list$station[i]
  dat.old = get(load(file = paste0(path_results,"attribution0/FGLS-full/", name.i, "fgls.RData")))
  dat.new = get(load(file = paste0(path_results,"attribution0/FGLS/", name.i, "fgls.RData")))
  l1 = length(dat.old$gps.gps$all.out$residual)
  l2 = length(dat.new$gps.gps$all.out$residual)
  l3 = length(dat.old$gps1.era1$all.out$residual)
  l4 = length(dat.new$gps1.era1$all.out$residual)
  data.vai[i,] = c(l1,l2,l3,l4)
}
colnames(data.vai)[1:4] = c("old.gps.gps","new.gps.gps", "old.gps1.era1","new.gps1.era1")
data.vai$diff.gps.gps = data.vai$new.gps.gps - data.vai$old.gps.gps
data.vai$diff.gps1.era1 = data.vai$new.gps1.era1 - data.vai$old.gps1.era1

save(data.vai, file = paste0(path_results,"attribution0/length_fgls.RData"))

Total.res1 = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))
Total.res2 = get(load(paste0(path_results,"attribution0/stats_test_real_data_new.RData")))

compar = cbind(Total.res1$config, Total.res2$config, data.vai, Total.res1$station)

