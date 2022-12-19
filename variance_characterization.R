# variance characterization from the moving variance 
sd.all= get(load( file = paste0(path_results,"attribution/sd.all_",  win.thres,"years_", nearby_ver,"screened.RData")))

range.var <- function(x){
  if(all(x==1)){
    NA
  }else{
    anu.mean = sapply(c(1:10), function(i) {mean(x[(one.year*(i-1)+1):(one.year*(i))],na.rm = TRUE)})
    anu.max = sapply(c(1:10), function(i) {max(x[(one.year*(i-1)+1):(one.year*(i))],na.rm = TRUE)})
    mean( (anu.max- anu.mean), na.rm = TRUE)
  }
}

range.all <- data.frame(matrix(NA, ncol = 6, nrow = length(sd.all)))
for (testi in c(1:6)) {
  name.var = list.test[testi]
  for (i in c(1:length(sd.all))) {
    sd.sta = sd.all[[i]]
    range.b = range.var(sd.sta$bef[[name.var]])
    range.a = range.var(sd.sta$aft[[name.var]])
    range.m = mean(range.b, range.a)
    range.all[i, testi] <- range.m
  }
}
