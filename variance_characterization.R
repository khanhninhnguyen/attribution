# variance characterization from the moving variance 
sd.all= get(load( file = paste0(path_results,"attribution/sd.all_",  win.thres=10,"years_", nearby_ver,"screened.RData")))
one.year = 365
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
mean.all <- data.frame(matrix(NA, ncol = 6, nrow = length(sd.all)))

for (testi in c(1:6)) {
  name.var = list.test[testi]
  for (i in c(1:length(sd.all))) {
    sd.sta = sd.all[[i]]
    # range.b = range.var(sd.sta$bef[[name.var]])
    # range.a = range.var(sd.sta$aft[[name.var]])
    # if(is.na(range.a)== TRUE & is.na(range.b)== TRUE){
    #   range.m = NA
    # }else{
    #   range.m = mean(range.b, range.a)
    # }
    # range.all[i, testi] <- range.m
    mean.all[i, testi] <- mean( c(mean(sd.sta$bef[[name.var]], na.rm = TRUE), mean(sd.sta$aft[[name.var]], na.rm = TRUE)), na.rm = TRUE)
  }
}
colnames(range.all) <- list.test
colnames(mean.all) <- list.test

a = reshape2::melt(range.all)
b = reshape2::melt(mean.all)
d = rbind(a,b)
d$sta = rep(c("range", "mean"), each = nrow(a))
ggplot(d, aes(value, colour = variable, linetype= sta)) + theme_bw()+
  stat_ecdf()+
  ylab(" CDF ") +
  theme(axis.text = element_text(size = 16),legend.text=element_text(size=12),
        axis.title = element_text(size=16))




