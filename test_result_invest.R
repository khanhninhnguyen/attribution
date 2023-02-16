# ivestigate the test results
Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))

all.var = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_var.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
all.arima = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_autocorrelation.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)


# T values as fc of distance 
k = 2
t.val = Total.res[,c(paste0("t", list.name.test[k]), "distance", "ver.distance")]
# dat.p = reshape2::melt(t.val, id = "distance")
colnames(t.val)[1] = "value"
t.val$g = 1
t.val$g[which(abs(t.val$value)<1.96)] = 0
t.val$g = as.factor(t.val$g)
ggplot(data = t.val, aes(x = distance, y = (ver.distance), col = g))+theme_bw()+
  geom_point()+
  # ylim(0,35)+ 
  geom_hline(yintercept = 1.96)+
  geom_hline(yintercept = 2.58)

Total.res[which(t.val$g==0 & t.val$distance<25),]

