# ivestigate the test results
Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))

all.var = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_var.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
all.arima = read.csv(file = paste0(path_results, "attribution0/FGLS_on_real_data_autocorrelation.txt"), sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)


# T values as fc of distance 
k = 6
t.val = Total.res[,c(paste0("t", list.name.test[k]), "distance", "ver.distance")]
# dat.p = reshape2::melt(t.val, id = "distance")
colnames(t.val)[1] = "value"
t.val$sig = 1
t.val$sig[which(abs(t.val$value)<1.96)] = 0
t.val$sig = as.factor(t.val$sig)
t.val$ver.dist = 100*floor(abs(t.val$ver.distance)/100)
t.val$ver.dist = as.factor(paste0(t.val$ver.dist, "-" , (t.val$ver.dist +100)))

p <- ggplot(data = t.val, aes(x = distance, y = abs(value), col = ver.dist))+theme_bw()+
  geom_point(size = 0.7)+
  ylab("abs t-values")+
  geom_hline(yintercept = 1.96, alpha = 0.5, lwd = 0.5)+
  geom_hline(yintercept = 2.58, alpha = 0.5, lwd = 0.5)+
  geom_hline(yintercept = 3.098, alpha = 0.5, lwd = 0.5)+
  labs(subtitle = list.name.test[k])

ggsave(paste0(path_results,"attribution0/fc.dist",list.name.test[k], ".jpg" ), 
       plot = p, width = 14, height = 10, units = "cm", dpi = 600)



