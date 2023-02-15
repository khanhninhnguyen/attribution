# ivestigate the test results
Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))

# T values as fc of distance 
t.val = Total.res[,c(paste0("t", list.name.test[4]), "distance")]
dat.p = reshape2::melt(t.val, id = "distance")

ggplot(data = dat.p, aes(x = distance, y = abs(value), col = variable))+theme_bw()+
  geom_line()+
  ylim(0,35)+ 
  geom_hline(yintercept = 1.96)+
  geom_hline(yintercept = 2.58)