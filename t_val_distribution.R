# study the distribution of t-value

additional = read.table( file = paste0(path_results, 
                                      "attribution/attribution.segment.200km.500m.",version_name,
                                      nb_test.ref,nearby_ver,limit.type, screen_value, tolerance_noise, ".additional.txt"))


x = data$t
hist(x,breaks = 100)
lines(density(x))
a = data[which(abs(x)>10),]
additional[which(abs(x)>10),]

# read all results  -------------------------------------------------------
res.all = data.frame()
for (i in 1:6) {
  name.test = list.test[i]
  data = read.table(file = paste0(path_results, "attribution/", name.test,".", duration = thres_period, "days.level",
                                  sig.level = significant.level, ".scr", screen_value, "no.test.", nb_test.ref, 
                                  ".nearby.", nearby.ver = nearby_ver, check_correlation, ".txt"))
  res.all <- rbind(res.all, data)
}

t.val = data.frame(test = rep(factor(list.test, levels = list.test), each = (nrow(res.all)/6)), t = res.all$t)
ggplot(t.val, aes(t, col = test, facets = test)) +
  # geom_density() + 
  geom_histogram(bins = 100)+
  facet_grid( ~ test, scales = "free_x")+ theme_bw()

ind = which(t.val$test == "gps.gps" & abs(t.val$t) >10)
res.all[ind,]
additional[(ind ),]

