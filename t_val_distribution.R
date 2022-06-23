# study the distribution of t-value

name.test = list.test[2]
data = read.table(file = paste0(path_results, "attribution/", name.test,".", duration = thres_period, "days.level",
                         sig.level = significant.level, ".scr", screen_value, "no.test.", nb_test.ref, 
                         ".nearby.", nearby.ver = nearby_ver, check_correlation, ".txt"))
additional = read.table( file = paste0(path_results, 
                                      "attribution/attribution.segment.200km.500m.",version_name,
                                      nb_test.ref,nearby_ver,limit.type, screen_value, tolerance_noise, ".additional.txt"))


x = data$t
hist(x,breaks = 100)
lines(density(x))
a = data[which(abs(x)>10),]
additional[which(abs(x)>10),]