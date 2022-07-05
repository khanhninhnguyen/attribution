data = read.table(file = paste0(path_results, "attribution/", name.test = "gps.era",".", duration = thres_period, "days.level",
                                sig.level = significant.level, ".scr", screen_value, "no.test.", nb_test.ref, 
                                ".nearby.", nearby.ver = nearby_ver, check_correlation, ".txt"))
additional = read.table( file = paste0(path_results, 
                                       "attribution/attribution.segment.200km.500m.",version_name,
                                       nb_test.ref,nearby_ver,limit.type, screen_value, tolerance_noise, ".additional.txt"),
                         stringsAsFactors = FALSE)
additional <- additional[which(data$nb.before>=100),]
res <- data.frame(matrix(NA, nrow = 0, ncol = 2))
for (d in seq(50,201,5)) {
  a = additional[which(additional$hdist <= d),]
  b = aggregate(nearby.station~ref.station+detected, data = a, FUN = length)
  nb.case = nrow(a)
  nb.br = nrow(b)
  res <- rbind(res, c(nb.case, nb.br))
}
plot(x = seq(50,201,5), y = res[,2], type = "l", xlab = "distance(km)", ylab = "No. breakpoints")
plot(x = seq(50,201,5), y = res[,1], type = "l", xlab = "distance(km)", ylab = "No. cases")


a = read.table(file = paste0(path_data_support,"gps_sta_CODE_REPRO_2015_OPER_combi.txt"), header = TRUE)
colnames(a) = c("name","lat","lon","height","altitude")
b = unlist(list.nearby.ngl)
b = as.character(b)
d = b[which(is.na(b) == FALSE)]
e = d[which(d %in% name_main == FALSE)]
length(unique(d))

list.nearby.ngl = read.table(file = paste0(path_results, "attribution/list.name.nearby.NGL.txt"), stringsAsFactors = FALSE)
b = unlist(list.nearby.ngl[,c(-1)])


list.nearby.ngl = nearby_search_distance(coor.file = paste0(path_data_support,"gps_sta_CODE_REPRO_2015_OPER_combi.txt"),
                       list.name = name_main, list.gnss = name.full, horizontal = 200, vertical=500,
                       version_name = "test", nearby.ver = "test")
  







