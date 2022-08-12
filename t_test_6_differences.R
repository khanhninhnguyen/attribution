# main program for t test 
source(paste0(path_code_att,"support_attribution_t_test.R"))
source(paste0(path_code_att,"ESS_ARMA_param_est.R"))
name = name_main

# search for nearby -------------------------------------------------------

list.nearby.ngl = nearby_list_ngl(path_results,  
                                  nearby.ver = nearby_ver, 
                                  version_name, 
                                  name_main = name_main)

save(list.nearby.ngl,  file = paste0(path_results, "attribution/", "list_nearby", nearby_ver, ".RData"))

# form 6 series of differences --------------------------------------------

if (screen_value == 1){
  res.seg <-  get(load(file = paste0(path_results, "meta_compare1/",nb_test=34, criterion="BM_BJ", "position_break_screened.RData")))
  
  list.nearby.station.homo <- homogeneous.nearby.segment.multiple(name.station = name,
                                                                  valid.file.ref = validation.file.ref,
                                                                  valid.file.near = validation.file.near,
                                                                  last.list = list.nearby.ngl,
                                                                  correct.ver = vertical.correction,
                                                                  path_series_ref =  path_series_main,
                                                                  path_series_near = path_series_nearby,
                                                                  duration = thres_period,
                                                                  percentage.required = no.day.per,
                                                                  nb_test = nb_test.ref,
                                                                  version.name = version_name,
                                                                  nearby.ver = nearby_ver,
                                                                  limit.type = "superposed",
                                                                  screen_value = screen_value) 
} else if (screen_value == ""){
  res.seg <-  get(load(file = paste0(path_results, "meta_compare/",nb_test=34, criterion="BM_BJ", "position_break.RData")))
  
  list.nearby.station.homo <- homogeneous.nearby.segment.multiple.all(name.station = name,
                                                                      valid.file.ref = validation.file.ref,
                                                                      valid.file.near = validation.file.near,
                                                                      last.list = list.nearby.ngl,
                                                                      correct.ver = vertical.correction,
                                                                      path_series_ref =  path_series_main,
                                                                      path_series_near = path_series_nearby,
                                                                      duration = thres_period,
                                                                      percentage.required = no.day.per,
                                                                      nb_test = nb_test.ref,
                                                                      version.name = version_name,
                                                                      nearby.ver = nearby_ver,
                                                                      limit.type = "superposed",
                                                                      screen_value = screen_value) 
}
file.name.list.cases = paste0(path_results, "attribution/", "list.nearby.station.homo", nb_test = nb_test.near,"-",
                              screen_value, tolerance_noise)

write.table(list.nearby.station.homo,  file = paste0(file.name.list.cases, ".txt") ,  sep = "\t", quote = FALSE, col.names = TRUE)
save(list.nearby.station.homo,  file = paste0(file.name.list.cases, ".RData"))

# significant test the change in mean -------------------------------------

# list.nearby.station.homo <- list.nearby.station.homo[which(list.nearby.station.homo$np.bef<60),]
for (name.test in list.test) {
  list.t.test <-significance.test.weighted.indi(homo.nearby = list.nearby.station.homo,
                                                path_series_ref = path_series_main,
                                                path_series_near = path_series_nearby,
                                                duration = thres_period, sig.level = significant.level,
                                                name.test = name.test,
                                                screening.range = screening.data.gps,
                                                screening = screening1,  nb_test = nb_test.ref,
                                                correct.ver = vertical.correction,
                                                version.name = version_name, nearby.ver = nearby_ver,
                                                limit.type = "superposed",
                                                screen_value = screen_value,
                                                check_correlation = 5) 
}


# creat the additional information file -----------------------------------

################## Add the noise information to the t-test results -------------------
list.t.test$test.offset = list.t.test$wei.mean.be - list.t.test$wei.mean.af
list.t.test$seg.offset <- NA

for (i in 1:nrow(list.t.test)){
  ind.brp = which(res.seg$station == list.t.test$ref.station[i] & res.seg$day == list.t.test$detected[i])
  offset.seg = res.seg$mean[ind.brp] - res.seg$mean[ind.brp+1]
  list.t.test$seg.offset[i] = - offset.seg 
}
data.last = list.t.test[,c(1,2,5,28)]
# H&V distance 
distances <- get(load( paste0(path_results, "attribution/", version_name, nearby.ver = "NGL", "distances-pairs.RData")))

data.last$hdist <- c()
data.last$vdist <- c()
for (i in 1:nrow(data.last)) {
  station.ref = data.last$ref.station[i]
  station.near = data.last$nearby.station[i]
  n.main = which(distances$main == station.ref & distances$nearby ==station.near)
  data.last$hdist[i] = distances$distances[n.main]
  data.last$vdist[i] = distances$ver.dist[n.main]
}

################## Add the representativeness error in additional file ##########
# GPS
data.last$range.var.gps <- c()
data.last$r.range.var.gps <- c()
duration = 365
for (i in 1:nrow(data.last)) {
  
  station.ref = data.last$ref.station[i]
  station.near = data.last$nearby.station[i]
  detection = data.last$detected[i]
  min = detection - duration +1
  max = detection + duration
  
  repre.ref = read.series.segment(path_series = path_series_main, station = station.ref , nb_test = nb_test.ref, criterion = "BM_BJ")
  month.var = repre.ref$var[which(repre.ref$date >= min & repre.ref$date < max)]
  month.var.diff = max(month.var, na.rm = TRUE) - min(month.var, na.rm = TRUE)
  r.month.var.diff = (max(month.var) - min(month.var))/ max(month.var)
  
  data.last$range.var.gps[i] <- round(month.var.diff, digits = 3)
  data.last$r.range.var.gps[i] <- round(r.month.var.diff, digits = 3)
  
}

write.table( data.last, file = paste0(path_results, 
                                      "attribution/attribution.segment.200km.500m.",version_name,
                                      nb_test.ref,nearby_ver,limit.type, screen_value, tolerance_noise, ".additional.txt"), 
             sep = "\t", quote = FALSE, col.names = TRUE)













