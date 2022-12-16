# study the different veriecal correction 
distances <- get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData")))
new.ver = get(load(file = paste0(path_results,"attribution/four_main_series_",window.thres,"year_", nearby_ver,".RData")))
old.ver0 = get(load(file = paste0(path_results,"attribution/OLD_VERTICAL/four_main_series_",window.thres,"year_", nearby_ver,".RData")))

old.ver = old.ver0[names(new.ver)]

h = sapply(c(1:length(old.ver)), function(x){
  station.ref = substr(names(old.ver)[x] ,start = 1, stop = 4)
  station.seg = substr(names(old.ver)[x] ,start = 17, stop = 20)
  ind = which(distances$main == station.ref & distances$nearby == station.seg)
  distances$ver.dist[ind]
})
seg = c(1:3650)
gps.gps1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$GPS[seg], na.rm = TRUE) - mean(old.ver[[x]]$GPS1[seg], na.rm = TRUE))
gps.gps2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$GPS[seg], na.rm = TRUE) - mean(new.ver[[x]]$GPS1[seg], na.rm = TRUE))

gps.gps1[which(h<0)] = -gps.gps1[which(h<0)] 
gps.gps2[which(h<0)] = -gps.gps2[which(h<0)] 
plot(abs(h), gps.gps1)