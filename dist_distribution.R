# investigate the distance between detections 

segment <-get(load(paste0(path_results,"meta_compare",screen.value = "1","/",nb_test = 34,criterion = "BM_BJ","position_break",".RData")))
station = unique(segment$station)
# station = name
path_raw = path_CODE_aux_ERA5_v1b
res = data.frame(matrix(NA, ncol = 2, nrow = 0))
for(l in 1:length(station)){ 
  time_series = get(load(paste0(path_raw,station[l],".RData")))
  seg.sta = segment[which(segment$station == station[l]),]
  seg.sta$d = seg.sta$end - seg.sta$begin +1
  if(nrow(seg.sta) >1){
    distance <- data.frame( station = rep(station[l], nrow(seg.sta)), dist = seg.sta$d)
    res <- rbind(res, distance)
  }
  # for(i in 1:nrow(seg.sta)){
  #   seg = time_series[c(seg.sta$begin[i]:seg.sta$end[i]),]
  #   # print(which(is.na(seg$signal) == TRUE)) verified that there is no NA value in the raw data
  #   
  # }
}

hist(res$dist[which(res$dist<365)], breaks = 100,
     main = "Histogram of the segment's length (<365) for the full CODE data set",
     xlab = "days")
hist(log(res$dist), breaks = 500, main = "Logarithmic distribution of the segment's length for the full CODE data set",
     xlab = " ")