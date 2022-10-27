# investigate the distance between detections 

nb.test=44
segment <-get(load(paste0(path_results,"meta_compare",screen.value = "","/",nb_test = nb.test,criterion = "BM_BJ","position_break",".RData")))

# station = name_main
path_raw = path_NGL_v2_R
list.all =list.files(paste0(path_homo,nb.test,"/"))
station = substr(list.all,start = 6, stop = 9)
# station = station[station != "ouag"]
# station = name_main

res = data.frame(matrix(NA, ncol = 2, nrow = 0))
for(l in c(1:length(station))){ 
  time_series = read.series.segment(path_raw, station[l], nb_test = nb.test, criterion = "BM_BJ", add.muvar = 0)
  seg.sta = time_series$seg$BM_BJ 
  if(nrow(seg.sta) >=1){
    d = sapply(c(1:nrow(seg.sta)), function(x) length(na.omit(time_series$signal[seg.sta$begin[x]:seg.sta$end[x]])))
    distance <- data.frame( station = rep(station[l], nrow(seg.sta)), dist = d)
    res <- rbind(res, distance)
  }
  # for(i in 1:nrow(seg.sta)){
  #   seg = time_series[c(seg.sta$begin[i]:seg.sta$end[i]),]
  #   # print(which(is.na(seg$signal) == TRUE)) verified that there is no NA value in the raw data
  #   
  # }
}

hist(res$dist[res$dist<=365], breaks = 100,
     main = "Histogram of the segment's length (<365) for all stations in NGL",
     xlab = "days")
hist(log(res$dist), breaks = 500, main = "Logarithmic distribution of the segment's length for the full CODE data set",
     xlab = " ")

