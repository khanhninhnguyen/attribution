library("lubridate")       # Install & load lubridate package

# we use the breakpoints from segmentation directly, with flag validation or outlier (80 days)
meta.compare =  get(load(file = paste0(path_results,"validation/",nb_test.ref,"-",criterion,"metacompa",screen.value="",".RData")))
# remove data in cluster of breakpoint or breakpoint from the ne --------

data <- get(load(file = paste0(path_results,"attribution/six_diff_series_1year_", nearby_ver,".RData")))

list.cluster.removed <- data.frame(matrix(NA, ncol = 3, nrow = 0))
colnames(list.cluster.removed) = c("station","start", "end")
for (i in c(1:length(data))) {
  case.name = names(data)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  station.seg = meta.compare[which(meta.compare$name == station.ref),]
  ind = which(station.seg$detected == breakpoint)
  noise = station.seg$noise[ind]
  rest = c(station.seg$noise[c(ind:nrow(station.seg))],0)
  if (noise > 0){
    j = 1
    last.ind = j
    while(rest[j+1] > rest[j]){
      j = j+1
      last.ind = j
    }
    length.cluster = rest[last.ind] 
    start.cluster.ind = ind - 1 + last.ind - (length.cluster-1)
    end.cluster.ind = start.cluster.ind + length.cluster - 1
    start.cluster = station.seg$detected[start.cluster.ind]
    end.cluster = station.seg$detected[end.cluster.ind]
    list.cluster.removed[(nrow(list.cluster.removed)+1),] <- c(station.ref, as.character(start.cluster), as.character(end.cluster))
    for (m in 1:6) {
      name.test = list.test[m]
      station.data = as.data.frame(data[[i]])
      remove.ind = which(station.data$date <= end.cluster & station.data$date > start.cluster)
      station.data[remove.ind,-7] <- NA
      data[[i]] <- station.data
    }
  }
}


# list cluster removed: 

short.list = list.cluster.removed[!duplicated(list.cluster.removed),]

save(short.list, file = paste0(path_results,"attribution/removed_possible_crenel_", nearby_ver,".RData"))
save(data, file = paste0(path_results,"attribution/six_diff_series_1year_rm_crenel_", nearby_ver,".RData"))





# Remove all breakpoint in +/- 1 year excluding the crenel ----------------
meta.compare.near =  get(load(file = paste0(path_results,"validation/",nb_test.near,"-",criterion,"metacompa",screen.value="",".RData")))
data.cr <- get(load(file = paste0(path_results,"attribution/six_diff_series_1year_rm_crenel_", nearby_ver,".RData")))

list.cluster.removed <- data.frame(matrix(NA, ncol = 3, nrow = 0))
colnames(list.cluster.removed) = c("station","start", "end")
for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  station.near = substr(case.name ,start = 17, stop = 20)
  station.seg = meta.compare[which(meta.compare$name == station.ref),]
  nearby.seg = meta.compare.near[which(meta.compare.near$name == station.near),]
  
  # list breakpoints in the period
  begin = breakpoint %m+% years(-1)   
  fin = breakpoint %m+% years(1)
  
  list.brp.main = station.seg$detected[which(station.seg$detected > begin & station.seg$detected < fin & station.seg$noise ==0)]
  list.brp.near = station.near$detected[which(nearby.seg$detected > begin & nearby.seg$detected < fin)]
  
}
i= 141
