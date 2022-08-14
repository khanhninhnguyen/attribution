library("lubridate")       # Install & load lubridate package
window.thres <- 2 # window period
# we use the breakpoints from segmentation directly, with flag validation or outlier (80 days)
meta.compare =  get(load(file = paste0(path_results,"validation/",nb_test.ref,"-",criterion,"metacompa",screen.value="",".RData")))
# remove data in cluster of breakpoint or breakpoint from the ne --------

data <- get(load(file = paste0(path_results,"attribution/six_diff_series_",window.thres,"year_", nearby_ver,".RData")))

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

save(short.list, file = paste0(path_results,"attribution/removed_possible_crenel_",window.thres,"year_", nearby_ver,".RData"))
save(data, file = paste0(path_results,"attribution/six_diff_series_rm_crenel_", window.thres,"year_",nearby_ver,".RData"))





# Remove all breakpoint in +/- 1 year excluding the crenel ----------------
meta.compare.near =  get(load(file = paste0(path_results,"validation/",nb_test.near,"-",criterion,"metacompa",screen.value="",".RData")))
data.cr <- get(load(file = paste0(path_results,"attribution/six_diff_series_rm_crenel_",window.thres,"year_", nearby_ver,".RData")))

list.sim.brp <- data.frame(matrix(NA, ncol = 4, nrow = 0))
colnames(list.sim.brp) = c("case","sim.detected", "near.left", "near.right")

for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  station.near = substr(case.name ,start = 17, stop = 20)
  station.seg = meta.compare[which(meta.compare$name == station.ref),]
  nearby.seg = meta.compare.near[which(meta.compare.near$name == station.near),]
  
  # list breakpoints in the period
  begin = breakpoint %m+% years(-window.thres)   
  fin = breakpoint %m+% years(window.thres)
  
  list.brp.main = station.seg$detected[which(station.seg$detected > begin & station.seg$detected < fin & station.seg$noise ==0)]
  list.brp.near = nearby.seg$detected[which(nearby.seg$detected > begin & nearby.seg$detected < fin)]
  
  list.all.brp <- c(list.brp.main, list.brp.near)
  station.data = as.data.frame(data.cr[[i]])
  # check to remove the data in nearby breaks very closed
  con = 0
  if(length(list.brp.near) > 0){
    diff.nb = list.brp.near - breakpoint
    con = length(which(abs(diff.nb) < 10))
    if(con>0){
      ind.remove = which(abs(diff.nb) < 10)
      list.point.rm = list.brp.near[ind.remove]
      point.rm.ind = which.max(abs(list.point.rm - breakpoint))
      point.rm = list.point.rm[point.rm.ind]
      list.sim.brp[(nrow(list.sim.brp)+1),c(1:2)] <- c(i, as.character(point.rm))
      list.all.brp <- c(list.brp.main, list.brp.near[-ind.remove])
      period.rm <- c( breakpoint,point.rm)
      station.data[which(station.data$date < max(period.rm) & station.data$date > min(period.rm)), -7] <- NA
    }
  }
  # check to remove all other breaks 
  list.others = list.all.brp[which(list.all.brp > begin & list.all.brp < fin)]
  list.others.ord = sort(unique(list.others), decreasing = FALSE)
  if(length(list.others.ord )>1){ 
    ind.brp = which(list.others.ord == breakpoint)
    close.left = list.others.ord[ind.brp-1]
    close.right = list.others.ord[ind.brp+1]
    if(length(close.left) > 0){ # may create NA because the point does not exist
      station.data[which(station.data$date < close.left), -7] <- NA
      list.sim.brp[(nrow(list.sim.brp)+1),c(1,3)] <- c(i, as.character(close.left))
    }
    if(length(close.right) > 0){ # may create NA because the point does not exist
      station.data[which(station.data$date > close.right), -7] <- NA
      list.sim.brp[(nrow(list.sim.brp)+1),c(1,4)] <- c(i, as.character(close.right))
    }
  }
  data.cr[[i]] <- station.data
}

last.list.remove = list.sim.brp[rowSums(is.na(list.sim.brp)) != (ncol(list.sim.brp)-1), ]


save(last.list.remove, file = paste0(path_results,"attribution/restricted_by_closed_brp_",window.thres,"year_", nearby_ver,".RData"))
save(data.cr, file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",window.thres,"year_", nearby_ver,".RData"))





