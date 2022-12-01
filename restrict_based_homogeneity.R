library("lubridate")       # Install & load lubridate package
window.thres <- 10 # window period
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
      station.data[remove.ind,c(list.test[list.test != ("gps.era")])] <- NA
      data[[i]] <- station.data
    }
  }
}

### Remove the data set that have very small number of points (sum<100)

length.data = sapply(c(1:length(data)), function(x){ nrow(na.omit(data[[x]])) 
})
name.removed.stations = names(data)[which(length.data <100)]
data[name.removed.stations] <- NULL

# list cluster removed: 

short.list = list.cluster.removed[!duplicated(list.cluster.removed),]
all.cases.name = names(data)
all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
list.case.remove = paste(short.list$station, short.list$end, sep = ".")
case.remove.ind = which(all.cases.ind %in% list.case.remove)
data.last = data[-case.remove.ind]

save(short.list, file = paste0(path_results,"attribution/removed_possible_crenel_",window.thres,"year_", nearby_ver,".RData"))
save(data.last, file = paste0(path_results,"attribution/six_diff_series_rm_crenel_", window.thres,"year_",nearby_ver,".RData"))

#NOTE: data after restriction is lost of date in the crenel, in the next step, we need to complete data by date again 



# Remove all breakpoint in +/- 1 year excluding the crenel different restriction for different test----------------
list.cre = get(load(file = paste0(path_results,"attribution/removed_possible_crenel_",window.thres,"year_", nearby_ver,".RData")))
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
  begin = breakpoint - (window.thres*365-1)
  fin = breakpoint + window.thres*365
  
  list.brp.main = station.seg$detected[which(station.seg$detected >= begin & station.seg$detected < fin)]
  # remove all second points of crenels
  cre = which(list.cre$station == station.ref)
  if(length(cre) >0){
    station.cre = as.Date(list.cre$end[cre])
    list.brp.main =  list.brp.main[which(list.brp.main %in% station.cre == FALSE)]
  }
  list.brp.near = nearby.seg$detected[which(nearby.seg$detected >= begin & nearby.seg$detected < fin)]
  
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
      station.data[which(station.data$date < max(period.rm) & station.data$date > min(period.rm)), list.test[c(2,3,4,6)]] <- NA
    }
  }
  # check to remove all other breaks in the four series combined the main and nearby 
  list.others = list.all.brp[which(list.all.brp >= begin & list.all.brp < fin)]
  list.others.ord = sort(unique(list.others), decreasing = FALSE)
  if(length(list.others.ord )>1){ 
    ind.brp = which(list.others.ord == breakpoint)
    close.left = list.others.ord[ind.brp-1]
    close.right = list.others.ord[ind.brp+1]
    if(length(close.left) > 0){ # may create NA because the point does not exist
      station.data[which(station.data$date < close.left),list.test[c(2,3,4,6)]] <- NA
      list.sim.brp[(nrow(list.sim.brp)+1),c(1,3)] <- c(i, as.character(close.left))
    }
    if(length(close.right) > 0){ # may create NA because the point does not exist
      station.data[which(station.data$date > close.right), list.test[c(2,3,4,6)]] <- NA
      list.sim.brp[(nrow(list.sim.brp)+1),c(1,4)] <- c(i, as.character(close.right))
    }
  }
  # check to remove all other breaks in the main series
  list.main.ord = sort(unique(list.brp.main), decreasing = FALSE)
  if(length( list.main.ord  )>1){ 
    ind.brp = which( list.main.ord  == breakpoint)
    close.left =  list.main.ord [ind.brp-1]
    close.right =  list.main.ord [ind.brp+1]
    if(length(close.left) > 0){ # may create NA because the point does not exist
      station.data[which(station.data$date < close.left),list.test[c(1)]] <- NA
    }
    if(length(close.right) > 0){ # may create NA because the point does not exist
      station.data[which(station.data$date > close.right), list.test[c(1)]] <- NA
    }
  }
  # check to remove all other breaks in the nearby series
  list.near.ord = sort(unique(c( breakpoint, list.brp.near)), decreasing = FALSE)
  if(length( list.near.ord  )>1){ 
    ind.brp = which( list.near.ord  == breakpoint)
    close.left =  list.near.ord [ind.brp-1]
    close.right =  list.near.ord [ind.brp+1]
    if(length(close.left) > 0){ # may create NA because the point does not exist
      station.data[which(station.data$date < close.left),list.test[c(5)]] <- NA
    }
    if(length(close.right) > 0){ # may create NA because the point does not exist
      station.data[which(station.data$date > close.right), list.test[c(5)]] <- NA
    }
  }
  data.cr[[i]] <- station.data
}

last.list.remove = list.sim.brp[rowSums(is.na(list.sim.brp)) != (ncol(list.sim.brp)-1), ]


save(last.list.remove, file = paste0(path_results,"attribution/restricted_by_closed_brp_",window.thres,"year_", nearby_ver,".RData"))
save(data.cr, file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",window.thres,"year_", nearby_ver,".RData"))


# checked results, all removed points is replaced by NA including the date column
#######   cCHECK WHY RESULT CONTAIN DATA FRAME AND LIST, WHY NOT CONSISTENT??????????????

