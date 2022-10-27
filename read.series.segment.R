# Read series and convert to date
read.series <- function(path_series,station, na.rm, add.full){
  series.dat = get(load(paste0(path_series,station,".RData")))
  DateConvert.near = as.Date(series.dat$date, format="%Y-%m-%d" )
  series.dat$date <- DateConvert.near
  if(na.rm == 1){
    series.dat <- na.omit(series.dat)
  }
  if(add.full == 1){
    series.dat <- tidyr::complete(series.dat, date = seq(min(series.dat$date), max(series.dat$date), by = "day"))
  }
  return(series.dat)
}

# function use to read the series and segmentation results -----------
read.series.segment <- function(path_series,station,nb_test, criterion, add.muvar){
  time_series = get(load(paste0(path_series,station,".RData")))
  DateConvert = as.Date(time_series$date, format="%Y-%m-%d" )
  time_series$date <- DateConvert
  segment.res <-get(load(paste0(path_homo,nb_test,"/homo_",station,nb_test,".RData")))
  mu <- rep(NA, nrow(time_series))
  # add the mean of segments
  if( add.muvar == 1){
    if (criterion!= ""){
      position.brp <- segment.res$seg[[criterion]]
      for (i in 1:nrow(position.brp)) {
        mu[position.brp$begin[i] : position.brp$end[i]] <- position.brp$mean[i]
      }
      mu[which(is.na(time_series$signal))] <- NA
      time_series$mu <-mu
    }
    # add the variance to each point
    time_series$var <- NA
    for (j in 1:12) {
      time_series$var[which(as.numeric(time_series$month)== j)] <- segment.res$variances[j]
    }
    time_series$var[which(is.na(time_series$signal))] <- NA
  }
  return(c(time_series,segment.res))
}
# read series and segmentation results after screening -----------
read.series.segment.screening <- function(path_series,station,nb_test, criterion){
  time_series = get(load(paste0(path_series,station,".RData")))
  DateConvert = as.Date(time_series$date, 'GMT')
  time_series$date <- DateConvert
  segment.res <-get(load(paste0(path_homo,nb_test,"/homo_",station,nb_test,".RData")))
  mu <- rep(NA, nrow(time_series))
  # add the mean of segments
  screened <- get(load( paste0(path_results,"meta_compare1/",nb_test,criterion,"position_break.RData")))
  screened.station <- screened[which(screened$station == station),]
  #updatee the number of breakpoints after screening
  segment.res$K[[criterion]] <- nrow(screened.station)
  segment.res$seg[[criterion]] <- screened.station[,1:3]
  if (criterion!= ""){
    position.brp <- segment.res$seg[[criterion]]
    for (i in 1:nrow(position.brp)) {
      mu[position.brp$begin[i] : position.brp$end[i]] <- position.brp$mean[i]
    }
    mu[which(is.na(time_series$signal))] <- NA
    time_series$mu <-mu
  }
  
  time_series$var <- NA
  for (j in 1:12) {
    time_series$var[which(as.numeric(time_series$month)== j)] <- segment.res$variances[j]
  }
  time_series$var[which(is.na(time_series$signal))] <- NA
  
  return(c(time_series,segment.res))
}
# add the results of ICL criterion
read.series.segment.ICL <- function(path_series, station, nb_test, cri){
  res <- read.series.segment(path_series,station,nb_test, criterion = "")
  ebs <- mget(load(file = paste0(path_results,"EBS/EBS",station,".RData")))
  res$K$ICL1 <-ebs$segment$K$ICL1
  res$K$ICL2 <-ebs$segment$K$ICL2
  list.criteria <- c("BM_BJ", "BM_slope", "Lav", "mBIC")
  
  #  if ((cri %in% list.criteria) == FALSE){
  #   K <- ebs$segment$K[[cri]]
  #   brp<- c()
  #   for (nb.brp in 1:(K-1)) {
  #     brp <- c(brp, which.max(ebs$segment[[cri]][[nb.brp]]))
  #   }
  #   # check is there 
  #   brp <- unique(brp)
  #   K <- length(brp)+1
  #   # find the corrected index because ebs run with the non na series 
  #   data.all <- setNames(data.frame(res$date,res$signal), c("date","signal"))
  #   setDT(data.all, keep.rownames = TRUE)[]
  #   end.date <- data.all$date[brp]
  #   list.brp.non.na <- subset(data.all, date %in%end.date)
  #   begin <-  as.numeric(c(1, list.brp.non.na$rn))
  #   end <- c()
  #   for (end.p in 1:(K-1)) {
  #     subdata <- data.all[which(data.all$date < list.brp.non.na$date[end.p])]
  #     end <- c(end,subdata$rn[nrow(subdata)])
  #   }
  #   end <- as.numeric(c(end, length(data.all$date)))
  #   add.na.mean <- unlist(res$Tot[[K]]$Tmu$mean, use.names=FALSE)
  #   add.na <- data.frame(begin,end,add.na.mean)
  #   colnames(add.na) <- c("begin","end","mean")
  #   res$seg[[cri]] <- add.na
  # }
  res$seg$ICL1 <- res$Tot[[res$K$ICL1]]$Tmu
  res$seg$ICL2 <- res$Tot[[res$K$ICL2]]$Tmu
  res$funct$ICL1 <- res$funct$BM_BJ
  res$funct$ICL2 <- res$funct$BM_BJ
  return(res)
}

read.series.segment.screening.ICL <- function(path_series,station,nb_test, criterion){
  segment.res <- read.series.segment.ICL(path_series,station,nb_test,criterion)
  screened <- get(load( paste0(path_results,"meta_compare1/",nb_test,criterion,"position_break.RData")))
  screened.station <- screened[which(screened$station == station),]
  #updatee the number of breakpoints after screening
  segment.res$K[[criterion]] <- nrow(screened.station)
  segment.res$seg[[criterion]] <- screened.station[,1:3]
  return(segment.res)
}


# merge 2 series ----------------------------------------------------------

merge_series <- function(path1, path2, station1, station2, duration, detection){
  series1 = read.series(path_series = path1, station = station1, na.rm = 0, add.full = 0)
  series2 = read.series(path_series = path2, station = station2, na.rm = 0, add.full = 0)
  
  both = inner_join(series1,series2, by = "date") 
  both$signal = both$GPS.y-both$GPS.x
  both$year <- both$year.x
  both$month <- both$month.x
  both <- tidyr::complete(both, date = seq(min(both$date), max(both$date), by = "day"))
  if (duration != 0){
    begin.day = detection - duration
    end.day = detection + duration
    cutoff <- both[which(both$date >= begin.day & both$date <= end.day),]
    cutoff <- na.omit(cutoff)
    out = cutoff
  } else{
    out = both
  }
  return(out)
}
  
# brp<- c()
# for (nb.brp in 1:length(ebs$segment$ICL1)) {
#   brp <- c(brp, which.max(ebs$segment$ICL1[[nb.brp]]))
# }
# data.all <- setNames(data.frame(res$date,res$signal), c("date","signal"))
# setDT(data.all, keep.rownames = TRUE)[]
# end.date <- data.all$date[brp]
# list.brp.non.na <- subset(data.all, date %in%end.date)
# begin <-  c(1, list.brp.non.na$rn)
# end <- c()
# for (end.p in 1:length(ebs$segment$ICL1)) {
#   subdata <- data.all[which(data.all$date < list.brp.non.na$date[end.p])]
#   end <- c(end,subdata$rn[nrow(subdata)])
# }
# end <- c(end, length(data.all$date))
# 
