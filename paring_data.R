# This function is used for pairing data 
window.thres = 2
data.cr = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
data.all <- data.cr
n0 = length(data.all)
data.out <- list()
for (i in c(1:n0)) {
  data.i = data.all[[i]]
  # data.i <- tidyr::complete(data.i, date = seq(min(data.i$date), max(data.i$date), by = "day"))
  
  case.name = names(data.all)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  
  day.sta = breakpoint %m+% years(-1)
  day.end = breakpoint %m+% years(+1)
    
  before = na.omit(data.i[which(data.i $date >= day.sta & data.i $date <= breakpoint),])
  after = na.omit(data.i[which(data.i $date > breakpoint & data.i $date <= day.end),])
  
  listday1 = after$date %m+% years(-1)
  listday2 = before$date[which(before$date %in% listday1 == TRUE )]
  listday3 = after$date[which(listday1 %in% listday2 == TRUE & duplicated(listday1) == FALSE)]
  
  if(length(listday2) > 100){
    before1 = before[which(before$date %in% listday2 == TRUE),]
    after1 = after[which(after$date %in% listday3 == TRUE),]
    # compute the pair data
    out <- data.frame( date = before1$date)
    for (j in c(1:6)) {
      out[list.test[j]] <- after1[list.test[j]] - before1[list.test[j]]
    }
   data.out [[case.name]] <- out 
  } 
  
}

save(data.out, file = paste0(path_results,"attribution/data.all_1year_", nearby_ver,"paired.RData"))

Y = data.out$`gope.2009-05-08.zdib`
save(Y, file = paste0(path_results,"attribution/data.all_1years_", nearby_ver,".gope.2009-05-08.zdib.paired.RData"))
