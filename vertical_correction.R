# Vertical correction for the nearby data in forming 6 series of differences 
four.series <- list()
six.series <- list()

# Read nearby stations, distance, list of breakpoints
distances <- get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData")))
nearby_list <- get(load(file = paste0("list_nearby", nearby_ver, ".RData")))
nearby_list1 <- nearby_list[which(is.na(nearby_list$Level1) == FALSE),]

# Read list nearby 
for (j in c(1:nrow(nearby_list1))){
  station.ref.j = nearby_list1$ref.sta[j]
  series.ref <- read.series(path_series = path_series_main, station = station.ref.j, na.rm = 1, add.full = 0)
  seg.ref <- list.brp[which(list.brp$name == station.ref.j),]
  nearby.list.j = na.omit(unlist(nearby_list1[which(nearby_list1$ref.sta == station.ref.j ),], use.names = FALSE))[-1]
  
  for (k in c(1:nrow(seg.ref))) {
    breakpoint <- seg.ref$detected[k]
    begin.point <- breakpoint %m+% years(-1)   
    end.point <- breakpoint %m+% years(1)
    
    for (l in c(1:length(nearby.list.j))) {
      series.near <- read.series(path_series = path_series_nearby, station = nearby.list.j[l], na.rm = 1, add.full = 0)
      both <- data.frame() 
      four_series_frame <- data.frame() 

      # vertical correction -----------------------------------------------------
      n.main = which(distances$main == station.ref.j & distances$nearby == nearby.list.j[l])
      ver.dist = distances$ver.dist[n.main]   #vertical disstance here is main-nearby
      series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
      series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
      # join 2 data frame
      both = inner_join(series.ref,series.near, by = "date")  # checked
      both <- both[which(both$date > begin.point & both$date <= end.point),]
      
      if(nrow(both) >0){
        # save 4 series  ----------------------------------------------------------
        four.series[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]] <- list(date = both$date,
                                                                                                          GPS = both$GPS.x,
                                                                                                          ERA = both$ERAI.x,
                                                                                                          GPS1 = both$GPS.y,
                                                                                                          ERA1 = both$ERAI.y)
        four_series_frame <- data.frame(date = both$date,
                                        GPS = both$GPS.x,
                                        ERA = both$ERAI.x,
                                        GPS1 = both$GPS.y,
                                        ERA1 = both$ERAI.y)
        write.table(four_series_frame, 
                    file = paste0(path_results, "attribution/four_series/", station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l], ".txt"),
                    sep = "\t", quote = FALSE)
        
        six.diff <- list()
        six_series_frame <- data.frame (date = both$date)
        
        for (m in 1:6) {
          name.test = list.test[m]
          var1 = diff.var(c(name.test))[1]
          var2 = diff.var(c(name.test))[2]
          six.diff[[name.test]] <- both[,c(var1)] - both[,c(var2)]
          six_series_frame[name.test] <- both[,c(var1)] - both[,c(var2)]
        }
        six.diff$date <- both$date
        six.series[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]] <- six.diff

        write.table(six_series_frame, 
                    file = paste0(path_results, "attribution/six_series/",station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l], ".txt"),
                    sep = "\t", quote = FALSE)
      }
    }
  }
}

save(four.series, file = paste0(path_results,"attribution/four_main_series_1year_", nearby_ver,".RData"))
save(six.series, file = paste0(path_results,"attribution/six_diff_series_1year_", nearby_ver,".RData"))


