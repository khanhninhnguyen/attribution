# Vertical correction for the nearby data in forming 6 series of differences 
# output include the four series of the raw without matching day main and nearby
window.thres <- 10

vertical_correction <- function(window.thres, path_results, version_name, nearby_ver, path_series_nearby, path_series_main){
  four.series <- list()
  six.series <- list()
  
  # Read nearby stations, distance, list of breakpoints
  distances <- get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData")))
  nearby_list <- get(load(file = paste0(path_results, "attribution/","list_nearby", nearby_ver, ".RData")))
  nearby_list1 <- nearby_list[which(is.na(nearby_list$Level1) == FALSE),]

  for (j in c(1:nrow(nearby_list1))){
    station.ref.j = nearby_list1$ref.sta[j]
    series.ref <- read.series(path_series = path_series_main, station = station.ref.j, na.rm = 1, add.full = 0)
    seg.ref <- list.brp[which(list.brp$name == station.ref.j),]
    nearby.list.j = na.omit(unlist(nearby_list1[which(nearby_list1$ref.sta == station.ref.j ),], use.names = FALSE))[-1]
    
    for (k in c(1:nrow(seg.ref))) {
      breakpoint <- seg.ref$detected[k]
      # begin.point <- breakpoint %m+% years(-window.thres)   
      # end.point <- breakpoint %m+% years(window.thres)
      begin.point <- breakpoint - (window.thres*365 -1)
      end.point <- breakpoint + window.thres*365
      
      for (l in c(1:length(nearby.list.j))) {
        series.near <- read.series(path_series = path_series_nearby, station = nearby.list.j[l], na.rm = 1, add.full = 0)
        both <- data.frame() 
        four_series_frame <- data.frame() 
        
        # vertical correction -----------------------------------------------------
        n.main = which(distances$main == station.ref.j & distances$nearby == nearby.list.j[l])
        ver.dist = distances$ver.dist[n.main]   #vertical disstance here is main-nearby
        # old version
        # series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
        # series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
        # new version
        ver.cor<-tryCatch(
          {
            read.table(
              file = paste0(path_data_support, "vertical-correction/", "cor_iwv_", nearby.list.j[l], "-", station.ref.j, ".txt"),
              header = TRUE
            )
          }, 
          error = function(e) {
            print(paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l]))
            NULL
          }
        )
        if(is.null(ver.cor) == FALSE){
          names(ver.cor)[1:4] <- c("year", "month", "a", "b")
          ver.cor.full = series.ref
          ver.cor.full[ , c(6,7)] <- as.data.frame(apply(ver.cor.full[ , c(6,7)], 2, as.numeric))
          ver.cor.full = left_join(ver.cor.full, ver.cor, by = c("month", "year"))
          ver.cor.full$GPS = ver.cor.full$GPS*ver.cor.full$a + ver.cor.full$b
          ver.cor.full$ERAI = ver.cor.full$ERAI*ver.cor.full$a + ver.cor.full$b
          series.ref0 = series.ref
          # series.ref0[c("GPS", "ERAI")] <- ver.cor.full[c("GPS", "ERAI")] # test the case where no vertical correction applied
          
          # join 2 data frame
          # both <- both[which(both$date > begin.point & both$date <= end.point),]
          # limit the raw series 
          series.ref1 <- series.ref0[which(series.ref0$date >= begin.point & series.ref0$date <= end.point),]
          series.near1 <- series.near[which(series.near$date >= begin.point & series.near$date <= end.point),]
          series.ref1 <- tidyr::complete(series.ref1, date = seq(begin.point, end.point, by = "day"))
          series.near1 <- tidyr::complete(series.near1, date = seq(begin.point, end.point, by = "day"))
          both = inner_join(series.ref1,series.near1, by = "date")  # checked

          condi1 = length(which(is.na(both$signal.y) == FALSE))
          if(condi1 >100){
            # save 4 series  ----------------------------------------------------------
            four_series_frame <- data.frame(date = both$date,
                                            GPS = both$GPS.x,
                                            ERA = both$ERAI.x,
                                            GPS1 = both$GPS.y,
                                            ERA1 = both$ERAI.y)
            four.series[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]] <- four_series_frame
            
            write.table(four_series_frame, 
                        file = paste0(path_results, "attribution/four_series/", station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l], window.thres, "yr.txt"),
                        sep = "\t", quote = FALSE)
            
            six_series_frame <- data.frame (date = both$date)
            
            for (m in 1:6) {
              name.test = list.test[m]
              var1 = diff.var(c(name.test))[1]
              var2 = diff.var(c(name.test))[2]
              six_series_frame[name.test] <- both[,c(var1)] - both[,c(var2)]
            }
            six.series[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]] <- six_series_frame
            #ADD 6 SERIES OF DIFFERENECES WITH FULL DATA INTO R DATA, IN CASE WE DON'T WANT TO LOSE INFORMATION DUE TO THE MATCHING BETWEEN MAIN AND NEARBY
            write.table(six_series_frame, 
                        file = paste0(path_results, "attribution/six_series/",station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l], window.thres, "yr.txt"),
                        sep = "\t", quote = FALSE)
            
          }
        }
      }
    }
  }
  
  save(four.series, file = paste0(path_results,"attribution/vertical_correction/four_main_series_",window.thres,"year_", nearby_ver,"0.RData"))
  save(six.series, file = paste0(path_results,"attribution/vertical_correction/six_diff_series_",window.thres,"year_", nearby_ver,"0.RData"))
  
}



