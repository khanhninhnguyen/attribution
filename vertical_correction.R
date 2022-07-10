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
  nearby.list.j = na.omit(unlist(nearby_list1[which(nearby_list1$ref.sta == ref.station.j ),], use.names = FALSE))[-1]
  
  for (k in c(1:nrow(seg.ref))) {
    brp <- seg.ref$detected[k]
    begin.point <- brp - 365
    end.point <- brp + 365
    for (l in c(1:length(nearby.list.j))) {
      series.near <- read.series(path_series = path_series_nearby, station = nearby.list.j[l], na.rm = 1, add.full = 0)

      # vertical correction -----------------------------------------------------
      n.main = which(distances$main == station.ref.j & distances$nearby == nearby.list.j[l])
      ver.dist = distances$ver.dist[n.main]   #vertical disstance here is main-nearby
      series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
      series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
    }
  }
  
}






