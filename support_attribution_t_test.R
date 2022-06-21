# ########## this file contains all functions used in the attribution ordered as in the attribution program 

# search for nearby stations ----------------------------------------------
# Reorder the name of station 
reorder_nearby = function(list.nearby.station.bef){
  
  nb.station.checked <- sapply(c(1:nrow(list.nearby.station.bef)), function(x) length(which(is.na(list.nearby.station.bef[x,])==F))) -1
  # Reoder the name of nearby station 
  for (i in 1:nrow(list.nearby.station.bef)){
    b <- rep(as.character(NA), (length(list.nearby.station.bef)-1))
    d <- sort(list.nearby.station.bef[i,2:length(list.nearby.station.bef)])
    if (length(d)>0){
      b[1:(length(d))] <- d
      list.nearby.station.bef[i,2:length(list.nearby.station.bef)] <-b
      
    }
  }
  res <- list.nearby.station.bef
  return(res)
}
# create list_nearby for NGL data 
nearby_list_ngl = function(path_results,  nearby.ver, version_name, name_main){
  colname = c("order", "ref.sta", "detect", "near1", "near2", "near3", "near4", "near5", "near6", "near7", "near8","near9", "near10")
  list.nearby.ngl = read.table(file = paste0(path_results, "attribution/list.name.nearby.NGL.txt"), col.names = colname, stringsAsFactors = FALSE)
  list.hdist.ngl = read.table(file = paste0(path_results, "attribution/list.hdist.nearby.200km.500m.NGL.txt"), col.names = colname, stringsAsFactors = FALSE)
  list.vdist.ngl = read.table(file = paste0(path_results, "attribution/list.vdist.nearby.200km.500m.NGL.txt"), col.names = colname, stringsAsFactors = FALSE)
  list.nearby.ngl = list.nearby.ngl[,c(-1, -3)]
  list.hdist.ngl = list.hdist.ngl[,c(-1, -3)]
  list.vdist.ngl = list.vdist.ngl[,c(-1, -3)]
  name = name_main
  
  # to have a consistent type if list.nearby
  list.nb = list()
  for (i in (1:81)) {
    nb.sta = list.nearby.ngl[which(list.nearby.ngl$ref.sta == name[i]),][,-1]
    nb = unique(as.vector(as.matrix(nb.sta)))
    nb[nb== "loz1"] <- NA
    list.nb[[name[i]]] <- na.omit(nb)
  }
  m1 <- max(lengths(list.nb))
  list.nearby <- as.data.frame(do.call(rbind, lapply(list.nb, `length<-`, m1)), stringsAsFactors= FALSE)
  names(list.nearby) <- paste0("Level", seq_along(list.nearby))
  list.nearby$ref.sta <- rownames(list.nearby)
  rownames(list.nearby) <- NULL
  list.nearby <- list.nearby[,c(29, 1:28)]
  list.distances <- list(main = c(), nearby = c(), distances = c(), ver.dist = c())
  for (p in 1:nrow(list.nearby)){
    ind = length(which(is.na(list.nearby[p, -1]) == FALSE))
    if (ind!=0){
      for (t in 1:ind) {
        list.distances$main <- c(list.distances$main, list.nearby$ref.sta[p])
        list.distances$nearby <- c( list.distances$nearby, list.nearby[p,t+1])
        list.nb.sta = list.nearby.ngl[which(list.nearby.ngl$ref.sta ==  list.nearby$ref.sta[p]),]
        position = which(list.nb.sta == list.nearby[p,t+1], arr.ind=TRUE)
        
        hdist.list = unique(as.numeric(list.hdist.ngl[which(list.nearby.ngl$ref.sta ==  list.nearby$ref.sta[p]),][position]))
        vdist.list = unique(as.numeric(list.vdist.ngl[which(list.nearby.ngl$ref.sta ==  list.nearby$ref.sta[p]),][position]))
        
        list.distances$distances <- c(list.distances$distances,  hdist.list)
        list.distances$ver.dist <- c(list.distances$ver.dist,  vdist.list)
      }
      
    }
  }
  save(list.distances, file = paste0(path_results, "attribution/", version_name,  nearby.ver, "distances-pairs.RData"))
  
  return(list.nearby)
}

# This function used to search for nearby stations based on the horizontal radius and vertical distances 
nearby_search_distance = function(coor.file, list.name, list.gnss, horizontal, vertical, version_name, nearby.ver){
  library(rgeos)
  library(geosphere)
  coor = read.table(coor.file, header = TRUE)
  colnames(coor) = c("name","lat","lon","height","altitude")
  last.list <- data.frame(matrix(NA, nrow = length(list.name) , ncol = 10))
  colnames(last.list) <- c("ref.sta","nearest","near1","near2","near3","near4","near5","near8","near9","near10")
  for (i in 1:nrow(coor)){
    if (coor$lon[i] >180){
      coor$lon[i] <- coor$lon[i]-360
    }
  }
  coor$name <- tolower(coor$name)
  lonlat = data.frame(coor$lon,coor$lat)
  lonlat = SpatialPoints(lonlat, proj4string=CRS("+proj=longlat +datum=WGS84"))
  count = 0
  length.nearby =c() # to check at last
  list.distances <- list(main = c(), nearby = c(), distances = c(), ver.dist = c())
  for (i in 1:length(list.name)){
    ind = which(coor$name == list.name[i])
    ref = lonlat[ind]
    coor.m =coor[-ind,]
    distance.ver = coor[ind,]$altitude - coor[-ind,]$altitude # vertical distance main - nearby
    lonlat1 = data.frame(coor.m$lon,coor.m$lat)
    spdf <- SpatialPointsDataFrame(coords = lonlat1, data = coor.m,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
    point = SpatialPoints(lonlat1, proj4string=CRS("+proj=longlat +datum=WGS84")) # checked
    distance = distGeo(point, ref)
    spdf$dist <- distance
    spdf$dist.ver <- distance.ver
    sort.fr <- spdf[order(distance),]
    sort.fr$dist.ver = as.numeric(sort.fr$dist.ver)*1000
    sort.fr$dist = as.numeric(sort.fr$dist)/1000
    if (vertical!=0){
      sort.fr = sort.fr[which(abs(sort.fr$dist.ver)<=vertical),]
    }else if (vertical ==0){
      sort.fr = sort.fr
    }
    posible = sort.fr[which(sort.fr$dist<=horizontal),]
    list.near = c()
    for (p in 1:nrow(posible)){
      era5.ind = length(which(list.gnss == posible$name[p]))
      if (era5.ind!=0){
        list.distances$main <- c(list.distances$main, list.name[i])
        list.distances$nearby <- c( list.distances$nearby, posible$name[p])
        list.distances$distances <- c(list.distances$distances, posible$dist[p])
        list.distances$ver.dist <- c(list.distances$ver.dist, posible$dist.ver[p])
        list.near = c(list.near, posible$name[p])
      }
    }
    last.list[i,1:(length(list.near)+1)] <- c(list.name[i],list.near)
    length.nearby = c(length.nearby, length(list.near))
    if(length(list.near)!=0){
      count = count+1
    }
  }
  save(list.distances, file = paste0(path_results, "attribution/", version_name, nearby.ver, "distances-pairs.RData"))
  return(last.list)
} # search for nearby station 

# check is the tested change point in a cluster -------------------------------------
check_cluster <- function(list.cluster, station, list.ind.cluster, min.d, max.d, breakpoint.c){
  cluster = list.cluster[which(list.cluster$name == station),]
  ind.sta = list.ind.cluster[which(list.ind.cluster$name == station),]
  if( nrow(ind.sta) > 0){
    cluster.ind = sapply(c(1:nrow(ind.sta)), function(x) which(cluster$end[ind.sta$deb[x]] < breakpoint.c & breakpoint.c < cluster$begin[ind.sta$fin[x]]), simplify = "array")
    cluster.pre = which(cluster.ind==1)
    if(length(cluster.pre) != 0){
      cls.begin = cluster$end[ind.sta$deb[cluster.pre]]
      cls.end = cluster$begin[ind.sta$fin[cluster.pre]]
    }else {
      cls.begin = max.d
      cls.end = min.d
    }
  } else {
    cls.begin = max.d
    cls.end = min.d
  }
  return(list(cls.begin = cls.begin, cls.end = cls.end))
}




# additional functions ----------------------------------------------------

mymax <- function(...,def= NA,na.rm=FALSE)
  if(!is.infinite(x<-suppressWarnings(max(...,na.rm=na.rm)))) x else def

mymin <- function(...,def=NA,na.rm=FALSE)
  if(!is.infinite(x<-suppressWarnings(min(...,na.rm=na.rm)))) x else def
diff.var <- function(name.test){
  if(name.test == "gps.gps"){
    varname = c("GPS.x", "GPS.y")
  }
  if(name.test == "gps.era"){
    varname = c("GPS.x", "ERAI.x")
  }
  if(name.test == "gps1.era"){
    varname = c("GPS.y", "ERAI.x")
  }
  if(name.test == "gps.era1"){
    varname = c("GPS.x", "ERAI.y")
  }
  if(name.test == "gps1.era1"){
    varname = c("GPS.y", "ERAI.y")
  }
  if(name.test == "era.era"){
    varname = c("ERAI.x", "ERAI.y")
  }
  return(varname)
}
check.sig <- function(x, sig.level){
  if (is.na(x) == TRUE){
    sig = NA
  }else{
    if (abs(x)<= sig.level){
      sig = 1
    } else{
      sig = 0
    } 
  }
  return(sig)
}
# form 6 series of differences --------------------------------------------

# find the homogeneous for all breakpoints withou screening 
homogeneous.nearby.segment.multiple.all = function( name.station, valid.file.ref, valid.file.near, last.list, path_series_ref, 
                                                    path_series_near, duration, percentage.required, correct.ver, nearby.ver,
                                                    version.name, nb_test, limit.type, screen_value  ){
  # Read the list of breakpoints in the main and the nearby 
  valid.ref = read.table(valid.file.ref, header = TRUE, sep = "\t",  colClasses = c("character", "Date","Date", "character","numeric",
                                                                                    "integer", "integer")
  )
  valid.near = read.table(valid.file.near, header = TRUE, sep = "\t",  colClasses = c("character", "Date","Date", "character","numeric",
                                                                                      "integer", "integer")
  )
  # List of the crenel
  list.cluster <-  read.table( file = paste0(path_results,"screening/lmin1.",criterion,".date_mean_weight_np.",version.name,nb_test.ref,".txt"),
                               sep = "\t", header = T,  colClasses = c("character", "Date","Date", "numeric","numeric",
                                                                       "integer", "integer")
  )
  list.ind.cluster = read.table(file = paste0(path_results,"screening/lmin1.",criterion,".screening.", version.name,nb_test.ref,".txt"),
                                sep="\t",  header = T, colClasses = c("character", "Date","integer", "integer",
                                                                      "integer", "integer", "numeric", "numeric",
                                                                      "numeric", "numeric", "numeric", "numeric",
                                                                      "numeric", "integer", "integer")
  )
  distances <- get(load( paste0(path_results, "attribution/", version_name,nearby.ver, "distances-pairs.RData")))
  
  validation = c()
  homo.code = c()
  homo = NA
  near.name = c()
  list.sta = c()
  possible = c()
  list.detec <-c(as.Date("00/00/00", "%y/%m/%d"))
  np.bef = c()
  np.aft = c()
  data.test <- list()
  for (j in name.station){
    ref.station = j # name staion consider
    ref.seg = valid.ref[which(valid.ref$name == j),]
    series.ref <- read.series(path_series = path_series_ref, station = ref.station, na.rm = 1, add.full = 0)
    
    nears.list = as.character(last.list[which(last.list$ref.sta == j),])
    
    if(nrow(ref.seg) >0 ){
      for (l in 1:nrow(ref.seg)){
        breakpoint = ref.seg$detected[l]
        other.brp = ref.seg$detected[-l]
        outlier.class =  ref.seg$noise[l]
        
        
        ## filter stations have no nearby 
        if (length(na.omit(nears.list)) == 1){
          homo = 0
          name.nearby = "nothave"
        } else if(length(na.omit(nears.list)) >= 2){ #consider the station have nearby <150km !!!! tiep tu day
          for (k in 2:length(na.omit(nears.list))){
            name.nearby = nears.list[k]
            near.seg =valid.near[which(valid.near$name==nears.list[k]),]
            series.near <- read.series(path_series = path_series_near, station = nears.list[k], na.rm = 1, add.full = 0)
            # correct the vertical dfference, data from here will be used in the weighted test 
            if(correct.ver == 1){
              n.main = which(distances$main == ref.station & distances$nearby ==name.nearby)
              ver.dist = distances$ver.dist[n.main]   #vertical disstance here is main-nearby
              h.dist = distances$distances[n.main] 
              # series.near$GPS <- series.near$GPS*(1+ver.dist*0.0004) ## old correction
              # series.near$ERAI <- series.near$ERAI*(1+ver.dist*0.0004)
              series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
              series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
            }else{
              series.near = series.near
            }
            
            # !!! add the point for ERA5
            data.joint = full_join(series.ref,series.near, by = "date")
            data.joint = data.joint[which(is.na(data.joint$signal.x) == FALSE | is.na(data.joint$signal.y) == FALSE),]
            list.nb1 = which(data.joint$date %in% breakpoint)
            list.nb2 = which(data.joint$date %in% near.seg$detected)
            # If not, we perform the test 
            
            both = inner_join(series.ref,series.near, by = "date")  # checked
            both$signal = both$GPS.y-both$GPS.x
            
            # Add the info about crenel---changed the min, max because the order is changed, print the begin and end date of crenel
            # check.cluster = check_cluster(list.cluster = list.cluster, station = ref.station, list.ind.cluster = list.ind.cluster,
            #                               min.d = breakpoint, max.d = breakpoint, breakpoint.c = breakpoint)
            # 
            # cls.begin = check.cluster$cls.begin
            # cls.end = check.cluster$cls.end
            
            #list.all.seg = c(other.brp)
            list.all.seg = c(other.brp, near.seg$detected)
            
            if (length(list.nb2) == 0){
              diff.brp = 100000
            }else{
              diff.brp = abs(c(sapply( list.nb1, "-", list.nb2 )))
            }
            if (length(diff.brp[which(diff.brp<= 10)]) >0 ){
              
              close.brealpoint = data.joint$date[list.nb2[which.min(diff.brp)]]
              list.2break = c(close.brealpoint, breakpoint)
              
              list.all.seg1 = list.all.seg[-which(list.all.seg== close.brealpoint)]
              begin.day = max(na.omit(c(breakpoint - duration +1, mymax(list.all.seg1[which(list.all.seg1 < breakpoint)]))))
              end.day = min(na.omit(c(breakpoint + duration, mymin(list.all.seg1[which(list.all.seg1 > breakpoint)]))))
              
              
              # before1 = both[which(both$date >= begin.day & both$date <= cls.begin & both$date <= min(list.2break)),]
              # after1 = both[which(both$date >  max(list.2break) & both$date >= cls.end & both$date <= end.day),]
              before1 = both[which(both$date >= begin.day & both$date <= min(list.2break)),]
              after1 = both[which(both$date >  max(list.2break) & both$date <= end.day),]
              listday1 = data.frame( date = before1$date, month = format(before1$date, "%m"), day = format(before1$date, "%d"))
              listday2 = data.frame( date = after1$date, month = format(after1$date, "%m"), day = format(after1$date, "%d"))
              list.day.joint = inner_join(listday1, listday2, by = c("month", "day"))
              
              before = before1[which(before1$date %in% list.day.joint$date.x == TRUE),]
              after = after1[which(after1$date %in% list.day.joint$date.y == TRUE),]
              
              required.nb.point = duration*percentage.required
              con1 = nrow(before) > required.nb.point  & nrow(after)>required.nb.point
              np.be = nrow(na.omit(before))
              np.af = nrow(na.omit(after))
              if (mean(con1)==1){
                possib = 2
                homo = 2
              }else{
                possib = 0
                homo=0
              }
            } else {
              
              if(length(list.all.seg) ==0 ){
                begin.day = breakpoint - duration +1
                end.day = breakpoint + duration
              } else if(length(list.all.seg) > 0  ){
                begin.day = max(na.omit(c(breakpoint - duration +1, mymax(list.all.seg[which(list.all.seg < breakpoint)]))))
                end.day = min(na.omit(c(breakpoint + duration, mymin(list.all.seg[which(list.all.seg > breakpoint)]))))
              }
              
              # superpose the 2 segment
              # before1 = both[which(both$date >= begin.day & both$date <= cls.begin & both$date <= breakpoint),]
              # after1 = both[which(both$date > breakpoint & both$date >= cls.end & both$date <= end.day),]
              # # superpose the 2 segment 
              before1 = both[which(both$date >= begin.day  & both$date <= breakpoint),]
              after1 = both[which(both$date > breakpoint & both$date <= end.day),]
              # 
              listday1 = data.frame( date = before1$date, month = format(before1$date, "%m"), day = format(before1$date, "%d"))
              listday2 = data.frame( date = after1$date, month = format(after1$date, "%m"), day = format(after1$date, "%d"))
              list.day.joint = inner_join(listday1, listday2, by = c("month", "day"))
              
              before = before1[which(before1$date %in% list.day.joint$date.x == TRUE),]
              after = after1[which(after1$date %in% list.day.joint$date.y == TRUE),]
              
              required.nb.point = duration*percentage.required
              if (outlier.class > 0){
                required.nb.point = 10
              }       
              
              con1 = nrow(before) > required.nb.point  & nrow(after)>required.nb.point
              np.be = nrow(na.omit(before))
              np.af = nrow(na.omit(after))
              if (mean(con1)==1){
                possib = 1
                homo=1 
              }else{
                possib = 0
                homo=0
              }
            }
            validation <- c( validation, ref.seg$valid[l])
            possible = c(possible,possib)
            homo.code = c(homo.code, homo)
            near.name = c(near.name, name.nearby)
            list.sta = c(list.sta, ref.station)
            list.detec = c(list.detec, breakpoint)
            np.bef = c(np.bef,  np.be)
            np.aft = c(np.aft,  np.af)
            if (possib > 0){
              data.test[[paste0(ref.station,".",as.character( breakpoint), ".", name.nearby)]] <- list(before = before, after = after)
            }
          }
        }
      }
    }
    
  }   
  save(data.test, file = paste0(path_results, "attribution/data_test_",version.name,nb_test,nearby.ver,limit.type,screen_value,"multiple.RData" ))
  out= data.frame(list.sta,list.detec[-1],validation,homo.code,near.name,possible, np.bef, np.aft)
  colnames(out) <- c("ref.station","detected","validated","homo", "nearby.station","possible", "np.bef", "np.aft" )
  out <- out[which( out$possible>0),]
  rownames(out) <- NULL
  return(out)
}
# find the homogeneous for all breakpoints after screening 
homogeneous.nearby.segment.multiple = function( name.station, valid.file.ref, valid.file.near, last.list, path_series_ref,
                                                path_series_near, duration, percentage.required, correct.ver, nearby.ver,
                                                version.name, nb_test, limit.type, screen_value  ){
  # Read the list of breakpoints in the main and the nearby
  valid.ref = read.table(valid.file.ref, header = TRUE, sep = "\t",  colClasses = c("character", "Date","Date", "character","numeric",
                                                                                    "integer", "integer")
  )
  valid.near = read.table(valid.file.near, header = TRUE, sep = "\t",  colClasses = c("character", "Date","Date", "character","numeric",
                                                                                      "integer", "integer")
  )
  # List of the crenel
  list.cluster <-  read.table( file = paste0(path_results,"screening/lmin1.",criterion,".date_mean_weight_np.",version.name,nb_test.ref,".txt"),
                               sep = "\t", header = T,  colClasses = c("character", "Date","Date", "numeric","numeric",
                                                                       "integer", "integer")
  )
  list.ind.cluster = read.table(file = paste0(path_results,"screening/lmin1.",criterion,".screening.", version.name,nb_test.ref,".txt"),
                                sep="\t",  header = T, colClasses = c("character", "Date","integer", "integer",
                                                                      "integer", "integer", "numeric", "numeric",
                                                                      "numeric", "numeric", "numeric", "numeric",
                                                                      "numeric", "integer", "integer")
  )
  
  distances <- get(load( paste0(path_results, "attribution/", version_name,nearby.ver, "distances-pairs.RData")))
  
  validation = c()
  homo.code = c()
  homo = NA
  near.name = c()
  list.sta = c()
  possible = c()
  list.detec <-c(as.Date("00/00/00", "%y/%m/%d"))
  np.bef = c()
  np.aft = c()
  data.test <- list()
  for (j in name.station){
    ref.station = j # name staion consider
    ref.seg = valid.ref[which(valid.ref$name == j),]
    series.ref <- read.series(path_series = path_series_ref, station = ref.station, na.rm = 1, add.full = 0)
    
    nears.list = as.character(last.list[which(last.list$ref.sta == j),])
    
    if(nrow(ref.seg) >0 ){
      for (l in 1:nrow(ref.seg)){
        breakpoint = ref.seg$detected[l]
        other.brp = ref.seg$detected[-l]
        
        ## filter stations have no nearby
        if (length(na.omit(nears.list)) == 1){
          homo = 0
          name.nearby = "nothave"
        } else if(length(na.omit(nears.list)) >= 2){ #consider the station have nearby <150km !!!! tiep tu day
          for (k in 2:length(na.omit(nears.list))){
            name.nearby = nears.list[k]
            near.seg =valid.near[which(valid.near$name==nears.list[k]),]
            series.near <- read.series(path_series = path_series_near, station = nears.list[k], na.rm = 1, add.full = 0)
            # correct the vertical dfference, data from here will be used in the weighted test
            if(correct.ver == 1){
              n.main = which(distances$main == ref.station & distances$nearby ==name.nearby)
              ver.dist = distances$ver.dist[n.main]   #vertical disstance here is main-nearby
              h.dist = distances$distances[n.main]
              # series.near$GPS <- series.near$GPS*(1+ver.dist*0.0004) ## old correction
              # series.near$ERAI <- series.near$ERAI*(1+ver.dist*0.0004)
              series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
              series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
            }else{
              series.near = series.near
            }
            
            # !!! add the point for ERA5
            data.joint = full_join(series.ref,series.near, by = "date")
            data.joint = data.joint[which(is.na(data.joint$signal.x) == FALSE | is.na(data.joint$signal.y) == FALSE),]
            list.nb1 = which(data.joint$date %in% breakpoint)
            list.nb2 = which(data.joint$date %in% near.seg$detected)
            # If not, we perform the test
            
            both = inner_join(series.ref,series.near, by = "date")  # checked
            both$signal = both$GPS.y-both$GPS.x
            
            # Add the info about crenel---changed the min, max because the order is changed, print the begin and end date of crenel
            check.cluster = check_cluster(list.cluster = list.cluster, station = ref.station, list.ind.cluster = list.ind.cluster,
                                          min.d = breakpoint, max.d = breakpoint, breakpoint.c = breakpoint)
            
            cls.begin = check.cluster$cls.begin
            cls.end = check.cluster$cls.end
            
            #list.all.seg = c(other.brp)
            list.all.seg = c(other.brp, near.seg$detected)
            
            if (length(list.nb2) == 0){
              diff.brp = 100000
            }else{
              diff.brp = abs(c(sapply( list.nb1, "-", list.nb2 )))
            }
            if (length(diff.brp[which(diff.brp<= 10)]) >0 ){
              
              close.brealpoint = data.joint$date[list.nb2[which.min(diff.brp)]]
              list.2break = c(close.brealpoint, breakpoint)
              
              list.all.seg1 = list.all.seg[-which(list.all.seg== close.brealpoint)]
              begin.day = max(na.omit(c(breakpoint - duration +1, mymax(list.all.seg1[which(list.all.seg1 < breakpoint)]))))
              end.day = min(na.omit(c(breakpoint + duration, mymin(list.all.seg1[which(list.all.seg1 > breakpoint)]))))
              
              
              before1 = both[which(both$date >= begin.day & both$date <= cls.begin & both$date <= min(list.2break)),]
              after1 = both[which(both$date >  max(list.2break) & both$date >= cls.end & both$date <= end.day),]
              listday1 = data.frame( date = before1$date, month = format(before1$date, "%m"), day = format(before1$date, "%d"))
              listday2 = data.frame( date = after1$date, month = format(after1$date, "%m"), day = format(after1$date, "%d"))
              list.day.joint = inner_join(listday1, listday2, by = c("month", "day"))
              
              before = before1[which(before1$date %in% list.day.joint$date.x == TRUE),]
              after = after1[which(after1$date %in% list.day.joint$date.y == TRUE),]
              
              required.nb.point = duration*percentage.required
              con1 = nrow(before) > required.nb.point  & nrow(after)>required.nb.point
              np.be = nrow(na.omit(before))
              np.af = nrow(na.omit(after))
              if (mean(con1)==1){
                possib = 2
                homo = 2
              }else{
                possib = 0
                homo=0
              }
            } else {
              
              if(length(list.all.seg) ==0 ){
                begin.day = breakpoint - duration +1
                end.day = breakpoint + duration
              } else if(length(list.all.seg) > 0  ){
                begin.day = max(na.omit(c(breakpoint - duration +1, mymax(list.all.seg[which(list.all.seg < breakpoint)]))))
                end.day = min(na.omit(c(breakpoint + duration, mymin(list.all.seg[which(list.all.seg > breakpoint)]))))
              }
              
              # superpose the 2 segment
              before1 = both[which(both$date >= begin.day & both$date <= cls.begin & both$date <= breakpoint),]
              after1 = both[which(both$date > breakpoint & both$date >= cls.end & both$date <= end.day),]
              
              listday1 = data.frame( date = before1$date, month = format(before1$date, "%m"), day = format(before1$date, "%d"))
              listday2 = data.frame( date = after1$date, month = format(after1$date, "%m"), day = format(after1$date, "%d"))
              list.day.joint = inner_join(listday1, listday2, by = c("month", "day"))
              
              before = before1[which(before1$date %in% list.day.joint$date.x == TRUE),]
              after = after1[which(after1$date %in% list.day.joint$date.y == TRUE),]
              
              required.nb.point = duration*percentage.required
              con1 = nrow(before) > required.nb.point  & nrow(after)>required.nb.point
              np.be = nrow(na.omit(before))
              np.af = nrow(na.omit(after))
              if (mean(con1)==1){
                possib = 1
                homo=1
              }else{
                possib = 0
                homo=0
              }
            }
            validation <- c( validation, ref.seg$valid[l])
            possible = c(possible,possib)
            homo.code = c(homo.code, homo)
            near.name = c(near.name, name.nearby)
            list.sta = c(list.sta, ref.station)
            list.detec = c(list.detec, breakpoint)
            np.bef = c(np.bef,  np.be)
            np.aft = c(np.aft,  np.af)
            data.test[[paste0(ref.station,".",as.character( breakpoint), ".", name.nearby)]] <- list(before = before, after = after)
          }
        }
      }
    }
    
  }
  save(data.test, file = paste0(path_results, "attribution/data_test_",version.name,nb_test,nearby.ver,limit.type,screen_value,"multiple.RData" ))
  out= data.frame(list.sta,list.detec[-1],validation,homo.code,near.name,possible, np.bef, np.aft)
  colnames(out) <- c("ref.station","detected","validated","homo", "nearby.station","possible", "np.bef", "np.aft" )
  out <- out[which( out$possible>0),]
  rownames(out) <- NULL
  return(out)
}

RobEstiMonthlyVariance.diff.S <- function(Y, name.var, alpha){
  Kmonth <- length(unique(Y$month))
  Y$month <- factor(Y$month, ordered = TRUE)
  z <- stats::aggregate(Y[name.var], by = Y[c("month", "year")], FUN = function(x) { re = diff(na.omit(x)); return(re)})
  sigma.est1 <- sapply(1:Kmonth,function(i) {
    e <- subset(z,z$month==levels(z$month)[i])
    ee <- unlist(e[name.var])
    if(length(ee) == 0){
      NA
    }else{
      robustbase::scaleTau2(ee)/sqrt(2-2*alpha) # change for AR(1), if white noise alspha = 0
    }
  })
  sigma.est <- rep(NA,12)
  sigma.est[as.numeric(levels(z$month))] <-  sigma.est1
  return(sigma.est)
}

# t test  -----------------------------------------------------------------
# Result of individual test, different than previous test, check correlation can be from 2 to 4 with 2 is the concorrected variance

significance.test.weighted.indi = function(homo.nearby, path_series_ref,  path_series_near, duration, sig.level, name.test,
                                           screening.range, screening ,version.name,nb_test, correct.ver, nearby.ver, limit.type,
                                           screen_value, check_correlation ){ 
  distances <- get(load( paste0(path_results, "attribution/", version_name, nearby.ver, "distances-pairs.RData")))
  data_test <-  get(load(paste0(path_results,"attribution/data_test_",version.name,nb_test,nearby.ver,limit.type,screen_value,"multiple.RData")))
  autocorr.order <- get(load( file = paste0(path_results, "attribution/subsampling.support.order.", 
                                            version.name = version.name,
                                            nb_test = nb_test,
                                            nearby.ver = nearby.ver,
                                            limit.type = "superposed",
                                            screen_value = screen_value, ".RData")))
  list.param = read.table(file =  paste0(path_results,"attribution/", "param.all_var_est_",  "_arima_alpha",alpha = significant.level,nb_test = nb_test.ref,"-",screen_value,tolerance_noise,".txt"), 
                          header = TRUE)
  nb.t.test <- nrow(homo.nearby)
  monthly.variance = list()
  
  var1 = diff.var(c(name.test))[1]
  var2 = diff.var(c(name.test))[2]
  rho <- rep(NA,nb.t.test)
  for (p in 1:nb.t.test){
    detection = homo.nearby$detected[p]
    
    station.ref = homo.nearby$ref.station[p]
    station.near = homo.nearby$nearby.station[p]
    # read two time series
    series.ref <- read.series(path_series = path_series_ref, station = station.ref, na.rm = 1, add.full = 0)
    series.near <- read.series(path_series = path_series_near, station =  station.near, na.rm = 1, add.full = 0)
    # vertical disstance here is main-nearby
    if(correct.ver == 1){
      n.main = which(distances$main == station.ref & distances$nearby ==station.near)
      ver.dist = distances$ver.dist[n.main]   
      series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
      series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
    }else{
      series.near = series.near
    }
    # Estimate the monthly variance used to test 
    both = inner_join(series.ref,series.near, by = "date")  # checked
    both[,name.test] = both[,c(var1)] - both[,c(var2)]
    
    both$year <- both$year.x
    both$month <- both$month.x
    
    both.pre = na.omit(both)
    both <- tidyr::complete(both, date = seq(min(both$date), max(both$date), by = "day"))
    both$year = droplevels( as.factor(format(both$date,format='%Y')))
    both$month = droplevels(as.factor(format(both$date,format='%m')))
    
    # sigma.est.month <- RobEstiMonthlyVariance.new(Y = (na.omit(both)), name.var = name.test)^2
    print(p)
    param = list.param[p,name.test]
    sigma.est.month <- RobEstiMonthlyVariance.diff.S(Y = (na.omit(both)), name.var = name.test, alpha = param)^2
    monthly.variance[[paste0(station.ref,"-",station.near )]] <- c(sigma = sigma.est.month)
    
    min = detection - duration
    max = detection + duration
    # Read data to test from the homogeneous function
    before =  data_test[[paste0( station.ref,".",as.character(detection), ".",  station.near)]]$before
    after = data_test[[paste0( station.ref,".",as.character(detection), ".",  station.near)]]$after
    
    before[,name.test] = before[,c(var1)] - before[,c(var2)] # compute the difference
    after[,name.test] = after[,c(var1)] - after[,c(var2)]
    # sort data due to autocorrelation
    # if (check_correlation == 1){
    #   order.auto = autocorr.order[which(autocorr.order$ref.station == station.ref &
    #                                       autocorr.order$detected == detection &
    #                                       autocorr.order$nearby.station ==  station.near),]
    #   max.order = order.auto$order
    #   list.sorted <- seq(1, nrow(before), max.order)
    #   before <- before[list.sorted,]
    #   after <- after[list.sorted,]
    # }
    if(check_correlation == 3){
      
      K.month = length(unique(before$month.x)) 
      list.month = sort(as.numeric(unique(before$month.x)), decreasing = FALSE) # the case of lacking months 
      
      # Estimate the order and ar coefficient
      list.order.ar <- lapply(c(1:K.month), function(x) {
        af = after[which(as.numeric(after$month.x) ==  list.month[x]),]
        be = before[which(as.numeric(before$month.x) ==  list.month[x]),]
        if (nrow(af) >2){
          rho_a = ar(af[,name.test], order.max = 1)
          rho_b = ar(be[,name.test], order.max = 1)
          y = c(rho_a$order, rho_b$order, rho_a$ar, rho_b$ar)
        } else{
          y = c(0,0)
        }
        return(y)
      })
      
      # Estimate the effective monthly variance 
      for (e.m in c(1:K.month)) {
        list.ar <- list.order.ar[[e.m]]
        af = after[which(as.numeric(after$month.x) == list.month[e.m]),]
        if (sum(list.ar[1:2]) == 1){
          mean.ar <- abs(list.ar[3])
        } else if(sum(list.ar[1:2]) == 2){
          mean.ar <- abs(mean( list.ar[3:4]))
        } else { mean.ar  <- 0 }
        n.eff <- (nrow(af))/(1 + 2*(1 - 1/ (nrow(af))) * mean.ar )
        sigma.est.month[list.month[e.m]] <- sigma.est.month[list.month[e.m]]*(nrow(af))/n.eff
      }
    }
    else if(check_correlation == 4){
      
      # Estimate the order and ar coefficient
      Y = rbind(before, after)
      rhoMaG = ar.ma.genton.2side(as.vector(unlist(Y[name.test])))
      
      # Estimate the effective monthly variance 
      n0 = nrow(Y)
      Neff.est = ESS(n = n0, rho = rhoMaG)
      var.eff2 = sigma.est.month*n0/Neff.est
      sigma.est.month = var.eff2
      rho[p] <- rhoMaG
    }else if(check_correlation == 5){
      
      list.ratio = read.table(file =  paste0(path_results,"attribution/", "ratio_correct_", "_arima_alpha",alpha = significant.level,nb_test = nb_test.ref,"-",screen_value,tolerance_noise,".txt"), 
                              header = TRUE)
      ratio_cor = list.ratio[p,name.test]
      var.eff2 = sigma.est.month*ratio_cor
      sigma.est.month = var.eff2
      rho[p] <- NA
    }
    
    # add the estimated the variance into the data and compute the weighted mean  
    
    before$month.var <- sapply(1:nrow(before), function(x) sigma.est.month[as.numeric(before$month.x[x])])
    before <- tidyr::complete(before, date = seq(min(before$date), max(before$date), by = "day"))
    before$month.var[which(is.na(before[,name.test]) ==T | before[,name.test] ==0)] <- NA
    weight.b <- (sum(1/before$month.var,na.rm = T))
    wei.mean.b<- (sum(before[,name.test]/before$month.var,na.rm = T))/weight.b
    
    after$month.var <- sapply(1:nrow(after), function(x) sigma.est.month[as.numeric(after$month.x[x])])
    after <- tidyr::complete(after, date = seq(min(after$date), max(after$date), by = "day"))
    after$month.var[which(is.na(after[,name.test]) ==T | after[,name.test] ==0)] <- NA
    weight.a <- (sum(1/after$month.var,na.rm = T))
    wei.mean.a<- (sum(after[,name.test]/after$month.var,na.rm = T))/weight.a
    
    homo.nearby$nb.before[p] = nrow(na.omit(before))
    homo.nearby$nb.after[p] = nrow(na.omit(after))
    
    # compute the statistics 
    t1 = before[[name.test]]
    t2 = after[[name.test]]
    t.stat <- (wei.mean.a - wei.mean.b)/(sqrt(1/weight.a + 1/weight.b))
    month.d = as.numeric(format(detection,"%m"))
    p.stat = round(pnorm(-abs(t.stat), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 3)
    homo.nearby$t[p] = round(t.stat, digits = 3)
    homo.nearby$p[p] = p.stat
    homo.nearby$sig[p] = check.sig(p.stat, sig.level)
    homo.nearby$wei.mean.be[p] <- round(wei.mean.b, digits = 3)
    homo.nearby$wei.mean.af[p] <- round(wei.mean.a, digits = 3)
    homo.nearby$offset[p] <- round((wei.mean.b - wei.mean.a), digits = 3) 
    homo.nearby$std.be[p] <- round(sd(t1, na.rm = T), digits = 3)
    homo.nearby$std.af[p] <- round(sd(t2, na.rm = T), digits = 3)
    homo.nearby$rel.std.be[p] <- round(sd(t1, na.rm = T)/(mean(before$GPS.x, na.rm = T)), digits = 3)
    homo.nearby$rel.std.af[p] <- round(sd(t2, na.rm = T)/(mean(after$GPS.x, na.rm = T)), digits = 3)
    homo.nearby$cor[p] <- round(cor(na.omit(both[, c(var1)]), na.omit(both[, c(var2)])), digits = 3)
    homo.nearby$weight.be[p] <- round(weight.b, digits = 3)
    homo.nearby$weight.af[p] <- round(weight.a, digits = 3)
    homo.nearby$corre1[p] <- round(cor(before[[var1]], before[[var2]], use = "complete.obs"), digits = 3)
    homo.nearby$corre2[p] <- round(cor(after[[var1]], after[[var2]], use = "complete.obs"), digits = 3)
    homo.nearby$SNR[p] <-  (abs(wei.mean.b - wei.mean.a))/(sqrt(1/weight.b + 1/weight.a))
    
  }
  test = homo.nearby[, -c(6:8)]
  # print the estimated monthlt variance
  write.table(test, 
              file = paste0(path_results, "attribution/", name.test,".", duration, "days.level", sig.level,
                            ".scr", screen_value, "no.test.", nb_test, ".nearby.", nearby.ver, check_correlation, ".txt"),
              sep = "\t", quote = FALSE)
  save(monthly.variance, file = paste0(path_results, "attribution/monthly-var-weigted.",name.test,".", 
                                       duration, "days.level", sig.level, ".scr", screen_value, "no.test.", nb_test,
                                       ".nearby.", nearby.ver, check_correlation, ".RData"))
  save(rho, file = paste0(path_results, "attribution/rho.",name.test,".", 
                          duration, "days.level", sig.level, ".scr", screen_value, "no.test.", nb_test,
                          ".nearby.", nearby.ver, check_correlation, ".RData"))
  
  homo.last = homo.nearby
  
  return(homo.last)
}

# plot results ------------------------------------------------------------
plot_attribution_fc_multi <- function(data1, data2, var.name, nb_test, station.ref, station.near, detection, nearby.ver,
                                      list.brp.ref, list.brp.near, type.meta.near, type.meta.ref, text1, text2, duration,
                                      ylab){
  begin.date = detection - duration
  end.date = detection + duration
  
  colors <- c(station.ref = "blue", station.near = "red", "diff" = "black","var" = "purple")
  
  names(colors) <- c( station.ref,  station.near, "diff", "var")
  t <- data2$date
  p <- ggplot(data1, aes(x = date)) + # need to change when plot for different pair, here I only consider the difference series p2 
    geom_line(aes(y = GPS.x , color = station.ref), size = 0.5) +
    geom_line(aes(y = GPS.y, color = station.near), size = 0.5)
  p<-p+ geom_vline(xintercept = list.brp.ref, size = 0.5, lty = 2, colour = "blue", alpha = 0.6)
  p<-p+ geom_vline(xintercept = detection, size = 1, lty = 1, colour = "blue", alpha = 1)
  if (length(list.brp.near) >0 ){
    p<-p + geom_vline(xintercept = list.brp.near, size = 0.5, lty = 2, colour = "red", alpha = 0.6)
  }
  p<- p + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                       breaks = function(x) seq.Date(from = begin.date,  to = end.date, by = "2 months"),
                       date_labels = "%Y-%m")
  p <- p+ scale_color_manual(values = colors)+
    theme_bw()+
    labs(x = "date",
         y = "IWV (kg/m2)",
         color = "Legend")+
    labs(subtitle = text1, size = 0.1) +
    theme(legend.text=element_text(size=2),
          axis.text = element_text(size=2),
          plot.subtitle = element_text(size = 0.1)) +
    theme_bw()
  # second figure 
  p1 <- ggplot(data2 , aes_string(x = "date")) +
    geom_path(data=data2,aes_string(y = var.name ), size = 0.5, alpha=0.5)+
    geom_vline(xintercept = detection, size = 1, lty = 2, colour = "black", alpha = 1)
  # geom_path(data=subset(data2, !is.na(data2[,var.name])),aes_string(y = var.name ), size = 0.5, alpha=0.5) 
  # geom_line(aes(y = mean), size = 1, color = "red")
  
  min.y <- ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]
  
  thres = max(data2$var)
  # p1 <- p1+ geom_line(data=subset(data2, !is.na(data2[,var.name])), aes(y = var- abs(min.y), color = "var"), size = 0.5, alpha=0.5)
  # p1 <- p1 + geom_hline(yintercept = min.y, color = "black", size = 0.3)
  #plot the known changes in the nearby --------
  type.known <- c("list.R","list.D","list.A","list.D","list.L", "list.E","list.U" )
  gps.meta <- c( "gps.gps", "gps.era1", "gps.era")
  if(var.name %in% gps.meta){
    # for (k in 1:7) {
    #   type.k = type.known[k]
    #   if(length(type.meta.near[[type.k]]) >0){
    #     list.brp.screened.y <- rep(min.y, length(type.meta.near[[type.k]]))
    #     p1 <- p1+ annotate("point", x = t[type.meta.near[[type.k]]], y = list.brp.screened.y , colour = "red", size = 1.5, alpha=0.5,shape = k-1)
    #   }
    # }
    # plot the known changes in the main station
    for (k in 1:4) {
      type.k = type.known[k]
      if(length(type.meta.ref[[type.k]]) >0){
        list.brp.screened.y <- rep(min.y, length(type.meta.ref[[type.k]]))
        p1 <- p1+ annotate("point", x = t[type.meta.ref[[type.k]]], y = list.brp.screened.y , colour = "blue", size = 1.5, alpha=0.5,shape = k+4)
      }
    }
  }
  
  
  min.x =  begin.date
  p1 <- p1+ scale_x_date(limits = as.Date(c(begin.date, end.date)),
                         breaks = function(x) seq.Date(from = begin.date,  to = end.date, by = "2 months"),
                         date_labels = "%Y-%m" )
  # p1 <- p1 + annotate("text", x = min.x+150, y = 1.8, label = text2, size = 6)
  p1 <- p1 + scale_color_manual(values = colors)+
    theme_bw()+  # it can overwrite some setting of theeme so put it at first 
    labs(x = "date",
         y = ylab,
         color = "Legend")+
    # scale_y_continuous(limits = c(-3, 2), breaks = seq(-3, 2, 1)) +
    # labs(subtitle = text2) +
    theme(legend.text=element_text(size=2),
          axis.text = element_text(size=20),
          axis.title.y =element_text(size=30),
          axis.title.x =element_text(size=18,  face = "bold")) 
  gA <- ggplotGrob(p)
  gB <- ggplotGrob(p1)
  gA$widths <- gB$widths
  jpeg(paste0(path_results,"figure/",nb_test, "/", nearby.ver, "/" , station.ref, detection, station.near,".jpeg"),width = 3000, height = 1800,res = 300) # change name
  # print(grid.arrange(gA, gB, nrow = 2))
  print(p1)
  dev.off()
  return(list(p, p1))
}









