# old code
check_cluster <- function(list.cluster, station, list.ind.cluster, breakpoint.c){
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

# new code use the meta_compare file, read only the cluster 
meeta.compare =  get(load(file = paste0(path_results,"validation/",nb_test.ref = 34,"-",criterion,"metacompa",screen.value="",".RData")))

## continue!!!!!
