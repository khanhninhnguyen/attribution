# study the different veriecal correction on the 4 differences: GPS-ERA', GPS'-ERA, GPS-GPS', ERA-ERA'

# draft -------------------------------------------------------------------
window.thres = 10
distances <- get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData")))
new.ver = get(load(file = paste0(path_results,"attribution/data.all_", window.thres,"years_", nearby_ver,"screened.RData")))
old.ver0 = get(load(file = paste0(path_results,"attribution/OLD_VERTICAL/six_diff_series_",window.thres,"year_", nearby_ver,".RData")))
ver0 = get(load(file = paste0(path_results,"attribution/six_diff_series_",window.thres,"year_", nearby_ver,"0.RData")))

old.ver = old.ver0[names(new.ver)]

h = sapply(c(1:length(old.ver)), function(x){
  station.ref = substr(names(old.ver)[x] ,start = 1, stop = 4)
  station.seg = substr(names(old.ver)[x] ,start = 17, stop = 20)
  ind = which(distances$main == station.ref & distances$nearby == station.seg)
  distances$ver.dist[ind]
})
seg = c(1:3650)
gps.gps1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$gps.gps[seg], na.rm = TRUE))
gps.gps2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$gps.gps[seg], na.rm = TRUE))
gps.gps0 = sapply(c(1:length(old.ver)), function(x)  mean(ver0[[x]]$gps.gps[seg], na.rm = TRUE))

gps.gps1[which(h<0)] = -gps.gps1[which(h<0)] 
gps.gps2[which(h<0)] = -gps.gps2[which(h<0)] 
gps.gps0[which(h<0)] = -gps.gps0[which(h<0)] 

plot(abs(h), gps.gps2)
plot(gps.gps1, gps.gps2)
plot(gps.gps0, gps.gps2)


gps.era1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$gps.era1[seg], na.rm = TRUE))
gps.era2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$gps.era1[seg], na.rm = TRUE))
gps.era0 = sapply(c(1:length(old.ver)), function(x)  mean(ver0[[x]]$gps.era1[seg], na.rm = TRUE))

gps1.era1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$gps1.era[seg], na.rm = TRUE))
gps1.era2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$gps1.era[seg], na.rm = TRUE))
gps1.era0 = sapply(c(1:length(old.ver)), function(x)  mean(ver0[[x]]$gps1.era[seg], na.rm = TRUE))

era.era1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$era.era[seg], na.rm = TRUE))
era.era2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$era.era[seg], na.rm = TRUE))
era.era0 = sapply(c(1:length(old.ver)), function(x)  mean(ver0[[x]]$era.era[seg], na.rm = TRUE))

dat = data.frame(
  GPS.GPS1 = c(mean(abs(gps.gps1), na.rm = TRUE), mean(abs(gps.gps2), na.rm = TRUE),  mean(abs(gps.gps0), na.rm = TRUE)),
  GPS.ERA1 = c(mean(abs(gps.era1), na.rm = TRUE), mean(abs(gps.era2), na.rm = TRUE),  mean(abs(gps.era0), na.rm = TRUE)),
  GPS1.ERA = c(mean(abs(gps1.era1), na.rm = TRUE), mean(abs(gps1.era2), na.rm = TRUE),  mean(abs(gps1.era0), na.rm = TRUE)),
  ERA.ERA1 = c(mean(abs(era.era1), na.rm = TRUE), mean(abs(era.era2), na.rm = TRUE),  mean(abs(era.era0), na.rm = TRUE))
)
dat$ver = as.factor(c( "Ver1", "Ver2", "Raw"))
a = reshape2::melt(dat, id = "ver")

ggplot(a, aes(x = variable, y = value, col = ver))+geom_point()+theme_bw()+ylab("Mean of the absolute mean differences")

dat = data.frame(
  GPS.GPS1 = c(sd((gps.gps1), na.rm = TRUE), sd((gps.gps2), na.rm = TRUE),  sd((gps.gps0), na.rm = TRUE)),
  GPS.ERA1 = c(sd((gps.era1), na.rm = TRUE), sd((gps.era2), na.rm = TRUE),  sd((gps.era0), na.rm = TRUE)),
  GPS1.ERA = c(sd((gps1.era1), na.rm = TRUE), sd((gps1.era2), na.rm = TRUE),  sd((gps1.era0), na.rm = TRUE)),
  ERA.ERA1 = c(sd((era.era1), na.rm = TRUE), sd((era.era2), na.rm = TRUE),  sd((era.era0), na.rm = TRUE))
)
dat$ver = as.factor(c( "Ver1", "Ver2", "Raw"))
a = reshape2::melt(dat, id = "ver")

ggplot(a, aes(x = variable, y = value, col = ver))+geom_point()+theme_bw()+ylab("Std of the mean differences")


dat = data.frame(raw = era.era0, ver1 = era.era1, ver2 = era.era2)
a = reshape2::melt(dat, id = "raw")

ggplot(a, aes(x = raw, y = value, col = variable))+geom_point()+theme_bw()+ylab("corrected")+
  geom_abline(slope = 1)

rmsd <- function(x,y){
  sd = (x-y)^2
  sqrt(mean(sd, na.rm = TRUE))
}

rmsd(era.era0, era.era1)
rmsd(era.era0, era.era2)





# correct  ----------------------------------------------------------------
# compare in separate segment based on the restrict length: four series are restricted with breakpoints by both main and nearby 
# compute the raw data without correction
window.thres = 10
data.res = get(load(file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",window.thres,"year_", nearby_ver,".RData")))
old.ver0 = get(load(file = paste0(path_results,"attribution/OLD_VERTICAL/six_diff_series_",window.thres,"year_", nearby_ver,".RData")))

list.break = data.frame(ref = substr(names(data.res), start = 1, stop = 4), 
                        brp = substr(names(data.res), start = 6, stop = 15),
                        nb = substr(names(data.res), start = 17, stop = 20))
list.break[] <- lapply(list.break, as.character)
list.break$brp = as.Date(list.break$brp , format = "%Y-%m-%d")

merge_series <- function(path1, path2, station1, station2, list.day){
  series1 = read.series(path_series = path1, station = station1, na.rm = 0, add.full = 0)
  series2 = read.series(path_series = path2, station = station2, na.rm = 0, add.full = 0)
  
  both = inner_join(series1,series2, by = "date") 
  both <- tidyr::complete(both, date = seq(min(both$date), max(both$date), by = "day"))
  both = both[which(both$date %in% list.day),]
  out <- list()
  for (m in 1:6) {
    name.test = list.test[m]
    var1 = diff.var(c(name.test))[1]
    var2 = diff.var(c(name.test))[2]
    out[name.test] <- both[,c(var1)] - both[,c(var2)]
  }
  return(out)
}

list.main = unique((list.break$ref))
l = sapply(c(1:length(list.main)), function(x) {
  list.s = list.break[which(list.break$ref == list.main[x]),]
  r = length(unique(list.s$brp))
  print(x)
  if(r ==1){
    y = rep(2, nrow(list.s))
  }else{
    y = rep(NA, nrow(list.s))
    for (i in c(1:nrow(list.s))) {
      case1 = list.s[i,]
      case2 = list.s[(which(list.s$nb == case1$nb & list.s$brp > case1$brp)[1]),]
      if(all(is.na(case2)) == FALSE){
        name1 = paste0(case1$ref,".",as.character(case1$brp), ".", case1$nb)
        name2 = paste0(case2$ref,".",as.character(case2$brp), ".", case2$nb)
        dat1 = data.res[[name1]]
        dat2 = data.res[[name2]]
        if(all(is.na(dat1$gps.gps))==TRUE){
          y[i] = 0
        }else if(all(is.na(dat1$gps.gps))==FALSE & all(is.na(dat2$gps.gps))==TRUE){
          y[i] = 2
        } else if(all(is.na(dat1$gps.gps))==FALSE & all(is.na(dat2$gps.gps))==FALSE){
          end.day = data.res[[name1]]$date[max(which(is.na(data.res[[name1]]$gps.gps) ==FALSE))]
          beg.day = data.res[[name2]]$date[min(which(is.na(data.res[[name2]]$gps.gps) ==FALSE))]
          if(end.day>beg.day){
            y[i] = 1
          }else{
            y[i] = 2
          }
        }
      }else{
        y[i] = 2
      }
    }
  }
  return(y)
})
list.break$side = unlist(l)
# sort out all segment will be used 
list.break = list.break[which(list.break$side!=0),]
list.break = list.break[-which(list.break$nb == "pama"),]

# compute the mean of differences for all segments with/without vertical correction
series1 = read.series(path_series = path_series_nearby, station = "pama", na.rm = 0, add.full = 0)
series2 = read.series(path_series = path_series_main, station = "medi", na.rm = 0, add.full = 0)

mean.all = data.frame(matrix(NA, ncol = 12, nrow = nrow(list.break)))
sd.all = data.frame(matrix(NA, ncol = 12, nrow = nrow(list.break)))

for (i in c(1:nrow(list.break))) {
  case.name = paste0(list.break$ref[i],".",as.character(list.break$brp[i]), ".", list.break$nb[i])
  # new vertical correction 
  new.ver = data.res[[case.name]]
  # old vertical correction 
  old.ver = old.ver0[[case.name]]
  # without vertical correction
  # choose the side (left or both) 
  if(list.break$side[i] == 1){
    list.day0 = new.ver$date[which(is.na(new.ver$gps.gps) == FALSE)]
    list.day0 = list.day0[which(list.day0<= list.break$brp[i])]
  }else if(list.break$side[i] == 2){
    list.day0 = new.ver$date[which(is.na(new.ver$gps.gps) == FALSE)]
  }
  if(length(list.day0) > 100) {
    old.ver =  old.ver[which(old.ver$date %in% list.day0 ),]
    raw.series = merge_series(path1 = path_series_main, path2 = path_series_nearby, 
                              station1 = list.break$ref[i], station2 = list.break$nb[i],
                              list.day = list.day0)
    ind.test = c(2,3,4,6)
    mean.all[i,] <- c( sapply(ind.test, function(x) mean( raw.series[[x]], na.rm = TRUE )),
                       colMeans(old.ver[list.test[ind.test]], na.rm = TRUE),
                       colMeans(new.ver[list.test[ind.test]], na.rm = TRUE))
    sd.all[i,] <- c( sapply(ind.test, function(x) sd( raw.series[[x]], na.rm = TRUE )),
                     sapply(ind.test, function(x) sd( unlist(old.ver[list.test[x]]), na.rm = TRUE )),
                     sapply(ind.test, function(x) sd( unlist(new.ver[list.test[x]]), na.rm = TRUE )))
    print(i)
  }
}

a = data.frame( value = colMeans(abs(sd.all), na.rm = TRUE))
a$ver = rep(c("Raw", "Ver1", "Ver2"), each =4)
a$variable =  rep(list.test[ind.test], 3)
ggplot(a, aes(x = variable, y = value, col = ver))+geom_point()+theme_bw()+ylab("Mean of the std of differences")

colnames(mean.all) <- paste( rep(list.test[ind.test], 3), rep(c("Raw", "Ver1", "Ver2"), each =4), sep = "")
dat = data.frame(raw = mean.all$gps1.eraRaw, ver1 = mean.all$gps1.eraVer1, ver2 = mean.all$gps1.eraVer2)
a = reshape2::melt(dat, id = "raw")

ggplot(a, aes(x = raw, y = value, col = variable))+geom_point()+theme_bw()+ylab("corrected")+
  geom_abline(slope = 1)+ylim(c(-4,7))+ xlim(c(-4,7))




