# distribution of real data 


# # Normalize data first:  ------------------------------------------------

# 1. Estimate the monthly variance in period +/- 2 years from the breakpoint: 

meta.compare.ref =  get(load(file = paste0(path_results,"validation/",nb_test.ref,"-",criterion,"metacompa",screen.value="",".RData")))

meta.compare.near =  get(load(file = paste0(path_results,"validation/",nb_test.near,"-",criterion,"metacompa",screen.value="",".RData")))

# Read nearby stations, distance, list of breakpoints
distances <- get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData")))
nearby_list <- get(load(file = paste0("list_nearby", nearby_ver, ".RData")))
nearby_list1 <- nearby_list[which(is.na(nearby_list$Level1) == FALSE),]

length.series <- c()
nb.brp <- c()
data.all.sd <- list()
data.all.iqr <- list()
dist.mean <- data.frame(matrix(NA, nrow = 0, ncol = (length(seq(-30,30,0.05))-1)))
# Read list nearby 
for (j in c(1:nrow(nearby_list1))){
  station.ref.j = nearby_list1$ref.sta[j]
  series.ref <- read.series(path_series = path_series_main, station = station.ref.j, na.rm = 1, add.full = 0)
  seg.ref <- list.brp[which(list.brp$name == station.ref.j),]
  nearby.list.j = na.omit(unlist(nearby_list1[which(nearby_list1$ref.sta == station.ref.j ),], use.names = FALSE))[-1]
  
  for (k in c(1:nrow(seg.ref))) {
    breakpoint <- seg.ref$detected[k]
    begin.point <- breakpoint %m+% years(-2)   
    end.point <- breakpoint %m+% years(2)
    
    for (l in c(1:length(nearby.list.j))) {
      series.near <- read.series(path_series = path_series_nearby, station = nearby.list.j[l], na.rm = 1, add.full = 0)
      both <- data.frame() 
      
      # vertical correction -----------------------------------------------------
      n.main = which(distances$main == station.ref.j & distances$nearby == nearby.list.j[l])
      ver.dist = distances$ver.dist[n.main]   #vertical disstance here is main-nearby
      series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
      series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
      # join 2 data frame
      both = inner_join(series.ref,series.near, by = "date")  # checked
      both <- both[which(both$date > begin.point & both$date <= end.point),]

      if(nrow(both) >30){
        length.series <- c(length.series, nrow(both))
        nearby.seg = meta.compare.near[which(meta.compare.near$name == nearby.list.j[l]),]
        if(nrow(nearby.seg) !=0){
          con = which( begin.point < nearby.seg$detected & end.point > nearby.seg$detected)
          con1 = which(abs(nearby.seg$detected[con] - breakpoint) <= 10)
          nb.brp <- c(nb.brp,length(con) - length(con1))
        }else{
          nb.brp <- c(nb.brp,0)
        }
        
        
        name.test = list.test[1]
        var1 = diff.var(c(name.test))[1]
        var2 = diff.var(c(name.test))[2]
        
        data = data.frame(month = both$month.x, year = both$month.x)
        data[,name.test] = both[,c(var1)] -  both[,c(var2)] # compute the difference
        
        month.std = RobEstiMonthlyVariance.diff.S(Y = (na.omit(data)), name.var = name.test, alpha = 0) # alpha maybe not need to correct here, just scale the data 
        
        std.t = rep(NA, nrow(data))
        for (m in 1:12) {
          if(month.std[m] != 0 & is.na(month.std[m]) == FALSE){
            std.t[which(as.numeric(data$month)==m)] <-  month.std[m]
          }
        }  
        if(sum(is.na(std.t)) != 0){
          data <- data[-which(is.na(std.t)),]
          data[,name.test] = data[,name.test]/(std.t[!is.na(std.t)])
        } else{
          data[,name.test] = data[,name.test]/(std.t)
        }
        data.all.sd[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]] <- data
      }
    }
  }
}
save(data.all.sd, file = paste0(path_results,"attribution/nor_series_4year_", nearby_ver,".RData"))


res = data.frame(size = length.series, nb.brp = nb.brp)

t1 = data.all$`albh.2015-12-29.sac4`
qqplot(t1)

# normalize by 2 steps median  ----------------------------------------------------

for (j in c(1:nrow(nearby_list1))){
  station.ref.j = nearby_list1$ref.sta[j]
  series.ref <- read.series(path_series = path_series_main, station = station.ref.j, na.rm = 1, add.full = 0)
  seg.ref <- list.brp[which(list.brp$name == station.ref.j),]
  nearby.list.j = na.omit(unlist(nearby_list1[which(nearby_list1$ref.sta == station.ref.j ),], use.names = FALSE))[-1]
  
  for (k in c(1:nrow(seg.ref))) {
    breakpoint <- seg.ref$detected[k]
    begin.point <- breakpoint %m+% years(-2)   
    end.point <- breakpoint %m+% years(2)
    
    for (l in c(1:length(nearby.list.j))) {
      series.near <- read.series(path_series = path_series_nearby, station = nearby.list.j[l], na.rm = 1, add.full = 0)
      both <- data.frame() 
      
      # vertical correction -----------------------------------------------------
      n.main = which(distances$main == station.ref.j & distances$nearby == nearby.list.j[l])
      ver.dist = distances$ver.dist[n.main]   #vertical disstance here is main-nearby
      series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
      series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
      # join 2 data frame
      both = inner_join(series.ref,series.near, by = "date")  # checked
      both <- both[which(both$date > begin.point & both$date <= end.point),]
      
      if(nrow(both) >30){
        length.series <- c(length.series, nrow(both))
        nearby.seg = meta.compare.near[which(meta.compare.near$name == nearby.list.j[l]),]
        if(nrow(nearby.seg) !=0){
          con = which( begin.point < nearby.seg$detected & end.point > nearby.seg$detected)
          con1 = which(abs(nearby.seg$detected[con] - breakpoint) <= 10)
          nb.brp <- c(nb.brp,length(con) - length(con1))
        }else{
          nb.brp <- c(nb.brp,0)
        }
        
        name.test = list.test[1]
        var1 = diff.var(c(name.test))[1]
        var2 = diff.var(c(name.test))[2]
        
        data = data.frame(month = both$month.x, year = both$month.x, date = both$date)
        data[,name.test] = both[,c(var1)] -  both[,c(var2)] # compute the difference
        before = data[which(data$date <= breakpoint),]
        after = data[which(data$date > breakpoint),]
        month.std = RobEstiMonthlyVariance.diff.S30(Y = (na.omit(data)), name.var = name.test, alpha = 0) # alpha maybe not need to correct here, just scale the data 
        
        std.tb = rep(NA, nrow(before))
        std.ta = rep(NA, nrow(after))
        
        for (m in 1:12) {
          if(month.std[m] != 0 & is.na(month.std[m]) == FALSE){
            std.ta[which(as.numeric(after$month)==m)] <-  month.std[m]
            std.tb[which(as.numeric(before$month)==m)] <-  month.std[m]
          }
        }  

        if(sum(is.na(std.tb)) != 0){
          before <- before[-which(is.na(std.tb)),]
          before[,name.test] = (before[,name.test] - median( before[,name.test])) /(std.tb[!is.na(std.tb)])
        }else if(sum(is.na(std.ta)) != 0){
          after <- after[-which(is.na(std.ta)),]
          after[,name.test] = (after[,name.test] - median(after[,name.test]))/(std.ta[!is.na(std.ta)])
        } else if (sum(is.na(std.tb)) == 0 & sum(is.na(std.ta)) == 0) {
          before[,name.test] = (before[,name.test] - median( before[,name.test])) /(std.tb[!is.na(std.tb)])
          after[,name.test] = (after[,name.test] - median(after[,name.test]))/(std.ta[!is.na(std.ta)])
        }
        data.all.sd[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]][1] <- before
        data.all.sd[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]][2] <- after
        if(nrow(before) >30){
          hist.b = hist(before[,name.test], breaks = seq(-30,30,0.05), plot = FALSE)
          dist.mean <- rbind(dist.mean, hist.b$counts)
        }
        if(nrow(after) > 30){
          hist.a = hist(after[,name.test], breaks = seq(-30,30,0.05), plot = FALSE)
          dist.mean <- rbind(dist.mean, hist.a$counts)
        }
        
        
      }
    }
  }
}


# error when the number of point for months is not limit (should be >30)
for (j in c(1:nrow(nearby_list1))){
  station.ref.j = nearby_list1$ref.sta[j]
  series.ref <- read.series(path_series = path_series_main, station = station.ref.j, na.rm = 1, add.full = 0)
  seg.ref <- list.brp[which(list.brp$name == station.ref.j),]
  nearby.list.j = na.omit(unlist(nearby_list1[which(nearby_list1$ref.sta == station.ref.j ),], use.names = FALSE))[-1]
  
  for (k in c(1:nrow(seg.ref))) {
    breakpoint <- seg.ref$detected[k]
    begin.point <- breakpoint %m+% years(-2)   
    end.point <- breakpoint %m+% years(2)
    
    for (l in c(1:length(nearby.list.j))) {
      series.near <- read.series(path_series = path_series_nearby, station = nearby.list.j[l], na.rm = 1, add.full = 0)
      both <- data.frame() 
      
      # vertical correction -----------------------------------------------------
      n.main = which(distances$main == station.ref.j & distances$nearby == nearby.list.j[l])
      ver.dist = distances$ver.dist[n.main]   #vertical disstance here is main-nearby
      series.near$GPS <- series.near$GPS*(exp(-5*(10**-4)*ver.dist)) ## change to the new correction: iwv(z1) = iwv(z2) exp(-0.0005/(z1-z2))
      series.near$ERAI <- series.near$ERAI*(exp(-5*(10**-4)*ver.dist))
      # join 2 data frame
      both = inner_join(series.ref,series.near, by = "date")  # checked
      both <- both[which(both$date > begin.point & both$date <= end.point),]
      
      if(nrow(both) >30){
        length.series <- c(length.series, nrow(both))
        nearby.seg = meta.compare.near[which(meta.compare.near$name == nearby.list.j[l]),]
        if(nrow(nearby.seg) !=0){
          con = which( begin.point < nearby.seg$detected & end.point > nearby.seg$detected)
          con1 = which(abs(nearby.seg$detected[con] - breakpoint) <= 10)
          nb.brp <- c(nb.brp,length(con) - length(con1))
        }else{
          nb.brp <- c(nb.brp,0)
        }
        
        name.test = list.test[1]
        var1 = diff.var(c(name.test))[1]
        var2 = diff.var(c(name.test))[2]
        
        data = data.frame(month = both$month.x, year = both$month.x, date = both$date)
        data[,name.test] = both[,c(var1)] -  both[,c(var2)] # compute the difference
        before = data[which(data$date <= breakpoint),]
        after = data[which(data$date > breakpoint),]
        month.std = RobEstiMonthlyVariance.diff.S30(Y = (na.omit(data)), name.var = name.test, alpha = 0) # alpha maybe not need to correct here, just scale the data 
        
        std.tb = rep(NA, nrow(before))
        std.ta = rep(NA, nrow(after))
        
        for (m in 1:12) {
          if(month.std[m] != 0 & is.na(month.std[m]) == FALSE){
            std.ta[which(as.numeric(after$month)==m)] <-  month.std[m]
            std.tb[which(as.numeric(before$month)==m)] <-  month.std[m]
          }
        }  
        
        if(sum(is.na(std.tb)) != 0){
          before <- before[-which(is.na(std.tb)),]
          sig.b = before[,name.test]/(std.tb[!is.na(std.tb)])
          before[,name.test] = (sig.b - median(sig.b))/IQR(sig.b)
        }else if(sum(is.na(std.ta)) != 0){
          after <- after[-which(is.na(std.ta)),]
          sig.a = after[,name.test]/(std.ta[!is.na(std.ta)])
          after[,name.test] = (sig.a - median(sig.a))/IQR(sig.a)
        } else if (sum(is.na(std.tb)) == 0 & sum(is.na(std.ta)) == 0) {
          sig.b = before[,name.test]/(std.tb[!is.na(std.tb)])
          before[,name.test] = (sig.b - median(sig.b))/IQR(sig.b)
          sig.a = after[,name.test]/(std.ta[!is.na(std.ta)])
          after[,name.test] = (sig.a - median(sig.a))/IQR(sig.a)
        }
        data.all.sd[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]][1] <- before
        data.all.sd[[paste0(station.ref.j,".",as.character( breakpoint), ".", nearby.list.j[l])]][2] <- after
        if(nrow(before) >30){
          hist.b = hist(before[,name.test], breaks = seq(-30,30,0.05), plot = FALSE)
          dist.mean <- rbind(dist.mean, hist.b$counts)
        }
        if(nrow(after) > 30){
          hist.a = hist(after[,name.test], breaks = seq(-30,30,0.05), plot = FALSE)
          dist.mean <- rbind(dist.mean, hist.a$counts)
        }
      }
    }
  }
}

ggqqplot(before, x = "gps.era",
         ggtheme = theme_pubclean())



