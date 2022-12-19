# study the different veriecal correction 
window.thres = 10
distances <- get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData")))
new.ver = get(load(file = paste0(path_results,"attribution/four_main_series_",window.thres,"year_", nearby_ver,".RData")))
old.ver0 = get(load(file = paste0(path_results,"attribution/OLD_VERTICAL/four_main_series_",window.thres,"year_", nearby_ver,".RData")))
ver0 = get(load(file = paste0(path_results,"attribution/four_main_series_",window.thres,"year_", nearby_ver,"0.RData")))

old.ver = old.ver0[names(new.ver)]

h = sapply(c(1:length(old.ver)), function(x){
  station.ref = substr(names(old.ver)[x] ,start = 1, stop = 4)
  station.seg = substr(names(old.ver)[x] ,start = 17, stop = 20)
  ind = which(distances$main == station.ref & distances$nearby == station.seg)
  distances$ver.dist[ind]
})
seg = c(1:3650)
gps.gps1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$GPS[seg], na.rm = TRUE) - mean(old.ver[[x]]$GPS1[seg], na.rm = TRUE))
gps.gps2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$GPS[seg], na.rm = TRUE) - mean(new.ver[[x]]$GPS1[seg], na.rm = TRUE))
gps.gps0 = sapply(c(1:length(old.ver)), function(x)  mean(ver0[[x]]$GPS[seg], na.rm = TRUE) - mean(ver0[[x]]$GPS1[seg], na.rm = TRUE))

gps.gps1[which(h<0)] = -gps.gps1[which(h<0)] 
gps.gps2[which(h<0)] = -gps.gps2[which(h<0)] 
gps.gps0[which(h<0)] = -gps.gps0[which(h<0)] 

plot(abs(h), gps.gps2)
plot(gps.gps1, gps.gps2)
plot(gps.gps0, gps.gps1)


gps.era1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$GPS[seg], na.rm = TRUE) - mean(old.ver[[x]]$ERA1[seg], na.rm = TRUE))
gps.era2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$GPS[seg], na.rm = TRUE) - mean(new.ver[[x]]$ERA1[seg], na.rm = TRUE))
gps.era0 = sapply(c(1:length(old.ver)), function(x)  mean(ver0[[x]]$GPS[seg], na.rm = TRUE) - mean(ver0[[x]]$ERA1[seg], na.rm = TRUE))

gps1.era1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$GPS1[seg], na.rm = TRUE) - mean(old.ver[[x]]$ERA[seg], na.rm = TRUE))
gps1.era2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$GPS1[seg], na.rm = TRUE) - mean(new.ver[[x]]$ERA[seg], na.rm = TRUE))
gps1.era0 = sapply(c(1:length(old.ver)), function(x)  mean(ver0[[x]]$GPS1[seg], na.rm = TRUE) - mean(ver0[[x]]$ERA[seg], na.rm = TRUE))

era.era1 = sapply(c(1:length(old.ver)), function(x)  mean(old.ver[[x]]$ERA[seg], na.rm = TRUE) - mean(old.ver[[x]]$ERA1[seg], na.rm = TRUE))
era.era2 = sapply(c(1:length(old.ver)), function(x)  mean(new.ver[[x]]$ERA[seg], na.rm = TRUE) - mean(new.ver[[x]]$ERA1[seg], na.rm = TRUE))
era.era0 = sapply(c(1:length(old.ver)), function(x)  mean(ver0[[x]]$ERA[seg], na.rm = TRUE) - mean(ver0[[x]]$ERA1[seg], na.rm = TRUE))

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



