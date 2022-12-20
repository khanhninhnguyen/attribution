# new characterization
# choose the longest segment from the screened data 

win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.gps"
one.year=365

list.break = data.frame(ref = substr(names(dat), start = 1, stop = 4), 
                        brp = substr(names(dat), start = 6, stop = 15),
                        nb = substr(names(dat), start = 17, stop = 20))
list.break[] <- lapply(list.break, as.character)
list.break$brp = as.Date(list.break$brp , format = "%Y-%m-%d")
list.main = unique((list.break$ref))
length.seg = matrix(NA, ncol = 2, nrow = 0)
for (i in c(1:length(list.main))) {
  list.s = list.break[which(list.break$ref == list.main[i]),]
  list.nb = split(list.s, list.s$nb)
  for (j in c(1:length(list.nb))) {
    list.ij = paste0(list.nb[[j]]$ref,".",as.character(list.nb[[j]]$brp), ".", list.nb[[j]]$nb)
    data.ij = dat[list.ij]
    length.all <- sapply(c(1:length(data.ij)), function(x){
      seg1 = data.ij[[x]][c(1:3650),]
      seg2 = data.ij[[x]][-c(1:3650),]
      y = c(nb.consecutive(list.day = seg1$date, x = seg1$gps.gps), nb.consecutive(list.day = seg2$date, x = seg2$gps.gps))
    })
    length.seg <- rbind(length.seg, t(length.all))
  }
}
list.break$len1 = length.seg[,1]
list.break$len2 = length.seg[,2]

res <- data.frame(seg = list.seg, side = list.side)
nb.consecutive <- function(list.day, x){
  a = list.day[which(is.na(x)== FALSE)]
  b = ts(a) 
  y = length(which(diff(b)==1))
  return(y)
}



