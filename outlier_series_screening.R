# this prog used to investigate the screening method 
data = get(load(file = paste0(path_results,"attribution/six_diff_series_1year_rm_crenel_restricted_closed_brp_", nearby_ver,".RData")))

which(names(data) == "auck.2016-07-22.whng")
which(names(data) == "brus.1997-05-14.eijs")
screen.out <- function(x) {
  y = na.omit(x)
  Q1 <- quantile(x, .25, na.rm = TRUE)
  Q3 <- quantile(x, .75, na.rm = TRUE)
  IQR <- IQR(x, na.rm = TRUE)
  up = Q3+IQR*1.5
  down = Q1-IQR*1.5
  removed <- which(x<down | x>up)
  if (length(removed) >0){
    x.screened = x[-removed]
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = removed, up = up, down = down))
}


# test ------

test = na.omit(data[[31]]$gps.era)


screen.out.sliding <- function(x,date) {
  dat = data.frame(x, date)
  dat.full = na.omit(dat)
  # solve the problem of gaps (period = 2 months)
  return(list(data = x.screened, point.rm = c(which(x<down | x>up))))
}
a = screen.out(test)
dat = data.frame(ind = c(1:length(test)), iwv = test)
ggplot( dat, aes(x=ind, y=iwv)) +
  geom_line() +
  geom_point() + 
  theme_bw()+ 
  geom_hline(yintercept = 1.225)+geom_hline(yintercept = -1.615)


for (i in c(1:(nrow(da)-59))) {
  samp = dat[c(i:(i+59)),c(1:2)]
  if (all(is.na(samp$x))){
    da$scr[i] <- NA
  }else{
    test = na.omit(samp)
    if(nrow(test) >= 30){
      # x.scr <- screen.out(test$x)$data
      # if(length(x.scr) != nrow(test)){
      #   da$scr[i] <- 1
      #   list.point.rm <- c(list.point.rm, test$date[screen.out(test$x)$point.rm])
      # }else{
      x.scr <- screen.out(diff(test$x))$data
      rm <- screen.out(diff(test$x))$point.rm
      rm.m = rm[-which(diff(rm) == 1)]
      if((length(x.scr)+1) != nrow(test)){
        da$scr[i] <- 1
        list.point.rm <- c(list.point.rm, test$date[rm.m])
      }else{
        da$scr[i] <- 0
      } 
    }else{
      da$scr[i] <- 2
    }
  }
}

list.rm = list.point.rm[!duplicated(list.point.rm)]
da$c <- 1
da$c[which(da$date %in% list.rm)] <- 2
da$c1 <- 1
da$c1[c(439,565,572,618)] <- 2

dat = da[-c(1:278),]
ggplot( dat, aes(x=date, y=x, col = c)) +
  geom_line() +
  geom_point() + 
  theme_bw()+ 
  geom_hline(yintercept = 2.732363 )+geom_hline(yintercept = -3.688095 )

a = diff(test)
screen.out(a)

two.steps.screening <- function(x,date){
  dat = data.frame(x = x, date = date)
  da <- tidyr::complete(dat, date = seq(min(dat$date), max(dat$date), by = "day"))
  
  
}



# list.rm = list.rm[-(which(diff(list.rm) == 1)+1)]




x = diff(da$x)

da$date[screen.out(da$x)$point.rm]



# sort all outliers 


for (l in 1:length(list.all)) {
  
}

# full -----

f = 59
dat = data.frame(x = data[[f]]$gps.gps, date = data[[f]]$date)
da <- tidyr::complete(dat, date = seq(min(dat$date), max(dat$date), by = "day"))
da$scr <- rep(NA, nrow(da))
da$up <- rep(NA, nrow(da))
da$down <- rep(NA, nrow(da))

list.point.rm <- as.Date(x = integer(0), origin = "1970-01-01")

for (i in c(1:(nrow(da)-59))) {
  ind.test = c(i:(i+59))
  samp = da[ind.test,2]
  if (all(is.na(samp$x))){
    da$scr[i] <- NA
  }else{
    y = unlist(samp)
    con = length(na.omit(y))
    if( con >= 30){
      up = screen.out(y)$up
      down =  screen.out(y)$down
      da$up[i+30] <- up
      da$down[i+30] <- down
      x = diff(y)
      x.scr <- screen.out(x)$data
      rm <- screen.out(x)$point.rm
      if(length(rm) > 0){
        da$scr[i] <- 1
        list.point.rm <- c(list.point.rm, da$date[rm+i])
      }else{
        da$scr[i] <- 0
      } 
    }else{
      da$scr[i] <- 2
    }
  }
}
list.rm = list.point.rm[!duplicated(list.point.rm)]
list.rm.n <- as.Date(x = integer(0), origin = "1970-01-01")
for (j in c(1:length(list.rm))){
  ind <- which(da$date == list.rm[j])
  subset <- da[c((ind-2): (ind+2)),]
  dist = abs(subset$x - mean(subset$x))
  list.rm.n <- c(list.rm.n, subset$date[which.max(dist)])
}
# keep only 1%
ind <- which(da$date %in% list.rm.n)
diff.x <- diff(da$x)  
dist.se <- abs(diff.x[ind])
nb.point = round(nrow(dat)/100)
thres = sort(dist.se, decreasing = TRUE)[nb.point]
list.rm.n <- da$date[ind[which(dist.se >= thres)]]

list.p <- da$date[screen.out(da$x)$point.rm]
list.pn <- as.Date(x = integer(0), origin = "1970-01-01")
for (k in c(1:length(list.p))) {
  ind <- which(da$date == list.p[k])
  subset <- da[c((ind-1): (ind+1)),]
  check.na = sum(is.na(subset$x))
  if(check.na >1){
    list.pn <- c(list.pn, list.p[k])
  }
}
list.all <- c(list.pn, list.rm.n)
da$c <- 1
da$c[which(da$date %in% list.all)] <- 2
ggplot( da, aes(x=date, y=x, col = c)) +
  geom_line(aes(x=date, y=up)) +
  geom_line(aes(x=date, y=down)) +
  # geom_line()+
  geom_point() + 
  scale_x_date(limits = as.Date(c(da$date[1], da$date[nrow(da)])),
                              breaks = function(x) seq.Date(from = da$date[1],  to = da$date[nrow(da)], by = "2 months"),
                              date_labels = "%Y-%m" )+
  ylab("GPS-GPS'")+
  theme_bw()


length(list.all)

