# this prog used to investigate the screening method 
data = get(load(file = paste0(path_results,"attribution/six_diff_series_1year_rm_crenel_restricted_closed_brp_", nearby_ver,".RData")))

which(names(data) == "auck.2016-07-22.whng")

test = na.omit(data[[31]]$gps.era)

screen.out <- function(x) {
  y = na.omit(x)
  Q1 <- quantile(y, .25)
  Q3 <- quantile(y, .75)
  IQR <- IQR(y)
  up = Q3+IQR*1.5
  down = Q1-IQR*1.5
  x.screened = x[which(x>down & x<up)]
  print(c(up, down))
  return(list(data = x.screened, point.rm = c(which(x<down | x>up))))
}

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

dat = data.frame(x = data[[31]]$gps.era, date = data[[31]]$date)
da <- tidyr::complete(dat, date = seq(min(dat$date), max(dat$date), by = "day"))
da$scr <- rep(NA, nrow(da))
list.point.rm <- as.Date(x = integer(0), origin = "1970-01-01")

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
























