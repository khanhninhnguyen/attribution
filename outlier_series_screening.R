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