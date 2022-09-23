# This function is used for pairing data,
# Minimum of points is set at 100
window.thres = 2
data.cr = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
data.all <- data.cr
n0 = length(data.all)
data.out <- list()
data.out1 <- list()

for (i in c(1:n0)) {
  data.i = data.all[[i]]
  # data.i <- tidyr::complete(data.i, date = seq(min(data.i$date), max(data.i$date), by = "day"))
  
  case.name = names(data.all)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  
  day.sta = breakpoint %m+% years(-1)
  day.end = breakpoint %m+% years(+1)
    
  before = na.omit(data.i[which(data.i $date >= day.sta & data.i $date <= breakpoint),])
  after = na.omit(data.i[which(data.i $date > breakpoint & data.i $date <= day.end),])
  
  listday1 = after$date %m+% years(-1)
  listday2 = before$date[which(before$date %in% listday1 == TRUE )]
  listday3 = after$date[which(listday1 %in% listday2 == TRUE & duplicated(listday1) == FALSE)]
  
  if(length(listday2) > 100){
    before1 = before[which(before$date %in% listday2 == TRUE),]
    after1 = after[which(after$date %in% listday3 == TRUE),]
    data.out1 [[case.name]] <- rbind(before1, after1)
    # compute the pair data
    out <- data.frame( date = before1$date)
    for (j in c(1:6)) {
      out[list.test[j]] <- after1[list.test[j]] - before1[list.test[j]]
    }
   data.out [[case.name]] <- out 
  } 
  
}

save(data.out, file = paste0(path_results,"attribution/data.all_1year_", nearby_ver,"paired.RData"))

Y = data.out$`gope.2009-05-08.zdib`
save(Y, file = paste0(path_results,"attribution/data.all_1years_", nearby_ver,".auck.2005-11-07.whng.RData"))


###############@----------------------------
##### TESTS on real data

# load("data.all_1years_NGL.gope.2009-05-08.zdib.paired.RData")
#data.all_1years_NGL.auck.2005-11-07.whng.paired.RData
head(Y)
y=Y$gps.era
dates=Y$date
n=length(y)
times=1:n

plot(dates,y,pch=16,cex=0.8)


residus=y-median(y)
plot(dates,residus^2)
lissage=loess(residus^2~times,degree=2) # revoir pour une meilleure estimation
lines(dates,lissage$fitted,col="red")

z=residus/sqrt(lissage$fitted)
plot(dates,z)
hist(z,breaks = 20)

prob.z=dnorm(z,0,1)
f=2
OS=log(prob.z)^(2*f)
plot(dates,OS)
ScOS=(OS-min(OS))/(max(OS)-min(OS))*10
plot(dates,ScOS)
thresh=7
abline(h=thresh,col="red")
rg=which(ScOS>thresh)

plot(dates,z,pch=16,cex=0.8)
points(dates[rg],z[rg],col="red")

plot(dates,y,pch=16,cex=0.8)
points(dates[rg],y[rg],col="red")






#mean(y,na.rm=TRUE)
#y[c(50,60,80)]=NA
reg=lm(z~1)
### OLS
print(coef(summary(reg)))
### Test with variance correction
library(sandwich)
V_HAC <- sandwich::NeweyWest(reg, lag = 1)
print(lmtest::coeftest(reg, sandwich::NeweyWest(reg, lag = 1))[, ])

