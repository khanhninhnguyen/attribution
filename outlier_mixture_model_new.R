source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))

# # Normalize data first:  ------------------------------------------------
window.thres = 2
data.cr = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
name.var = list.test[2]
dist.mean <- data.frame(matrix(NA, nrow = 0, ncol = (length(seq(-30,30,0.05))-1)))
data.all <- list()

for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  data.i = data.cr[[i]]
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  before =  data.i[which(data.i$date <= breakpoint),]
  after =  data.i[which(data.i$date > breakpoint),]
  if(nrow(na.omit(before)) > 30){
    bef.norm = one.step.norm(Y = before, list.test[1])
    if(length(bef.norm)>30){
      bef.norm.all <- list()
      for (k in c(1:6)) {
        norm.dat = one.step.norm(Y = before, list.test[k])
        names(norm.dat) <- NULL
        bef.norm.all[[list.test[k]]] <- norm.dat
      }
      bef.norm.all[["date"]] = before$date
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$bef <- bef.norm.all
    }
  }
  if(nrow(na.omit(after)) > 30){
    aft.norm = one.step.norm(Y = after, list.test[1])
    if(length(aft.norm) > 30){
      aft.norm.all <- list()
      for (k in c(1:6)) {
        norm.dat <- one.step.norm(Y = after, list.test[k])
        names(norm.dat) <- NULL
        aft.norm.all[[list.test[k]]] <- norm.dat
      }
      aft.norm.all[["date"]] = after$date
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$aft <- aft.norm.all
    }
  }
}
# save(dist.mean, file = paste0(path_results,"attribution/dist.mean_2years_", nearby_ver,".RData"))
save(data.all, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"onestep.RData"))
dat = get(load( file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"onestep.RData")))

# test for station GOPE --------------------------------------------
case.name = names(data.cr)[i]
d <- data.frame(Group = rep(c("raw","norm"), each = length(data.all[[case.name]]$bef$gps.era1)), 
                Sample =c(before$gps.era1, data.all[[case.name]]$bef$gps.era1))
d$Group <- as.factor(d$Group)
ggplot(d, aes(x = Sample, colour = Group)) + 
  geom_density()+scale_x_continuous(breaks=seq(-6,7,1))+  theme_bw()


plot(x, type = "l")
points(x = a$point.rm, y = x[a$point.rm], col = "red", cex=2)
points(x = b$point.rm, y = x[b$point.rm], col = "green")


d <- data.frame(Group = rep(c("norm", "s1","s2"), each = length(x)), 
                Sample = c(x, x1, x2))
                
qplot(sample=Sample, data=d, color=as.factor(Group))+theme_bw()

                
                
                


# see the qqplot of all cases ---------------------------
name.var = list.test[2]

res.x = list()
res.y = list()
for (i in c(1:length(dat))) {
  x = dat[[i]]$bef
  if(is.null(x)==FALSE & sum(is.na(x))!= length(x) ){
    a = qqnorm(x, plot.it = FALSE)
    res.x[[i]] <- a$x
    res.y[[i]] <- a$y
  }
}
all.x= unlist(res.x, use.names = FALSE)
lim = seq(-3.2,3.2,0.01)
limy = rep(NA, length(lim))
for (r in c(2:length(lim))) {
  mini = lim[r-1]
  maxi = lim[r]
  y <- c()
  for (l in c(1: length(res.y))) {
    ind = which(res.x[[i]] > mini & res.x[[i]] < maxi)
    if(length(ind)!=0){ y = c(y,res.y[[i]][ind])}
  }
  limy[r] <- mean(y, na.rm = TRUE)
}
plot( lim, limy, main = "Mean of QQ plots", xlab= "Theoretical Quantiles", ylab = "Sample Quantiles")



# screening all data ------------------------------------------------------

n1 = length(dat)
res <- data.frame(matrix(NA, nrow = n1, ncol = 3))
length.rm = data.frame(matrix(NA, nrow = n1, ncol = 2))
for (i in seq(1,n1,1)) {
  x = dat[[i]]$bef$gps.gps
  date1 = dat[[i]]$bef$date
  if(is.null(x)==FALSE & sum(is.na(x))!= length(x) ){
    a = screen.O(x)
    b = screen.diff.o(x, dif = 1, thres1 = 3/(sqrt(2)), thres2 = 3,sdt = 1)
    s1 = rep(0, length(x))
    s1[a$point.rm] <- 1
    s2 = rep(0, length(x))
    s2[b$point.rm] <- 2
    x1 <- x
    x1[a$point.rm] <- NA
    x2 <- x
    x2[b$point.rm] <- NA
    res[i,] <- c(shapiro.test(x)$p.value, shapiro.test(x1)$p.value, shapiro.test(x2)$p.value)
    length.rm[i,] <- c(length(a$point.rm), length(b$point.rm))
    Y = data.frame( date = date1, val = x, s = as.factor(s1+s2))
    Y$s <- sapply(Y$s, function(x) {
      if(x==0){y = "n"
      }else if (x ==1) {y="s1"
      }else if(x==2){y="s2"
      }else{y="both"}} )
    
    Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
    p <- ggplot(data = Y1, aes(x = date, y= val))+
      geom_line()+theme_bw()+geom_point( aes(col = s), size = 0.7)+
      scale_color_manual(values=c("n"= "black", "s1"= " green", "s2" = "red", "both" ="blue"))+
      labs(subtitle = paste0(names(data.all)[i], "    P.val = ", paste(round( res[i,], digits = 3), collapse  = ",  ")))+
      theme(axis.text = element_text(face="bold",size=14))
    jpeg(paste0(path_results,"attribution/test.scr/",names(data.all)[i], "b.jpeg"),
         width = 4000, height = 1800,res = 300) # change name
    print(p)
    dev.off()
  }
}

colnames(res) <- c("raw", "scr1", "scr2")
boxplot(res, main = "boxplot of p-value")
save(res, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"p.val.RData"))

r = which(res$scr1<0.01)
which(names(dat) == "medi.2006-05-19.rovi")

d <- data.frame(Group= rep(c("raw","scr1","scr2"), each = length(bef.norm)), Sample=c(before$gps.era1, bef.norm))
qplot(sample=Sample, data=d, color=as.factor(Group))

# test screening method on simulation 

res <- data.frame(matrix(NA, ncol = 2, nrow = 1000))
list.rm <- c()
for (j in c(1:1000)) {
  set.seed(j)
  x = rnorm(1000,0,1)
  x[seq(100, 1000, 100)] <- rep(5, 10)
  r = screen.O(x)
  list.rm <- c(list.rm, r$point.rm)
  res[j,] <- c(shapiro.test(x)$p.value, shapiro.test(r$data)$p.value)
}
list.rm[! duplicated(list.rm)]
summary(res)

which(names(dat) == "auck.2005-11-07.hamt")


# see the percentage of normalized data out of 3 sigma in +/-1 year only but sliding var and sliding mean is computed 
# from the day before to avoid the edge effect

window.thres = 2
data.cr = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
name.var = list.test[2]
dist.mean <- data.frame(matrix(NA, nrow = 0, ncol = (length(seq(-30,30,0.05))-1)))
data.all <- list()

for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  data.i = data.cr[[i]]
  data.i <- tidyr::complete(data.i, date = seq(min(data.i$date), max(data.i$date), by = "day"))
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  before =  data.i[which(data.i$date <= breakpoint),]
  after =  data.i[which(data.i$date > breakpoint),]
  if(nrow(na.omit(before)) > 30){
    ind.sta = which(before$date == max(breakpoint %m+% years(-1), min(before$date)))
    if(length(ind.sta)>0){
      bef.norm = one.step.norm(before, name.var = name.var, estimator = "Sca") 
      bef.norm <- bef.norm[ind.sta:length(bef.norm)]
      data.all <- c(data.all, length(which(abs(bef.norm)>3))/length(bef.norm))
    }
  }
  if(nrow(na.omit(after)) > 30){
    ind.end = which(after$date == min(breakpoint %m+% years(1), max(after$date)))
    if(length(ind.end) > 0){
      aft.norm = one.step.norm(after, name.var = name.var, estimator = "Sca") 
      aft.norm = aft.norm[1:ind.end]
      data.all <- c(data.all, length(which(abs(aft.norm)>3))/ length(aft.norm))
    }
  }
}






# normalized stations 
da = data.cr$`auck.2005-11-07.whng`

y = one.step.norm(da, name.var = "gps.gps", estimator = "Sca")
save(Y, file = paste0(path_results,"auck.2005-11-07.whng.RData"))
Y = data.frame(date = da$date, gps.gps = y)
rownames(Y) <- NULL
a = get(load( file = paste0(path_results,"auck.2005-11-07.whng.RData")))


# investigate the concatenated data  --------------------------------------
# took normalized data 
dat = get(load( file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"normalized1.RData")))
# FIND WHY WE EXCESS AT 0
all.dat <- c()
zero <- data.frame(matrix(NA, ncol = 2, nrow = length(dat)))
for (r in 1:length(dat)) {
  if(length(dat[[r]]$bef$gps.gps)>30){
    ind1 = length(dat[[r]]$bef$gps.gps)-30
    all.dat = c(all.dat, dat[[r]]$bef$gps.gps[30:ind1])
  } 
  if(length(dat[[r]]$aft$gps.gps)>30){
    ind2 = length(dat[[r]]$aft$gps.gps)-30
    all.dat = c(all.dat, dat[[r]]$aft$gps.gps[30:ind2])
  }
  if (length(dat[[r]]$bef$gps.gp) >0){
    zero[r,1] <- length(which(dat[[r]]$bef$gps.gps<=0 & dat[[r]]$bef$gps.gps>-0.02))
  }
  if (length(dat[[r]]$aft$gps.gp) >0){
    zero[r,2] <- length(which(dat[[r]]$aft$gps.gps<=0 & dat[[r]]$aft$gps.gps>-0.02))
  }
}
save(all.dat, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"normalized.concat.RData"))
x = get(load(file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"normalized.concat.RData")))
x = na.omit(all.dat)
a = hist(x, 
         main = "Histogram of normalized data", 
         xlab = "", 
         breaks = 100,
         xlim=c(-5, 5),
         prob = TRUE)

# qqnorm(x)
x2 <- a$breaks

fun1 <- dnorm(x2, mean = mean(all.dat, na.rm = T), sd = sd(all.dat, na.rm = T))
fun2 <- dnorm(x2, mean = 0, sd = 1)

lines(x2, fun1, col = 2, lwd = 2)
lines(density(na.omit(x)), col = 4, lwd = 2)
lines(x2, trueCDF, col = 2, lwd = 2)

case.name = names(data.cr)[i]
station.ref = substr(case.name ,start = 1, stop = 4)
station.near = substr(case.name ,start = 17, stop = 20)
data.i = data.cr[[i]]
breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
before =  data.i[which(data.i$date <= breakpoint),]
Y1 <- tidyr::complete(before, date = seq(min(before$date), max(before$date), by = "day"))
bef = Y1[which(Y1$date > (breakpoint-366)),]
plot(dat$`albh.2015-12-29.pgc5`$bef$gps.gps, type="l")
abline(h = 0)
abline(h = -0.05)

lines(bef$gps.gps, col="red")

