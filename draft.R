# draft -------------------------------------------------------------------

# test the scaling factor between the variance of the difference a --------

phi = 0.8
theta = -0.5
var0 = 1
var_x = (1 + 2*phi*theta+theta**2)*var0/(1-phi**2)
print(var_x)
scale.diff = 2*(1-phi)*(1+theta**2-theta+phi*theta)/((1 + 2*phi*theta+theta**2))
scale.ma = 2*(1+ theta**2-theta)/(1+theta**2)
nb.sim = 1000
var.est = rep(NA, nb.sim)
var.diff.est = rep(NA, nb.sim)

for (i in 1:nb.sim) {
  sim.series <- arima.sim(model = list(ar = phi, ma = theta), n = 1000, sd = sqrt(var0))
  sim.diff <- diff(sim.series)
  var.est[i] <- var(sim.series)
  var.diff.est[i] <- var(sim.diff)/scale.diff
}
summary(var.est)
summary(var.diff.est)


# test the impact of estimate 2n and the mean -----------------------------
n = 200
nb.sim = 1000
var.2n = rep(NA, nb.sim)
var.mean = rep(NA, nb.sim)

for (i in 1:nb.sim) {
  x = rnorm(n, mean = 0, sd = 1)
  var.2n[i] = var(x)
  var.mean[i] = (var(x[1:100]) + var(x[101:200]))/2
}
summary(var.2n)
summary(var.mean)
var(var.2n)
var(var.mean)


nearby_search_distance(coor.file, list.name, list.gnss, horizontal, vertical, version_name, nearby.ver)
  
a = read.table(file = validation.file.ref, header = TRUE)

list.all =list.files(paste0(path_homo,"34/"))
name.full = substr(list.all,start = 6, stop = 9)

name.full %in% list.nearby.station.homo





# check the parameters from ARMA ------------------------------------------

res <- data.frame(matrix(NA, ncol = 3, nrow = 100))
for (i in c(1:100)) {
  x = arima.sim(model = list(ar=0.5), n = 1000, sd=1) + 1
  fit = arima(x, order = c(1,0,0), include.mean = TRUE)
  res[i,] <- c(mean(x), fit$coef)
}

summary(res$X2)
summary(res$X3)

yt <- arima.sim(list(order=c(1,0,0), ar=.5), n=500)
xt <- yt + 10   


# impact of gaps on the param ---------------------------------------------

nb.sim = 1000
gaps.list = seq(10, 80, 20)
n0 = 360
sdv = 1
Res1 <- list()
for (k in 1:length(gaps.list)) {
  # ar.order <- data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
  # ma.order <- data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
  arma.order <- data.frame(matrix(NA, ncol = 3, nrow = nb.sim))
  for (i in 1:nb.sim) {
    # sim.ar <- arima.sim(model = list(ar = 0.3), n = n0, sd = sdv)
    # sim.ma <- arima.sim(model = list(ma = 0.3), n = n0, sd = sdv)
    sim.arma <- arima.sim(model = list(ar = 0.3, ma = 0), n = n0, sd = sdv)
    gaps =  rbinom(n0, 1, gaps.list[k]/100)
    ind.gaps = which(gaps != 0)
    # sim.ar[ind.gaps] <- NA
    # sim.ma[ind.gaps] <- NA 
    sim.arma[ind.gaps] <- NA
    
    # arfit2 = fit.arima(sim.ar, select.sig = 1)
    # mafit2 = fit.arima(sim.ma, select.sig = 1)
    # armafit2 = fit.arima(sim.arma, select.sig = 1)
  
    # ar.order[i,] <- 
    # ma.order[i,] <- c(unlist(mafit2))
    arma.order[i,c(1,2)] <- arima(sim.arma, order = c(1,0,0), method = "ML")$coef
    arma.order[i,3] <- classic(sim.arma)
  }
  # ar.2 = ar.order[which(ar.order$X3==1 & ar.order$X4==0),]
  # ma.2 = ma.order[which(ma.order$X3==0 & ma.order$X4==1),]
  # arma.2 = arma.order[which(arma.order$X3==1 & arma.order$X4==1),]
  # 
  # res = c(nrow(ar.2), nrow(ma.2), nrow(arma.2))
  # names(res) <- c("ar1", "ma1", "arma11")
  Res1[[k]] <-arma.order
} 

Res1.df = bind_rows(Res1)
b <- data.frame( value = Res1.df$X3, gaps = rep(as.factor(gaps.list), each = nb.sim ))

jpeg(paste0(path_results,"figure/ATTRIBUTION-TEST/" , "per-atuoarima-gaps-phi0.5.jpeg"),
     width = 3000, height = 2000,res = 300) # change name
ggplot(b, aes(x = gaps, col = gaps, y =value))+
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)  +
# scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0,1))+
  geom_hline(yintercept =0.3) +
  labs(y = "phi", x = "Gaps percentage(%)")+ 
  theme_bw()
dev.off()


# distribution of real data after restricted homogeneity
source(paste0(path_code_att,"sliding_variance.R"))

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
  if(nrow(na.omit(before)) > 31){
    norm.dat = one.step.norm(before, name.var = name.var, estimator = "Sca") 
    if(length(bef.norm)>30){
      hist.b = hist(bef.norm, breaks = seq(-30,30,0.05), plot = FALSE)
      dist.mean <- rbind(dist.mean, hist.b$counts/length(bef.norm))
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$bef<- bef.norm
    }
  }
  if(nrow(na.omit(after)) > 31){
    aft.norm = one.step.norm(after, name.var = name.var, estimator = "Sca") 
    if(length(aft.norm) > 30){
      hist.a = hist(aft.norm, breaks = seq(-30,30,0.05), plot = FALSE)
      dist.mean <- rbind(dist.mean, hist.a$counts/length(aft.norm))
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$aft <- aft.norm
    }
  }
}
save(dist.mean, file = paste0(path_results,"attribution/dist.mean_2years_", nearby_ver,".RData"))
save(data.all, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,".RData"))



res <- data.frame(t = seq(-29.95,30,0.05), count = colMeans(dist.mean))
res1 <- data.frame(t = seq(-29.95,30,0.1), count = e)

d = res[c(400:800),]

ggplot(d, aes(x = t, y = count))+geom_col() + theme_bw()+xlim(-10,10)
fitG = function(x,y,mu,sig,scale){
  
  f = function(p){
    d = p[3]*dnorm(x,mean=p[1],sd=p[2])
    sum((d-y)^2)
  }
  
  optim(c(mu,sig,scale),f)
}
p1 <- fitG(x = d$t, y = d$count, mu = 0, sig = 1, scale =1)
p = p1$par
plot(d$t, d$count)
lines(d$t,p[3]*dnorm(d$t,p[1],p[2]))

# MEAN VALUE INDICATE THE PURE GAUSSIAN --> LOOK INTO A SMALLER REGION WHERE NB OF POINT >0 IN OUTLIER PERIOD-----
a = colSums(dist.mean)

res <- data.frame(t = seq(-29.95,30,0.05), count = colMeans(dist.mean))

d <- tibble(Group=rep(1:2694, each=1200), Sample=as.vector(t(dist.mean)))

e = unlist(sapply(c(21), function(x) lengths(data.all[[x]])))
# d <- tibble(Group=rep(1:(length(e)), times=e), Sample=unlist(data.all))
d <- data.frame(Group=rep((21), times=e[21]), Sample=unlist(data.all[c(21)]))
qplot(sample=Sample, data=d, color=as.factor(Group))

plot( before$gps.era1, type = "l")
lines(bef.norm, col="red")
abline(h = -2.5)

d <- data.frame(Group= rep(c("raw","norm"), each = length(bef.norm)), Sample=c(before$gps.era1, bef.norm))
qplot(sample=Sample, data=d, color=as.factor(Group))+theme_bw()

d$Group <- as.factor(d$Group)
ggplot(d, aes(x = Sample, colour = Group, fill = Group)) + 
  geom_histogram(alpha = 0.5, position = "identity", bins = 200)+
  theme_bw()

# choose only europe stations ----------------------------------

list.pos = read.table(file = paste0(path_data_support,"gps_sta_CODE_REPRO_2015_OPER_combi.txt"),header = TRUE)
colnames(list.pos) <- c("name", "lat", "lon", "hei","alt")
list.f = list.pos[which(list.pos$lat>30 & list.pos$lat<50 & list.pos$lon >-10 & list.pos$lon <30),]
name.f = tolower(list.f$name)
list.name.eu = name_main[name_main%in%name.f]

name.var = list.test[3]
dist.mean <- data.frame(matrix(NA, nrow = 0, ncol = (length(seq(-30,30,0.05))-1)))
data.all <- list()
for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  if(station.ref %in% list.name.eu){
    station.near = substr(case.name ,start = 17, stop = 20)
    data.i = data.cr[[i]]
    breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
    before =  data.i[which(data.i$date <= breakpoint),]
    after =  data.i[which(data.i$date > breakpoint),]
    if(nrow(na.omit(before)) > 31){
      bef.norm = two.step.norm(Y = before, name.var)
      if(length(bef.norm)>30){
        # hist.b = hist(bef.norm, breaks = seq(-30,30,0.05), plot = FALSE)
        # dist.mean <- rbind(dist.mean, hist.b$counts)
        data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$bef<- bef.norm
      }
    }
    if(nrow(na.omit(after)) > 31){
      aft.norm = two.step.norm(Y = after, name.var)
      if(length(aft.norm) > 30){
        # hist.a = hist(aft.norm, breaks = seq(-30,30,0.05), plot = FALSE)
        # dist.mean <- rbind(dist.mean, hist.a$counts)
        data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$aft <- aft.norm
      }
    }
  }
}

e = unlist(sapply(c(217:220), function(x) lengths(data.all[[x]])))
# d <- tibble(Group=rep(1:(length(e)), times=e), Sample=unlist(data.all))
d <- data.frame(Group=rep((1:length(e)), times=e), Sample=unlist(data.all[217:220]))
qplot(sample=Sample, data=d, color=as.factor(Group))+theme_bw()


d$Group <- as.factor(d$Group)
ggplot(d, aes(x = Sample, colour = Group)) + 
  geom_density()+scale_x_continuous(breaks=seq(-6,7,1))+
  # geom_histogram(alpha = 0.5, position = "identity", bins = 200)+
  theme_bw()+  theme(legend.position="none")


coun = 0
for (i in c(217)) {
  case.name = names(data.all)[i]
  l = which( names(data.cr) == case.name)
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  data.i = data.cr[[l]]
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  before =  data.i[which(data.i$date <= breakpoint),]
  after =  data.i[which(data.i$date > breakpoint),]
  if(nrow(na.omit(before)) > 31){
    coun = coun+1
    bef.norm = two.step.norm(Y = before, name.var)
    data.b = data.frame(date = before$date, val = bef.norm)
    Y1 <- tidyr::complete(data.b, date = seq(min(data.b$date), max(data.b$date), by = "day"))
    p <- ggplot(data = Y1, aes(x = date, y= val))+
      geom_line()+theme_bw()+labs(subtitle = paste0("case", coun))
    jpeg(paste0(path_results,"attribution/",station.ref,".",as.character( breakpoint), ".", station.near, "b.jpeg"), width = 3000, height = 1800,res = 300) # change name
    print(p)
    dev.off()
    
  }
  if(nrow(na.omit(after)) > 31){
    coun = coun+1
    aft.norm = two.step.norm(Y = after, name.var)
    data.a = data.frame(date = after$date, val = aft.norm)
    Y1 <- tidyr::complete(data.a, date = seq(min(data.a$date), max(data.a$date), by = "day"))
    p <- ggplot(data = Y1, aes(x = date, y = val))+
      geom_line()+theme_bw()+labs(subtitle = paste0("case", coun))
    jpeg(paste0(path_results,"attribution/",station.ref,".",as.character( breakpoint), ".", station.near, "a.jpeg"), width = 3000, height = 1800,res = 300) # change name
    print(p)
    dev.off()
  }
}

data.all = get(load(file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,".RData")))
##################################################################
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
      bef.norm = sliding.median(before, name.var = "gps.gps") 
      bef.norm.all <- data.frame(l = bef.norm, date = before$date)
      bef.norm.all <- bef.norm.all[c(ind.sta:length(bef.norm)),]
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$bef <- bef.norm.all$l
    }
  }
  if(nrow(na.omit(after)) > 30){
    ind.end = which(after$date == min(breakpoint %m+% years(1), max(after$date)))
    if(length(ind.end) > 0){
      aft.norm = sliding.median(after, name.var = "gps.gps") 
      aft.norm.all <- data.frame(l = aft.norm, date = after$date)
      aft.norm.all <- aft.norm.all[c(ind.sta:length(aft.norm)),]
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$aft <- aft.norm.all$l
    }
  }
}

save(data.all, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"normalized1.RData"))



