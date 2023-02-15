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







# simulation test impact of sliding normalization -------------------------
t0 = as.Date("2003-01-01", format = "%Y-%m-%d")
raw = c()
norm = c()
for(l in c(1:1000)){
  set.seed(l)
  # x = arima.sim(model = list(ar = 0.3), n = 1000, sd = 1)
  x = rnorm(1000, 0 ,1)
  raw = c(raw, x)
  y = one.step.norm(Y = data.frame(date = seq(t0, t0+999,1), x = x), "x", estimator = "Sca", length.wind = 500 )
  norm = c(norm, y)
}

a = hist(x, 
         main = "Histogram of normalized data with trimmed mean", 
         xlab = "", 
         breaks = 500,
         xlim=c(-5, 5),
         prob = TRUE)
x2 <- a$breaks
fun2 <- dnorm(x2, mean = 0, sd = 1)

lines(x2, fun2, col = 2, lwd = 2)




# test loess function in estimate the mean 

tot <- data.frame(matrix(NA, ncol = 2, nrow = 1000))
for (l in c(1:1000)) {
  n0 = 60
  prob.outliers = 0.1
  times = c(1:n0)
  tot1 <- list()
  for (i in 1:nb.sim) {
    set.seed(i)
    x = SimulatedSeries(n = n0, P = 5, prob.outliers = 0.1, size.outliers = 3, rho = 0, theta = 0)$Y
    tot[i,1] <- median(x, na.rm = TRUE)
    a <- loess(x~times,degree=1, span = 60/n0, normalize = FALSE, na.action = na.exclude)
    tot[i,2] <- a$fitted[30]
  }
}
nb.sim=1000


set.seed(30)
N <- 300
combined <- data.frame(futureChange = c(rnorm(N, mean = 0, sd = 1), rnorm(N, mean = 0, sd = 5)),
                       direction = rep(c("long", "not long"), each = N))

lower.limit <- min(combined$futureChange)
upper.limit <- max(combined$futureChange)
long.density <- density(subset(combined, direction == "long")$futureChange, from = lower.limit, to = upper.limit, n = 2^10)
not.long.density <- density(subset(combined, direction == "not long")$futureChange, from = lower.limit, to = upper.limit, n = 2^10)

density.difference <- long.density$y - not.long.density$y
intersection.point <- long.density$x[which(diff(density.difference > 0) != 0) + 1]

ggplot(combined, aes(futureChange, fill = direction)) + geom_density(alpha = 0.2) + 
  geom_vline(xintercept = intersection.point, color = "red")

library(tidyverse)
load("data.all_1years_NGL.auck.2005-11-07.whng.RData")
head(Y)
dim(Y)

# Adding the missing data corresponding to missing date
YY <- Y
Y <- pad(Y)
na.Y <- which(is.na(Y$gps.era))

# We create the data adapted to the model: one year before the break and one year after
date.detected.break="2005-11-07"
rg.break <- which(Y$date==date.detected.break)
one.year=365
time.series <- seq(1-one.year/2,one.year/2,1)

Data.mod <- Y %>% dplyr::select(gps.era,date) %>%
  mutate(signal=gps.era) %>% dplyr::select(-gps.era) %>%
  slice((rg.break-one.year+1):(rg.break+one.year)) %>%
  mutate(Group=c(rep("Left",one.year),rep("Right",one.year))) %>%
  mutate(Xt=rep(time.series,2)) %>% dplyr::select(-date)

for (i in 1:4){
  eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=rep(cos(i*time.series*(2*pi)/one.year),2),sin",i,"=rep(sin(i*time.series*(2*pi)/one.year),2))")))
}
head(Data.mod)

#model
fit.signal=lm(signal~Group+Xt+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=Data.mod)
summary(fit.signal)
fit.signal.hac=lmtest::coeftest(fit.signal, sandwich::NeweyWest(fit.signal, lag = 1))[, ]
fit.signal.hac=as.data.frame(fit.signal.hac)
pval=fit.signal.hac$`Pr(>|t|)`

# Variable selection by hand
fit.signal=lm(signal~Group+Xt+sin2,data=Data.mod)
fit.signal.hac=lmtest::coeftest(fit.signal, sandwich::NeweyWest(fit.signal, lag = 1))[, ]
fit.signal.hac=as.data.frame(fit.signal.hac)
fit.signal.hac
a=which.max(fit.signal.hac[2:dim(fit.signal.hac)[1],]$`Pr(>|t|)`)
max(fit.signal.hac[2:dim(fit.signal.hac)[1],]$`Pr(>|t|)`)
rownames(fit.signal.hac)[a+1]

# final model-> sin2+constante+rupture (GroupR)+tendance (Xt)

signal.predicted <- rep(NA,length(Data.mod$signal))
signal.predicted[!((1:length(Data.mod$signal)) %in% na.Y)] <- fit.signal$fitted.values

Left <- 1:one.year
plot(Data.mod$Xt[Left],Data.mod$signal[Left],col="red", type="l",xlab="time",ylab="signal")
lines(Data.mod$Xt[Left],signal.predicted[Left],col="red")

Right <- (one.year+1):(2*one.year)
lines(Data.mod$Xt[Right],Data.mod$signal[Right],col="blue")
lines(Data.mod$Xt[Right],signal.predicted[Right],col="blue")

if(length(a)){
  print("1")
}else{ print ('2')}

vf1Fixed <- varFixed( ~ sin1+cos1 )
vf1 <- Initialize(vf1Fixed, data = dataMod)
varWeights( vf1)


x = as.matrix(data.frame(i = rep(1, 200), j = rep(c(0,1), each = 100), xt = c(1:200)))
x = as.matrix(data.frame(i = rep(1, 200), j = rep(c(0,1), each = 100)))

M =  t(x) %*% x
det(M)
solve(M)
a = rep(NA, 1000)
for (i in c(1:1000)) {
  set.seed(i)
  y = rnorm(200, 0,1)
  dat = data.frame(y = y, j = rep(c(0,1), eahc = 100), x = c(1:200))
  fit.i = lm(y~., data = dat)
  hac = sandwich::kernHAC(fit.i, prewhite = 1,approx = c("ARMA(1,1"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
  fit.hac=lmtest::coeftest(fit.i,df=fit.i$df.residual,vcov.= hac)[, ] %>% as.data.frame()
  
  a[i] = fit.hac$`Pr(>|t|)`[2]
}

res = data.frame(matrix(NA, ncol = 3, nrow = 1000))
res.bic = data.frame(matrix(NA, ncol = 4, nrow = 1000))
for (i in c(1:1000)) {
  set.seed(i)
  x = rnorm(1000, 0, 1)
  fit.b = forecast::auto.arima(x , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean = FALSE,lambda = NULL,
                               max.p = 2, max.q = 2, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
  arma1= Arima(x, order = c(1,0,1), include.mean = FALSE)
  whitenoise = Arima(x, order = c(0,0,0), include.mean = FALSE)
  # mod_capt <- capture.output(forecast::auto.arima(x , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean = FALSE,lambda = NULL,
  #                                                max.p = 2, max.q = 2, start.p = 0, trace = TRUE, allowdrift = FALSE,  approximation=FALSE))
  # ind.c = sapply(c(1:length(mod_capt)), function(x)  grepl( c("ARIMA"), mod_capt[x]))
  # ind.c = which(ind.c == TRUE)
  # ind.c =   ind.c[-c((length(ind.c)-1),length(ind.c))]
  # mod.list = mod_capt[ind.c]
  # a = sapply(c(1:length(mod.list)), function(x) as.numeric(strsplit(mod.list[x], split = ":")[[1]][2]))
  # b = sapply(c(1:length(mod.list)), function(x) strsplit(strsplit(mod.list[x], split = ":")[[1]][1], split = " ")[[1]][2])
  
  res.bic[i,] = c( whitenoise$bic, arma1$bic, whitenoise$loglik, arma1$loglik)
  res[i,] = arimaorder(fit.b)
}


a = sapply(c(1:1000), function(x) model.iden(as.numeric(unlist(res[x,]))))
a = read.series(path_series = path_series_nearby, station = "pama", na.rm = F, add.full = 0 )
b = read.series(path_series = path_series_main, station = "medi", na.rm = F, add.full = 0 )



# comparison bw 1 year vs 10 years ----------------------------------------

all.1 = get(load(file = paste0(path_results, "attribution/range_mean_var", win.thres=1,".RData")))
all.10 = get(load(file = paste0(path_results, "attribution/range_mean_var", win.thres=10,".RData")))

a = as.data.frame(all.1[[2]])
a$name = substr(rownames(a), start = 1, stop = 4)
b = as.data.frame(all.10[[2]])
b$name = substr(rownames(b), start = 1, stop = 4)
d = inner_join(a,b, by = c("name"))

res = data.frame(matrix(NA, ncol = 3, nrow = 1000))
res.coef = data.frame(matrix(NA, ncol = 2, nrow = 1000))
for (i in c(1:1000)) {
  set.seed(i)
  x = arima.sim(model = list(ar=0.6, ma = -0.4), n =1000, sd=1) 
  # set.seed(i+1000)
  # y= arima.sim(model = list(ar=0.44, ma = 0.22), n =1000, sd=1) 
  fit.b = forecast::auto.arima(x, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean = FALSE,lambda = NULL,
                               max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
  res[i,] = arimaorder(fit.b)
  coef = rep(NA,2)
  coef[which(fit.b$coef!=0)] <- fit.b$coef
  res.coef[i,] = coef
}


a = sapply(c(1:1000), function(x) model.iden(as.numeric(unlist(res[x,]))))
table(a)

a = get(load(file = paste0(path_results,"attribution/ver1/list.segments.selected", win.thres,".RData")))
a1 = na.omit(a)
a2 = get(load(file = paste0(path_results,"attribution/ver7/list.segments.selected", win.thres=10,".RData")))
a2 = na.omit(a2)

b = get(load(file = paste0(path_results,"attribution/ver3/six.models", win.thres,".RData")))
b1 = na.omit(b)
b2 = get(load(file = paste0(path_results,"attribution/ver2/six.models", win.thres,".RData")))
b2 = na.omit(b2)



d =  get(load(file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData")))

var.list = c(1:5)
res.tot = rep(NA, 5)
for (j in c(1:5)) {
  res.j = rep(NA, 1000)
  for(i in c(1:1000)){
    set.seed(i)
    x1 = rnorm(300, 0,var.list[j])
    set.seed(i+1000)
    x2 = rnorm(300, 0,var.list[j])
    res.i = t.test(x1, x2)
    res.j[i] = res.i$p.value
  }
  res.tot[j] = length(which(res.j<0.05))
}



