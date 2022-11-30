# this prog is used to investigate the multicolinearity problem 
# inspect by the autocorrelation between regressors

# inspect the variance inflation index 
# this prog used to apply the ancova/fgls on the real data 
library(tidyverse)   
library(attempt)
library(nlme)

source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"newUsed_functions.R"))
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365

all.cases.name = names(dat)
all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
unique.ind = match(unique(all.cases.ind), all.cases.ind )
gps.era.dat = dat[unique.ind]

data.test = gps.era.dat
list.ind = c(1:length(data.test))
tot.res <- data.frame(matrix(NA, ncol = 4, nrow = length(list.ind)))
Res <- list()
for (k in list.ind) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")

  # Contruction of the dataset 
  Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>% 
    mutate(complete.time=1:(2*one.year*win.thres)) %>% 
    mutate(Xt=complete.time-one.year*win.thres/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  res.hac.1step <- Test_OLS_vcovhac_1step(Data.mod)
  res.hac <- Test_OLS_vcovhac(Data.mod)
  tot.res[k,4] = car::vif(res.hac.1step$fit.ols)[1]
  ind.jump = which(rownames(res.hac$fit.hac) == "Jump")
  ind.Xt = which(rownames(res.hac$fit.hac) == "Xt")
  
  if(length(rownames(res.hac$fit.hac)) > 2){
    if( length(ind.jump) != 0){
      if(length(ind.Xt) != 0){
        ind.j = which(colnames(res.hac$vcov.para) == "Jump")
        ind.i = which(rownames(res.hac$vcov.para) == "Xt") 
        tot.res[k,1] =  res.hac$vcov.para[ind.i,ind.j]/(sqrt(res.hac$vcov.para[ind.i,ind.i]*res.hac$vcov.para[ind.j,ind.j]))
      }
      tot.res[k,3] <- car::vif(res.hac$fit.ols)[ind.jump]
    }else{tot.res[k,3] =  NA }
  }else{
    tot.res[k,1] =  -2 
    tot.res[k,3] =  -1 
    }
  
  # res.hac <- Test_OLS_vcovhac(Data.mod)
  ########################
  # compute the correlation between regressors: jump and trend 
  tot.res[k,2] = res.hac.1step$vcov.para[3,2]/(sqrt(res.hac.1step$vcov.para[2,2]*res.hac.1step$vcov.para[3,3]))
  # tot.res[[name.dataset]] <- res.hac
  Res[[name.dataset]] <- list(full = res.hac.1step, selec = res.hac)
  print(k)
}
save(Res, file = paste0(path_results, "attribution/all.hac.", win.thres, "years.RData"))

colnames(tot.res) <- c("corr.selected", "corr.full","VIF.selected", "VIF.full")
# tot.res[,1] <- if_else(is.na(tot.res[,1]), -2, tot.res[,1])
hist(tot.res[,1] , breaks =50, main = "Histogram of the covariance after variable selection ", xlab = "")
hist(tot.res[,2] , breaks =50, main = "Histogram of the covariance before variable selection ", xlab = "")
hist(tot.res[,3] , breaks =50, main = "Histogram of the VIF after variable selection", xlab = "")
hist(tot.res[,4] , breaks =50, main = "Histogram of the VIF before variable selection", xlab = "")
save(tot.res, file = paste0(path_results, "attribution/multicolinear.", win.thres, "years.RData"))


# analyze results ---------------------------------------------------------
# plot individual cases
win.thres = 10
res = get(load(file = paste0(path_results, "attribution/all.hac.", win.thres, "years.RData")))

sig.com <- function(x, ver, vari.name, feature){
  out = rep(NA, length(x))
  for (i in c(1:length(x))) {
    res.i = x[[i]][[ver]]$fit.hac
    ind = which(rownames(res.i) == vari.name)
    if(length(ind) == 0){
      out[i] = NA
    }else{
      out[i] = round(unlist(res.i[ind, feature]), digits = 5)
    }
  }
  return(out)
}
hac.full <- sig.com(res, ver = "full", vari.name = "Jump", feature = "Pr(>|z|)")
hac.sel <- sig.com(res, ver = "selec", vari.name = "Jump", feature = "Pr(>|z|)")
hac.sel[which(is.na(hac.sel)==TRUE)] <- 1
hac.full.x <- sig.com(res, ver = "full", vari.name = "Xt", feature = "Pr(>|z|)")
hac.sel.x <- sig.com(res, ver = "selec", vari.name = "Xt", feature = "Pr(>|z|)")

ind.plot = which(hac.sel.x<0.01)
ind.plot = ind.plot[-which(ind.plot %in% c(69, 124, 125, 138))]
for (j in ind.plot) {
  vif.f = round(car::vif( res[[j]][["full"]]$fit.ols)[1], digits = 1)
  # vif.s = round(car::vif( res[[j]][["selec"]]$fit.ols)[1], digits = 1)
  plot_HAC(case.name = names(res)[j], res.i = res[[j]][["full"]], data.in = dat[[names(res)[j]]], name.var = "gps.era", ver = "10yf", add.subtitle = vif.f)
  plot_HAC(case.name = names(res)[j], res.i = res[[j]][["selec"]], data.in = dat[[names(res)[j]]], name.var = "gps.era", ver = "10ys", add.subtitle = "")
}

a = names(data.test)[ind.plot]
# histogram of VIF

tot.res = get(load( file = paste0(path_results, "attribution/multicolinear.", win.thres=1, "years.RData")))
tot.res$corr.selected[which(is.na(tot.res$corr.selected) == TRUE)]<- -2
tot.res[which(tot.res$VIF.selected>3),]

a = names(data.test)[which(tot.res$VIF.selected>3)]
tot.res[which(names(data.test) %in% a),]


# check the VIF of jump and trend -----------------------------------------

Res <- data.frame(matrix(NA, nrow = length(list.ind), ncol = 10))
for (k in list.ind) {
  name.dataset = names(data.test)[k]
  Y.with.NA = data.test[[k]]
  date.detected.break = as.Date(substr(name.dataset,start = 6, stop = 15) , format = "%Y-%m-%d")
  
  # Contruction of the dataset 
  Data.mod <- Y.with.NA %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>% 
    mutate(complete.time=1:(2*one.year*win.thres)) %>% 
    mutate(Xt=complete.time-one.year*win.thres/2) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  res.hac.1step <- Test_OLS_vcovhac_1step(Data.mod)
  Res[k,] = car::vif( res.hac.1step$fit.ols)
}
tot.res$VIF.selected[which(is.na(tot.res$VIF.selected)==TRUE)] <- -1


# check the theoretical VIF 
vif.HAC <- function(mod, vcov.matrix, ...) {
  if (any(is.na(coef(mod)))) 
    stop ("there are aliased coefficients in the model")
  v <- vcov.matrix
  assign <- attr(model.matrix(mod), "assign")
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  }
  else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("model contains fewer than 2 terms")
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) *
      det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) result <- result[, 1]
  else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  return(list(VIF = result, R = R))
}

vif.HAC1 <- function(mod, vcov.matrix, ...) {
  if (any(is.na(coef(mod)))) 
    stop ("there are aliased coefficients in the model")
  v <- vcov.matrix
  assign <- attr(model.matrix(mod), "assign")
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  }
  else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("model contains fewer than 2 terms")
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign != term)
    for (sub in subs) {
      result[sub, 1] <- det(as.matrix(R[1, 1])) *
        det(as.matrix(R[sub, sub])) / detR
      result[sub, 2] <- length(sub)
    }
  }
  if (all(result[, 2] == 1)) result <- result[, 1]
  else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  return(list(VIF = result, R = R))
}
one.year = 365
win.thres = 10
y = rnorm(n = 730, 0, 1)
Data.mod <- data.frame( signal = rep(1, one.year*win.thres*2)) %>%
  mutate(Jump=c(rep(0,one.year*win.thres),rep(1,one.year*win.thres))) %>% 
  mutate(complete.time=1:(2*one.year*win.thres)) %>% 
  mutate(Xt=complete.time-one.year*win.thres/2)
for (i in 1:4){
  eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
}
Data.mod <- Data.mod %>% dplyr::select(-complete.time)

a = as.matrix(Data.mod)
cov.m = solve(t(a)%*%a) # OLS case
R = cov2cor(cov.m)
R.inv = solve(R[-1,-1])
diag(R.inv)

r= rep(NA, length(res))
for ( i in c(1:length(res))){ 
  a = res[[i]]$full$fit.hac$Estimate
  r[i] = a[2]*a[3]
}



# check probability of confusion situation --------------------------------
source(paste0(path_code_att, "simulate_time_series.R"))
nb.sim = 10000
n = 2860 # mean length of long time series 
t = c(1:n) - n/4
res <- data.frame(matrix(NA, ncol = 3, nrow = nb.sim))
res.coef <- data.frame(matrix(NA, ncol = 3, nrow = nb.sim))
trend.0 = 0.0001
mu0 = -0.36
ar = 0
sd0 = 0.7

tot.res <- data.frame(matrix(NA, ncol= 4, nrow = 5))
for (k in c(1:6)) {
  for (j in c(1:nb.sim)) {
    set.seed(j)
    noise.j = simulate.general(burn.in = 10000,
                               arma.model = c(ar,0),
                               hetero = 0,
                               monthly.var = 0,
                               sigma = sd0,
                               N = n,
                               gaps = 0,
                               outlier = 0,
                               prob.outliers = 0,
                               size.outliers = 0)
    signal = t*trend.0 + noise.j
    # signal = rnorm(n, 0, 1)
    signal[(n/2+1):n] = signal[(n/2+1):n] + mu0[k]
    Data.mod <- data.frame(y = signal, trend = t, mu = rep(c(0,1), each = (n/2)))
    lr <- lm(y~., data = Data.mod)
    summa = coeftest(lr)[, ] %>% as.data.frame()
    res[j,] <- summa$`Pr(>|t|)`
    res.coef[j,] <- lr$coefficients
  }
  tot.res[k,1] = length(which(res$X3 <0.05 & res$X2>0.05))
  tot.res[k,2] = length(which(res$X3 >0.05 & res$X2>0.05))
  tot.res[k,3] = length(which(res$X3 <0.05 & res$X2<0.05))
  tot.res[k,4] = length(which(res$X3 >0.05 & res$X2<0.05))
}
tot.res$jump.val = mu0
a = reshape2::melt(tot.res, id = "jump.val")
ggplot(data = a, aes(x = jump.val, y = value/nb.sim, col=variable))+geom_point()+theme_bw()+ylab("Proportion of significant estimates")+xlab("Jump")

a = which(res$X3 <0.05 & res$X2<0.05)
b = res.coef[a,]
table(b$X2*b$X3<0)


# check the confusion in the real data ------------------------------------
ver = "full"
jump.est <- sig.com(res, ver = ver, vari.name = "Jump", feature = "Estimate")
trend.est <- sig.com(res, ver = ver, vari.name = "Xt", feature = "Estimate")
jump.p <- sig.com(res, ver = ver, vari.name = "Jump", feature = "Pr(>|z|)")
trend.p <- sig.com(res, ver = ver, vari.name = "Xt", feature = "Pr(>|z|)")
length.data = sapply(c(1:length(data.test)), function(x) length(na.omit(data.test[[x]]$gps.era)))
length.data[which(jump.p < 0.05 & trend.p<0.05)]
length.data[which(jump.p < 0.05 & trend.p>0.05)]
length.data[which(jump.p > 0.05 & trend.p<0.05)]
length.data[which(jump.p > 0.05 & trend.p>0.05)]

length(which(jump.p < 0.05 & trend.p<0.05))
length(which(jump.p < 0.05 & trend.p>0.05))
length(which(jump.p > 0.05 & trend.p<0.05))
length(which(jump.p > 0.05 & trend.p>0.05))


# inspect the impact of multicollinearity ---------------------------------

sim.collinear <- function(nb.sim, off.set, trend, n, ar, burn.in, hetero, sigma.sim){
  coef.all = data.frame(matrix(NA, ncol = 3, nrow = nb.sim))
  var.all = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
  t = c((-n/2):(n/2-1))
  for (i in c(1:nb.sim)) {
    set.seed(i)
    y = simulate.general(N = n, arma.model = c(ar,0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim),
                         monthly.var = 0)
    y[(n/2):n] = y[(n/2):n]+off.set
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = sigma.sim, t = t)
    
    ols.fit = lm(signal~jump, data = Data.mod)
    
    ols.fit.t = lm(signal~jump+t, data = Data.mod)
    
    coef.all[i,] = c(ols.fit$coefficients[2],ols.fit.t$coefficients[2:3] )
    var.all[i,] = c(vcov(ols.fit)[2,2],vcov(ols.fit.t )[2,2])
  }
  return(list(coef = coef.all, var = var.all))
}
L = seq(100, 6000, 500)
off.set=0
trend = 0
ar=0
hetero=0
sigma.sim=0.4
res = list()
for (j in c(1:length(L))) {
  Lj= L[j]
  a = sim.collinear(nb.sim=1000, off.set=off.set, trend=trend, n=Lj, ar=ar, burn.in=0, hetero=hetero, sigma.sim=sigma.sim)
  res[[j]] = a
}
save(res, file = paste0(path_results, "attribution/sim.collinear.o",off.set,"t",trend,"ar",ar,"h",hetero,".RData"))
res = get(load(file = paste0(path_results, "attribution/sim.collinear.o",off.set,"t",trend,"ar",ar,"h",hetero,".RData")))

wt = sapply(c(1:length(res)), function(x) mean(res[[x]]$coef$X1))
w = sapply(c(1:length(res)), function(x) mean(res[[x]]$coef$X2))

wt = sapply(c(1:length(res)), function(x) sd(res[[x]]$coef$X1))
w = sapply(c(1:length(res)), function(x) sd(res[[x]]$coef$X2))

wt = sapply(c(1:length(res)), function(x) mean(res[[x]]$var$X1))
w = sapply(c(1:length(res)), function(x) mean(res[[x]]$var$X2))


###CHECK WHY THE vcov RETURN DIFFERENT VALUE FOR VAR(BETA) IN THE OLS CASE???/
# individual case-----
nb.sim = 1000
off.set = 0.3
trend = 0
n=2000
T1 = n/2
a = cos(2*pi*(c(1:n)/T1))
var.m = 1
var.t = var.m - 0.9*a
# var.all = seq(0, 0.5, 0.1)
ar = 0.3
coef.all = data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
var.all = data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
t = c((-n/2):(n/2-1))
t1 = c(1:n)
for (i in c(1:nb.sim)) {
  set.seed(i)
  y =  simulate.general(N = n, arma.model = c(0,0), burn.in = 0, hetero = 0, sigma = sqrt(1),
                            monthly.var = 0)
  y[(n/2):n] = y[(n/2):n]+off.set
  Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), t1=t1, t=t)
  for (j in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",j,"=cos(j*t1*(2*pi)/T1),sin",j,"=sin(j*t1*(2*pi)/T1))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-t1)
  Data.mod <- Data.mod[-1,]
  ols.fit = lm(signal~t, data = Data.mod)
  
  ols.fit.t = lm(signal~., data = Data.mod)
  ols.fit.cos = lm(signal~t+cos1+cos2+cos3+cos4, data = Data.mod)
  ols.fit.sin = lm(signal~t+sin1+sin2+sin3+sin4, data = Data.mod)
  
  coef.all[i,] = c(ols.fit$coefficients[2],ols.fit.t$coefficients[2], ols.fit.cos$coefficients[2], ols.fit.sin$coefficients[2])
  var.all[i,] = c(vcov(ols.fit)[2,2],vcov(ols.fit.t )[2,2], vcov(ols.fit.cos )[2,2], vcov(ols.fit.sin )[2,2])
  # var.all[i,] = c(vcov.para[2,2],vcov.para.t[2,2])
}
# plot results ----- 
# histo of the variance
start.c = round(min(var.all), digits = 3) - 0.001
end.c = round(max(var.all), digits = 3) +0.001
hist(var.all$X1, col=rgb(0,0,1,0.2), seq(start.c, end.c, 0.0001), xaxt='n',
     xlim = c(start.c, end.c), main = "Histogram of var(jump)", xlab = "")
hist(var.all$X3, col=rgb(1,0,0,0.2), seq(start.c, end.c, 0.0001), add=TRUE, xlab = "")
legend('topright', c('jump', 'jump+trend'),
       fill=c(rgb(0,0,1,0.2), rgb(1,0,0,0.2)))
axis(side=1, at=seq(start.c, end.c, 0.001), labels=seq(start.c, end.c, 0.001))

# histo of the estimates
start.c = round(min(coef.all), digits = 2) - 0.01
end.c = round(max(coef.all), digits = 2) +0.01
hist(coef.all$X1, col=rgb(0,0,1,0.2), breaks = seq(start.c, end.c, 0.01), xlim = c(start.c, end.c), main = "Histogram of jump estimates", xlab = "")
hist(coef.all$X2, col=rgb(1,0,0,0.2), breaks = seq(start.c, end.c, 0.01), add=TRUE, xlab = "")
legend('topright', c('jump', 'jump+trend'),
       fill=c(rgb(0,0,1,0.2), rgb(1,0,0,0.2)))

hist(coef.all$X3, breaks = 100)


round(summary(coef.all$X1), digits = 5)
round(summary(coef.all$X2), digits = 5)
round(summary(var.all$X1), digits = 5)
round(summary(var.all$X2), digits = 5)
# the case wihout  the jump, check the multicollinearity between the trend and sin + cos 

# without jump ----------
for (i in c(1:nb.sim)) {
  set.seed(i)
  y =  simulate.general(N = n, arma.model = c(0,0), burn.in = 0, hetero = 0, sigma = sqrt(1),
                        monthly.var = 0)
  Data.mod = data.frame(signal = y, t1=t1, t=t)
  for (j in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",j,"=cos(j*t1*(2*pi)/T1),sin",j,"=sin(j*t1*(2*pi)/T1))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-t1)
  ols.fit.t = lm(signal~t, data = Data.mod)
  
  ols.fit = lm(signal~., data = Data.mod)
  
  coef.all[i,] = c(ols.fit$coefficients[2],ols.fit.t$coefficients[2:3] )
  var.all[i,] = c(vcov(ols.fit)[2,2],vcov(ols.fit.t )[2,2])
}
start.c = round(min(var.all), digits = 6) - 0.0000001
end.c = round(max(var.all), digits = 6) +0.0000001
hist(var.all$X1, col=rgb(0,0,1,0.2), seq(start.c, end.c, 0.00000001), xaxt='n',
     xlim = c(start.c, end.c), main = "Histogram of var(jump)", xlab = "")
hist(var.all$X2, col=rgb(1,0,0,0.2), seq(start.c, end.c, 0.00000001), add=TRUE, xlab = "")
legend('topright', c('jump', 'jump+trend'),
       fill=c(rgb(0,0,1,0.2), rgb(1,0,0,0.2)))
axis(side=1, at=seq(start.c, end.c, 0.0000001), labels=seq(start.c, end.c,  0.0000001))

# histo of the estimates
start.c = round(min(coef.all), digits = 2) - 0.01
end.c = round(max(coef.all), digits = 2) +0.01
hist(coef.all$X1, col=rgb(0,0,1,0.2), breaks = seq(start.c, end.c, 0.01), xlim = c(start.c, end.c), main = "Histogram of jump estimates", xlab = "")
hist(coef.all$X2, col=rgb(1,0,0,0.2), breaks = seq(start.c, end.c, 0.01), add=TRUE, xlab = "")
legend('topright', c('jump', 'jump+trend'),
       fill=c(rgb(0,0,1,0.2), rgb(1,0,0,0.2)))

a = data.frame(rep(1,(n-1)), Data.mod[,c(3,5)])
a = as.matrix(a)
t(a)%*%a
solve(t(a)%*%a)

