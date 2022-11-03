# this prog is used to investigate the multicolinearity problem 
# inspect by the autocorrelation between regressors

# inspect the variance inflation index 
# this prog used to apply the ancova/fgls on the real data 
library(tidyverse)   
library(attempt)
library(nlme)

source(paste0(path_code_att, "newUsed_functions.R"))
win.thres = 1
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
win.thres = 1
res = get(load(file = paste0(path_results, "attribution/all.hac.", win.thres, "years.RData")))

sig.com <- function(x, ver, vari.name){
  out = rep(NA, length(x))
  for (i in c(1:length(x))) {
    res.i = x[[i]][[ver]]$fit.hac
    ind = which(rownames(res.i) == vari.name)
    if(length(ind) == 0){
      out[i] = NA
    }else{
      out[i] = unlist(res.i$`Pr(>|z|)`[ind])
    }
  }
  return(out)
}
hac.full <- sig.com(res, ver = "full", vari.name = "Jump")
hac.sel <- sig.com(res, ver = "selec", vari.name = "Jump")
hac.sel[which(is.na(hac.sel)==TRUE)] <- 1
hac.full.x <- sig.com(res, ver = "full", vari.name = "Xt")
hac.sel.x <- sig.com(res, ver = "selec", vari.name = "Xt")

ind.plot = which(hac.sel.x<0.01)
ind.plot = ind.plot[-which(ind.plot %in% c(69, 124, 125, 138))]
for (j in ind.plot) {
  vif.f = round(car::vif( res[[j]][["full"]]$fit.ols)[1], digits = 1)
  # vif.s = round(car::vif( res[[j]][["selec"]]$fit.ols)[1], digits = 1)
  plot_HAC(case.name = names(res)[j], res.i = res[[j]][["full"]], data.in = dat[[names(res)[j]]], name.var = "gps.era", ver = "10yf", add.subtitle = vif.f)
  plot_HAC(case.name = names(res)[j], res.i = res[[j]][["selec"]], data.in = dat[[names(res)[j]]], name.var = "gps.era", ver = "10ys", add.subtitle = vif.s)
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


