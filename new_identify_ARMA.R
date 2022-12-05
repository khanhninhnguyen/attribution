# new function to identify the model of noise from the OLS residual and normalization
source(paste0(path_code_att, "newUsed_functions.R"))
source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))
library(nlme)
# Heteroskedastic ---------------------------------------------------------

# from OLS resdiual 
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365
last_signif <- function(signal, pq, alpha){  
  nb.or <- sum(pq)
  pq1 = rep(NA,3)
  while ( identical(as.numeric(pq1), pq) == FALSE) { # iteratively identify the model, stop when the model are the same after the significant check
    pq1 = pq
    if(nb.or==0){
      pandcoef <- list(p.value = rep(-1,4),coef = rep(0,4))
    }else{
      fitARIMA = try(arima( signal, pq, method="ML"), TRUE)
      if (class(fitARIMA) == "try-error"){
        fitARIMA = fit.b
      }
      pandcoef <- p.and.coef(fitARIMA, pq, nb.or)
    }
    pq = check_sig(p.val = pandcoef$p.value, alpha = alpha)
    nb.or <- sum(pq)
  }
  return(list( pq = pq, pandcoef = pandcoef))
}


# if we want to not dupplicate gps-era case
# all.cases.name = names(dat)
# all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
# unique.ind = match(unique(all.cases.ind), all.cases.ind )
# gps.era.dat = dat[unique.ind]
# data.test = gps.era.dat 
# n = length(data.test)
# sd.all <- rep(NA, n)


# OLS fit compute the normalized residual -----------------------------------------------------------------

one.year=365

compute.norm.resi <- function(dat.i, name.test){
  one.year=365
  scr = screen.O(Y = dat.i , name.var = name.test, method = 'sigma', global.mu = 0, iter = 1, estimator = "Sca", fix.thres = 0, loes = 0, loes.method = 0)
  sd0 = unlist(scr$sd.est[[1]])
  Xt = c(1:nrow(dat.i))- nrow(dat.i)/2
  
  Data.mod = data.frame( signal = dat.i[,name.test],Xt = Xt, complete.time = c(1:nrow(dat.i)))
  for (k in c(1:4)){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",k,"=cos(k*complete.time*(2*pi)/one.year),sin",k,"=sin(k*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  ols.fit = lm(signal~., data = Data.mod)
  ols.resi = rep(NA, nrow(dat.i))
  ols.resi[which(is.na(Data.mod$signal)==FALSE)] <- ols.fit$residuals 
  norm.res = ols.resi/sd0
  return(norm.res)
}

residus = list()
for (testi in c(1:6)) {
  name.test = list.test[testi]
  res.testi = list()
  for (i in c(1:length(dat))) {
    dat.all = dat[[i]]
    dat.bef = dat.all[c(1: (one.year*10)),]
    dat.aft = dat.all[-c(1: (one.year*10)),]
    res.bef = compute.norm.resi(dat.i = dat.bef, name.test)
    res.aft = compute.norm.resi(dat.i = dat.aft, name.test)
    res.testi[[names(dat)[i]]] <- c(res.bef, res.aft)
  }
  residus[[name.test]] <- res.testi
}

# fit the ARIMA 

fit.arima <- function(signal.test){
  fit.b = forecast::auto.arima(signal.test , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                               max.p = 2, max.q = 2, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
  
  pq <- arimaorder(fit.b)
  # order.init[k, c((testi*3-2): (testi*3))] <- pq
  options(warn = 1)
  
  refit0 = last_signif(signal = signal.test, pq, alpha = significant.level)
  pq = refit0$pq
  return(list(pq = pq, coef = refit0$pandcoef$coef, p = refit0$pandcoef$p.value))
}

order.arma <- list()
for (testi in c(1:6)) {
  name.test = list.test[testi]
  res.testi = residus[[name.test]]
  order.arma.bef = data.frame(matrix(NA, ncol = 3, nrow = length(dat)))
  order.arma.aft = data.frame(matrix(NA, ncol = 3, nrow = length(dat)))
  for (i in c(1:length(dat))) {
    dat.i = res.testi[[i]]
    bef.residus = dat.i[1:(one.year*10)]
    aft.residus = dat.i[-c(1:(one.year*10))]
    order.arma.bef[i,] = fit.arima(bef.residus)$pq
    order.arma.aft[i,] = fit.arima(aft.residus)$pq
  }
  order.arma[[name.test]] <- list(order.arma.bef, order.arma.aft)
}

