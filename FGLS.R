# FGLS 
library(nlme)
# Function are used 
# construct design matrix 
construct.design <- function(data.df, name.series, break.ind, one.year){
  Data.mod <- data.df %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(complete.time=1:nrow(data.df)) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  n0 = nrow(data.df)
  Data.mod$right = c(rep(0, break.ind), rep(1, (n0-break.ind)))
  Data.mod$left = rep(1, n0)
  
  return(Data.mod)
}
# construct the correlation matrix 
cor_matrix <- function(phi, theta, n.full, n.true){
  if(phi ==0 & theta ==0){
    y = rep(1, n.full)
  }else{
    y = ARMAacf(ar = phi, ma = theta, lag.max = n.true-1, pacf = FALSE)
  }
  return(toeplitz(y))
}

GLS <- function(phi, theta, var.t, design.matrix){
  # compute the variance covariance matrix 
  ind1 = which(is.na(design.matrix$signal)==FALSE)
  var.t.na = var.t[ind1]
  var.matrix = diag(sqrt(var.t.na))
  cor.matrix = cor_matrix(phi, theta, n.full = nrow(design.matrix), n.true = nrow(na.omit(design.matrix)))
  if(phi==0&theta==0){
    cov.var= diag(var.t.na)
  }else{
    cov.var0 = var.matrix  %*%  cor.matrix 
    cov.var = cov.var0 %*%  var.matrix
  }
  
  # estimate
  X = as.matrix(design.matrix%>% dplyr::select(-signal))
  X = X[ind1,]
  term2 = solve(cov.var)
  term1 = t(X) %*% term2 %*% X
  term3 = solve(term1)
  beta = term3 %*% t(X) %*% term2  %*% (as.matrix(design.matrix$signal[ind1]))
  var.beta = term3 
  
  # form the frame of result as in the ols 
  residual = design.matrix$signal[ind1] - X%*%beta
  t.val = beta/sqrt((diag(var.beta)))
  p.val = round(pnorm(-abs(t.val), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 4)
  fit.gls = data.frame(Estimate = beta)
  fitted.gls = X%*%beta
  fit.gls$`Std. Error` = sqrt(diag(var.beta))
  fit.gls$`t value` = t.val
  fit.gls$`Pr(>|t|)` = p.val
  
  return(list(Coefficients = beta, t.table = fit.gls, vcov = var.beta, residual = residual, fitted = fitted.gls))
}
WLS <- function(var.t, design.matrix){
  # compute the variance covariance matrix 
  cov.var = diag(sqrt(var.t))
  
  # estimate
  X = as.matrix(design.matrix%>% dplyr::select(-signal))
  
  term2 = solve(cov.var)
  term1 = t(X) %*% term2 %*% X
  term3 = solve(term1)
  beta = term3 %*% t(X) %*% term2  %*% (as.matrix(design.matrix$signal))
  var.beta = term3 
  
  # form the frame of result as in the ols 
  residual = design.matrix$signal - X%*%beta
  t.val = beta/sqrt((diag(var.beta)))
  p.val = round(pnorm(-abs(t.val), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 4)
  fit.gls = data.frame(Estimate = beta)
  fitted.gls = X%*%beta
  fit.gls$`Std. Error` = NA
  fit.gls$`t value` = t.val
  fit.gls$`Pr(>|t|)` = p.val
  
  return(list(Coefficients = beta, fit.wls = fit.gls, vcov = var.beta, residual = residual, fitted = fitted.gls))
}
chooseparam <- function(noise.model, arima.fit){
  if(noise.model[1]==1 & noise.model[3]==1){
    phi = arima.fit$coef[1]
    theta = arima.fit$coef[2]
  }else if(noise.model[1]==1 & noise.model[3]==0){
    phi = arima.fit$coef[1]
    theta = 0
  }else if(noise.model[1]==0 & noise.model[3]==1){
    phi = 0
    theta = arima.fit$coef[1]
  }else if(noise.model[1]==0 & noise.model[3]==0){
    phi = 0
    theta = 0
  }
  return(list(phi = phi, theta = theta))
}
RobEstiSlidingVariance.S <- function(Y, name.var, alpha, estimator, length.wind){# require date in the dataset, return std at time t
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = nrow(Y1)
  sigma.est1 = rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-(length.wind-1),1)
    end = min(n, i+length.wind)
    x.i = x[begin:end]
    x.idiff = (x.i)
    thre = 30
    if(i < 30|i>(n-30)){thre = 16}
    if(length(na.omit(x.idiff)) <= thre){
      sd <- NA
    }else{
      sd <- my.estimator(estimator = estimator, x.idiff)
    }
    sigma.est1[i] <- sd
  }
  # linear regression of the variance for gaps  MAYBE REPLACE BY INTERPOLATION FUNCTION
  s = sigma.est1
  if (sum(is.na(s)) != 0 & sum(is.na(s)) != length(s)){
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    if(na.ind[1] == 1){
      ind.stop = which(is.na(s)==FALSE)[1]-1
      na.ind <- na.ind[-c(1:ind.stop)]
    }else if (is.na(s[n]) == 1){
      m = which(is.na(s)==FALSE)
      ind.start = m[length(m)]
      na.ind <- na.ind[-which(na.ind %in% c(ind.start:n))]
    }
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  sigma.est = s[which(Y1$date %in% Y$date)]
  return(sigma.est^2)
}

FGLS <- function(design.m, tol, day.list, noise.model){
  resi0 = rep(NA, nrow(design.m))
  # call expression
  list.para <- colnames(design.m)[2:dim(design.m)[2]]
  mod.X <-  list.para %>% stringr::str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
  # ols
  ols.fit = lm( mod.expression, data = design.m)
  resi0[which(is.na(design.m$signal)==FALSE)] <- ols.fit$residuals
  old.coef = ols.fit$coefficients
  # estimate initial moving variance 
  Y0 = data.frame(date = day.list, residus = resi0)
  w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
  change1 = 10
  i=0
  while (change1 > tol) {
    # design.m$w0 = w0
    # estimate phi and theta 
    fit.wls = WLS(var.t = w0, design.matrix = design.m)
    normalized.res = fit.wls$residual/sqrt(w0)
    # noise.model = fit.arima(normalized.res )
    arima.fit = arima(x =normalized.res, order = noise.model, include.mean = FALSE)
    # # GLS
    coef.arma = chooseparam(noise.model = noise.model, arima.fit = arima.fit)
    fit.gls = GLS(phi = coef.arma$phi, theta = coef.arma$theta, var.t = w0, design.matrix = design.m)
    # use th gls function 
    # fit.gls <- eval(parse(text=paste0("gls(",mod.expression,",data=design.m,correlation =", cor.struct, "na.action=na.omit,weights=varFixed(value = ~w0)",")")))
    # update moving variance 
    fit.val = fit.gls$fitted
    resi0 = fit.gls$residual
    t.table = fit.gls$fit.gls
    Y0 = data.frame(date = day.list, residus = resi0)
    w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
    change1 = sum((fit.gls$Coefficients - old.coef)^2)
    print(old.coef)
    old.coef = fit.gls$Coefficients
    i=1+i
    print(coef.arma)
    print(change1)
  }
  return(list( coefficients = fit.gls$Coefficients, var = w0^2, residual = resi0, fit = fit.val, t.table=t.table, i=i))
}
my.estimator <- function(estimator,x){
  x1 = na.omit(x)
  if(estimator == "mad"){
    n0 = length(x1)
    f <- qnorm(3/4)*(1-0.7612/n0 - 1.123/(n0^2))
    y = mad(x1, constant = 1/f)
  }else if(estimator == "Qn"){
    y = robustbase::Qn(x1)
  }else if(estimator == "Sca"){
    y = robustbase::scaleTau2(x1, consistency = "finiteSample")
  }
  return(y)
}

p.and.coef <- function(fitARIMA, pq1, nb.or){
  test.sig = coeftest(fitARIMA)
  ord = pq1[c(1,3)]
  orde = c(rbind(ord,ord-1))
  orde[which(orde <0)] <- 0
  ind.param = which(orde >0)
  orde[ind.param] <- fitARIMA$coef[1:nb.or]
  p.value <- rep(-1, 4)
  p.value[ ind.param] <- test.sig[,4][1:nb.or]
  return(list(p.value = p.value , coef = orde))
}
# return significant order
check_sig <- function(p.val, alpha){
  ar.or = length(which(p.val[1:2] >= 0 & p.val[1:2] <= alpha))
  ma.or = length(which(p.val[3:4] >= 0 & p.val[3:4] <= alpha))
  return(c(ar.or, 0, ma.or))
}

fit.arima <- function(signal.test){
  fit.b = forecast::auto.arima(signal.test , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                               max.p = 2, max.q = 2, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
  
  pq <- forecast::arimaorder(fit.b)
  # order.init[k, c((testi*3-2): (testi*3))] <- pq
  # options(warn = 2)
  
  refit0 = last_signif(signal = signal.test, pq, alpha = significant.level, fit.b = fit.b)
  pq = refit0$pq
  
  if( any(pq > 1)){
    fit.b = forecast::auto.arima( signal.test, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                  max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
    pq = forecast::arimaorder(fit.b)
  }
  test.pq = pq
  refit1 = last_signif(signal = signal.test, pq, alpha = significant.level, fit.b = fit.b)
  
  pq = refit1$pq
  if(identical(as.numeric(test.pq),pq) == FALSE){print(c(test.pq,pq))}
  return(list(pq = pq, coef = refit1$pandcoef$coef, p = refit1$pandcoef$p.value))
}
last_signif <- function(signal, pq, alpha, fit.b){  
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

FGLS1 <- function(design.m, tol, day.list, noise.model, length.wind0){
  start_time <- Sys.time()
  resi0 = rep(NA, nrow(design.m))
  ind1 = which(is.na(design.m$signal)==FALSE)
  # call expression
  list.para <- colnames(design.m)[2:dim(design.m)[2]]
  mod.X <-  list.para %>% stringr::str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
  # ols
  ols.fit = lm( mod.expression, data = design.m)
  resi0[which(is.na(design.m$signal)==FALSE)] <- ols.fit$residuals
  old.coef = ols.fit$coefficients
  # estimate initial moving variance 
  Y0 = data.frame(date = day.list, residus = resi0)
  w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = length.wind0)
  change1 = 10
  change2 = 10
  change3 =10
  coef.old= list(phi = 0, theta =0)
  i=0
  # print(i)
  while (change3 > tol) {
    # design.m$w0 = w0
    # estimate phi and theta 
    # fit.wls = WLS(var.t = w0, design.matrix = design.m)
    design.m$signal[which(is.na(w0)==TRUE)] = NA
    design.mWLS = design.m
    design.mWLS$w0 = w0
    fit.wls = lm(mod.expression, data = design.mWLS, weights=w0)
    # print(i*10)
    resi0 <- design.m$signal -as.matrix(design.m[,-1]) %*% as.matrix(fit.wls$coefficients)
    normalized.res = resi0/sqrt(w0)
    # noise.model = fit.arima(normalized.res )
    arima.fit = arima(x =normalized.res, order = noise.model, include.mean = FALSE)
    # # GLS
    coef.arma = chooseparam(noise.model = noise.model, arima.fit = arima.fit)
    fit.gls = GLS(phi = coef.arma$phi, theta = coef.arma$theta, var.t = w0, design.matrix = design.m)
    # use th gls function 
    # fit.gls <- eval(parse(text=paste0("gls(",mod.expression,",data=design.m,correlation =", cor.struct, "na.action=na.omit,weights=varFixed(value = ~w0)",")")))
    # update moving variance 
    fit.val = fit.gls$fitted
    resi0 = rep(NA, nrow(design.m))
    ind1 = which(is.na(design.m$signal)==FALSE)
    resi0[ind1] = fit.gls$residual
    t.table = fit.gls$t.table
    Y0 = data.frame(date = day.list, residus = resi0)
    w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = length.wind0)
    normal.beta = (fit.gls$Coefficients - old.coef)/(sqrt(diag(fit.gls$vcov)))
    change2 = change1-sum((normal.beta)^2)
    change1 = sum((normal.beta)^2)
    # change3
    old.coef = fit.gls$Coefficients
    i=1+i
    j=0
    if(i>10){
      if(change2<(tol/100)){j=1}
      break
    }
    if(identical(noise.model,c(0,0,1))){
      change3 = abs(coef.arma$theta - coef.old$theta)
    }else{
      change3 = abs(coef.arma$phi - coef.old$phi)
    }
    coef.old = coef.arma                
  }
  print(i)
  end_time <- Sys.time()
  design.m$residual = resi0
  design.m$norm.res = normalized.res
  design.m$date = day.list
  return(list( coefficients = fit.gls$Coefficients, var = w0, residual = resi0, fit = fit.val, t.table=t.table, coef.arma = coef.arma,  
               i=i, j = j, change1= change1, all.out = fit.gls, t = (end_time - start_time), design.matrix = design.m))
}
cor.struct <- function(noise.model){
  if(identical(noise.model,c(0,0,0))){
    y = "NULL"
  }else if (identical(noise.model,c(1,0,0))){
    y = "corAR1(form = ~ 1)"
  }else if(identical(noise.model,c(0,0,1))){
    y = "corARMA(form = ~ 1, q=1)"
  }else if(identical(noise.model,c(1,0,1))){
    y = "corARMA(form = ~ 1, p =1, q=1)"
  }
}
cor.struct.hac <- function(noise.model){
  if(identical(noise.model,c(0,0,0))){
    y = "AR(1)"
  }else if (identical(noise.model,c(1,0,0))){
    y = "AR(1)"
  }else if(identical(noise.model,c(0,0,1))){
    y = "ARMA(1,1)"
  }else if(identical(noise.model,c(1,0,1))){
    y = "ARMA(1,1)"
  }
}

FGLS2 <- function(design.m, tol, day.list, noise.model){
  start_time <- Sys.time()
  ind1 = which(is.na(design.m$signal)==FALSE)
  resi0 = rep(NA, nrow(design.m))
  # call expression
  list.para <- colnames(design.m)[2:dim(design.m)[2]]
  mod.X <-  list.para %>% stringr::str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
  # ols
  ols.fit = lm( mod.expression, data = design.m)
  resi0[which(is.na(design.m$signal)==FALSE)] <- ols.fit$residuals
  old.coef = ols.fit$coefficients
  n = length(ols.fit$residuals)
  # estimate initial moving variance 
  Y0 = data.frame(date = day.list, residus = resi0)
  w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
  change1 = 10
  change2 = 10
  cor.structure = cor.struct(noise.model)
  i=0
  while (change1 > tol) {
    design.m$w0 = w0
    # estimate phi and theta 
    # fit.wls = WLS(var.t = w0, design.matrix = design.m)
    # normalized.res = fit.wls$residual/sqrt(w0)
    # # noise.model = fit.arima(normalized.res )
    # arima.fit = arima(x =normalized.res, order = noise.model, include.mean = FALSE)
    # # # GLS
    # coef.arma = chooseparam(noise.model = noise.model, arima.fit = arima.fit)
    # fit.gls = GLS(phi = coef.arma$phi, theta = coef.arma$theta, var.t = w0, design.matrix = design.m)
    # use th gls function 
    fit.gls <- eval(parse(text=paste0("gls(",mod.expression,",data=design.m,correlation =", cor.structure, ",na.action=na.omit,weights=varFixed(value = ~w0)",")")))
    # update moving variance 
    fit.val = fit.gls$fitted
    resi0 = rep(NA, nrow(design.m))
    resi0[ind1] = fit.gls$residuals
    t.table =  lmtest::coeftest(fit.gls,df=(n-10))[, ] %>% as.data.frame()
    Y0 = data.frame(date = day.list, residus = resi0)
    w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
    normal.beta = (fit.gls$coefficients - old.coef)/(sqrt(diag(fit.gls$varBeta)))
    change2 = change1-sum((normal.beta)^2)
    change1 = sum((normal.beta)^2)
    # print(old.coef)
    old.coef = fit.gls$Coefficients
    i=1+i
    j=0
    if(i>10){
      if(change2<(tol/100)){j=1}
      break
    }
    # print(coef.arma)
    # print(change1)
  }
  print(i)
  end_time <- Sys.time()
  return(list( coefficients = fit.gls$Coefficients, var = w0, residual = resi0, fit = fit.val, all.out = fit.gls,
               t.table=t.table,i=i, j =j, change1= change1, phi = fit.gls$modelStruct$corStruct, t =  (end_time - start_time)))
}
FGLS3 <- function(design.m, tol, day.list, noise.model){
  resi0 = rep(NA, nrow(design.m))
  # call expression
  list.para <- colnames(design.m)[2:dim(design.m)[2]]
  mod.X <-  list.para %>% stringr::str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
  # ols
  ols.fit = lm( mod.expression, data = design.m)
  resi0[which(is.na(design.m$signal)==FALSE)] <- ols.fit$residuals
  old.coef = ols.fit$coefficients
  # estimate initial moving variance 
  Y0 = data.frame(date = day.list, residus = resi0)
  w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
  change1 = 10
  change2 = 10
  change3 =10
  coef.old= list(phi = 0, theta =0)
  i=0
  while (change1 > tol) {
    # design.m$w0 = w0
    # estimate phi and theta 
    fit.wls = WLS(var.t = w0, design.matrix = design.m)
    normalized.res = fit.wls$residual/sqrt(w0)
    # noise.model = fit.arima(normalized.res )
    arima.fit = arima(x =normalized.res, order = noise.model, include.mean = FALSE)
    # # GLS
    coef.arma = chooseparam(noise.model = noise.model, arima.fit = arima.fit)
    fit.gls = GLS(phi = coef.arma$phi, theta = coef.arma$theta, var.t = w0, design.matrix = design.m)
    # use th gls function 
    # fit.gls <- eval(parse(text=paste0("gls(",mod.expression,",data=design.m,correlation =", cor.struct, "na.action=na.omit,weights=varFixed(value = ~w0)",")")))
    # update moving variance 
    fit.val = fit.gls$fitted
    resi0 = fit.gls$residual
    t.table = fit.gls$t.table
    Y0 = data.frame(date = day.list, residus = resi0)
    w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
    normal.beta = (fit.gls$Coefficients - old.coef)/(sqrt(diag(fit.gls$vcov)))
    change2 = change1-sum((normal.beta)^2)
    change1 = sum((normal.beta)^2)
    # change3
    old.coef = fit.gls$Coefficients
    i=1+i
    j=0
    if(i>10){
      if(change2<(tol/100)){j=1}
      break
    }
    change3 = abs(coef.arma$phi - coef.old$phi)
    coef.old = coef.arma                
  }
  print(i)
  return(list( coefficients = fit.gls$Coefficients, var = w0, residual = resi0, fit = fit.val, t.table=t.table,i=i, j = j, change1= change1))
}
remove_na_2sides <- function(df, name.series){
  a = which(is.na(df[name.series])== FALSE)
  df = df[c(min(a):(max(a))), ]
  return(df)
}
model.iden <- function(order){
  model = c()
  if (identical(order, c(1,0,1))){ model = "ARMA(1,1)"}
  else if (identical(order, c(1,0,0))){ model = "AR(1)"}
  else if (identical(order, c(0,0,1))){ model = "MA(1)"}
  else if (identical(order, c(0,0,0))){ model = "White"}
  else if (identical(order, c(2,0,0))){ model = "AR(2)"}
  else if (identical(order, c(2,0,1))){ model = "ARMA(2,1)"}
  else if (identical(order, c(1,0,2))){ model = "ARMA(1,2)"}
  else if (identical(order, c(0,0,2))){ model = "MA(2)"}
  else if (identical(order, c(2,0,2))){ model = "ARMA(2,2)"}
  
  return(model)
}


