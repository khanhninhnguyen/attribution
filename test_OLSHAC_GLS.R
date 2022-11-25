# this prog is used for comparison: OLS+HAC vs GLS
source(paste0(path_code_att,"simulate_time_series.R"))
library(sandwich)
# gls function new version
var_cova <- function(Sig, phi, theta, n){
  M = 9*n
  # compute the variance of et, 
  Sig.k = Sig[2:(M-n)] 
  Sig.k1 = Sig[1:(M-n-1)] 
  k = c(0:(M-n-2))
  
  # to compute the variance 
  c1 = (phi^(2*k))
  c2 = c1*(theta^2)
  c3 = (phi^(2*k+1))*theta
  c4 = (phi^(2*k-1))*theta
  
  term1 = rev(c2) * Sig.k1
  term2 = rev(c1) * Sig.k
  term3 = rev(c3) * Sig.k1
  term4 = rev(c4) * Sig.k
  # variance of point t0-1
  e = term4[length(term4)]
  et = sum(term1+ term2+ term3 + term4) - e
  
  var.t = rep(NA, n)
  for (i in c(1:n)) {
    s = Sig[(M-n+i)]
    s1 = Sig[(M-n+i-1)]
    err = (theta/phi)*s1
    c5 = (s1)*(theta**2) + s + theta*phi*(s1)
    ei = c5 + (phi**2)*(et + s1*theta/phi)
    var.t[i] = ei
    et = ei
  }
  
  # to commpute the covariance 
  var.cov = matrix(NA, ncol = n, nrow = n)
  for (l in c(1:(n-1))) {
    var.et = var.t[l]
    sigma.et = Sig[(M-n+l)]
    for (l1 in c((l+1):(n))) {
      var.cov[l, l1] = (phi^(l1-l))*var.et + (phi^(l1-l-1))*theta*sigma.et
    }
  }
  var.cov[is.na(var.cov)] = 0
  var.cov = var.cov + t(var.cov)
  diag(var.cov) = var.t
  
  return(var.cov)
}
gls.true <- function(var.t, phi, theta, design.matrix, trend){
  if(phi ==0 & theta ==0){
    cov.var = diag(var.t)
  }else{
    Sig = rep(var.t,9)
    cov.var = var_cova(Sig = Sig, phi = phi, theta = theta, n = nrow(design.matrix))
  }
  X = data.frame(intercept = rep(1,n))
  if(trend ==0){
    X = cbind(X, design.matrix[c("jump")])
  }else{
    X = cbind(X,design.matrix[c("jump", "Xt")])
  }
  X = as.matrix(X)
  term1 = t(X) %*% (solve(cov.var)) %*% X
  beta = solve(term1) %*% t(X) %*% (solve(cov.var)) %*% (as.matrix(design.matrix$signal))
  var.beta = solve(term1)
  return(list(beta = beta, var.beta = var.beta))
}
# initial condition ---------------------------------------------------------
nb.sim = 1
n = 200
sig.m = 1
sig.v = 0.9
T1 = n/2
a = cos(2*pi*(c(1:n)/T1)-pi)
var.t = (sig.m - sig.v*a)
res = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
res.var = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
t = c(1:n)-n/2
# off.set = c(0.1, 0.2, 0.3, 0.4, 0.5)

# comparison with heteroskedasticity --------------------------------------
#### IID case with the estimated variance in GLS
offset = seq(0, 0.5, 0.1)
param.test = offset
Res.fin = data.frame(matrix(NA, ncol = 3, nrow = length(param.test)))
for (l in c(1:length(param.test))) {
  off.set = param.test[l]
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  # plot(var.t)
  for (i in c(1:nb.sim)) {
    set.seed(i)
    y = simulate.general(N = n, arma.model = c(0,0), burn.in =0, hetero = 0, sigma = sqrt(sig.m),
                                   monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + off.set
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = var.t, t = t)
    
    # Test with vcov (HAC)
    ols.fit = lm(signal~jump, data = Data.mod)
    vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE,kernel = "Quadratic Spectral", sandwich = TRUE)
    # vcov.para=sandwich::vcovHC(ols.fit, type = "HC3")
    fit.hac=lmtest::coeftest(ols.fit,df=(n-1),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.ols=lmtest::coeftest(ols.fit,df=(n-1))[, ] %>% as.data.frame()

    # Test with gls
    # et = ols.fit$residuals^2
    # Data.mod.e = data.frame(signal = et, sin1 = sin(2*pi*t/T1), cos1 = cos(2*pi*t/T1))
    # lm.et = lm(signal~., data = Data.mod.e)
    # var.e = fitted(lm.et)
    # gls.fit = gls(signal~jump, data = Data.mod, correlation = corAR1( form = ~t), weights = varFixed(value = ~var.e ))
    # fit.gls=lmtest::coeftest(gls.fit,df=(n-1))[, ] %>% as.data.frame()
    gls.fit.true = gls(signal~jump, data = Data.mod,  correlation =  NULL, weights = NULL)
    fit.gls.true =lmtest::coeftest(gls.fit.true,df=(n-1))[, ] %>% as.data.frame()
    
    tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fit.gls.true = fit.gls.true)
    coef.res[[i]] = list( ols = ols.fit$coefficients, gls = gls.fit.true$coefficients)
    var.res[[i]] = list( ols = vcov(ols.fit), gls = gls.fit.true$varBeta, hac = vcov.para)

  }
  # significance level
  pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
  pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
  pval.gls.true <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls.true[2,4]))
  p.val = c(length(which(pval.ols>0.05)),
            length(which(pval.gls.true>0.05)), length(which(pval.hac>0.05)))
  Res.fin[l,] = p.val
  # r[[l]] = length(which( pval.gls>0.05))
}

# save(Res.fin, file = paste0(path_results, "res.RData"))
colnames(Res.fin) <- c("ols", "gls.t", "hac")
res = (nb.sim - Res.fin)/nb.sim
res$delta = param.test
dat.plot =reshape2::melt(res, id = "delta")
ggplot(dat.plot, aes(x = delta, y = value, col = variable))+
  geom_point() + theme_bw()+
  ylab("Significant jumps predicted rate")
# var (beta)
var.df = data.frame(ols = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$ols[1,1])),
                    hac = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$hac[1,1])),
                    gls = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$gls[1,1])))
summary(var.df)

# comparison in the AR(1) -------------------------------------------------
ar.val = seq(0, 0.8, 0.2)
# offset = seq(0.1, 0.5, 0.1)
param.test = ar.val
off.set = 0
Res.fin = data.frame(matrix(NA, ncol = 3, nrow = length(param.test)))
# ar = 0.3
for (l in c(1:length(param.test))) {
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  ar = param.test[l]
  for (i in c(1:nb.sim)) {
    set.seed(i)
    y = simulate.general(N = n, arma.model = c(ar,0), burn.in = 1000, hetero = 0, sigma = sqrt(sig.m),
                         monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + off.set
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = var.t, t = t)
    
    # Test with vcov (HAC)
    ols.fit = lm(signal~jump, data = Data.mod)
    vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE,kernel = "Quadratic Spectral", sandwich = TRUE)
    fit.hac=lmtest::coeftest(ols.fit,df=(n-1),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.ols=lmtest::coeftest(ols.fit,df=(n-1))[, ] %>% as.data.frame()
  
    gls.fit.true = gls(signal~jump, data = Data.mod,  correlation =  corAR1(form = ~t), weights = NULL)
    fit.gls.true =lmtest::coeftest(gls.fit.true,df=(n-1))[, ] %>% as.data.frame()
    
    tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fit.gls.true = fit.gls.true)
    coef.res[[i]] = list( ols = ols.fit$coefficients, gls = gls.fit.true$coefficients)
    var.res[[i]] = list( ols = vcov(ols.fit), gls = gls.fit.true$varBeta, hac = vcov.para)
    
  }
  # significance level
  pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
  pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
  pval.gls.true <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls.true[2,4]))
  p.val = c(length(which(pval.ols>0.05)),
            length(which(pval.gls.true>0.05)), length(which(pval.hac>0.05)))
  Res.fin[l,] = p.val
}

save(Res.fin, file = paste0(path_results, "res.2000.RData"))
colnames(Res.fin) <- c("ols", "gls.t", "hac")
res = (nb.sim - Res.fin)/nb.sim
name.x = "jump"
# res = Res.fin/nb.sim
res[name.x] = param.test
dat.plot =reshape2::melt(res, id = name.x)
ggplot(dat.plot, aes(x = jump, y = value, col = variable))+
  geom_point() + theme_bw()+
  ylab("True positive rate")
# var (beta)
var.df = data.frame(ols = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$ols[1,1])),
                    hac = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$hac[1,1])),
                    gls = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$gls[1,1])))
summary(var.df)


# comparison in ARMA(1,1) -------------------------------------------------

for (l in c(1:length(param.test))) {
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  
  for (i in c(1:nb.sim)) {
    set.seed(i)
    # y = trend*t + simulate.general(N = n, arma.model = c(ar,0), burn.in = 1000, hetero = 1, sigma = sqrt(var.t),
    #                      monthly.var = 0)
    y = simulate.general(N = n, arma.model = c(ar.val,0), burn.in = 1000, hetero = 1, sigma = sqrt(var.t),
                         monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + off.set
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = var.t, t = t)
    
    # Test with vcov (HAC)
    ols.fit = lm(signal~jump, data = Data.mod)
    vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE,kernel = "Quadratic Spectral", sandwich = TRUE)
    # vcov.para=sandwich::vcovHC(ols.fit, type = "HC3")
    fit.hac=lmtest::coeftest(ols.fit,df=(n-1),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.ols=lmtest::coeftest(ols.fit,df=(n-1))[, ] %>% as.data.frame()
    
    # Test with gls
    # et = ols.fit$residuals^2
    # Data.mod.e = data.frame(signal = et, sin1 = sin(2*pi*t/T1), cos1 = cos(2*pi*t/T1))
    # lm.et = lm(signal~., data = Data.mod.e)
    # var.e = fitted(lm.et)
    # gls.fit = gls(signal~jump, data = Data.mod, correlation = corAR1( form = ~t), weights = varFixed(value = ~var.e ))
    # fit.gls=lmtest::coeftest(gls.fit,df=(n-1))[, ] %>% as.data.frame()
    gls.fit.true = gls(signal~jump, data = Data.mod,  correlation =  corAR1(form = ~t), weights = varFixed(value = ~var.t ))
    fit.gls.true =lmtest::coeftest(gls.fit.true,df=(n-1))[, ] %>% as.data.frame()
    
    tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fit.gls.true = fit.gls.true)
    coef.res[[i]] = list( ols = ols.fit$coefficients, gls = gls.fit.true$coefficients)
    var.res[[i]] = list( ols = vcov(ols.fit), gls = gls.fit.true$varBeta, hac = vcov.para)
    
  }
  # significance level
  pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
  pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
  pval.gls.true <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls.true[2,4]))
  p.val = c(length(which(pval.ols>0.05)),
            length(which(pval.gls.true>0.05)), length(which(pval.hac>0.05)))
  Res.fin[l,] = p.val
  # r[[l]] = length(which( pval.gls>0.05))
}


# heteroskedastic noise ---------------------------------------------------
# offset = seq(0.1, 0.5, 0.1)
off.set = 0.3
var.all = seq(0, 0.8, 0.2)
param.test = var.all
Res.fin = data.frame(matrix(NA, ncol = 3, nrow = length(param.test)))
for (l in c(1:length(param.test))) {
  sig.v = param.test[l]
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  var.t = (sig.m - sig.v*a)
  plot(var.t)
  for (i in c(1:nb.sim)) {
    set.seed(i)

    y = simulate.general(N = n, arma.model = c(0.3,0), burn.in = 1000, hetero = 1, sigma = sqrt(var.t),
                         monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + off.set
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = var.t, t = t)
    
    # Test with vcov (HAC)
    ols.fit = lm(signal~jump, data = Data.mod)
    vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE,kernel = "Quadratic Spectral", sandwich = TRUE)
    fit.hac=lmtest::coeftest(ols.fit,df=(n-1),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.ols=lmtest::coeftest(ols.fit,df=(n-1))[, ] %>% as.data.frame()
    
    gls.fit.true = gls(signal~jump, data = Data.mod,  correlation =  corAR1(form = ~t), weights = varFixed(value = ~var.t ))
    fit.gls.true =lmtest::coeftest(gls.fit.true,df=(n-1))[, ] %>% as.data.frame()
    gls.fit.true =  gls.true(var.t, phi = 0, theta =0, design.matrix = Data.mod, trend = 0)
    fit.gls.true.t = gls.fit.true$beta/sqrt((diag(gls.fit.true$var.beta)))
    fit.gls.true.p = round(pnorm(-abs(fit.gls.true.t), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 4)
    fit.gls.true = fit.hac
    fit.gls.true$Estimate = gls.fit.true$beta
    fit.gls.true$`t value` = fit.gls.true.t
    fit.gls.true$`Pr(>|t|)` = fit.gls.true.p
    # 
    # Sig = rep(var.t,9)
    # cov.var = var_cova(Sig = Sig, phi = 0.3, theta = 0, n = n)
    # w = (sqrt(diag(cov.var)))
    # w = w - (mean(w)-1)
    # 
    # Data.mod$signal = y/w
    # ols.nor = lm(signal~jump, data = Data.mod)
    # vcov.para.nor=sandwich::kernHAC(ols.nor,prewhite = FALSE,kernel = "Quadratic Spectral", sandwich = TRUE)
    # fit.hac.norm=lmtest::coeftest(ols.nor,df=(n-1),vcov.=vcov.para.nor)[, ] %>% as.data.frame()
    # tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fit.gls.true = fit.hac.norm)
    
    tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fit.gls.true = fit.gls.true)
    coef.res[[i]] = list( ols = ols.fit$coefficients, gls = gls.fit.true$beta)
    var.res[[i]] = list( ols = vcov(ols.fit), gls = gls.fit.true$var.beta, hac = vcov.para)
    
  }
  # significance level
  pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
  pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
  pval.gls.true <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls.true[2,4]))
  p.val = c(length(which(pval.ols>0.05)),
            length(which(pval.gls.true>0.05)), length(which(pval.hac>0.05)))
  Res.fin[l,] = p.val
  # r[[l]] = length(which( pval.gls>0.05))
}
colnames(Res.fin) <- c("ols", "gls.t", "hac")
res = (nb.sim- Res.fin)/nb.sim
name.x = "delta"
# res = Res.fin/nb.sim
res[name.x] = param.test
dat.plot =reshape2::melt(res, id = name.x)
ggplot(dat.plot, aes(x = delta, y = value, col = variable))+
  geom_point() + theme_bw()+
  ylab("True positive rate") 


# heteroskedastic + AR(1) -----------------------------------------------------------------------

off.set = 0.3
# var.all = seq(0, 0.5, 0.1)
ar.val = seq( 0, 0.8, 0.2)
param.test =ar.val
Res.fin = data.frame(matrix(NA, ncol = 3, nrow = length(param.test)))
trend = -0.01
for (l in c(1:length(param.test))) {
  # sig.v = param.test[l]
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  ar = param.test[l]
  # var.t = (sig.m - sig.v*a)
  # plot(var.t)
  for (i in c(1:nb.sim)) {
    set.seed(i)
    
    y = trend * t + simulate.general(N = n, arma.model = c(ar,0), burn.in = 1000, hetero = 1, sigma = sqrt(var.t),
                         monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + off.set
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = var.t, t = t)
    
    # Test with vcov (HAC)
    ols.fit = lm(signal~jump+t, data = Data.mod)
    vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE,kernel = "Quadratic Spectral", sandwich = TRUE)
    fit.hac=lmtest::coeftest(ols.fit,df=(n-1),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.ols=lmtest::coeftest(ols.fit,df=(n-1))[, ] %>% as.data.frame()
    
    gls.fit.true = gls(signal~jump+t, data = Data.mod,  correlation =  corAR1(form = ~t), weights = varFixed(value = ~var.t ))
    fit.gls.true =lmtest::coeftest(gls.fit.true,df=(n-1))[, ] %>% as.data.frame()
    
    tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fit.gls.true = fit.gls.true)
    coef.res[[i]] = list( ols = ols.fit$coefficients, gls = gls.fit.true$coefficients)
    var.res[[i]] = list( ols = vcov(ols.fit), gls = gls.fit.true$varBeta, hac = vcov.para)
    
  }
  # significance level
  pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
  pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
  pval.gls.true <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls.true[2,4]))
  p.val = c(length(which(pval.ols>0.05)),
            length(which(pval.gls.true>0.05)), length(which(pval.hac>0.05)))
  Res.fin[l,] = p.val
  # r[[l]] = length(which( pval.gls>0.05))
}
colnames(Res.fin) <- c("ols", "gls.t", "hac")
res = (nb.sim - Res.fin)/nb.sim
name.x = "rho"
# res = Res.fin/nb.sim
res[name.x] = param.test
dat.plot =reshape2::melt(res, id = name.x)
ggplot(dat.plot, aes(x = rho, y = value, col = variable))+
  geom_point() + theme_bw()+
  ylab("True Positive rate") 

# test the variance estimation --------------------------------------------
tot.res <- data.frame(matrix(NA, ncol = 200, nrow = nb.sim))
for (i in c(1:nb.sim)) {
  set.seed(i)
  # y = simulate.general(N = n, arma.model = c(0,0), burn.in = 1000, hetero = 0, sigma = sqrt(var.t),
  #                      monthly.var = 0)
  # y[(n/2):n] <- y[(n/2):n] + offset
  y = rnorm(n, mean = 0, sd = sig.m )
  Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = var.t, t = t)
  
  # Test with vcov (HAC)
  ols.fit = lm(signal~jump, data = Data.mod)
  
  # Test with gls
  et = ols.fit$residuals^2
  Data.mod.e = data.frame(signal = et, sin1 = sin(2*pi*t/T1), cos1 = cos(2*pi*t/T1))
  lm.et = lm(signal~., data = Data.mod.e)
  var.e = fitted(lm.et)

  tot.res[i,] =var.e
  
}

a. = sapply(c(1:n), function(x) sd(tot.res[,x]))
dat = data.frame(x = c(1:200), value = a, sd = a.)
ggplot(dat, aes(x=x,y=value)) +
  geom_pointrange(aes(ymin=value-sd, ymax=value+sd))+ 
theme_bw() + 
  ylim(c(0,1.2))


