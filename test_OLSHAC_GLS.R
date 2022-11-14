# this prog is used for comparison: OLS+HAC vs GLS
source(paste0(path_code_att,"simulate_time_series.R"))
library(sandwich)
# initial condition ---------------------------------------------------------
nb.sim = 10000
n = 200
sig.m = 0.6
sig.v = 0.1
T1 = n/2
a = cos(2*pi*(c(1:n)/T1))
var.t = (sig.m - sig.v*a)
res = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
res.var = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
t = c(1:n)-n/2
# off.set = c(0.1, 0.2, 0.3, 0.4, 0.5)
ar.val = 0.3
offset = c(0.1, 0.2, 0.3, 0.4, 0.5)
# comparison with heteroskedasticity --------------------------------------
#### IID case with the estimated variance in GLS
param.test = offset
Res.fin = data.frame(matrix(NA, ncol = 4, nrow = length(param.test)))
for (l in c(1:length(param.test))) {
  offset  = param.test[l]
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  for (i in c(1:nb.sim)) {
    set.seed(i)
    y = simulate.general(N = n, arma.model = c(ar.val,0), burn.in = 1000, hetero = 1, sigma = sqrt(var.t),
                         monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + offset
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = var.t, t = t)
    
    # Test with vcov (HAC)
    ols.fit = lm(signal~jump, data = Data.mod)
    vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE,kernel = "Quadratic Spectral", sandwich = TRUE)
    # vcov.para=sandwich::vcovHC(ols.fit, type = "HC3")
    fit.hac=lmtest::coeftest(ols.fit,df=(n-1),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.ols=lmtest::coeftest(ols.fit,df=(n-1))[, ] %>% as.data.frame()

    # Test with gls
    et = ols.fit$residuals^2
    Data.mod.e = data.frame(signal = et, sin1 = sin(2*pi*t/T1), cos1 = cos(2*pi*t/T1))
    lm.et = lm(signal~., data = Data.mod.e)
    var.e = fitted(lm.et)
    gls.fit = gls(signal~jump, data = Data.mod, correlation = corAR1(form = ~t), weights = varFixed(value = ~var.e ))
    fit.gls=lmtest::coeftest(gls.fit,df=(n-1))[, ] %>% as.data.frame()
    gls.fit.true = gls(signal~jump, data = Data.mod,  correlation = corAR1(form = ~t), weights = varFixed(value = ~var.t ))
    fit.gls.true =lmtest::coeftest(gls.fit.true,df=(n-1))[, ] %>% as.data.frame()
    
    tot.res[[i]] = list(fit.hac =fit.hac, fit.gls= fit.gls, fit.ols = fit.ols, fit.gls.true = fit.gls.true)
    coef.res[[i]] = list( ols = ols.fit$coefficients, gls = gls.fit$coefficients)
    var.res[[i]] = list( ols = vcov(ols.fit), gls = gls.fit$varBeta, hac = vcov.para)
    # tot.res[[i]] = list(fit.gls.true = fit.gls.true)
    
  }
  # significance level
  pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
  pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
  pval.gls <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls[2,4]))
  pval.gls.true <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls.true[2,4]))
  p.val = c(length(which(pval.ols>0.05)), length(which(pval.gls>0.05)), 
            length(which(pval.gls.true>0.05)), length(which(pval.hac>0.05)))
  Res.fin[l,] = p.val
  # r[[l]] = length(which( pval.gls>0.05))
}
colnames(Res.fin) <- c("ols", "gls", "gls.t", "hac")
res = (nb.sim - Res.fin)/nb.sim
Res.fin$offset = param.test
dat.plot =reshape2::melt(res, id = "offset")
ggplot(dat.plot, aes(x = offset, y = (1-value/nb.sim), col = variable))+
  geom_point() + theme_bw()+
  ylab("TPR")
# var (beta)
var.df = data.frame(ols = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$ols[1,1])),
                    hac = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$hac[1,1])),
                    gls = unlist(sapply(c(1:nb.sim), function(x) var.res[[x]]$gls[1,1])))
summary(var.df)


x = as.matrix(rep(1,n))
ols.v = var(y)/(t(x) %*% x)
xtx_inv <- solve(t(x) %*% x)

xtx_inv
hac.v = xtx_inv %*% t(x) %*% diag(var.t) %*% x %*% xtx_inv
hac.v

# comparison in the AR(1) -------------------------------------------------

# comparison in ARMA(1,1) -------------------------------------------------


