# this prog is used for comparison: OLS+HAC vs GLS
source(paste0(path_code_att,"simulate_time_series.R"))
library(sandwich)
# initial condition ---------------------------------------------------------
nb.sim = 10000
n = 100
sig.m = 1
sig.v = 0.8
T1 = n/2
a = cos(2*pi*(c(1:n)/T1))
var.t = (sig.m - sig.v*a)
ar = c(0, 0)
res = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
res.var = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
t = c(1:n)-n/2
# comparison with heteroskedasticity --------------------------------------
#### IID case 
tot.res <- list()
coef.res <- list()
var.res <- list()
for (i in c(1:nb.sim)) {
  set.seed(i)
  y = simulate.general(N = n, arma.model = ar, burn.in = 0, hetero = 1, sigma = sqrt(var.t),
                       monthly.var = 0)
  y[(n/2):n] <- y[(n/2):n] + 0.3
  Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = var.t)
 
  # Test with vcov (HAC)
  ols.fit = lm(signal~jump, data = Data.mod)
  vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE,kernel = "Quadratic Spectral", sandwich = TRUE)
  # vcov.para=sandwich::vcovHC(ols.fit, type = "HC3")
  fit.hac=lmtest::coeftest(ols.fit,df=(n-1),vcov.=vcov.para)[, ] %>% as.data.frame()
  fit.ols=lmtest::coeftest(ols.fit,df=(n-1))[, ] %>% as.data.frame()
  
   # Test with gls
  gls.fit = gls(signal~jump, data = Data.mod, correlation = NULL, weights = varFixed(value = ~var.t))
  fit.gls=lmtest::coeftest(gls.fit,df=(n-1))[, ] %>% as.data.frame()
  tot.res[[i]] = list(fit.hac =fit.hac, fit.gls= fit.gls, fit.ols = fit.ols)
  coef.res[[i]] = list( ols = ols.fit$coefficients, gls = gls.fit$coefficients)
  var.res[[i]] = list( ols = vcov(ols.fit), gls = gls.fit$varBeta, hac = vcov.para)
}
# significance level
pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
pval.gls <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls[2,4]))
table(pval.gls<0.05)
table(pval.ols<0.05)
table(pval.hac<0.05)
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


