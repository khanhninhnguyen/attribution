# this prog is used to compare the FPR and TPR of the test of the change in mean between the OLS-HAC and the GLS

# input: what do you want to test. Ex: TPR of test when data is AR(1) with different rho
x.axis = "rho"
y.axis = ifelse(TPR == 1, "TPR", "FPR")

# specify simulation model
trend.sim = 0
# specify the regression model 
trend.reg = 0
mod.expression = ifelse(trend.reg == 1, "signal~jump+Xt", "signal~jump")

# specify the condition
mod.sim <- function(heteroscedastic, autocorr, var.inno, list.param.sig, list.param.ar, x.axis, individual, n, T1){
  a = cos(2*pi*(c(1:n)/T1)-pi)
  if(heteroscedastic == 1 & autocorr == 0){
    hetero = 1
    burn.in = 0
    sigma.t = sapply(c(1:length(list.param.sig)), function(x) var.inno - a*list.param.sig[x])
    ar = rep(0, length(list.param.sig))
  } else if(heteroscedastic == 0 & autocorr == 1){
    hetero = 0
    burn.in = 1000
    sigma.t = rep(var.inno, length(list.param.ar))
    ar = list.param.ar
  } else if (heteroscedastic == 1 & autocorr == 1){
    hetero = 1
    burn.in = 1000
    if(x.axis == "rho"){
      sigma.0 = var.inno - a*list.param.sig[individual]
      sigma.t = matrix( rep(sigma.0,  length(list.param.ar)) , ncol =  length(list.param.ar), byrow = FALSE )
      ar = list.param.ar
    } else if(x.axis == "sig.v"){
      sigma.t = sapply(c(1:length(list.param.sig)), function(x) var.inno - a*list.param.sig[x])
      ar = list.param.ar[individual]
    }
  }
  return(list(hetero = hetero, burn.in = burn.in, sigma.t = sigma.t, ar = ar))
}
gen.test = mod.sim(heteroscedastic = 1, autocorr = 1, var.inno = 1, list.param.ar = c(0.1, 0.3, 0.9), list.param.sig = c(0.1, 0.3), x.axis = "rho", individual = 2, n = 200, T1 = 100)

nb.sim = 10
n = 200
off.set = 0.3

res.var = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
t = c(1:n)-n/2
Res.fin = data.frame(matrix(NA, ncol = 3, nrow = length(gen.test$ar)))

hetero = gen.test$hetero
burn.in = gen.test$burn.in
for (l in c(1:length(gen.test$ar))) {
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  ar = gen.test$ar[l]
  sigma.sim = gen.test$sigma.t[,l]
  for (i in c(1:nb.sim)) {
    set.seed(i)
    y = simulate.general(N = n, arma.model = c(ar,0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim),
                         monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + off.set
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), var.t = sigma.sim, t = t)
    
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




