# this prog is used to compare the FPR and TPR of the test of the change in mean between the OLS-HAC and the GLS
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"FGLS.R"))
# source(paste0(path_code_att,"newUsed_functions.R"))

# input: what do you want to test. Ex: TPR of test when data is AR(1) with different rho------------------
off.set = 0.3
heteroscedast = 1
autocor = 0
x.axis = "rho"
one.year = 365

y.axis = ifelse(off.set !=0, "TPR", "FPR")
nb.sim = 100
n = 200
t = c(1:n)-n/2
list.param.ar = seq(0, 0.9, 0.15)
list.param.sig = seq(0.1, 0.35, 0.05)

# specify simulation model
trend.sim = 0
# specify the regression model 
trend.reg = 10

# specify the condition
mod.sim <- function(heteroscedastic, autocorr, var.inno, list.param.sig, list.param.ar, x.axis, individual, n, T1){
  a = cos(2*pi*(c(1:n)/T1) -pi)
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
      ar = rep(list.param.ar[individual],length(list.param.sig))
    }
  }
  return(list(hetero = hetero, burn.in = burn.in, sigma.t = sigma.t, ar = ar))
}
gen.test = mod.sim(heteroscedastic = heteroscedast, autocorr = autocor, var.inno = 0.4, list.param.ar = list.param.ar, 
                   list.param.sig = list.param.sig, x.axis = x.axis, individual = 6, n = n, T1 = n/2)



res.var = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
Res.fin = data.frame(matrix(NA, ncol = 4, nrow = length(gen.test$ar)))

hetero = gen.test$hetero
burn.in = gen.test$burn.in
sigma.sim = gen.test$sigma.t[,3]
ar=0.45
tot.res <- list()
coef.res <- list()
var.res <- list()
# test simulation FGLS
for (i in c(1:nb.sim)) {
  set.seed(i)
  y = simulate.general(N = n, arma.model = c(ar,0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim),
                       monthly.var = 0)
  y[(n/2):n] <- y[(n/2):n] + off.set
  df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2014-07-31"), by="days"))
  Data.mod = construct.design(data.df = df, name.series = "y", break.ind = 100)
  
  # Test with vcov (HAC)
  list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
  mod.X <-  list.para %>% stringr::str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
  # ols
  ols.fit = lm(mod.expression, data = Data.mod)
  vcov.para=sandwich::kernHAC(ols.fit, prewhite = TRUE, approx = c("AR(1)"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
  
  
  fit.hac=lmtest::coeftest(ols.fit,df=(n-2-trend.reg),vcov.=vcov.para)[, ] %>% as.data.frame()
  fit.ols=lmtest::coeftest(ols.fit,df=(n-2-trend.reg))[, ] %>% as.data.frame()
  
  #GLS with true covariance matrix
  # gls.fit.true = gls.true(var.t = Data.mod$weight1, phi = ar, theta = 0, design.matrix = Data.mod, trend = trend.reg)
  
  # FGLS
  fit.gls = FGLS(design.m = Data.mod, tol=0.00001, day.list = df$date)
  
  tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fit.gls = fit.gls)
  coef.res[[i]] = list( ols = ols.fit$coefficients, gls = fit.gls$coefficients)
  var.res[[i]] = list( ols = vcov(ols.fit), gls = fit.gls$varBeta, hac = vcov.para)
  
}
# significance level
pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
pval.gls.true <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls.true[2,4]))
pval.gls <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls$`Pr(>|t|)`[2,4]))

p.val = c(length(which(pval.ols>0.05)),
          # length(which(pval.gls.true>0.05)),
          length(which(pval.gls>0.05)), 
          length(which(pval.hac>0.05)))