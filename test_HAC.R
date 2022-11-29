# test hac
# effectiveness of the HC and HAC in different case 
source(paste0(path_code_att,"simulate_time_series.R"))

nb.sim = 1000
off.set = 0.3
trend = 0.6
n=200
T1 = n/6
a = cos(2*pi*(c(1:n)/T1))
var.m = 0.4
var.t = var.m - 0.35*a
# var.all = seq(0, 0.5, 0.1)
ar = 0.3
coef.all = data.frame(matrix(NA, ncol = 5, nrow = nb.sim))
var.all = data.frame(matrix(NA, ncol = 5, nrow = nb.sim))
t = c((-n/2):(n/2-1))
t1 = c(1:n)
for (i in c(1:nb.sim)) {
  set.seed(i)
  y =  simulate.general(N = n, arma.model = c(0.3,0), burn.in = 0, hetero = 0, sigma = sqrt(var.t),
                        monthly.var = 0)
  y[(n/2):n] = y[(n/2):n]+off.set
  Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2))
 
  ols.fit = lm(signal~jump, data = Data.mod)
  
  # HAC kernel 
  # vcov.ker0=sandwich::kernHAC(ols.fit,prewhite = 0, adjust = TRUE, kernel = "Quadratic Spectral", sandwich = TRUE)
  # vcov.ker1=sandwich::kernHAC(ols.fit,prewhite = 1, adjust = TRUE,  kernel = "Quadratic Spectral", sandwich = TRUE)
  # vcov.hc3=sandwich::vcovHC(ols.fit, type = "HC3", sandwich = TRUE)
  # vcov.hac0=sandwich::vcovHAC(ols.fit, prewhite = FALSE)
  # vcov.hac1=sandwich::vcovHAC(ols.fit, prewhite = TRUE)
  
  p.k0 = lmtest::coeftest(ols.fit,df=(n-2),vcov.=vcov.ker0)[, ] %>% as.data.frame()
  p.k1 = lmtest::coeftest(ols.fit,df=(n-2),vcov.=vcov.ker1)[, ] %>% as.data.frame()
  p.hc3 = lmtest::coeftest(ols.fit,df=(n-2),vcov.=vcov.hc3)[, ] %>% as.data.frame()
  p.hac0 = lmtest::coeftest(ols.fit,df=(n-2))[, ] %>% as.data.frame()
  p.hac1 = lmtest::coeftest(ols.fit,df=(n-2),vcov.=vcov.hac1)[, ] %>% as.data.frame()
  coef.all[i,] = c(p.k0[2,4], p.k1[2,4], p.hc3[2,4], p.hac0[2,4], p.hac1[2,4] )
  var.all[i,] = c(vcov.ker0[2,2], vcov.ker1[2,2], vcov.hc3[2,2], vcov.hac0[2,2], vcov.hac1[2,2])
  
}
1 - sapply(c(1:5), function(x) length(which(coef.all[,x] > 0.05)))/nb.sim

sapply(c(1:5), function(x) length(which(coef.all[,x] < 0.05)))/nb.sim

