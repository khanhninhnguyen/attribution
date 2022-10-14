
# this prog used to test the HC and HAC variance estimator  ---------------
source(paste0(path_code_att,"simulate_time_series.R"))
library(sandwich)
gls.func <- function(Y, X, sigma.matrix){
  sigma.inv = solve(sigma.matrix)
  beta = solve(t(X)%*%sigma.inv%*%X)%*%t(X)%*%sigma.inv%*%Y
  return(beta)
}
ols.func <- function(Y, X){
  beta = solve(t(X)%*%X)%*%t(X)%*%Y
  return(beta)
}

nb.sim = 10000
n = 10
sig.m = 1
sig.v = 0.8
T1 = n*2
true.var = (sig.m - sig.v*cos(2*pi*(c(1:n)/T1)))

res = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
res.beta = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))

x = as.matrix(rep(1, n))
Sig = matrix(0, nrow = n, ncol = n)
diag(Sig) <- true.var

var.beta.OLS = solve(t(x)%*%x)%*%t(x)%*%Sig%*%x%*%solve(t(x)%*%x)
var.beta.GLS = solve(t(x)%*%solve(Sig)%*%x)


for (i in c(1:nb.sim)) {
  set.seed(i)
  y = simulate.general(N = n, auto = 0, arma.model = c(0,0), burn.in = 0, hetero = 1, sigma = sqrt(true.var),
                       monthly.var = 0)
  y = y + 1
  
  fit1=lm(y~x)
  # Covariance matrix of beta
  hc0=vcovHC(fit1,type="HC0")
  hc1=vcovHC(fit1,type="HC1")
  hc2=vcovHC(fit1,type="HC2")
  hc3=vcovHC(fit1,type="HC3")
  sandwich.1 = sandwich(fit1, bread. = bread(fit1), meat. = meat(fit1))
  V_HAC <- vcovHAC(fit1)
  
  res[i,] <- c(as.numeric(hc0), as.numeric(hc1), as.numeric(hc2), as.numeric(hc3), as.numeric(sandwich.1), as.numeric(V_HAC))
  # fit ols and gls 
  beta.ols = ols.func(X = x, Y = y)
  beta.gls = gls.func(Y = y, X = x, sigma.matrix = Sig)
  res.beta[i,] <- c(beta.ols, beta.gls)
}

summary(res)
apply(res, 2, var)
summary(res.beta)
apply(res.beta, 2, var)

a = 10+0.125*(50+60*(c(1:100)))
b = 0.25*exp(0.04*(50+60*(c(1:10))))




n1 = 100000
r1 = data.frame(matrix(NA, nrow = n1, ncol = 10))
r2 = data.frame(matrix(NA, nrow = n1, ncol = 10))

for (i in c(1:n1)) {
  x = rnorm(10, 0 , 1)
  z = x-mean(x)
  r1[i,] <- x
  r2[i,] <- z
}

a1 = colMeans(r1)
a2 = colMeans(r2)
mean(a1)
mean(a2)

apply(r1, 2, var)

mean(apply(r2, 2, sd))









