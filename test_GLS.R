# this prog used to test how GLS work 
source(paste0(path_code_att,"simulate_time_series.R"))
library(sandwich)
# inverse matrix : solve(mat)
# transpose matrix: t(matrix) 
gls.func <- function(Y, X, sigma.matrix){
  sigma.inv = solve(sigma.matrix)
  beta = solve(t(X)%*%sigma.inv%*%X)%*%t(X)%*%sigma.inv%*%Y
  return(beta)
}
ols.func <- function(Y, X){
  beta = solve(t(X)%*%X)%*%t(X)%*%Y
  return(beta)
}

# ols.func(Y = y, X= x)
# Sig = diag(length(y))
# a = gls.func(Y = y, X = x, sigma.matrix = Sig)

# homoskedatic 
nb.sim = 10000
n = 10
res = data.frame(matrix(NA, ncol = 3, nrow = nb.sim))
res.var = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
mu = rep(NA, nb.sim)
alpha = 0.5

x = as.matrix(rep(1, n))
s0 = c(1, alpha^(c(1:(n-1))))
Sig = (1/(1-alpha**2)) * toeplitz(s0)

for (i in c(1:nb.sim)) {
  set.seed(i)
  # y = as.matrix(rnorm(n, 0, 1)) 
  # y = y * b
  y = simulate.general(N = n, auto = 1, arma.model = c(alpha,0), burn.in = 1000, hetero = 0, sigma = 1,
                       monthly.var = 0)
  y = y + 1
  # mu[i] = sum(y/(b**2))/(sum(1/(b**2)))
 
  beta.ols = ols.func(X = x, Y = y)
  # Sig = b*diag(length(y)) heteroskedastic 
  beta.gls = gls.func(Y = y, X = x, sigma.matrix = Sig)
  
  fit1=lm(y~x)
  # e1=fit1$resid
  hc=vcovHC(fit1,type="HC0")
  # a = lmtest::coeftest(fit1, vcov. = hc)
  # res.hc = as.data.frame(a[, ])
  
  res[i,] <- c(beta.gls, beta.ols, as.numeric(hc))
  # res.var[i,] <- c(sum((y - mean(y))^2)/(n-1), sum((y - mean(y))^2)/(n-1))
}
summary(res)
summary(sqrt(res.var))
apply(res, 2, sd)
apply(sqrt(res.var), 2, sd)


sig.m = 1
sig.v = 0.8
T1 = n
a = cos(2*pi*(c(1:n)/T1))
b = sig.m - sig.v*a
sum(y/(b**2))/(sum(1/(b**2)))


# compute the true covariance matrix of ARMA 
inv.s = solve(Sig)
t1 = solve(t(x)%*%x)%*%t(x)
d.ols = sqrt(t1%*%Sig%*%t(t1))
d.gls = sqrt(solve(t(x)%*%solve(Sig)%*%x))
   # (1/(1-0.9**2)) * solve(t(x)%*% inv.s %*% x)
d1 = solve(t(x)%*%x)


fit1=lm(y_diff~x_diff)
print(summary(fit1))
e1=fit1$resid
# Heteroskedasticity-Consistent Covariance Matrix Estimation
#hc0=vcovHC(fit1,type="const")
#print(sqrt(diag(hc0)))
# type=c("const","HC","HC0","HC1","HC2","HC3","HC4")

# HC0 is the White estimator
hc1=vcovHC(fit1,type="HC0")
print(sqrt(diag(hc1)))
#Heteroskedasticity and autocorrelation consistent (HAC) estimation
#of the covariance matrix of the coefficient estimates in a
#(generalized) linear regression model.
hac1=vcovHAC(fit1,sandwich=T)
print(sqrt(diag(hac1)))

lmtest::coeftest(fit.signal, sandwich::NeweyWest(fit.signal, lag = 1))[, ]


print(lmtest::coeftest(fit1, vcov. = hc1))

