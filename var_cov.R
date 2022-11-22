#this function is used to compute the variance covariance matrix of the ARMA(1,1) with the heteroskedastic noise
rm(list = ls())

var_cova <- function(Sig, phi, theta, n){
  M = 9*n
  k = c(0:(M-n-2))
  
  Sig.k = Sig[2:(M-n)] 
  Sig.k1 = Sig[1:(M-n-1)] 
  
  c1 = (phi^(2*k))
  c2 = c1*(theta^2)
  c3 = (phi^(2*k+1))*theta
  c4 = (phi^(2*k-1))*theta
  
  term1 = rev(c2) * Sig.k1
  term2 = rev(c1) * Sig.k
  term3 = rev(c3) * Sig.k1
  term4 = rev(c4) * Sig.k
  
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
  return(var.t)
}

Sig = c(rep(c(1,2,3), 300),1)
# Sig = rep(1,901)
a = var_cova(Sig, phi = 0.5, theta = 0, n =100)
phi = 0.6
theta = -0.2
(1+ 2*theta*phi + theta**2)/(1-phi**2)
va = rep(NA, 100)
v.all = data.frame(matrix(NA, ncol = 2, nrow = 100))
for (t in c(802:901)) {
  v = 0
  for (k in c(0:800)) {
    print(t-k)
    a1 = (phi**(2*k)) * (theta**2) * Sig[t-k-1]
    a2 = (phi**(2*k)) * Sig[t-k]
    a3 = (phi**(2*k+1)) * theta * Sig[t-k-1]
    a4 = (phi**(2*k-1)) * theta * Sig[t-k]
    v.k = a1+a2+a3+a4
    print(v.k)
    v = v.k+v
    v.all[(k+1),2] = v.k
  }
  v = v - (phi**(-1) * theta * Sig[t])
  ind = t -801
  va[ind] = v
}

# test in simulated series
nb.sim = 10000
n = 100
sig.m = 0.6
sig.v = 0.4
T1 = n/2
a = cos(2*pi*(c(1:n)/T1) )
var.t = (sig.m - sig.v*a)
dat.all <- data.frame(matrix(NA, ncol = n, nrow = nb.sim))
for (i in c(1:nb.sim)) {
  set.seed(i)
  y = simulate.general(N = n, arma.model = c(0.9,0), burn.in =1000, hetero = 1, sigma = sqrt(var.t),
                       monthly.var = 0)
  dat.all[i,] <- y
}

b = sapply(c(1:n), function(x) var(dat.all[,x]))
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
cov1 = ar1_cor(n, 0.9)
e = diag(sqrt(var.t)) %*% cov1 %*% diag(sqrt(var.t))
f= diag(e)/(1-phi**2)
phi = 0.9


# comparison 
# olivier res
Sig =  rep(var.t, 9)
Sig = c(Sig[1], Sig)
var.O = var_cova(Sig, phi = 0.9, theta = 0, n =100)
# decomposition
cov1 = ar1_cor(n, 0.9)
e = diag(sqrt(var.t)) %*% cov1 %*% diag(sqrt(var.t))
f= diag(e)/(1-phi**2)
# effective sample size 
# s = (n + 2*(phi^(n+1)) - n*(phi**2) - 2*phi)/(n*((1-phi)**2))
# f = var.t*s





