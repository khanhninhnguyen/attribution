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
a = var_cova(Sig, phi = 0.6, theta = -0.2, n =100)
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



