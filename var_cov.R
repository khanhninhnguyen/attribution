#this function is used to compute the variance covariance matrix of the ARMA(1,1) with the heteroskedastic noise

n = 100
M = 9*n
phi = 0.5
theta = -0.4
k = c(1:M)

(phi^(2*k))*(theta^2)