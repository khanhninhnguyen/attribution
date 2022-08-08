rm(list = ls())

nb.sim = 10000
n0 = 60
sigma = data.frame(matrix(NA, ncol =10, nrow = nb.sim))
alpha = 0.5
for (i in 1:nb.sim) {
  # x =  rnorm(n = n0, mean = 0, sd = 1)
  x = arima.sim(model = list(ar = alpha), n = n0, sd = 1)
  x[(floor(n0/2)+1):n0] <- x[(floor(n0/2)+1):n0] +1
  a = diff(x)
  sigma[i,1] = sd(x)
  sigma[i,2] <- mad(x)
  sigma[i,3] <- robustbase::scaleTau2(x)
  sigma[i,4] <- robustbase::Qn(x)
  sigma[i,5] <- robustbase::scaleTau2(x)
  sigma[i,6] = sd(a)/sqrt(2-2*alpha)
  sigma[i,7] <- mad(a)/sqrt(2-2*alpha)
  sigma[i,8] <- robustbase::scaleTau2(a)/sqrt(2-2*alpha)
  sigma[i,9] <- robustbase::Qn(a)/sqrt(2-2*alpha)
  sigma[i,10] <- robustbase::scaleTau2(a)/sqrt(2-2*alpha)
}
colnames(sigma) <- c("classic", "MAD", "Sn", "Qn","ScaleTAu", "classic.diff","MAD.diff", "Sn.diff", "Qn.diff","ScaleTAu.diff")


boxplot(sigma[,c(6:10)], ylab = "Std.est", main = "n = 60, offset = 1, std = 1, phi = 0.5")
#abline(h = 1)
abline(h = 1/sqrt((1-alpha**2)))

