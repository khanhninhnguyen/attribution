# draft -------------------------------------------------------------------

# test the scaling factor between the variance of the difference a --------

phi = 0.8
theta = -0.5
var0 = 1
var_x = (1 + 2*phi*theta+theta**2)*var0/(1-phi**2)
print(var_x)
scale.diff = 2*(1-phi)*(1+theta**2-theta+phi*theta)/((1 + 2*phi*theta+theta**2))
scale.ma = 2*(1+ theta**2-theta)/(1+theta**2)
nb.sim = 1000
var.est = rep(NA, nb.sim)
var.diff.est = rep(NA, nb.sim)

for (i in 1:nb.sim) {
  sim.series <- arima.sim(model = list(ar = phi, ma = theta), n = 1000, sd = sqrt(var0))
  sim.diff <- diff(sim.series)
  var.est[i] <- var(sim.series)
  var.diff.est[i] <- var(sim.diff)/scale.diff
}
summary(var.est)
summary(var.diff.est)


# test the impact of estimate 2n and the mean -----------------------------
n = 200
nb.sim = 1000
var.2n = rep(NA, nb.sim)
var.mean = rep(NA, nb.sim)

for (i in 1:nb.sim) {
  x = rnorm(n, mean = 0, sd = 1)
  var.2n[i] = var(x)
  var.mean[i] = (var(x[1:100]) + var(x[101:200]))/2
}
summary(var.2n)
summary(var.mean)
var(var.2n)
var(var.mean)


library(gitcreds)








