# this prog used to find the good approch for the normalization of data 
source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"simulate_time_series.R"))

# compare moving window vs loess:
nb.sim = 10000
list.sd = seq(0.2,0.8, 0.2)
length.month1 = c(31,28,31,30,31,30,31,31,30,31,30,31)
L = 365
n = 365
t = c(1:n)
case = 4
std.t = list.sd[case]*2*sin(2*pi*t/L-pi/2)/2 +1
std = std.t[seq(1,365,length.out = 12)]
res1 <- data.frame(matrix(NA, nrow = nb.sim, ncol = 365))
res2 <- data.frame(matrix(NA, nrow = nb.sim, ncol = 365))
res3 <- data.frame(matrix(NA, nrow = nb.sim, ncol = 365))

t0 <- as.Date("2021-01-01", format = "%Y-%m-%d" )
t.year <- seq(t0, t0+364,1)
for (i in c(1:nb.sim)) {
  set.seed(i)
  # sim.ar <- simulate.series.2(mean.1 = 0, sigma.1 = std,
  #                             N = n,  arma.model = c(0, 0), length.month = length.month1)
  sim.ar = simulate.general(auto = 0,
                       burn.in = 0,
                       arma.model = c(0,0),
                       hetero = 1,
                       monthly.var = 0,
                       sigma = std.t,
                       N = n,
                       gaps = 0,
                       outlier = 1,
                       prob.outliers = 0.1,
                       size.outliers = 3)
  # sim.ar = rnorm(n, 0, 1)

  std.est <- RobEstiSlidingVariance.S(Y = data.frame(date = t.year[1:n], y = sim.ar),
                                      name.var = "y",
                                      alpha = 0, estimator = "Sca",
                                      length.wind = 60)

  std.est1 <- RobEstiSlidingVariance.WLS (Y = data.frame(date = t.year[1:n], y = sim.ar),
                                          name.var = "y",
                                          alpha = 0,
                                          length.wind = 60, loes.method = "gaussian")

  std.est2 <- RobEstiSlidingVariance.WLS (Y = data.frame(date = t.year[1:n], y = sim.ar),
                                          name.var = "y",
                                          alpha = 0,
                                          length.wind = 60, loes.method = "symmetric")

  res1[i,] <- std.est 
  res2[i,] <- std.est1 
  res3[i,] <- std.est2

}

s1 <- colMeans(res1)

s2 <- colMeans(res2)

s3 <- colMeans(res3)

sd1 = apply(res1,2,sd)
sd2 = apply(res1,2,sd)
sd3 = apply(res1,2,sd)

# std.t = rep(1, n)
MSE1 = (colMeans(res1) -std.t  )^2 +  apply(res1,2,var)
MSE2 = (colMeans(res2) -std.t  )^2 +  apply(res2,2,var)
MSE3 = (colMeans(res3) -std.t  )^2 +  apply(res2,2,var)

# plot(t.year, s3)
# plot(t.year, s2, col = "red")
# plot(t.year, s1, col = "blue")

dat = data.frame( t = t.year, true = std.t, MW.ScaleTau = s1, loess.g = s2, loess.s = s3,
                  bias.MW.ScaleTau = abs(s1 - std.t), bias.loess.g = abs(s2 - std.t),  bias.loess.s = abs(s3 - std.t))
a = reshape2::melt(dat, id = "t")
jpeg(paste0(path_results,"attribution/variances/bias.model2",list.sd[case], ".Loess.0.jpeg" ),
     width = 3000, height = 1500,res = 300)
ggplot(data = a, aes(x = t, y = value, col = variable)) + 
  geom_line()+
  theme_bw()+
  ylab("Mean of estimated scale")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()
a = as.numeric(res1 %>% summarise_if(is.numeric, sd))
b = as.numeric(res2 %>% summarise_if(is.numeric, sd))
d = as.numeric(res3 %>% summarise_if(is.numeric, sd))


dat = data.frame( t = t.year, MW.ScaleTau = a, loess.g = b, loess.s = d)
a = reshape2::melt(dat, id = "t")
jpeg(paste0(path_results,"attribution/variances/SD.model2",list.sd[case], ".Loess.0.jpeg" ),
     width = 3000, height = 1500,res = 300)
ggplot(data = a, aes(x = t, y = value, col = variable)) + 
  geom_line()+
  theme_bw()+
  ylab("Sd of estimated scale")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()

dat = data.frame( t = t.year, MW.ScaleTau = MSE1, loess.g = MSE2, loess.s = MSE3)
a = reshape2::melt(dat, id = "t")
jpeg(paste0(path_results,"attribution/variances/MSE2",list.sd[case], ".Loess.0.jpeg" ),
     width = 3000, height = 1500,res = 300)
ggplot(data = a, aes(x = t, y = value, col = variable)) + 
  geom_line()+
  theme_bw()+
  ylab("MSE")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()














