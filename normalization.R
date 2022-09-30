# this prog used to find the good approch for the normalization of data 
source(paste0(path_code_att,"sliding_variance.R"))

# compare moving window vs loess:
nb.sim = 1000
list.sd = seq(0.2,0.8, 0.2)
length.month1 = c(31,28,31,30,31,30,31,31,30,31,30,31)
L = 365
n = 365
t = c(1:n)
std.t = list.sd[1]*2*sin(2*pi*t/L-pi/2)/2 +1
std = std.t[seq(1,365,length.out = 12)]
res1 <- data.frame(matrix(NA, nrow = nb.sim, ncol = 365))
res2 <- data.frame(matrix(NA, nrow = nb.sim, ncol = 365))
t0 <- as.Date("2021-01-01", format = "%Y-%m-%d" )
t.year <- seq(t0, t0+364,1)
for (i in c(1:nb.sim)) {
  set.seed(i)
  sim.ar <- simulate.series.2(mean.1 = 0, sigma.1 = std,
                              N = n,  arma.model = c(0, 0), length.month = length.month1)
  # SimData = SimulatedSeries(n = n, 
  #                           P = choose_model(1), 
  #                           prob.outliers = 0.1,
  #                           size.outliers = 5, 
  #                           rho = 0, theta = 0)
  
  # 
  std.est <- RobEstiSlidingVariance.S(Y = data.frame(date = t.year, y = sim.ar), 
                                      name.var = "y", 
                                      alpha = 0, estimator = "Sca",
                                      length.wind = 60)
  
  std.est1 <- RobEstiSlidingVariance.WLS (Y = data.frame(date = t.year, y = sim.ar), 
                                          name.var = "y", 
                                          alpha = 0, 
                                          length.wind = 60)
  res1[i,] <- std.est
  res2[i,] <- std.est1
  
}

s1 <- colMeans(res1)

s2 <- colMeans(res2)

lines(t.year, s1)
plot(t.year, s2)
dat = data.frame( t = t.year, true = std.t, MW.ScaleTau = s1, loess = s2, 
                  bias1 = abs(s1 - std.t), bias2 = abs(s2 - std.t))
a = reshape2::melt(dat, id = "t")
jpeg(paste0(path_results,"attribution/variances/bias.model1.0.2.Loess60.jpeg" ),
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

dat = data.frame( t = t.year, MW.ScaleTau = a, loess = b)
a = reshape2::melt(dat, id = "t")
jpeg(paste0(path_results,"attribution/variances/bias.model1.0.2.Loess.sd60.jpeg" ),
     width = 3000, height = 1500,res = 300)
ggplot(data = a, aes(x = t, y = value, col = variable)) + 
  geom_line()+
  theme_bw()+
  ylab("Mean of estimated scale")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()