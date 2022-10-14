# test the screening methods 
source(paste0(path_code_att,"support_test_screening_methods.R"))
source(paste0(path_code_att,"support_screening.R"))
source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"UsedFunctions.R"))

nb.sim = 10000
res <- list()
length.list = c(10, 20, 30, 60, 120, 360, 720)
outlier.list = seq(0.00, 0.1, 0.01)
rho.list = seq(0.0, 0.9, 0.15)
var.list <- outlier.list
size.outliers = 3

choose_model <- function(x){
  if(x == 1){ P = 5 } # fixed value of outlier
  else if(x == 2){ P = 3 } # normal outlier 
  else if(x == 3){ P = 4 } # normal overlayed outlier 
  else if(x == 4){ P = 7 } # normal overlayed outlier + arma data
  else if(x == 5){ P = 8 } # normal overlayed outlier + replaced at a threshold
  else if(x == 6){ P = 1 } # no outlier
  return(P)
}
model.out = 4
P = choose_model(model.out)
for (l in c(1:length(var.list))) {
  n0 = 10
  # prob.outliers = var.list[l]
  # print(prob.outliers)
  sigma = data.frame(matrix(NA, ncol =5, nrow = nb.sim))
  rho = 0
  theta = 0
  off.set = 0
  ratio = sqrt((1 + 2*theta*rho + theta**2)/(1 - rho**2))
  for (i in 1:nb.sim) {
    set.seed(i)
    # x = SimulatedSeries(n = n0, P = P, prob.outliers = prob.outliers, size.outliers = size.outliers, rho = rho, theta = theta)$Y
    x <- rnorm(n = n0, mean = 0, sd =1)
    # x <- arima.sim(model = list(ar = rho), n = n0, sd = 1, n.start = 10000)
    # x[(floor(n0/2)+1):n0] <- x[(floor(n0/2)+1):n0]+off.set
    # x[(floor(n0*3/4)+1):n0] <- x[(floor(n0*3/4)+1):n0]+off.set
    # alpha.e = arima(x, order = c(1,0,0), method = "ML")$coef[1]
    # a = diff(x)
    # if it is AR(1) add this factor on the estimator : *sqrt(1-rho^2)
    sigma[i,1] <- sd(x)
    sigma[i,2] <- mad(x)
    sigma[i,3] <- robustbase::Sn(x)
    sigma[i,4] <- robustbase::Qn(x)
    sigma[i,5] <- robustbase::scaleTau2(x)
    # sigma[i,1] = sd(a)/sqrt(2-2*alpha.e)
    # sigma[i,2] <- mad(a)/sqrt(2-2*alpha.e)
    # sigma[i,3] <- robustbase::Sn(a)/sqrt(2-2*alpha.e)
    # sigma[i,4] <- robustbase::Qn(a)/sqrt(2-2*alpha.e)
    # sigma[i,5] <- robustbase::scaleTau2(a)/sqrt(2-2*alpha.e)
  }
  colnames(sigma) <- c("SD", "MAD", "Sn", "Qn","ScaleTAu")
  std <- sapply(c(1:5), function(x) IQR(sigma[,x]))
  bias <- sapply(c(1:5), function(x) mean(sigma[,x])-1)
  res$sdt[[l]] <- std
  res$bias[[l]] <- bias
  # b = data.frame(estimator = rep(colnames(sigma), each = nb.sim), std = unlist(sigma))
  # p <- ggplot(b, aes(x = estimator, y = std))+ 
  #   geom_boxplot()+theme_bw()+
  #   labs(subtitle = paste0("phi = ", alpha, ", n = ", n0, ", offset = ", off.set))+
  #   geom_hline(yintercept = 1)+
  #   theme(axis.text = element_text(size = 15))
  # if(n0==30){p <- p +  ylim(c(0,2))
  # }else{p <- p +  ylim(c(0.5,1.5))}
  # 
  # jpeg(paste0(path_results,"attribution/variances/sim.e", "phi = ", alpha, ", n = ", n0, ", offset = ", off.set, "o2.jpeg"),
  #      width = 3000, height = 1800,res = 300) # change name
  # print(p)
  # dev.off()
  
}


# Box plot to easily compare different estimators -----------------
boxplot(sigma, ylab = "Std.est", main = "n = 60, offset = 1, std = 1, phi = 0.5")
abline(h = 1)
abline(h = 1/sqrt((1-alpha**2)))

sigma = data.frame(matrix(NA, ncol =5, nrow = nb.sim))
d <- c()
for (k in c(1:nb.sim)){
  set.seed(k)
  x = arima.sim(model = list(  ar = alpha ), n = n0, sd = sqrt(1-alpha**2))
  x[(floor(n0/2)+1):n0] <- x[(floor(n0/2)+1):n0]+2
  a = diff(x)
  d <- c(d,a[15])
}

hist(d)
summary(d)

# Scatter plot ----------------
# IQR
lab.x = "Outlier percentage"

iqr = as.data.frame(res$sdt)
colnames(iqr) <- as.character(var.list)
# ratio.l = sqrt((1 + 2*theta*rho.list + theta**2)/(1 - rho.list**2))
# true.bias <- iqr[1,1]*ratio.l
# iqr <- rbind(iqr, c(true.bias))
iqr$est <- c("SD", "MAD", "Sn", "Qn","ScaleTAu")
a = reshape2::melt(iqr, id ="est")
data.plot <- data.frame()
jpeg(paste0(path_results,"attribution/variances/IQR.model",model.out, "rho = ", rho, " theta = ", theta, "out.size = ",  size.outliers, "offset = ", off.set, prob.outliers,"60.jpeg" ),
     width = 3000, height = 1500,res = 300)
ggplot(data = a, aes(x = variable, y = value, col = est)) + 
  geom_point()+
  theme_bw()+
  xlab(lab.x)+
  scale_y_continuous( name = "IQR", limits = c(0, ceiling(max(a$value)*10)/10), breaks = seq(0.2, ceiling(max(a$value)*10)/10, 0.1))+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()
# BIAS
iqr = as.data.frame(res$bias)
colnames(iqr) <- as.character(var.list)
# true.bias <- sample.b 
# iqr <- rbind(iqr, c(true.bias))
iqr$est <- c("SD", "MAD", "Sn", "Qn","ScaleTAu")
a = reshape2::melt(iqr, id ="est")
lim1 = floor(min(a$value)*10)/10
lim2 = ceiling(max(a$value)*10)/10
data.plot <- data.frame()
jpeg(paste0(path_results,"attribution/variances/bias.model", model.out, "rho = ", rho, " theta = ", theta, "out.size = ",  size.outliers, "offset = ", off.set, prob.outliers, "60.jpeg" ),
     width = 3000, height = 1500,res = 300)
ggplot(data = a, aes(x = variable, y = value, col = est)) + 
  geom_point()+
  theme_bw()+
  xlab(lab.x)+
  scale_y_continuous( name = "Bias", limits = c(lim1, lim2), breaks = round(seq(lim1, lim2, 0.1), digits=2)) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()

d <- rep(NA, 10000)
r = 1/(1-0.9**2)
e = ARMAacf(ar=0.5, lag.max = 100)
s=0
for (j in (1:99)) {
  s = s+ (1-j/100)*e[j]
}
e = (1-2*s/99)*(1/(1-0.5**2))
for (i in c(1:10000)) {
  # x <- rnorm(20, 0 ,1)
  x <- arima.sim(model = list(ar = 0.5), n =100, sd = 1, n.start =1000 )
  d[i] <- var(x)
}

summary(d)

sample.b <- c()
for (rho in rho.list) {
  n=60
  ac = ARMAacf(ar=rho , lag.max = 60)
  s=0
  for (j in (1:59)) {
    s = s+ (1-j/60)*ac[j]
  }
  true.v = 1
  a = (1-(2*s)/(n-1))
  e = (a*true.v/n0)*1.35
  sample.b <- c(sample.b,e)  
}
sav <- sample.b

sample.b - sav


# compare with weighted mean estimation -----------------------------------

for (l in c(1:length(var.list))) {
  n0 = 60
  prob.outliers = 0.1
  sigma <- data.frame(matrix(NA, nrow = nb.sim, ncol = 3))
  sigma.loes <- list()
  for (i in 1:nb.sim) {
    set.seed(i)
    # x = SimulatedSeries(n = n0, P = 5, prob.outliers = 0.1, size.outliers = 3, rho = 0, theta = 0)$Y
    x = SimulatedSeries(n = n0, P = 2, prob.outliers = 0.1, size.outliers = 3, rho = 0, theta = 0)$Y
    sigma[i, 1] <- robustbase::scaleTau2(x)
    sigma[i, 2] <- loess.sd(x, alpha = 0.5)[30]
    sigma[i, 3] <- sd(x)
    
  }
 
}


# compare moving window vs loess:---------------------------
nb.sim = 10000
list.sd = seq(0.2,0.8, 0.2)
length.month1 = c(31,28,31,30,31,30,31,31,30,31,30,31)
L = 365
n = 365
t = c(1:n)
std.t = list.sd[4]*2*sin(2*pi*t/L-pi/2)/2 +1
std = std.t[seq(1,365,length.out = 12)]
res1 <- data.frame(matrix(NA, nrow = nb.sim, ncol = 365))
res2 <- data.frame(matrix(NA, nrow = nb.sim, ncol = 365))
t0 <- as.Date("2021-01-01", format = "%Y-%m-%d" )
t.year <- seq(t0, t0+364,1)
for (i in c(1:nb.sim)) {
  set.seed(i)
  sim.ar <- rnorm(n,0,1)
  # sim.ar <- simulate.series.2(mean.1 = 0, sigma.1 = std,
  #                             N = n,  arma.model = c(0, 0), length.month = length.month1)
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
                                          length.wind = 60, loes.method = "gaussian")
  res1[i,] <- std.est
  res2[i,] <- std.est1
  
}

s1 <- colMeans(res1)

s2 <- colMeans(res2)

lines(t.year, s1)
plot(t.year, s2)

sd.c = (1-1/n)
dat = data.frame( t = t.year, true = rep(sd.c,n), MW.ScaleTau = s1, loess = s2)
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


