rm(list = ls())
path_results=paste0("/home/knguyen/Documents/PhD/","Results/")
library("ggplot2")

nb.sim = 10000
res <- list()
length.list = c(10,20,30,60,90,120,240,480,720)
r <- c()
for (l in c(1:length(length.list))) {
  n0 = length.list[l]
  sigma = data.frame(matrix(NA, ncol =5, nrow = nb.sim))
  alpha = 0
  off.set = 0
  for (i in 1:nb.sim) {
    set.seed(i)
    if(alpha !=0){
      x = arima.sim(model = list(  ar = alpha ), n = n0, sd = sqrt(1-alpha**2))
    }else{
      x = rnorm(n = n0, mean = 0, sd = 1)
      # adding outlier 
      set.seed(i)
      out = 5*sample(c(-1,1,0),n0,replace=T, prob = c(0.005, 0.005, 0.99))
      out.ind = which(out!=0)
      x[out.ind] <- x[out.ind]+out[out.ind]
      # set.seed(2)
      # ind.out = rbinom(n0, 1, 0.01)
      # x[ind.out] <- -(max(x)+2)
      # x = arima.sim(model = list( order = c(0, 0, 0)), n = n0, sd = sqrt(1-alpha**2))
    }  
    r <- c(r,length(out.ind))
    x[(floor(n0/2)+1):n0] <- x[(floor(n0/2)+1):n0]+off.set
    x[(floor(n0*3/4)+1):n0] <- x[(floor(n0*3/4)+1):n0]+off.set
    alpha.e = arima(x, order = c(1,0,0), method = "ML")$coef[1]
    a = diff(x)
    sigma[i,1] = sd(x)
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
  bias <- sapply(c(1:5), function(x) median(sigma[,x])-1)
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


# Box plot to easily compare different estimators ------
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

# Scatter plot ------
iqr = as.data.frame(res$sdt)
colnames(iqr) <- as.character(length.list)
iqr$est <- c("SD", "MAD", "Sn", "Qn","ScaleTAu")
a = reshape2::melt(iqr, id ="est")
data.plot <- data.frame()

ggplot(data = a, aes(x = variable, y = value, col = est)) + 
  geom_point()+
  theme_bw()+
  xlab("Sample size")+
  ylab("IQR")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))

iqr = as.data.frame(res$bias)
colnames(iqr) <- as.character(length.list)
iqr$est <- c("SD", "MAD", "Sn", "Qn","ScaleTAu")
a = reshape2::melt(iqr, id ="est")
data.plot <- data.frame()

ggplot(data = a, aes(x = variable, y = value, col = est)) + 
  geom_point()+
  theme_bw()+
  xlab("Sample size")+
  ylab("Bias")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))






# compare different ways to simulate outliers 
n0 = 100
alpha = 0
nb.sim = 10000
off.set = 0
sigma = data.frame(matrix(NA, ncol =10, nrow = nb.sim))
for (i in 1:nb.sim) {
  set.seed(i)
  if(alpha !=0){
    x = arima.sim(model = list(  ar = alpha ), n = n0, sd = sqrt(1-alpha**2))
  }else{
    x = rnorm(n = n0, mean = 0, sd = 1)
    # adding outlier 1
    set.seed(i)
    out = 5*sample(c(-1,1,0),n0,replace=T, prob = c(0.005, 0.005, 0.99))
    out.ind = which(out!=0)
    x1 <- x
    x1[out.ind] <- out[out.ind]
    # adding outlier 2
    set.seed(i)
    ind.out = rbinom(n0, 1, 0.01)
    o = sample(c(5,-5), length(which(ind.out!=0)), replace = T)
    x2 <- x
    x2[ind.out] <- o
  }  
  x[(floor(n0/2)+1):n0] <- x[(floor(n0/2)+1):n0]+off.set
  x[(floor(n0*3/4)+1):n0] <- x[(floor(n0*3/4)+1):n0]+off.set
  sigma[i,1] = sd(x1)
  sigma[i,2] <- mad(x1)
  sigma[i,3] <- robustbase::Sn(x1)
  sigma[i,4] <- robustbase::Qn(x1)
  sigma[i,5] <- robustbase::scaleTau2(x1)
  sigma[i,6] = sd(x2)
  sigma[i,7] <- mad(x2)
  sigma[i,8] <- robustbase::Sn(x2)
  sigma[i,9] <- robustbase::Qn(x2)
  sigma[i,10] <- robustbase::scaleTau2(x2)
}

list.est = c("sd", "mad", "Sn", "Qn", "ScaleTau")
colnames(sigma) <- c(paste0(list.est, ".1"), paste0(list.est, ".2"))
sigma$iter <- c(1:10000)
a = reshape2::melt(sigma, id = "iter")
ggplot(data = a, aes(x = iter, y = value, col = variable)) + 
  geom_point()+
  theme_bw()+
  xlab("Sample size")+
  ylab("Bias")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))





