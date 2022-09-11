rm(list = ls())
path_results=paste0("/home/knguyen/Documents/PhD/","Results/")
library("ggplot2")
source(paste0("/home/knguyen/Documents/PhD/","Code/attribution/","UsedFunctions.R"))

nb.sim = 10000
res <- list()
length.list = c(10, 20, 30, 60, 120, 360, 720)
outlier.list = seq(0.00, 0.1, 0.01)
outlier.list = seq(0.00, 0.1, 0.01)
rho.list = seq(0.0, 0.9, 0.15)
var.list <- rho.list
size.outliers = 3

choose_model <- function(x){
  if(x == 1){ P = 5 } # fixed value of outlier
  else if(x == 2){ P = 3 } # normal outlier 
  else if(x == 3){ P = 4 } # normal overlayed outlier 
  else if(x == 4){ P = 7 } # normal overlayed outlier + arma data
  return(P)
}
model.out = 4
P = choose_model(model.out)
for (l in c(1:length(var.list))) {
  n0 = 1000
  prob.outliers = 0
  sigma = data.frame(matrix(NA, ncol =5, nrow = nb.sim))
  rho = var.list[l]
  theta = 0
  off.set = 0
  ratio = sqrt((1 + 2*theta*rho + theta**2)/(1 - rho**2))
  for (i in 1:nb.sim) {
    set.seed(i)
    x = SimulatedSeries(n = n0, P = P, prob.outliers = prob.outliers, size.outliers = size.outliers, rho = rho, theta = theta)$Y
    # x <- rnorm(n = n0, mean = 0, sd =1)
    x[(floor(n0/2)+1):n0] <- x[(floor(n0/2)+1):n0]+off.set
    x[(floor(n0*3/4)+1):n0] <- x[(floor(n0*3/4)+1):n0]+off.set
    # alpha.e = arima(x, order = c(1,0,0), method = "ML")$coef[1]
    # a = diff(x)
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
lab.x = "Rho"

iqr = as.data.frame(res$sdt)
colnames(iqr) <- as.character(var.list)
true.bias <- iqr[1,1]*ratio
iqr <- rbind(iqr, c(true.bias))
iqr$est <- c("SD", "MAD", "Sn", "Qn","ScaleTAu", "true")
a = reshape2::melt(iqr, id ="est")
data.plot <- data.frame()
jpeg(paste0(path_results,"attribution/variances/IQR.model",model.out, "rho = ", rho, " theta = ", theta, "out.size = ",  size.outliers, "offset = ", off.set, prob.outliers,"1000.jpeg" ),
     width = 3000, height = 1500,res = 300)
ggplot(data = a, aes(x = variable, y = value, col = est)) + 
  geom_point()+
  theme_bw()+
  xlab(lab.x)+
  scale_y_continuous( name = "IQR", limits = c(0, ceiling(max(a$value)*10)/10), breaks = seq(0, ceiling(max(a$value)*10)/10, 0.1))+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()
# BIAS
iqr = as.data.frame(res$bias)
colnames(iqr) <- as.character(var.list)
true.bias <- iqr[1,1]*ratio
iqr <- rbind(iqr, c(true.bias))
iqr$est <- c("SD", "MAD", "Sn", "Qn","ScaleTAu", "true")
a = reshape2::melt(iqr, id ="est")
lim1 = floor(min(a$value)*10)/10
lim2 = ceiling(max(a$value)*10)/10
data.plot <- data.frame()
jpeg(paste0(path_results,"attribution/variances/bias.model", model.out, "rho = ", rho, " theta = ", theta, "out.size = ",  size.outliers, "offset = ", off.set, prob.outliers, "1000.jpeg" ),
     width = 3000, height = 1500,res = 300)
ggplot(data = a, aes(x = variable, y = value, col = est)) + 
  geom_point()+
  theme_bw()+
  xlab(lab.x)+
  scale_y_continuous( name = "Bias", limits = c(lim1, lim2), breaks = seq(lim1, lim2, 0.1))+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15,face="bold"))
dev.off()




