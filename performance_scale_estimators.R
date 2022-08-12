rm(list = ls())
path_results=paste0("/home/knguyen/Documents/PhD/","Results/")
library("ggplot2")

nb.sim = 10000
for (l in c(30,60,90,120)) {
  n0 = l
  sigma = data.frame(matrix(NA, ncol =5, nrow = nb.sim))
  alpha = 0.5
  off.set = 0
  for (i in 1:nb.sim) {
    set.seed(i)
    if(alpha !=0){
      x = arima.sim(model = list(  ar = alpha ), n = n0, sd = sqrt(1-alpha**2))
    }else{
      x = rnorm(n = n0, mean = 0, sd = 1)
      # x = arima.sim(model = list( order = c(0, 0, 0)), n = n0, sd = sqrt(1-alpha**2))
    }  
  
    x[(floor(n0/2)+1):n0] <- x[(floor(n0/2)+1):n0]+off.set
    x[(floor(n0*3/4)+1):n0] <- x[(floor(n0*3/4)+1):n0]+off.set
    alpha.e = arima(x, order = c(1,0,0), method = "ML")$coef[1]
    a = diff(x)
    # sigma[i,1] = sd(x)
    # sigma[i,2] <- mad(x)
    # sigma[i,3] <- robustbase::Sn(x)
    # sigma[i,4] <- robustbase::Qn(x)
    # sigma[i,5] <- robustbase::scaleTau2(x)
    sigma[i,1] = sd(a)/sqrt(2-2*alpha.e)
    sigma[i,2] <- mad(a)/sqrt(2-2*alpha.e)
    sigma[i,3] <- robustbase::Sn(a)/sqrt(2-2*alpha.e)
    sigma[i,4] <- robustbase::Qn(a)/sqrt(2-2*alpha.e)
    sigma[i,5] <- robustbase::scaleTau2(a)/sqrt(2-2*alpha.e)
  }
  colnames(sigma) <- c("classic", "MAD", "Sn", "Qn","ScaleTAu")
  b = data.frame(estimator = rep(colnames(sigma), each = nb.sim), std = unlist(sigma))
  p <- ggplot(b, aes(x = estimator, y = std))+ 
    geom_boxplot()+theme_bw()+
    labs(subtitle = paste0("phi = ", alpha, ", n = ", n0, ", offset = ", off.set))+
    geom_hline(yintercept = 1)+
    theme(axis.text = element_text(size = 15))
  if(n0==30){p <- p +  ylim(c(0,2))
  }else{p <- p +  ylim(c(0.5,1.5))}

  jpeg(paste0(path_results,"attribution/variances/sim.e", "phi = ", alpha, ", n = ", n0, ", offset = ", off.set, ".jpeg"),
       width = 3000, height = 1800,res = 300) # change name
  print(p)
  dev.off()
}



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

