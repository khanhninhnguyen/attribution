# test the screening methods 
source(paste0(path_code_att,"support_test_screening_methods.R"))
source(paste0(path_code_att,"support_screening.R"))

gaps.list <- c(5,10,15)
k = 1
n0 = 60
nb.sim = 10000

# AR(1) -------------------------------------------------------------------
res = data.frame(matrix(NA, ncol = 4, nrow = nb.sim))

for (i in c(1:nb.sim)) {
  sim.ar <- arima.sim(model = list(ar = 0.3), n = n0, sd = 1)
  ind.out = 3*rbinom(n0, 1, gaps.list[k]/100)
  sim.ar <- sim.ar + ind.out*sign(sim.ar)
  scr1 = screen.out(sim.ar)
  scr2 = screen.out(diff(sim.ar))
  
  nb.rm1 = length(scr1$point.rm)
  nb.rm2 = length(scr2$point.rm)
  
  list.out = which(ind.out!=0)
  list.scr1 = scr1$point.rm
  list.scr2 = scr2$point.rm
  for (j in c(1:length(list.scr2))){
    ind = list.scr2[j]
    if(ind<10){ 
      subset <-  sim.ar[c(1 : (ind+10))]
    }else if (ind>(n0-10)){
      sim.ar[c((ind-10): (ind))]
    }else{
      subset <- sim.ar[c((ind-10): (ind+10))]
    }
    dist = abs(sim.ar[c(ind:(ind+1))] - mean(subset, na.rm = TRUE))
    ind.n = which.max(dist)
    if(ind.n != 1){
      list.scr2[j] <- list.scr2[j]+1
    }
  }
  res[i,] <- c(nb.rm1, nb.rm2, length(which(list.scr1 %in% list.out == TRUE)), length(which(list.scr2 %in% list.out == TRUE))) 
}
colnames(res) <- c("Normal", "Difference", "Normal.t", "Difference.t")
boxplot(res[,c(1:2)]/n0, main = "Percentage", xlab="Percentage of outlier")
abline(h = 0.05)

per.true = res[,c(3:4)]/res[,c(1:2)]
boxplot(per.true, main = "Percentage", xlab="Percentage of corrected outlier")
abline(h = 0.95)


# Seasonal bias -----------------------------------------------------------
res = data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
L = 365
t = c(1:720)
A = 2
s.bias = A*cos(2*pi*t/L)
for (i in c(1:nb.sim)) {
  sim.ar <- arima.sim(model = list(ar = 0.5), n = n0, sd = 1) + s.bias
  ind.out = 3*rbinom(n0, 1, gaps.list[k]/100)
  sim.ar <- sim.ar + ind.out*sign(sim.ar)
  scr1 = screen.out(sim.ar)
  scr2 = screen.out(diff(sim.ar))
  
  nb.rm1 = length(scr1$point.rm)
  nb.rm2 = length(scr2$point.rm)
  
  list.out = which(ind.out!=0)
  list.scr1 = scr1$point.rm
  list.scr2 = scr2$point.rm
  for (j in c(1:length(list.scr2))){
    ind = list.scr2[j]
    if(ind<10){ 
      subset <-  sim.ar[c(1 : (ind+10))]
    }else if (ind>(n0-10)){
      sim.ar[c((ind-10): (ind))]
    }else{
      subset <- sim.ar[c((ind-10): (ind+10))]
    }
    dist = abs(sim.ar[c(ind:(ind+1))] - mean(subset, na.rm = TRUE))
    ind.n = which.max(dist)
    if(ind.n != 1){
      list.scr2[j] <- list.scr2[j]+1
    }
  }
  res[i,] <- c(nb.rm1, nb.rm2, length(which(list.scr1 %in% list.out == TRUE)), length(which(list.scr2 %in% list.out == TRUE))) 
}
colnames(res) <- c("Normal", "Difference", "Normal.t", "Difference.t")
boxplot(res[,c(1:2)]/n0, main = "Percentage", xlab="Percentage of detected outlier")
abline(h = 0.05)

per.true = res[,c(3:4)]/res[,c(1:2)]
boxplot(per.true, main = "Percentage", xlab="Percentage of correct detected outlier")
abline(h = 0.95)

# Monthly variance --------------------------------------------------------

res = data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
L = 365
t = c(1:730)
std.t = 1.2*sin(2*pi*t/L-pi/2)/2 +1
std = std.t[seq(1,365,length.out = 12)]
length.month1 = c(31,28,31,30,31,30,31,31,30,31,30,31)

for (i in c(1:nb.sim)) {
  sim.ar <- simulate.series.2(mean.1 = 0, sigma.1 = std,
                              N = n0,  arma.model = c(0.5, 0), length.month = length.month1)
  ind.out = 3*rbinom(n0, 1, gaps.list[k]/100)
  sim.ar <- sim.ar + ind.out*sign(sim.ar)
  scr1 = screen.out(sim.ar)
  scr2 = screen.out(diff(sim.ar))
  
  nb.rm1 = length(scr1$point.rm)
  nb.rm2 = length(scr2$point.rm)
  
  list.out = which(ind.out!=0)
  list.scr1 = scr1$point.rm
  list.scr2 = scr2$point.rm
  for (j in c(1:length(list.scr2))){
    ind = list.scr2[j]
    if(ind<10){ 
      subset <-  sim.ar[c(1 : (ind+10))]
    }else if (ind>(n0-10)){
      sim.ar[c((ind-10): (ind))]
    }else{
      subset <- sim.ar[c((ind-10): (ind+10))]
    }
    dist = abs(sim.ar[c(ind:(ind+1))] - mean(subset, na.rm = TRUE))
    ind.n = which.max(dist)
    if(ind.n != 1){
      list.scr2[j] <- list.scr2[j]+1
    }
  }
  res[i,] <- c(nb.rm1, nb.rm2, length(which(list.scr1 %in% list.out == TRUE)), length(which(list.scr2 %in% list.out == TRUE))) 
}
colnames(res) <- c("Normal", "Difference", "Normal.t", "Difference.t")
boxplot(res[,c(1:2)]/n0, main = "Percentage", xlab="Percentage of detected outlier")
abline(h = 0.05)

per.true = res[,c(3:4)]/res[,c(1:2)]
boxplot(per.true, main = "Percentage", xlab="Percentage of correct detected outlier")
abline(h = 0.95)


# sliding window ----------------------------------------------------------

for (i in c(1:(nrow(da)-60))) {
  ind.test = c(i:(i+60))
  samp = da[ind.test,2]
  list.point.rm <- c()
  if (all(is.na(samp$x))){
    da$scr[i] <- NA
  }else{
    y = unlist(samp)
    con = length(na.omit(y))
    if( con >= 30){
      up = screen.out(y)$up
      down =  screen.out(y)$down
      da$up[i+30] <- up
      da$down[i+30] <- down
      x = diff(y)
      x.scr <- screen.out(x)$data
      rm <- screen.out(x)$point.rm
      if(length(rm) > 0){
        da$scr[i] <- 1
        list.point.rm <- c(list.point.rm, (rm+i))
      }else{
        da$scr[i] <- 0
      } 
    }else{
      da$scr[i] <- 2
    }
  }
}
list.rm = list.point.rm[!duplicated(list.point.rm)]
list.rm.n <- as.Date(x = integer(0), origin = "1970-01-01")
for (j in c(1:length(list.rm))){
  ind <- which(da$date == list.rm[j])
  subset <- da[c((ind-2): (ind+2)),]
  dist = abs(subset$x - mean(subset$x))
  list.rm.n <- c(list.rm.n, subset$date[which.max(dist)])
}
# keep only 1%
ind <- which(da$date %in% list.rm.n)
diff.x <- diff(da$x)  
dist.se <- abs(diff.x[ind])
nb.point = round(nrow(dat)/100)
thres = sort(dist.se, decreasing = TRUE)[nb.point]
list.rm.n <- da$date[ind[which(dist.se >= thres)]]

list.p <- da$date[screen.out(da$x)$point.rm]
list.pn <- as.Date(x = integer(0), origin = "1970-01-01")
for (k in c(1:length(list.p))) {
  ind <- which(da$date == list.p[k])
  subset <- da[c((ind-1): (ind+1)),]
  check.na = sum(is.na(subset$x))
  if(check.na >1){
    list.pn <- c(list.pn, list.p[k])
  }
}

sliding <- function(series){
  res <- data.frame(matrix(NA, ncol = 3, nrow = length(series)))
  colnames(res) <- c("scr", "up", "down")
  list.point.rm <- c()
  n0 = length(series)
  
  for (i in c(1:(n0-59))) {
    ind.test = c(i:(i+59))
    samp = series[ind.test]
    if (all(is.na(samp))){
      res$scr[i+30] <- NA
    }else{
      y = samp
      con = length(na.omit(y))
      if( con >= 30){
        x = diff(y)
        rm <- screen.out(x)$point.rm
        if(length(rm) > 0){
          res$scr[i+30] <- 1
          list.point.rm <- c(list.point.rm, (rm+i))
          print(i)
        }else{
          res$scr[i+30] <- 0
        } 
      }else{
        res$scr[i+30] <- 2
      }
    }
  }
  list.rm = list.point.rm[!duplicated(list.point.rm)]
  # this step is to distinguish the first or second point in the difference is outlier - need to be studied more 
  for (j in c(1:length(list.rm))){
    ind = list.rm[j]
    if(ind<10){ 
      subset <- series[c(1 : (ind+10))]
    }else if (ind>(n0-10)){
      subset <- series[c((ind-10): (ind))]
    }else{
      subset <- series[c((ind-10): (ind+10))]
    }
    dist = abs(series[c(ind:(ind+1))] - mean(subset, na.rm = TRUE))
    ind.n = which.max(dist)
    if(ind.n != 1){
      list.rm[j] <- list.rm[j]+1
    }
  }
}

series =  simulate.series.2(mean.1 = 0, sigma.1 = std,
                       N = n0,  arma.model = c(0.5, 0), length.month = length.month1)
plot(x)




# compare with other screening method -------------------------------------
# replace the difference by the mean of absolute difference bw the tested points with the point before and after ---------------
res = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
# Box plot

for (i in c(1:nb.sim)) {
  sim.ar <- arima.sim(model = list(ar = 0.5), n = n0, sd = 1)
  ind.out = 3*rbinom(n0, 1, gaps.list[k]/100)
  sim.ar <- sim.ar + ind.out*sign(sim.ar)
  scr1 = screen.iqr(sim.ar)
  scr2 = screen.diff(sim.ar, dif = 1)
  scr3 = screen.mad(sim.ar)
  
  nb.rm1 = length(scr1$point.rm)
  nb.rm2 = length(scr2$point.rm)
  nb.rm3 = length(scr3$point.rm)
  
  list.out = which(ind.out!=0)
  list.scr1 = scr1$point.rm
  list.scr2 = scr2$point.rm
  list.scr3 = scr3$point.rm
  
  res[i,] <- c(nb.rm1, nb.rm2, nb.rm3,
               length(which(list.scr1 %in% list.out == TRUE)), 
               length(which(list.scr2 %in% list.out == TRUE)),
               length(which(list.scr3 %in% list.out == TRUE))) 
}
colnames(res) <- c("Normal", "Difference","MAD", "Normal.t", "Difference.t", "MAD.t")
boxplot(res[,c(1:3)]/n0, main = "Percentage", xlab="Percentage of outlier")
abline(h = 0.05)

per.true = res[,c(4:6)]/res[,c(1:3)]
boxplot(per.true, main = "Percentage", xlab="Percentage of corrected outlier")
abline(h = 0.95)


# scatter plot: percentage of corrected detection as a function of rho --------
list.rho = seq(0.1,0.9, 0.2)
L = 365
t = c(1:730)
A = 1
s.bias = A*cos(2*pi*t/L)
res.tol <- list()
for (l in c(1:length(list.rho))) {
  rho = list.rho[l]
  res = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
  colnames(res) <- c("Normal", "Difference","MAD", "Normal.t", "Difference.t", "MAD.t")
  for (i in c(1:nb.sim)) {
    sim.ar <- arima.sim(model = list(ar = rho), n = n0, sd = 1) 
    ind.out = 3*rbinom(n0, 1, gaps.list[k]/100)
    sim.ar <- sim.ar + ind.out*sign(sim.ar)
    scr1 = screen.iqr(sim.ar)
    scr2 = screen.diff(sim.ar, dif = 1)
    scr3 = screen.mad(sim.ar)
    
    nb.rm1 = length(scr1$point.rm)
    nb.rm2 = length(scr2$point.rm)
    nb.rm3 = length(scr3$point.rm)
    
    list.out = which(ind.out!=0)
    list.scr1 = scr1$point.rm
    list.scr2 = scr2$point.rm
    list.scr3 = scr3$point.rm
    
    res.i <- c(nb.rm1*100/nb.sim, nb.rm2*100/nb.sim, nb.rm3*100/nb.sim,
                 length(which(list.scr1 %in% list.out == TRUE))/nb.rm1, 
                 length(which(list.scr2 %in% list.out == TRUE))/nb.rm2,
                 length(which(list.scr3 %in% list.out == TRUE))/nb.rm3) 
    res.i[is.na(res.i)] <- 0
    res[i,] <- res.i
  }
  res.tol[[l]] <- res
}

std.all = as.data.frame(t(sapply(1:5, function(x){ sapply(as.data.frame(res.tol[[x]]), sd)})))
mean.all = as.data.frame(t(sapply(1:5, function(x){ sapply(as.data.frame(res.tol[[x]]), mean)})))

plot.d = data.frame(methods = rep(colnames(mean.all[,c(1:3)]), each = 5),
                    mean = unlist(mean.all[,c(1:3)]), 
                    std = unlist(std.all[,c(1:3)]),
                    phi = rep(list.rho, 3))
ggplot(data = plot.d , aes( x = phi, y = mean, col = methods))+
   theme_bw()+geom_line()+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.04)+
  ylab("Per. outlier detection(%)")+
  labs(subtitle = "imposed 5% outliers")


plot.d1 = data.frame(methods = rep(colnames(mean.all[,c(4:6)]), each = 5),
                    mean = unlist(mean.all[,c(4:6)]), 
                    std = unlist(std.all[,c(4:6)]),
                    phi = rep(list.rho, 3))
ggplot(data = plot.d1 , aes( x = phi, y = mean, col = methods))+
  theme_bw()+geom_line()+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.04)+
  ylab("Per. outlier corrected detection(%)")+
  labs(subtitle = "imposed 5% outliers")


# scatter plot bias  ------------------------------------------------------

list.rho = seq(0.1,0.9, 0.2)
L = 365
t = c(1:730)
A = 1
list.bias = seq(0.2,1.6,0.2)
res.tol <- list()
for (l in c(1:length(list.bias))) {
  s.bias = A*cos(2*pi*t/L)
  res = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
  colnames(res) <- c("Normal", "Difference","MAD", "Normal.t", "Difference.t", "MAD.t")
  for (i in c(1:nb.sim)) {
    sim.ar <- arima.sim(model = list(ar = 0.5), n = n0, sd = 1) + s.bias
    ind.out = 3*rbinom(n0, 1, gaps.list[k]/100)
    sim.ar <- sim.ar + ind.out*sign(sim.ar)
    scr1 = screen.iqr(sim.ar)
    scr2 = screen.diff(sim.ar, dif = 1)
    scr3 = screen.mad(sim.ar)
    
    nb.rm1 = length(scr1$point.rm)
    nb.rm2 = length(scr2$point.rm)
    nb.rm3 = length(scr3$point.rm)
    
    list.out = which(ind.out!=0)
    list.scr1 = scr1$point.rm
    list.scr2 = scr2$point.rm
    list.scr3 = scr3$point.rm
    
    res.i <- c(nb.rm1*100/nb.sim, nb.rm2*100/nb.sim, nb.rm3*100/nb.sim,
               length(which(list.scr1 %in% list.out == TRUE))/nb.rm1, 
               length(which(list.scr2 %in% list.out == TRUE))/nb.rm2,
               length(which(list.scr3 %in% list.out == TRUE))/nb.rm3) 
    res.i[is.na(res.i)] <- 0
    res[i,] <- res.i
  }
  res.tol[[l]] <- res
}

std.all = as.data.frame(t(sapply(1:length(list.bias), function(x){ sapply(as.data.frame(res.tol[[x]]), sd)})))
mean.all = as.data.frame(t(sapply(1:length(list.bias), function(x){ sapply(as.data.frame(res.tol[[x]]), mean)})))

plot.d = data.frame(methods = rep(colnames(mean.all[,c(1:3)]), each = length(list.bias)),
                    mean = unlist(mean.all[,c(1:3)]), 
                    std = unlist(std.all[,c(1:3)]),
                    amp = rep(list.bias, 3))
ggplot(data = plot.d , aes( x = amp, y = mean, col = methods))+
  theme_bw()+geom_line()+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.04)+
  ylab("Per. outlier detection(%)")+
  labs(subtitle = "imposed 5% outliers")


plot.d1 = data.frame(methods = rep(colnames(mean.all[,c(4:6)]), each = length(list.bias)),
                     mean = unlist(mean.all[,c(4:6)]), 
                     std = unlist(std.all[,c(4:6)]),
                     amp = rep(list.bias, 3))
ggplot(data = plot.d1 , aes( x = amp, y = mean, col = methods))+
  theme_bw()+geom_line()+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.04)+
  ylab("Per. outlier corrected detection(%)")+
  labs(subtitle = "imposed 5% outliers")

library("car")
qqPlot(sim.ar, id = FALSE)

# when data is monthly variance + AR --------------------------------------

list.rho = seq(0.1,0.9, 0.2)
L = 365
t = c(1:730)
list.sd = seq(0.2,0.8, 0.2)
res.tol <- list()
n0=730
length.month1 = c(31,28,31,30,31,30,31,31,30,31,30,31)
for (l in c(1:length(list.sd))) {
  std.t = list.sd[l]*2*sin(2*pi*t/L-pi/2)/2 +1
  std = std.t[seq(1,365,length.out = 12)]
  res = data.frame(matrix(NA, ncol = 6, nrow = nb.sim))
  colnames(res) <- c("Normal", "Difference","MAD", "Normal.t", "Difference.t", "MAD.t")
  for (i in c(1:nb.sim)) {
    sim.ar <- simulate.series.2(mean.1 = 0, sigma.1 = std,
                                N = n0,  arma.model = c(0.5, 0), length.month = length.month1)
    ind.out = 3*rbinom(n0, 1, gaps.list[k]/100)
    sim.ar <- sim.ar + ind.out*sign(sim.ar)
    scr1 = screen.iqr(sim.ar)
    scr2 = screen.diff(sim.ar, dif = 1)
    scr3 = screen.mad(sim.ar)
    
    nb.rm1 = length(scr1$point.rm)
    nb.rm2 = length(scr2$point.rm)
    nb.rm3 = length(scr3$point.rm)
    
    list.out = which(ind.out!=0)
    list.scr1 = scr1$point.rm
    list.scr2 = scr2$point.rm
    list.scr3 = scr3$point.rm
    
    res.i <- c(nb.rm1*100/nb.sim, nb.rm2*100/nb.sim, nb.rm3*100/nb.sim,
               length(which(list.scr1 %in% list.out == TRUE))/nb.rm1, 
               length(which(list.scr2 %in% list.out == TRUE))/nb.rm2,
               length(which(list.scr3 %in% list.out == TRUE))/nb.rm3) 
    res.i[is.na(res.i)] <- 0
    res[i,] <- res.i
  }
  res.tol[[l]] <- res
}

std.all = as.data.frame(t(sapply(1:length(list.sd), function(x){ sapply(as.data.frame(res.tol[[x]]), sd)})))
mean.all = as.data.frame(t(sapply(1:length(list.sd), function(x){ sapply(as.data.frame(res.tol[[x]]), mean)})))

plot.d = data.frame(methods = rep(colnames(mean.all[,c(1:3)]), each = length(list.sd)),
                    mean = unlist(mean.all[,c(1:3)]), 
                    std = unlist(std.all[,c(1:3)]),
                    amp = rep(list.sd, 3))
ggplot(data = plot.d , aes( x = amp, y = mean, col = methods))+
  theme_bw()+geom_line()+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.04)+
  ylab("Per. outlier detection(%)")+
  labs(subtitle = "imposed 5% outliers")


plot.d1 = data.frame(methods = rep(colnames(mean.all[,c(4:6)]), each = length(list.sd)),
                     mean = unlist(mean.all[,c(4:6)]), 
                     std = unlist(std.all[,c(4:6)]),
                     amp = rep(list.sd, 3))
ggplot(data = plot.d1 , aes( x = amp, y = mean, col = methods))+
  theme_bw()+geom_line()+
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.04)+
  ylab("Per. outlier corrected detection(%)")+
  labs(subtitle = "imposed 5% outliers")

# simplest simulation -----------------------------------------------------
res = data.frame(matrix(NA, ncol = 5, nrow = nb.sim))
for (i in c(1:nb.sim)) {
  # sim.ar <- rnorm(n0, sd = 1, mean = 0)
  sim.ar <- arima.sim(model = list(ar = 0.8), n = n0, sd = 1, n.start = 1000)
  
  # ind.out = 10*rbinom(n0, 1, gaps.list[k]/100)
  # sim.ar <- sim.ar + ind.out*sign(sim.ar)
  scr1 = screen.iqr(sim.ar)
  up1 = 2.69792/(sqrt(1-0.5**2))
  down1 = -2.69792/(sqrt(1-0.5**2))
  removed1 <- which(sim.ar<down1 | sim.ar>up1)
  
  # scr2 = screen.sca(sim.ar1)
  # up2 = 2.69792/(sqrt(1-0.5**2))
  # down2 = -2.69792/(sqrt(1-0.5**2))
  # removed2 <- which(sim.ar1<down2 | sim.ar1>up2)
  
  scr2 = screen.mad(sim.ar)
  scr3 = screen.sca(sim.ar)
  scr4 = screen.qn(sim.ar)
  
  nb.rm1 = length(removed1 )
  nb.rm2 = length(scr1$point.rm )
  nb.rm3 = length(scr2$point.rm  )
  nb.rm4 = length(scr3$point.rm )
  nb.rm5 = length(scr4$point.rm )
  
  # list.out = which(ind.out!=0)
  res[i,] <- c(true = nb.rm1, IQR = nb.rm2, MAD = nb.rm3, Sca = nb.rm4, Qn = nb.rm5)
}
summary(res/n0)
boxplot(res/n0)
abline(h = 0.007)

plot.d = data.frame(method = rep(c("true", "IQR", "MAD", "Sca", "Qn"), each = nb.sim), variable = unlist(res)/n0)

plot <- ggplot(data = plot.d, aes(y = variable, x = method)) + 
  theme_bw()+
  geom_boxplot(fatten = NULL) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  geom_hline(yintercept = 0.007)+
  ylab("outlier percentage")+
  labs(subtitle = "Normal distribution, mean = 0, sd =1, n = 100, phi = 0.5")
print(plot)

library(tseries)
adf.test(sim.ar)
set.seed(123)
x <- arima.sim(model = list(ar = 0.8, ma = -0.2), n = 1000)
set.seed(122)
x1 <- arima.sim(model = list(ar = 0.8, ma = -0.2), n = 1000)
set.seed(122)
sim.ar2 <-  arima.sim(model = list(ar = 0.8, ma = -0.2), n = 1000)














# normal data, iterative screening, no normalization *

Y = as.data.frame(dat[[30]]$bef)
Y <- data.frame(date = Y1$date, gps.gps = rnorm(1000, 0,1))
Y1 <- tidyr::complete(Y, date = seq(min(Y$date), min(Y$date)+999, by = "day"))
a <- screen.O(Y, name.var = "gps.gps", method = 'sigma')

res <- rep(NA, 1000)
for(r in c(1:1000)){
  set.seed(r)
  Y <- data.frame(date = Y1$date, gps.gps = rnorm(1000, 0,1))
  a <- screen.O(Y, name.var = "gps.gps", method = 'sigma', iter = 1)
  # a <- screen.O1(Y$gps.gps, method = "def")
  res[r] <- length( a$point.rm)
  print(r)
}
summary(res)
hist(res)

# looking at a specific case, where MAD detect more and more modulate 











