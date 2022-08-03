# test the screening methods 
source(paste0(path_code_att,"support_test_screening_methods.R"))

gaps.list <- c(5,10,15)
k = 1
n0 = 720
screen.out <- function(x) {
  y = na.omit(x)
  Q1 <- quantile(x, .25, na.rm = TRUE)
  Q3 <- quantile(x, .75, na.rm = TRUE)
  IQR <- IQR(x, na.rm = TRUE)
  up = Q3+IQR*1.5
  down = Q1-IQR*1.5
  removed <- which(x<down | x>up)
  if (length(removed) >0){
    x.screened = x[-removed]
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = removed, up = up, down = down))
}
nb.sim = 1000
res = data.frame(matrix(NA, ncol = 4, nrow = nb.sim))

# AR(1) -------------------------------------------------------------------

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

x = c(1:365)
list.std = list()
list.range = seq(0, 0.9, 0.1)*2
for (k in 1:length(list.range)) {
  first.mean = list.range[k] * sin(2*pi*x/L-pi/2)/2 +1
  second.mean = list.range[k] * sin(2*pi*x/L+pi/2)/2 +1 
  list.std[[2*k-1]] <- first.mean[c(seq(1,length(first.mean), 35), 365)]
  list.std[[2*k]] <- second.mean[c(seq(1,length(second.mean), 35), 365)]
}

res = data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
L = 365
t = c(1:720)
std.t = 1.2*sin(2*pi*t/L-pi/2)/2 +1
std = std.t[seq(1,365,length.out = 12)]
length.month1 = c(31,28,31,30,31,30,31,31,30,31,30,31)

for (i in c(1:nb.sim)) {
  sim.ar <- simulate.series.2(mean.1 = 0, sigma.1 = std,
                              N = n0, c(0.5, 0), length.month = length.month1)
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

for (i in c(1:(nrow(da)-59))) {
  ind.test = c(i:(i+59))
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

sliding <- function(x){
  res <- data.frame(matrix(NA, ncol = 3, nrow = length(x)))
  colnames(res) <- c("scr", "up", "down")
  list.point.rm <- c()
  
  for (i in c(1:(length(x)-59))) {
    ind.test = c(i:(i+59))
    samp = x[ind.test]
    if (all(is.na(samp))){
      res$scr[i+30] <- NA
    }else{
      y = samp
      con = length(na.omit(y))
      if( con >= 30){
        up = screen.out(y)$up
        down =  screen.out(y)$down
        res$up[i+30] <- up
        res$down[i+30] <- down
        x = diff(y)
        x.scr <- screen.out(x)$data
        rm <- screen.out(x)$point.rm
        if(length(rm) > 0){
          res$scr[i] <- 1
          list.point.rm <- c(list.point.rm, (rm+i))
        }else{
          res$scr[i+30] <- 0
        } 
      }else{
        res$scr[i+30] <- 2
      }
    }
  }
  
  ind <- which(da$date %in% list.rm.n)
  list.rm = list.point.rm[!duplicated(list.point.rm)]
  
}






