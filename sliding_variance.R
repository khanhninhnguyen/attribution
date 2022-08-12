# sliding window variance 

RobEstiMonthlyVariance.diff.S30 <- function(Y, name.var, alpha){
  Kmonth <- length(unique(Y$month))
  Y$month <- factor(Y$month, ordered = TRUE)
  z <- stats::aggregate(Y[name.var], by = Y[c("month", "year")], FUN = function(x) { re = diff(na.omit(x)); return(re)})
  sigma.est1 <- sapply(1:Kmonth,function(i) {
    e <- subset(z,z$month==levels(z$month)[i])
    ee <- unlist(e[name.var])
    if(length(ee) <= 29){
      NA
    }else{
      robustbase::scaleTau2(ee)/sqrt(2-2*alpha) # change for AR(1), if white noise alspha = 0
    }
  })
  sigma.est <- rep(NA,12)
  sigma.est[as.numeric(levels(z$month))] <-  sigma.est1
  return(sigma.est)
}
 
RobEstiSlidingVariance.S<- function(Y, name.var, alpha){# require date in the dataset 
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = nrow(Y1)
  sigma.est1 = rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-30,1)
    end = min(n, i+30)
    x.i = x[begin:end]
    x.idiff = diff(x.i)
    thre = 30
    if(i < 30|i>(n-30)){thre = 15}
    if(length(na.omit(x.idiff)) <= thre){
      sd <- NA
    }else{
      sd <- robustbase::scaleTau2(x.idiff, na.rm = TRUE )/sqrt(2-2*alpha) # change for AR(1), if white noise alspha = 0
    }
    sigma.est1[i] <- sd
  }
  # linear regression of the variance for gaps  MAYBE REPLACE BY INTERPOLATION FUNCTION
  s = sigma.est1
  if (sum(is.na(s)) != 0){
    na.x = as.numeric(is.na(s[-length(s)]))
    diff.nax = diff(na.x)
    start.na = which(diff.nax==1)
    for (k in c(1:length(start.na))) {
      n.na = min(which(diff.nax[start.na[k]:(n-1)] == -1))
      x1 = s[start.na[k]]
      x2 = s[start.na[k] +  n.na]
      slope = (x2-x1)/n.na
      s[c((start.na[k]+1):(start.na[k]+n.na-1))] <- s[start.na[k]]+c(1:(n.na-1))*slope
    }
  }
  
  
 
  sigma.est = sigma.est1[which(Y1$date %in% Y$date)]
  return(sigma.est)
}

