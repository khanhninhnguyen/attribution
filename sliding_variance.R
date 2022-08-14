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
 
RobEstiSlidingVariance.S <- function(Y, name.var, alpha){# require date in the dataset, return std at time t
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
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  sigma.est = s[which(Y1$date %in% Y$date)]
  return(sigma.est)
}

sliding.median <- function(Y, name.var){
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = length(x)
  slide.med <- rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-30,1)
    end = min(n, i+30)
    x.i = x[begin:end]
    slide.med[i] <- median(x.i, na.rm = TRUE)
  }
  slide.med = slide.med[which(Y1$date %in% Y$date)]
  return(slide.med)
}
  
sliding.IQR <- function(Y, name.var){
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = length(x)
  slide.med <- rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-30,1)
    end = min(n, i+30)
    x.i = x[begin:end]
    slide.med[i] <- median(x.i, na.rm = TRUE)
  }
  slide.med = slide.med[which(Y1$date %in% Y$date)]
  return(slide.med)
}










