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
 
RobEstiSlidingVariance.S.diff <- function(Y, name.var, alpha){# require date in the dataset, return std at time t
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
  if (sum(is.na(s)) != 0 & sum(is.na(s)) != length(s)){
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    if(na.ind[1] == 1){
      ind.stop = which(is.na(s)==FALSE)[1]-1
      na.ind <- na.ind[-c(1:ind.stop)]
    }else if (is.na(s[n]) == 1){
      m = which(is.na(s)==FALSE)
      ind.start = m[length(m)]
      na.ind <- na.ind[-which(na.ind %in% c(ind.start:n))]
    }
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  sigma.est = s[which(Y1$date %in% Y$date)]
  return(sigma.est)
}

RobEstiSlidingVariance.S <- function(Y, name.var, alpha, estimator, length.wind){# require date in the dataset, return std at time t
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = nrow(Y1)
  sigma.est1 = rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-(length.wind-1),1)
    end = min(n, i+length.wind)
    x.i = x[begin:end]
    x.idiff = (x.i)
    thre = 30
    if(i < 30|i>(n-30)){thre = 16}
    if(length(na.omit(x.idiff)) <= thre){
      sd <- NA
    }else{
      sd <- my.estimator(estimator = estimator, x.idiff)
    }
    sigma.est1[i] <- sd
  }
  # linear regression of the variance for gaps  MAYBE REPLACE BY INTERPOLATION FUNCTION
  s = sigma.est1
  if (sum(is.na(s)) != 0 & sum(is.na(s)) != length(s)){
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    if(na.ind[1] == 1){
      ind.stop = which(is.na(s)==FALSE)[1]-1
      na.ind <- na.ind[-c(1:ind.stop)]
    }else if (is.na(s[n]) == 1){
      m = which(is.na(s)==FALSE)
      ind.start = m[length(m)]
      na.ind <- na.ind[-which(na.ind %in% c(ind.start:n))]
    }
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  sigma.est = s[which(Y1$date %in% Y$date)]
  return(sigma.est)
}

RobEstiSlidingVariance.WLS <- function(Y, name.var, alpha, length.wind, loes.method){# require date in the dataset, return std at time t
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = nrow(Y1)
  slide.var = rep(NA, n)
  times=1:n
  residus = x - median(x, na.rm = TRUE)
  lissage=loess(residus^2~times,degree=1, span = length.wind/n, normalize = FALSE, na.action = na.exclude, family = loes.method)
  slide.var[which(is.na(x) == FALSE)] <- (na.omit(lissage$fitted))
  # linear regression of the variance for gaps  MAYBE REPLACE BY INTERPOLATION FUNCTION
  s = slide.var
  if (sum(is.na(s)) != 0 & sum(is.na(s)) != length(s)){
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    if(na.ind[1] == 1){
      ind.stop = which(is.na(s)==FALSE)[1]-1
      na.ind <- na.ind[-c(1:ind.stop)]
    }else if (is.na(s[n]) == 1){
      m = which(is.na(s)==FALSE)
      ind.start = m[length(m)]
      na.ind <- na.ind[-which(na.ind %in% c(ind.start:n))]
    }
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  var.t = sapply(s, function(x) ifelse(x>0, x, abs(x)))
  std.t <- sqrt(var.t)
  sigma.est = std.t[which(Y1$date %in% Y$date)]
  return(sigma.est)
}

RobEstiSlidingVariance.log.WLS <- function(Y, name.var, alpha, length.wind){# require date in the dataset, return std at time t
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = nrow(Y1)
  slide.var = rep(NA, n)
  times=1:n
  residus = x - median(x, na.rm = TRUE)
  adding = 0
  y = log(residus^2 + adding)
  y[which(y == -Inf)] <- NA
  lissage=loess(y~times,degree=1, span = length.wind/n, normalize = FALSE, na.action = na.exclude)
  slide.var[which(is.na(y ) == FALSE)] <- (na.omit(lissage$fitted))
  # linear regression of the variance for gaps  MAYBE REPLACE BY INTERPOLATION FUNCTION
  s = slide.var
  if (sum(is.na(s)) != 0 & sum(is.na(s)) != length(s)){
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    if(na.ind[1] == 1){
      ind.stop = which(is.na(s)==FALSE)[1]-1
      na.ind <- na.ind[-c(1:ind.stop)]
    }else if (is.na(s[n]) == 1){
      m = which(is.na(s)==FALSE)
      ind.start = m[length(m)]
      na.ind <- na.ind[-which(na.ind %in% c(ind.start:n))]
    }
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  
  # var.t = sapply(s, function(x) ifelse(x>1, x, 1))
  std.t <- sqrt(exp(s)-adding)
  sigma.est = std.t[which(Y1$date %in% Y$date)]
  return( sigma.est)
}

sliding.median <- function(Y, name.var, length.wind){
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = length(x)
  slide.med <- rep(NA, n)
  l = rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-(length.wind-1),1)
    end = min(n, i+length.wind)
    x.i = x[begin:end]
    thre = 30
    if(i < 30|i>(n-30)){thre = 16}
    if(length(na.omit(x.i))> thre){
      slide.med[i] <- median(x.i, na.rm = TRUE)
      l[i] <- length(na.omit(x.i))
    }
  }
  slide.med = slide.med[which(Y1$date %in% Y$date)]
  # l = l[which(Y1$date %in% Y$date)]
  return(slide.med)
  # return(l)
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

two.step.norm <- function(Y, name.var){
  norm2 = rep(NA, nrow(Y))
  std.t <- RobEstiSlidingVariance.S.diff(Y, name.var, alpha = 0)
  slide.med <- sliding.median(Y, name.var)
  # step 1: normalize by std and median
  norm1 <- unlist((Y[name.var] - slide.med)/std.t)
  # step 2: avoid impact of constant in AR(1) by taking the IQR
  ind.na = which(is.na(norm1))
  norm1 <- na.omit(unlist(norm1,use.names = FALSE))
  if(length(ind.na) ==0){
    norm2<- (norm1 - median(norm1))/IQR(norm1)
  }else{  norm2[-ind.na] <- (norm1 - median(norm1))/IQR(norm1)}
  return(norm2)  
}

sliding.mean <- function(Y, name.var, length.wind){
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = length(x)
  slide.med <- rep(NA, n)
  l = rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-(length.wind-1),1)
    end = min(n, i+length.wind)
    x.i = x[begin:end]
    thre = 30
    if(i < 30|i>(n-30)){thre = 16}
    if(length(na.omit(x.i))> thre){
      slide.med[i] <- mean(x.i, na.rm = TRUE, trim = 0.1)
      l[i] <- length(na.omit(x.i))
    }
  }
  slide.med = slide.med[which(Y1$date %in% Y$date)]
  # l = l[which(Y1$date %in% Y$date)]
  return(slide.med)
  # return(l)
}


one.step.norm <- function(Y, name.var, estimator, length.wind, global.mu, loes, loes.method){
  norm2 = rep(NA, nrow(Y))
  # slide.mean <- sliding.median(Y, name.var, length.wind)
  # slide.mean <- sliding.mean(Y, name.var, length.wind)
  if(global.mu == 1){
    slide.mean = rep(median(unlist(Y[name.var]), na.rm = TRUE), nrow(Y))
  }else{
    slide.mean <- sliding.median(Y, name.var, length.wind = length.wind)
  }
  if (loes == 1){
    # w = 1/(std.t**2) to compute the weighted mean 
    # slide.mean = sum(y*w, na.rm = TRUE)/sum(w, na.rm = TRUE)
    std.t <- RobEstiSlidingVariance.WLS( Y, name.var, alpha = 0, length.wind = length.wind, loes.method = loes.method)
  } else{
    std.t <- RobEstiSlidingVariance.S(Y, name.var, alpha = 0, estimator = estimator, length.wind= length.wind)
  }
  # step 1: normalize by std and median
  norm1 <- unlist((Y[name.var] - slide.mean)/std.t)
  return(list(norm.sig = norm1, mu = slide.mean, sd = std.t))
}

loess.sd <- function(x, alpha){
  n=length(x)
  times=1:n
  residus = x - median(x, na.rm = TRUE)
  lissage=loess(residus^2~times,degree=1, span = alpha, normalize = TRUE, na.action = na.exclude)
  slide.var <- rep(NA, n)
  slide.var[which(is.na(x) == TRUE)] <- approx(times, predict(lissage), xout= which(is.na(x)==TRUE))$y
  slide.var[which(is.na(x) == FALSE)] <- na.omit(lissage$fitted)
  std.t <- sqrt(sapply(slide.var, function(y) ifelse(y>0, y, abs(y))))
  return(std.t)
}

my.estimator <- function(estimator,x){
  x1 = na.omit(x)
  if(estimator == "mad"){
    y = mad(x1)
  }else if(estimator == "Qn"){
    y = robustbase::Qn(x1)
  }else if(estimator == "Sca"){
    y = robustbase::scaleTau2(x1)
  }
  return(y)
}
# check
