# list of screening methods
screen.iqr <- function(x) {
  y = na.omit(x)
  Q1 <- quantile(x, .25, na.rm = TRUE)
  Q3 <- quantile(x, .75, na.rm = TRUE)
  IQR <- IQR(x, na.rm = TRUE)
  up = Q3+IQR*1.723869
  down = Q1-IQR*1.723869
  removed <- which(x<down | x>up)
  if (length(removed) >0){
    x.screened = x[-removed]
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = removed, up = up, down = down))
}
screen.diff <- function(x, dif){
  n0 = length(x)
  if(dif == 1){
    y = diff(x)
  }
  Q1 <- quantile(y, .25, na.rm = TRUE)
  Q3 <- quantile(y, .75, na.rm = TRUE)
  IQR <- IQR(y, na.rm = TRUE)
  up = Q3+IQR*1.67
  down = Q1-IQR*1.67
  removed <- which(y<down | y>up)
  list.rm <- removed
  if (length(list.rm) >0){
    for (j in c(1:length(removed))){
      ind = removed[j]
      ind.n = 11
      if(ind<10){ 
        subset <-  x[c(1 : (ind+10))]
        ind.n = ind
      }else if (ind>(n0-10)){
        subset <-  x[c((ind-10): n0)]
      }else{
        subset <- x[c((ind-10): (ind+10))]
      }
      std.2 <- c()
      for (k in c(1:2)) {
        std.2[k] <- sd(subset[-(ind.n+k-1)])
      }
      ind.new = which.min(std.2)
      if(ind.new != 1){
        list.rm[j] <- list.rm[j]+1
      }
    }
    x.screened = x[-list.rm]
    list.rm = list.rm[!duplicated(list.rm)]
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = list.rm, up = up, down = down))
}
screen.mad <- function(x) {
  y = abs(x - median(x))/mad(x)
  removed <- which(y>2.7)
  if (length(removed) >0){
    x.screened = x[-removed]
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = removed))
}
screen.sca <- function(x) {
  y = abs(x - median(x))/(robustbase::scaleTau2	(x))
  removed <- which(y>2.7)
  if (length(removed) >0){
    x.screened = x[-removed]
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = removed))
}
screen.qn <- function(x, thres) {
  y = abs(x - median(x))/(robustbase::Qn(x))
  removed <- which(y>thres)
  if (length(removed) >0){
    x.screened = x[-removed]
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = removed))
}
screen.diff.o <- function(x, dif, thres1, thres2, sdt){
  n0 = length(na.omit(x))
  if(dif == 1){
    y = diff(x)
  }
  # Q1 <- quantile(y, .25, na.rm = TRUE)
  # Q3 <- quantile(y, .75, na.rm = TRUE)
  # IQR <- IQR(y, na.rm = TRUE)
  up = thres1
  down =-thres1
  removed <- which(y<down | y>up)
  list.rm <- removed
  if (length(list.rm) >0){
    for (j in c(1:length(removed))){
      ind = removed[j]
      std.2 <- c()
      for (k in c(1:2)) {
        std.2[k] <- x[ind+k-1] - median(x, na.rm = TRUE)
      }
      ind.new = which.max(abs(std.2))
      if(ind.new != 1){
        list.rm[j] <- list.rm[j]+1
      }
    }
    # x.screened = x[-list.rm]
    list.rm = list.rm[!duplicated(list.rm)]
    list.rm  <- list.rm[which(x[list.rm]<(-thres2)| x[list.rm]>thres2)]
    if(length(list.rm) !=0){
      # check probability 
      detect.val = abs(x[list.rm ])
      removed <- removed[order(detect.val, decreasing = TRUE)]
      detect.val <- abs(x[list.rm])
      limit = max(detect.val)
      last.rm <- c()
      for (i in seq(max(limit,thres2), min(limit,thres2),-0.1)) {
        detect = which(detect.val > i)
        E = n0*2*pnorm(-i, mean = 0, sd = sdt)
        npj = round(length(detect)-E)
        if(npj >0){
          detect.rm = detect[1:npj]
          last.rm <- c(last.rm, list.rm[detect.rm])
          removed <- removed[-detect.rm]
          detect.val <- detect.val[-detect.rm]
        }
      }
      removed = last.rm[!duplicated(last.rm)]
    }else {
      removed <- c()
    }
  }
  if(length(removed) != 0){  x.screened = x[-removed]}
  else{removed <- c()
  x.screened = x}
  return(list(data = x.screened, point.rm = removed, up = up, down = down))
}
screen.qn.o <- function(x, thres, sdt) {
  n = length(na.omit(x))
  # y = abs(x - median(x))/(robustbase::Qn(x))
  removed <- which(abs(x)>thres)
  last.rm <- c()
  if (length(removed) >0){
    # reorder ind of removal 
    detect.val = abs(x[removed])
    removed <- removed[order(detect.val, decreasing = TRUE)]
    detect.val <- abs(x[removed])
    limit = max(detect.val)
    for (i in seq(limit,thres,-0.1)) {
      detect = which(detect.val > i)
      E = n*2*pnorm(-i, mean = 0, sd = sdt)
      npj = round(length(detect)-E)
      if(npj >0){
        detect.rm = detect[1:npj]
        last.rm <- c(last.rm, removed[detect.rm])
        removed <- removed[-detect.rm]
        detect.val <- detect.val[-detect.rm]
      }
    }
    removed = unique(last.rm)
    if(is.null(removed) == TRUE){
      removed <- c()
      x.screened = x
    }else{
      x.screened = x[-removed]
    }
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = removed))
}

myround <- function(y){
  if(y == 0.5){z <- 1}
  else{z <- round(y)}
  return(z)
}
scr.O <- function(x, method, estimator, fix.thres){
  n = length(x)
  mean0 = median(x, na.rm = TRUE)
  if(method == "IQR"){
    Q1 <- quantile(x, .25, na.rm = TRUE)
    Q3 <- quantile(x, .75, na.rm = TRUE)
    IQR <- IQR(x, na.rm = TRUE)
    up = Q3+IQR*1.723869
    down = Q1-IQR*1.723869
  }else if(method == "sigma"){
    sdt <- my.estimator(estimator = estimator, x)
    up = 3*sdt
    down = -3*sdt
  }else if(method == "def"){
    sdt = 1
    up = fix.thres
    down = -fix.thres
    mean0 = 0
  }
  candidate <- which(x<down | x>up)
  if (length(candidate) >0){
    # reorder ind of removal 
    detect.val = abs(x[candidate])
    removed <- candidate[order(detect.val, decreasing = TRUE)]
    detect.val <- abs(x[removed])
    limit = floor(max(detect.val)*10)/10
    last.rm <- c()
    thres <- floor(min(abs(c(up,down)))*10)/10
    for (i in seq(limit,thres,-0.1)) {
      detect = which(detect.val > i)
      E = n*2*pnorm(-i, mean = mean0, sd = sdt)
      npj = myround(length(detect)-E)
      if(npj >0){
        detect.rm = detect[1:npj]
        last.rm <- c(last.rm, removed[detect.rm])
        removed <- removed[-detect.rm]
        detect.val <- detect.val[-detect.rm]
      }
    }
    candidate.out = unique(last.rm)
  }else{ candidate.out <- c()}
  if( length(candidate.out) >0 ){ 
    x.out <- x
    x.out[candidate.out] <- NA
  }else{
    x.out <- x
  }
  return(list(x = x.out, point.rm = candidate.out))
}
screen.O <- function(Y, name.var, method, iter, estimator, fix.thres,  global.mu, loes, loes.method){
  last.rm <- list()
  removed <- 1
  y <- unlist(Y[name.var])
  plo <- list()
  i = 0
  x <- y
  normalized <- list()
  mu.est <- list() 
  sd.est <- list()
  while(length(removed) > 0){
    i <- i+1
    if(iter ==1){
      normalise <- one.step.norm(Y, name.var = name.var, estimator = estimator, length.wind = 60, global.mu = global.mu,
                                 loes = loes, loes.method = loes.method)
      x <- normalise$norm.sig
      names(x) <- NULL
    }
    normalized[[i]] <- x
    mu.est[[i]] <- normalise$mu
    sd.est[[i]] <- normalise$sd
    re = scr.O(x, method = method, estimator = estimator, fix.thres = 0)
    plo[[i]] <- x
    Y[re$point.rm, name.var] <- NA
    x <- unlist(Y[name.var])
    last.rm[[i]] <- re$point.rm
    removed <- re$point.rm
  }
  
  if(length(last.rm) != 0){
    x.screened = y
    x.screened[unlist(last.rm)] <- NA
  } else{
    x.screened = y
  }
  return(list(data = x.screened, point.rm = last.rm, normalized = normalized, mu.est = mu.est, sd.est = sd.est))
}
# make some iteration with this block 
# NEED TO SEE THE PLOT AFTER EACH ITERATION, WHY THEY REMOVE TOO MANY POINTS?
# IF NEEDED, DO SIMULATIONS
screen.O1 <- function(y, method,estimator){
  last.rm <- c()
  removed <- 1
  plo <- list()
  i = 0
  x <- y
  while(length(removed) > 0){
    i <- i+1
    # if(iter ==1){
    #   x <- one.step.norm(Y, name.var = name.var)
    #   names(x) <- NULL
    # }
    re = scr.O(x, method = method,estimator)
    plo[[i]] <- x
    x[re$point.rm] <- NA
    last.rm <- c(last.rm,re$point.rm)
    removed <- re$point.rm
  }
  print(last.rm)
  if(length(last.rm) != 0){
    x.screened = y[-last.rm]
  } else{
    x.screened = y
  }
  return(list(data = x.screened, point.rm = last.rm))
}

my.estimator <- function(estimator,x){
  x1 = na.omit(x)
  if(estimator == "mad"){
    n0 = length(x1)
    f <- qnorm(3/4)*(1-0.7612/n0 - 1.123/(n0^2))
    y = mad(x1, constant = 1/f)
  }else if(estimator == "Qn"){
    y = robustbase::Qn(x1)
  }else if(estimator == "Sca"){
    y = robustbase::scaleTau2(x1, consistency = "finiteSample")
  }
  return(y)
}

scr.O1 <- function(x, method, estimator, fix.thres){ # make it 2 sides 
  n = length(x)
  mean0 = median(x, na.rm = TRUE)
  if(method == "IQR"){
    Q1 <- quantile(x, .25, na.rm = TRUE)
    Q3 <- quantile(x, .75, na.rm = TRUE)
    IQR <- IQR(x, na.rm = TRUE)
    up = Q3+IQR*1.723869
    down = Q1-IQR*1.723869
  }else if(method == "sigma"){
    sdt <- my.estimator(estimator = estimator, x)
    up = 3*sdt
    down = -3*sdt
  }else if(method == "def"){
    sdt = 1
    up = fix.thres
    down = -fix.thres
    mean0 = 0
  }
  candidate <- which(x<down | x>up)
  if (length(candidate) >0){
    # reorder ind of removal 
    detect.val = abs(x[candidate])
    removed <- candidate[order(detect.val, decreasing = TRUE)]
    detect.val <- abs(x[removed])
    limit = floor(max(detect.val)*10)/10
    last.rm <- c()
    thres <- floor(min(abs(c(up,down)))*10)/10
    for (i in seq(limit,thres,-0.1)) {
      detect = which(detect.val > i)
      E = n*2*pnorm(-i, mean = mean0, sd = sdt)
      npj = myround(length(detect)-E)
      if(npj >0){
        detect.rm = detect[1:npj]
        last.rm <- c(last.rm, removed[detect.rm])
        removed <- removed[-detect.rm]
        detect.val <- detect.val[-detect.rm]
      }
    }
    candidate.out = unique(last.rm)
  }else{ candidate.out <- c()}
  if( length(candidate.out) >0 ){ 
    x.out <- x
    x.out[candidate.out] <- NA
  }else{
    x.out <- x
  }
  return(list(x = x.out, point.rm = candidate.out))
}

# plot function

plot_screening <- function(case.name, data.in, var.name, side){
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  
  Y = data.in[[case.name]]
  if(side == "bef"){
    y = subset(Y, Y$date <= breakpoint)
    y <- tidyr::complete(y, date = seq((breakpoint-364), breakpoint, by = "day"))
  }else{
    y = subset(Y, Y$date > breakpoint)
    y <- tidyr::complete(y, date = seq((breakpoint+1), (breakpoint+365), by = "day"))
  }
  var.name =  var.name
  
  method <- screen.O(Y = y, name.var = var.name, method = 'sigma', iter = 1, estimator = "Sca", fix.thres = 0, 
                     loes = 0, loes.method ="", global.mu = 0)
  mu = as.data.frame(method$mu.est, col.names = paste("mu", c(1:length(method$mu.est))))
  sd1 = as.data.frame(method$sd.est, col.names = paste("sd", c(1:length(method$mu.est))))
  data.i = cbind(date = y$date, mu)
  for(x in c(1:length(method$normalized))) {
    data.i[paste("lower",x, sep = "")] = mu[,x] - 3* sd1[,x]
    data.i[paste("upper",x, sep = "")] = mu[,x] + 3* sd1[,x]
    series.rm <- rep(NA, length(mu[,x]))
    if(x < length(method$normalized)){
      series.rm[unlist(method$point.rm[[x]])] <- unlist(y[method$point.rm[[x]], var.name])
    }
    data.i[paste("rm",x, sep = "")] = series.rm 
  }
  data.i$y = unlist(y[var.name])
  
  Y.i <- data.i
  a = reshape2::melt(Y.i, id = "date")
  gr = c(c(1:length(method$normalized)), rep(c(1:length(method$normalized)), each = 3), 0)
  gr[c(1:length(method$normalized))*3+length(method$normalized)] <- 10
  a$col <- as.factor(rep(gr, each = 365 ))
  b = a[which(a$col == 10),]
  d = a[which(a$col != 10),]
  b$col <- as.factor(rep(c(1:length(method$normalized)), each = 365 ))
  p <- ggplot(data = d, aes(x = date, y = value, group = variable))+
    geom_line(aes(colour= col))+ theme_bw()
  p <- p + geom_point( data = b, aes(y = value, colour= col), size=2)
  p <- p + theme(axis.text = element_text(size = 15),
                 axis.title = element_text(size=15,face="bold"),
                 plot.subtitle = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
  jpeg(paste0(path_results,"attribution/" , case.name,side, ".screened.", ".jpeg"),
       width = 3000, height = 1500,res = 300)
  
  print(p)
  dev.off()
}
