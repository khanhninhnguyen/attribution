# list of screening methods
screen.iqr <- function(x) {
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
screen.qn.o <- function(x, thres, sdt) {
  n = length(x)
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

screen.diff.o <- function(x, dif, thres, sdt){
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
        std.2[k] <- sd(subset[-(ind.n+k-1)], na.rm = TRUE)
      }
      ind.new = which.min(std.2)
      if(ind.new != 1){
        list.rm[j] <- list.rm[j]+1
      }
    }
    # x.screened = x[-list.rm]
    list.rm = list.rm[!duplicated(list.rm)]
    # check probability 
    detect.val = abs(x[list.rm ])
    removed <- removed[order(detect.val, decreasing = TRUE)]
    detect.val <- abs(x[list.rm])
    limit = max(detect.val)
    last.rm <- c()
    for (i in seq(max(limit,thres), min(limit,thres),-0.1)) {
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
    removed =  last.rm[!duplicated(last.rm)]
    if(is.null(removed) == TRUE){
      removed <- c()
      x.screened = x
    }else{
      x.screened = x[-removed]
    }
    
  }else{ x.screened = x}
  return(list(data = x.screened, point.rm = removed, up = up, down = down))
}
