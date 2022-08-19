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

scr.O <- function(x, method ){
  n = length(x)
  if(method == "IQR"){
    Q1 <- quantile(x, .25, na.rm = TRUE)
    Q3 <- quantile(x, .75, na.rm = TRUE)
    IQR <- IQR(x, na.rm = TRUE)
    up = Q3+IQR*1.723869
    down = Q1-IQR*1.723869
  }else if(method == "sigma"){
    sdt <- robustbase::scaleTau2(x, na.rm = TRUE)
    up = 3*sdt
    down = -3*sdt
  }
  candidate <- which(x<down | x>up)
  if (length(candidate) >0){
    # reorder ind of removal 
    detect.val = abs(x[candidate])
    removed <- candidate[order(detect.val, decreasing = TRUE)]
    detect.val <- abs(x[removed])
    limit = max(detect.val)
    last.rm <- c()
    thres <- min(abs(c(up,down)))
    for (i in seq(limit,thres,-0.1)) {
      detect = which(detect.val > i)
      E = n*2*pnorm(-i, mean = median(x, na.rm = TRUE), sd = robustbase::scaleTau2(x, na.rm = TRUE))
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
    }
    candidate.out <- removed
  }else{ candidate.out <- c()}
  if( length(candidate.out) >0 ){ 
    x.out <- x
    x.out[candidate.out] <- NA
  }else{
    x.out <- x
  }
  return(list(x.out = x.out, candidate.out = candidate.out))
}
screen.O <- function(Y, name.var, method){
  last.rm <- c()
  removed <- 1
  y <- Y[name.var] 
  while(length(removed) > 0){
    x <- one.step.norm(Y, name.var = name.var)
    names(x) <- NULL
    r = scr.O(x, method = method)
    x <- r$x.out
    Y[name.var] <- x
    last.rm <- c(last.rm,r$candidate.out)
    removed <- r$candidate.out
    print(removed)
  }
  if(length(last.rm) != 0){
    x.screened = y[-last.rm]
  } else{
    x.screened = y
  }
  return(list(data = x.screened, point.rm = last.rm))
}
# make some iteration with this block 
# NEED TO SEE THE PLOT AFTER EACH ITERATION, WHY THEY REMOVE TOO MANY POINTS?
# IF NEEDED, DO SIMULATIONS



