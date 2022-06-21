# Model identification
# function ----------------------------------------------------------------
# return pvalue and order 
p.and.coef <- function(fitARIMA, pq1, nb.or){
  test.sig = coeftest(fitARIMA)
  ord = pq1[c(1,3)]
  orde = c(rbind(ord,ord-1))
  orde[which(orde <0)] <- 0
  ind.param = which(orde >0)
  orde[ind.param] <- fitARIMA$coef[1:nb.or]
  p.value <- rep(-1, 4)
  p.value[ ind.param] <- test.sig[,4][1:nb.or]
  return(list(p.value = p.value , coef = orde))
}
# return significant order
check_sig <- function(p.val, alpha){
  ar.or = length(which(p.val[1:2] >0 & p.val[1:2] <= alpha))
  ma.or = length(which(p.val[3:4] >0 & p.val[3:4] <= alpha))
  return(c(ar.or, 0, ma.or))
}
# return last significant order when re-estimate the parameter with significant order
last_signif <- function(signal, pq, alpha){  
  nb.or <- sum(pq)
  pq1 = rep(NA,3)
  while ( identical(as.numeric(pq1), pq) == FALSE) { # iteratively identify the model, stop when the model are the same after the significant check
    pq1 = pq
    if(nb.or==0){
      pandcoef <- list(p.value = rep(-1,4),coef = rep(0,4))
    }else{
      fitARIMA = try(arima( signal, pq, method="ML"), TRUE)
      if (class(fitARIMA) == "try-error"){
        fitARIMA = fit.b
      }
      pandcoef <- p.and.coef(fitARIMA, pq, nb.or)
    }
    pq = check_sig(p.val = pandcoef$p.value, alpha = alpha)
    nb.or <- sum(pq)
  }
  return(list( pq = pq, pandcoef = pandcoef))
}
diff.var <- function(name.test){
  if(name.test == "gps.gps"){
    varname = c("GPS.x", "GPS.y")
  }
  if(name.test == "gps.era"){
    varname = c("GPS.x", "ERAI.x")
  }
  if(name.test == "gps1.era"){
    varname = c("GPS.y", "ERAI.x")
  }
  if(name.test == "gps.era1"){
    varname = c("GPS.x", "ERAI.y")
  }
  if(name.test == "gps1.era1"){
    varname = c("GPS.y", "ERAI.y")
  }
  if(name.test == "era.era"){
    varname = c("ERAI.x", "ERAI.y")
  }
  return(varname)
}

data_test <- get(load(paste0(path_results,"attribution/data_test_",version_name,nb_test = nb_test.ref,
                             nearby.ver = nearby_ver,limit.type = limit.type,screen_value,"multiple.RData")))

list.print <- get(load(file = paste0(path_results, "attribution/", "list.nearby.station.homo", nb_test = nb_test.near,"-",
                      screen_value, tolerance_noise, ".RData")))
# model identification ----------------------------------------------------
segment.or = "before" # need to change in the sigmal when compute the difference too in the loop
order.init <- data.frame(matrix(NA, ncol = 18, nrow = nrow(list.print)))
coef.init <- data.frame(matrix(NA, ncol = 24, nrow = nrow(list.print)))
pval.init <- data.frame(matrix(NA, ncol = 24, nrow = nrow(list.print)))
order.sig <- data.frame(matrix(NA, ncol = 18, nrow = nrow(list.print)))
coef.sig <- data.frame(matrix(NA, ncol = 24, nrow = nrow(list.print)))
pval.sig <- data.frame(matrix(NA, ncol = 24, nrow = nrow(list.print)))

for (testi in 1:6) {
  
  name.test = list.test[testi]
  var1 = diff.var(c(name.test))[1]
  var2 = diff.var(c(name.test))[2]
  # list.var <- get(load(file = paste0(path_results, "attribution/monthly-var-weigted.",name.test,".", 
  #                                    duration, "days.level", sig.level, ".scr", screen_value, "no.test.", nb_test = nb_test.ref,
  #                                    ".nearby.", nearby.ver = nearby_ver, check_correlation, ".RData")))
  for (k in 1:nrow(list.print)) {
    # read data for each case
    station.ref = list.print$ref.station[k]
    station.near = list.print$nearby.station[k]
    detection = list.print$detected[k]
    data = as.data.frame(data_test[[paste0( station.ref,".",as.character(detection), ".",  station.near)]][segment.or])
    colnames(data) <- gsub(paste0(segment.or,'.'), '', colnames(data))
    
    data[,name.test] = data[,c(var1)] -  data[,c(var2)] # compute the difference
    data$month = data$month.x
    data$year = data$year.x
    
    month.std = RobEstiMonthlyVariance.diff.S(Y = (na.omit(data)), name.var = name.test, alpha = 0)
    
    std.t = rep(NA, nrow(data))
    for (m in 1:12) {
      if(month.std[m] != 0 & is.na(month.std[m]) == FALSE){
        std.t[which(as.numeric(data$month.x)==m)] <-  month.std[m]
      }
    }  
    if(sum(is.na(std.t)) != 0){
      data <- data[-which(is.na(std.t)),]
      signal.test = data[,name.test]/(std.t[!is.na(std.t)])
    } else{
      signal.test = data[,name.test]/(std.t)
    }
    
    # 1. fit initial auto.arima order
    
    fit.b = forecast::auto.arima(signal.test , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                 max.p = 2, max.q = 2, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
    
    pq <- arimaorder(fit.b)
    order.init[k, c((testi*3-2): (testi*3))] <- pq
    options(warn = 1)
    
    # 2. 3. fit initial arima param and its significance level
    
    refit0 = last_signif(signal = signal.test, pq, alpha = significant.level)
    
    coef.init[k, c((testi*4-3): (testi*4))] <-  refit0$pandcoef$coef
    pval.init[k, c((testi*4-3): (testi*4))] <-  refit0$pandcoef$p.value
    pq = refit0$pq
    
    # 4. Choose significant param (order) and refit with order >1
    
    if( any(pq > 1)){
      fit.b = forecast::auto.arima( signal.test, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                    max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
      pq = arimaorder(fit.b)
    }
    
    refit1 = last_signif(signal = signal.test, pq, alpha = significant.level)
    
    pq = refit1$pq
    order.sig[k, c((testi*3-2): (testi*3))] <- pq
    coef.sig[k, c((testi*4-3): (testi*4))] <-  refit1$pandcoef$coef
    pval.sig[k, c((testi*4-3): (testi*4))] <-  refit1$pandcoef$p.value
    
  }
  print(c((testi*3-2):(testi*3)))
}

col.names <- c(rbind(paste0(list.test, ".p1"),paste0(list.test, ".p2"),paste0(list.test, ".q1"),paste0(list.test, ".q2")))
col.names1 <- c(rbind(paste0(list.test, ".p"),paste0(list.test, ".d"),paste0(list.test, ".q")))

colnames(coef.init) <- col.names
colnames(pval.init) <- col.names
colnames(coef.sig) <- col.names
colnames(pval.sig) <- col.names
colnames(order.init) <- col.names1
colnames(order.sig) <- col.names1

file.name.model = paste0(segment.or, "_arima_alpha",alpha = significant.level, nb_test = nb_test.ref, "-", screen_value, ".txt")


write.table(order.init,  file = paste0(path_results,"attribution/", "order_init_", file.name.model),
            sep = "\t", quote = FALSE, col.names = TRUE)
write.table(order.sig,  file = paste0(path_results,"attribution/", "order_sig_pmax1_", file.name.model),
            sep = "\t", quote = FALSE, col.names = TRUE)
write.table(coef.init,  file = paste0(path_results,"attribution/", "coef_init_", file.name.model),
            sep = "\t", quote = FALSE, col.names = TRUE)
write.table(coef.sig,  file = paste0(path_results,"attribution/", "coef_sig_pmax1_", file.name.model),
            sep = "\t", quote = FALSE, col.names = TRUE)
write.table(pval.init,  file = paste0(path_results,"attribution/", "p.value_init_", file.name.model),
            sep = "\t", quote = FALSE, col.names = TRUE)
write.table(pval.sig,  file = paste0(path_results,"attribution/", "p.value_sig_pmax1_", file.name.model),
            sep = "\t", quote = FALSE, col.names = TRUE)

# decide the last model to correct t-test ---------------------------------

order.sig.bef = read.table(file = paste0(path_results,"attribution/", "order_sig_pmax1_before_arima_alpha",alpha = significant.level, nb_test = nb_test.ref,"-",screen_value,".txt"))
order.sig.aft = read.table(file = paste0(path_results,"attribution/", "order_sig_pmax1_after_arima_alpha",alpha = significant.level, nb_test = nb_test.ref,"-",screen_value,".txt"))

model.iden <- function(order){
  model = c()
  if (identical(order, c(1,0,1))){ model = "arma11"}
  else if (identical(order, c(1,0,0))){ model = "ar1"}
  else if (identical(order, c(0,0,1))){ model = "ma1"}
  else if (identical(order, c(0,0,0))){ model = "white"}
  return(model)
}

list.order.bef = data.frame(matrix(NA, nrow = nrow(list.print), ncol = 6))
list.order.aft = data.frame(matrix(NA, nrow = nrow(list.print), ncol = 6))

colnames(list.order.bef) <- list.test
colnames(list.order.aft ) <- list.test

for (x in 1:6) {
  for (k in 1:nrow(list.print)) {
    list.order.bef[k,x] <- model.iden(as.numeric(order.sig.bef[k,c((x*3-2):(x*3))]))
    list.order.aft[k,x] <- model.iden(as.numeric(order.sig.aft[k,c((x*3-2):(x*3))]))
    
  }
}

# choose only one model for 2 sides ----------------------------------------

sort.fc <- function(x,y){
  res <- c()
  if(x == y){
    res = x
  }else{
    if (x == "white" & y %in% c("ar1", "ma1") ){
      res = y
    }else if(x %in% c("ar1", "ma1") & y == "white"){
      res = x
    }else if(x == "arma11" & y == "white"){
      res = y
    }else if(y == "arma11" & x == "white"){
      res = x
    }else if(x == "ar1" & y != "white" ){
      res = x
    }else if(x != "white" & y == "ar1" ){
      res = y
    }else if(x == "ma1" & y == "arma11" ){
      res = x
    }else if(x == "arma11" & y == "ma1" ){
      res = y
    }
  } 
  return(res)
}
last.all.model <- data.frame(matrix(NA, ncol = 6, nrow = nrow(list.order.bef)))

for (i in 1:6) {
  last.model <- unlist(sapply(c(1:nrow(list.order.bef)), function(x) sort.fc(list.order.bef[x,i], list.order.aft[x,i])))
  last.all.model[,i] <- unlist(last.model)
}

write.table(last.all.model,  file = paste0(path_results,"attribution/", "last.model", "_arima_alpha", alpha = significant.level, nb_test = nb_test.ref,"-",screen_value,tolerance_noise,".txt"),
            sep = "\t", quote = FALSE, col.names = TRUE)

# from the last model, estimate the parameter and effective ratio ----------------------------
model.order <- function(name.model){
  if(name.model == "ar1"){ order = c(1,0,0)}
  if(name.model == "ma1"){ order = c(0,0,1)}
  if(name.model == "arma11"){ order = c(1,0,1)}
  if(name.model == "white"){ order = c(0,0,0)}
  return(order)
}
model.all = read.table(file = paste0(path_results,"attribution/", "last.model", "_arima_alpha",alpha = significant.level, nb_test = nb_test.ref,"-",screen_value,tolerance_noise,".txt"))
ratio.cor <- data.frame(matrix(NA, ncol = 6, nrow = nrow(list.order.bef)))
param.all <- data.frame(matrix(NA, ncol = 6, nrow = nrow(list.order.bef)))
for (testi in 1:6) {
  
  name.test = list.test[testi]
  var1 = diff.var(c(name.test))[1]
  var2 = diff.var(c(name.test))[2]
  
  for (k in 1:nrow(list.print)) {
    # read data for each case
    station.ref = list.print$ref.station[k]
    station.near = list.print$nearby.station[k]
    detection = list.print$detected[k]
    data = data_test[[paste0( station.ref,".",as.character(detection), ".",  station.near)]]
    bef = data$before
    aft = data$after
    L = nrow(bef)
    bef[,name.test] = bef[,c(var1)] -  bef[,c(var2)] # compute the difference
    aft[,name.test] = aft[,c(var1)] -  aft[,c(var2)] # compute the difference
    
    model.spe = model.all[k,testi]
    if (model.spe == "white"){
      r = 1 # n/ne ratio
      param = 0
    }else if (model.spe == "ar1"){
      bef.p = arima( bef[,name.test], order = c(model.order("ar1")), method = "ML")
      aft.p = arima( aft[,name.test], order = c(model.order("ar1")), method = "ML")
      ar = (max(0,bef.p$coef[1]) + max(0,aft.p$coef[1]))/2
      param = ar
      r = (L + 2*(ar**(L+1)) -L*(ar**2) - 2*ar)/(L*((1-ar)**2))
    }else if (model.spe == "ma1"){
      bef.p = arima( bef[,name.test], order = c(model.order("ma1")), method = "ML")
      aft.p = arima( aft[,name.test], order = c(model.order("ma1")), method = "ML")
      ma =(max(0,bef.p$coef[1]) + max(0,aft.p$coef[1]))/2
      param = ma/(1 + ma**2)
      r = (1 + 2 *(1-1/L)*(ma/(ma**2+1)))
    }else if (model.spe == "arma11"){
      bef.p = arima( bef[,name.test], order = c(model.order("arma11")), method = "ML")
      aft.p = arima( aft[,name.test], order = c(model.order("arma11")), method = "ML")
      ar = (max(0,bef.p$coef[1]) + max(0,aft.p$coef[1]))/2
      ma = (min(0,bef.p$coef[2]) + min(0,aft.p$coef[2]))/2
      rho1 = (1+ar*ma)*(ar+ma)/(1 + 2*ar*ma + ma**2)
      param = rho1
      r = 1 + (2*rho1*(ar**L - L*ar+L - 1)/(L*((1-ar)**2)))
    }
    ratio.cor[k, testi] <- r
    param.all[k, testi] <- param
  }
}
colnames(ratio.cor) <- list.test
colnames(param.all) <- list.test

write.table(ratio.cor,  file = paste0(path_results,"attribution/", "ratio_correct_",  "_arima_alpha",alpha = significant.level,nb_test = nb_test.ref,"-",screen_value,tolerance_noise,".txt"),
            sep = "\t", quote = FALSE, col.names = TRUE)
write.table(param.all,  file = paste0(path_results,"attribution/", "param.all_var_est_",  "_arima_alpha",alpha = significant.level,nb_test = nb_test.ref,"-",screen_value,tolerance_noise,".txt"),
            sep = "\t", quote = FALSE, col.names = TRUE)

