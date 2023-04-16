# Final prog to study the auto.arima to put into the thesis
source(paste0(path_code_att,"FGLS.R"))
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"support_characterization.R"))

# SIMULATION ---------------------------------------------------------------
## AR(1) model -------------------------------------------------------------
length.list = seq(200, 2000, 200)
coef.list = seq(0, 0.8, 0.1)

nb.sim = 1000
tot.res = data.frame(n = length.list, 
                     TPR.white = rep(NA, length(length.list)),
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1

set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = fit.arima(y.ar)
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_AR1.RData"))

## MA(1) model -------------------------------------------------------------
length.list = seq(200, 2000, 200)
coef.list = seq(0, 0.8, 0.1)

nb.sim = 1000
tot.res = data.frame(n = length.list, 
                     TPR.white = rep(NA, length(length.list)),
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1

set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=0,ma=ar0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = fit.arima(y.ar)
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_MA1.RData"))
## ARMA(1,1) model -------------------------------------------------------------

length.list = seq(200, 2000, 200)
coef.list = seq(0, 0.8, 0.1)

nb.sim = 1000
tot.res = data.frame(n = length.list, 
                     TPR.white = rep(NA, length(length.list)),
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1

set.seed(1)
tot.res = list()
for (l in c(1:length(length.list))) {
  n = length.list[l]
  tot.fit = list()
  for (i in c(1:length(coef.list))) {
    fit.i = list()
    ar0 = coef.list[i]
    ma0 = 0.3-ar0
    for (j in c(1:nb.sim)) {
      y.ar = simulate.general1(N = n, arma.model = c(ar=ar0,ma=ma0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      # fit
      fit.ar = fit.arima(y.ar)
      fit.i[[j]] = fit.ar
    }
    tot.fit[[i]] = fit.i
  }
  tot.res[[l]] = tot.fit
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_ARMA1.RData"))




# PLOT
arma = get(load(file = paste0(path_results, "attribution0/performance_autoarima_ARMA1.RData")))
ar = get(load(file = paste0(path_results, "attribution0/performance_autoarima_AR1.RData")))
ma = get(load(file = paste0(path_results, "attribution0/performance_autoarima_MA1.RData")))

get_data <- function(list.ini, param.val, true.model, details){
  nb.sim = length(list.ini[[1]][[1]])
  # choose details case at specific length 
  if(details!=0){
    tot.df = data.frame(matrix(NA, ncol = length(list.ini[[1]]), nrow = nb.sim))
    for (j in c(1:length(list.ini[[1]]))) {
      model.est = sapply(c(1:nb.sim), function(x) model.iden(list.ini[[5]][[j]][[x]]$pq))
      tot.df[,j] = model.est
    }
  }else{
    
    if(param.val == 0){
      tot.df = data.frame(matrix(NA, ncol = length(list.ini), nrow = length(list.ini[[1]])))
      for (i in c(1:length(list.ini))) {
        for (j in c(1:length(list.ini[[1]]))) {
          model.est = sapply(c(1:nb.sim), function(x) model.iden(list.ini[[i]][[j]][[x]]$pq))
          y = length(which(model.est == true.model))
          tot.df[j,i] = y
        }
      }
    }else{ # could be filtered the case of true model 
      tot.df = data.frame(matrix(NA, ncol = 2*length(list.ini), nrow = length(list.ini[[1]])))
      for (i in c(1:length(list.ini))) {
        for (j in c(1:length(list.ini[[1]]))) {
          phi = sapply(c(1:nb.sim), function(x) list.ini[[i]][[j]][[x]]$coef[1])
          theta = sapply(c(1:nb.sim), function(x) list.ini[[i]][[j]][[x]]$coef[3])
          tot.df[j,(2*i-1)] = mean(phi, na.rm = TRUE)
          tot.df[j,(2*i)] = mean(theta, na.rm = TRUE)
        }
      }
    }
  }
  return(tot.df)
}

a = get_data(list.ini = ar, param.val = 0, true.model = "AR(1)", details = 3)



