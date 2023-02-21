# study the performance of auto.arima 
source(paste0(path_code_att,"FGLS.R"))
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"support_characterization.R"))

# As a function of length -------------------------------------------------
length.list = seq(200, 2000, 200)
nb.sim = 1000
tot.res = data.frame(n = length.list, 
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1
## As a function of length -------------------------------------------------
set.seed(1)
for (i in c(1:length(length.list))) {a
  n = length.list[i]
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  for (j in c(1:nb.sim)) {
    y.ar = simulate.general1(N = n, arma.model = c(ar=0.3,ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.ma = simulate.general1(N = n, arma.model = c(ar=0,ma=0.3), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.arma = simulate.general1(N = n, arma.model = c(ar=0.7,ma=-0.4), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.ar = fit.arima(y.ar)
    fit.ma = fit.arima(y.ma)
    fit.arma = fit.arima(y.arma)
    TPR[j,] = c( model.iden(fit.ar$pq), model.iden(fit.ma$pq), model.iden(fit.arma$pq))
  }
  colnames(TPR) = c("ar", "ma", "arma")
  tot.res[i,c(2:4)] = c(length(which(TPR$ar == "AR(1)")), 
                  length(which(TPR$ma == "MA(1)")),
                  length(which(TPR$arma == "ARMA(1,1)")))
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_length.RData"))
a = get(load( file = paste0(path_results, "attribution0/performance_autoarima_length.rds")))






