# study the performance of auto.arima 
source(paste0(path_code_att,"FGLS.R"))
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"support_characterization.R"))

# As a function of length -------------------------------------------------
length.list = seq(200, 2000, 200)
nb.sim = 10000
tot.res = data.frame(n = length.list, 
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1
## As a function of length -------------------------------------------------
set.seed(1)
for (i in c(1:length(length.list))) {
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
a = get(load( file = paste0(path_results, "attribution0/performance_autoarima_length.RData")))
colnames(a) = c("N", "AR(1)", "MAR(1)", "ARMA(1,1)")
dat.p = reshape2::melt(a, id="N")
p = ggplot(data = dat.p, aes(x = N, y = value/nb.sim, col = variable))+
  theme_bw()+
  geom_point(size = 0.3)+
  geom_hline(yintercept = 0.95, lwd = 0.3)+
  ylab("TPR")+
  theme(axis.text.x = element_text(size = 5), 
      axis.text.y = element_text(size = 5),
      legend.text=element_text(size=4),
      axis.title = element_text(size = 5), 
      legend.key.size = unit(0.3, "cm"), 
      plot.tag = element_text(size = 5), 
      plot.subtitle = element_text(size = 5),
      legend.title=element_blank())

ggsave(paste0(path_results,"attribution0/TPR.auto.arima.identification.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)


# coefficient dependence AR(1) --------------------------------------------------
coef.list = seq(0, 0.8, 0.1)
set.seed(1)
TPR.ar <- rep(NA, nb.sim)
tot.res <- list()
for (i in c(1:length(coef.list))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  for (j in c(1:nb.sim)) {
    y.ar = simulate.general1(N = n, arma.model = c(ar=coef.list[i],ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.ar = fit.arima(y.ar)
    TPR.ar[j] = model.iden(fit.ar$pq)
  }
  tot.res[[i]] = TPR.ar
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_ar.RData"))


# coefficient dependence MA(1)  --------------------------------------------------
coef.list = seq(0, 0.8, 0.1)
set.seed(1)
TPR.ma <- rep(NA, nb.sim)
tot.res <- list()
for (i in c(1:length(coef.list))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  for (j in c(1:nb.sim)) {
    y.ma = simulate.general1(N = n, arma.model = c(ar=0,ma=coef.list[i]), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.ma = fit.arima(y.ma)
    TPR.ma[j] = model.iden(fit.ma$pq)
  }
  tot.res[[i]] = TPR.ma
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_ma.RData"))

# coefficient dependence ARMA(1,1) --------------------------------------------------
coef.list.ar = seq(0, 0.8, 0.1)
coef.list.ma = 0.3-coef.list.ar

set.seed(1)
TPR.ar <- rep(NA, nb.sim)
tot.res <- list()
for (i in c(1:length(coef.list.ar))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  for (j in c(1:nb.sim)) {
    y.ar = simulate.general1(N = n, arma.model = c(ar=coef.list.ar[i],ma=coef.list.ma[i]), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.ar = fit.arima(y.ar)
    TPR.ar[j] = model.iden(fit.ar$pq)
  }
  tot.res[[i]] = TPR.ar
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_arma.RData"))

