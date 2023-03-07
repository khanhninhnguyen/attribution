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
tot.fit = list()
for (i in c(1:length(length.list))) {
  n = length.list[i]
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.ar = simulate.general1(N = n, arma.model = c(ar=0.3,ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.ma = simulate.general1(N = n, arma.model = c(ar=0,ma=0.3), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.arma = simulate.general1(N = n, arma.model = c(ar=0.7,ma=-0.4), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit
    fit.ar = fit.arima(y.ar)
    fit.ma = fit.arima(y.ma)
    fit.arma = fit.arima(y.arma)
    TPR[j,] = c( model.iden(fit.ar$pq), model.iden(fit.ma$pq), model.iden(fit.arma$pq))
    fit.i[[j]] = list( ar = fit.ar, ma = fit.ma, arma = fit.arma)
  }
  colnames(TPR) = c("ar", "ma", "arma")
  tot.res[i,c(2:4)] = c(length(which(TPR$ar == "AR(1)")), 
                  length(which(TPR$ma == "MA(1)")),
                  length(which(TPR$arma == "ARMA(1,1)")))
  tot.fit[[i]] = fit.i
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_length.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_length_all.RData"))

# a = get(load( file = paste0(path_results, "attribution0/performance_autoarima_length.RData")))
# colnames(a) = c("N", "AR(1)", "MAR(1)", "ARMA(1,1)")
# dat.p = reshape2::melt(a, id="N")
# p = ggplot(data = dat.p, aes(x = N, y = value/nb.sim, col = variable))+
#   theme_bw()+
#   geom_point(size = 0.3)+
#   geom_hline(yintercept = 0.95, lwd = 0.3)+
#   ylab("TPR")+
#   theme(axis.text.x = element_text(size = 5), 
#       axis.text.y = element_text(size = 5),
#       legend.text=element_text(size=4),
#       axis.title = element_text(size = 5), 
#       legend.key.size = unit(0.3, "cm"), 
#       plot.tag = element_text(size = 5), 
#       plot.subtitle = element_text(size = 5),
#       legend.title=element_blank())
# 
# ggsave(paste0(path_results,"attribution0/TPR.auto.arima.identification_length.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)


# coefficient dependence AR(1) --------------------------------------------------
coef.list = seq(0, 0.8, 0.1)
set.seed(1)
TPR.ar <- rep(NA, nb.sim)
tot.res <- list()
tot.fit = list()

for (i in c(1:length(coef.list))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.ar = simulate.general1(N = n, arma.model = c(ar=coef.list[i],ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.ar = fit.arima(y.ar)
    TPR.ar[j] = model.iden(fit.ar$pq)
    fit.i[[j]] = fit.ar
  }
  tot.res[[i]] = TPR.ar
  tot.fit[[i]] = fit.i
}


save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_ar.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_coef_ar_all.RData"))

# b = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_ar.RData")))
# coef.list = seq(0.1, 0.8, 0.1)
# 
# tpr = sapply(c(1:length(b)), function(x) length(which(b[[x]]=="AR(1)")))[-1]/nb.sim
# dat.p = data.frame(ar = coef.list, tpr = tpr)
# p = ggplot(data = dat.p, aes(x = ar, y = tpr))+
#   theme_bw()+
#   geom_point(size = 0.3)+
#   geom_hline(yintercept = 0.95, lwd = 0.3)+
#   ylab("TPR")+
#   theme(axis.text.x = element_text(size = 5), 
#         axis.text.y = element_text(size = 5),
#         legend.text=element_text(size=4),
#         axis.title = element_text(size = 5), 
#         legend.key.size = unit(0.3, "cm"), 
#         plot.tag = element_text(size = 5), 
#         plot.subtitle = element_text(size = 5),
#         legend.title=element_blank())
# ggsave(paste0(path_results,"attribution0/TPR.auto.arima.identification_ar.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

# coefficient dependence MA(1)  --------------------------------------------------
coef.list = seq(0, 0.8, 0.1)
set.seed(1)
TPR.ma <- rep(NA, nb.sim)
tot.res <- list()
tot.fit = list()
for (i in c(1:length(coef.list))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.ma = simulate.general1(N = n, arma.model = c(ar=0,ma=coef.list[i]), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.ma = fit.arima(y.ma)
    TPR.ma[j] = model.iden(fit.ma$pq)
    fit.i[[j]] = fit.ma
  }
  tot.res[[i]] = TPR.ma
  tot.fit[[i]] = fit.i
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_ma.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_coef_ma_all.RData"))

# coefficient dependence ARMA(1,1) --------------------------------------------------
coef.list.ar = seq(0, 0.8, 0.1)
coef.list.ma = 0.3-coef.list.ar
tot.fit = list()
set.seed(1)
TPR.ar <- rep(NA, nb.sim)
tot.res <- list()
for (i in c(1:length(coef.list.ar))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.arma = simulate.general1(N = n, arma.model = c(ar=coef.list.ar[i],ma=coef.list.ma[i]), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.arma = fit.arima(y.arma)
    TPR.ar[j] = model.iden(fit.arma$pq)
    fit.i[[j]] = fit.arma
  }
  tot.res[[i]] = TPR.ar
  tot.fit[[i]] = fit.i
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_arma.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_coef_arma_all.RData"))


# plot --------------------------------------------------------------------
## length
all.length = get(load(file = paste0(path_results, "attribution0/performance_autoarima_length_all.RData")))
#NEED TO BE CHECKED ONLY WHEN THEY IDENTIRY THE TRUE MODEL
phi = sapply(c(1:length(all.length)), function(y) sapply(c(1:length(all.length[[y]])), function(x)all.length[[1]][[x]]$arma$coef[1])) 
phi.rms = sapply(c(1:length(all.length)), function(x) sqrt( sum((phi[,x] - 0.3)^2)))
tpr.ar1 = sapply(c(1:length(ar1)), function(x) which(ar1[[x]]=="AR(1)"))

boxplot(phi)
ar1 = get(load( file = paste0(path_results, "attribution0/performance_autoarima_length.RData")))

ar1 = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_ar_all.RData")))
ma1 = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_ma_all.RData")))
arma11 = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_arma_all.RData")))

tpr.ar1 = sapply(c(1:length(ar1)), function(x) length(which(ar1[[x]]=="AR(1)")))[-1]/nb.sim
tpr.ma1 = sapply(c(1:length(ma1)), function(x) length(which(ma1[[x]]=="MA(1)")))[-1]/nb.sim
tpr.arma11 = sapply(c(1:length(arma11)), function(x) length(which(arma11[[x]]=="ARMA(1,1)")))[-1]/nb.sim

ar1 = sapply(c(1:length(all.length[[2]])), function(x) model.iden(all.length[[2]][[x]]$arma$pq))
ind.ar = which(ar1=="AR(1)")

