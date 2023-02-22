# use to check the test result by simulation
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"FGLS.R"))
nb.sim = 1000
name.case = "pots.2016-06-05.d020"
fgls.res = get(load(file = paste0(path_results,"attribution0/FGLS-full/",name.case, "fgls.RData")))
fgls.res1 = get(load(file = paste0(path_results,"attribution0/FGLS/",name.case, "fgls.RData")))

sigma.sim = fgls.res$gps.gps$var
hetero = 1
burn.in = 1000
ar = fgls.res$gps.gps$coef.arma$phi
ma = fgls.res$gps.gps$coef.arma$theta
n = length(sigma.sim)
set.seed(1)
design.m = fgls.res1$gps.gps$design.matrix
Data.mod = design.m[,c(1:11)]
noise.model = c(ifelse(ar!=0, 1, 0), 0, ifelse(ma!=0, 1, 0))
names(noise.model) =NULL
ind.na = which(is.na(Data.mod$signal)==TRUE)
TPR = rep(NA, nb.sim)
for (i in c(1:nb.sim)) {
  y.arma = simulate.general1(N = n, arma.model = c(ar=ar,ma=ma), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
  y.arma[c((n/2):n)] = y.arma[c((n/2):n)]+0.2
  y.arma[ind.na] = NA
  Data.mod$signal = y.arma
  fit.fgls = FGLS1(design.m = Data.mod, tol=0.01, day.list = design.m$date, noise.model = noise.model, length.wind0 = 60)
  TPR[i] = fit.fgls$t.table$`t value`[9]
  print(i)
}
table(abs(TPR)<1.96)
save(TPR, file = paste0(path_results, "attribution0/TPR.RData"))

