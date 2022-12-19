# new function to identify the model of noise from the OLS residual and normalization
source(paste0(path_code_att, "newUsed_functions.R"))
source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))
library(nlme)
# Heteroskedastic ---------------------------------------------------------

# from OLS resdiual 
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.era"
one.year=365
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

# if we want to not dupplicate gps-era case
# all.cases.name = names(dat)
# all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
# unique.ind = match(unique(all.cases.ind), all.cases.ind )
# gps.era.dat = dat[unique.ind]
# data.test = gps.era.dat 
# n = length(data.test)
# sd.all <- rep(NA, n)


# OLS fit compute the normalized residual -----------------------------------------------------------------

one.year=365

compute.norm.resi <- function(dat.i, name.test){
  one.year=365
  scr = screen.O(Y = dat.i , name.var = name.test, method = 'sigma', global.mu = 0, iter = 1, estimator = "Sca", fix.thres = 0, loes = 0, loes.method = 0)
  sd0 = unlist(scr$sd.est[[1]])
  Xt = c(1:nrow(dat.i))- nrow(dat.i)/2
  
  Data.mod = data.frame( signal = dat.i[,name.test],Xt = Xt, complete.time = c(1:nrow(dat.i)))
  for (k in c(1:4)){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",k,"=cos(k*complete.time*(2*pi)/one.year),sin",k,"=sin(k*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  ols.fit = lm(signal~., data = Data.mod)
  ols.resi = rep(NA, nrow(dat.i))
  ols.resi[which(is.na(Data.mod$signal)==FALSE)] <- ols.fit$residuals 
  norm.res = ols.resi/sd0
  return(norm.res)
}

residus = list()
for (testi in c(1:6)) {
  name.test = list.test[testi]
  res.testi = list()
  for (i in c(1:length(dat))) {
    dat.all = dat[[i]]
    dat.bef = dat.all[c(1: (one.year*10)),]
    dat.aft = dat.all[-c(1: (one.year*10)),]
    res.bef = compute.norm.resi(dat.i = dat.bef, name.test)
    res.aft = compute.norm.resi(dat.i = dat.aft, name.test)
    res.testi[[names(dat)[i]]] <- c(res.bef, res.aft)
  }
  residus[[name.test]] <- res.testi
}

save(residus, file = paste0(path_results,"attribution/norm.residual.ols.RData"))
# fit the ARIMA 
residus = get(load(file = paste0(path_results,"attribution/norm.residual.ols.RData")))
fit.arima <- function(signal.test){
  fit.b = forecast::auto.arima(signal.test , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                               max.p = 2, max.q = 2, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
  
  pq <- arimaorder(fit.b)
  # order.init[k, c((testi*3-2): (testi*3))] <- pq
  options(warn = 1)
  
  refit0 = last_signif(signal = signal.test, pq, alpha = significant.level)
  pq = refit0$pq
  
  if( any(pq > 1)){
    fit.b = forecast::auto.arima( signal.test, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                  max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
    pq = arimaorder(fit.b)
  }

  refit1 = last_signif(signal = signal.test, pq, alpha = significant.level)
  
  pq = refit1$pq
  return(list(pq = pq, coef = refit1$pandcoef$coef, p = refit1$pandcoef$p.value))
}

order.arma <- list()
coef.arma <- list()
for (testi in c(1:6)) {
  name.test = list.test[testi]
  res.testi = residus[[name.test]]
  order.arma.bef = data.frame(matrix(NA, ncol = 3, nrow = length(dat)))
  order.arma.aft = data.frame(matrix(NA, ncol = 3, nrow = length(dat)))
  coef.arma.bef = data.frame(matrix(NA, ncol = 4, nrow = length(dat)))
  coef.arma.aft = data.frame(matrix(NA, ncol = 4, nrow = length(dat)))
  for (i in c(1:length(dat))) {
    dat.i = res.testi[[i]]
    bef.residus = dat.i[1:(one.year*10)]
    aft.residus = dat.i[-c(1:(one.year*10))]
    bef.fit = fit.arima(bef.residus)
    aft.fit = fit.arima(aft.residus)
    order.arma.bef[i,] = bef.fit$pq
    order.arma.aft[i,] = aft.fit$pq
    coef.arma.bef[i,] = bef.fit$coef
    coef.arma.aft[i,] = aft.fit$coef
  }
  order.arma[[name.test]] <- list(order.arma.bef, order.arma.aft)
  coef.arma[[name.test]] <- list(coef.arma.bef, coef.arma.aft)
}
save(order.arma, file = paste0(path_results,"attribution/order.model.arma.restrict.RData"))
save(coef.arma, file = paste0(path_results,"attribution/coef.model.arma.restrict.RData"))

a = get(load(file = paste0(path_results,"attribution/order.model.arma.restrict.RData")))
list.model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)", "AR(2)", "MA(2)", "ARMA(1,2)", "ARMA(2,1)", "ARMA(2,2)")
# list.model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)")
# barplot of models------
model.iden <- function(order){
  model = c()
  if (identical(order, c(1,0,1))){ model = "ARMA(1,1)"}
  else if (identical(order, c(1,0,0))){ model = "AR(1)"}
  else if (identical(order, c(0,0,1))){ model = "MA(1)"}
  else if (identical(order, c(0,0,0))){ model = "White"}
  else if (identical(order, c(2,0,0))){ model = "AR(2)"}
  else if (identical(order, c(2,0,1))){ model = "ARMA(2,1)"}
  else if (identical(order, c(1,0,2))){ model = "ARMA(1,2)"}
  else if (identical(order, c(0,0,2))){ model = "MA(2)"}
  else if (identical(order, c(2,0,2))){ model = "ARMA(2,2)"}
  
  return(model)
}
length.data = nrow(a$gps.era[[1]])
six.model = data.frame(matrix(NA, ncol = 6, nrow = length.data))
for (i in 1:length(list.test)) {
  name.test = list.test[i]
  six.model[,i] = sapply(c(1:length.data), function(x) model.iden(as.numeric(unlist(a[[name.test]][[2]][x,]))))
}

# remove the duplicated 
name.ref = substr(names(dat) ,start = 1, stop = 15)

list.dup = which(duplicated(name.ref))
six.model[list.dup,1] <- NA
six.values = c()
for (i in 1:length(list.test)) {
  value.count = sapply(c(list.model), function(x) length(which(six.model[,i] == x)))
  six.values <- c( six.values, value.count)
}
res.plot = data.frame(series = rep(list.test, each = 9), mod = rep(list.model, 6), value = six.values*100/784)
res.plot$value[which(res.plot$series == "gps.era")] <- res.plot$value[which(res.plot$series == "gps.era")]*784/156
res.plot$series = factor(res.plot$series, 
                         levels=list.test)
res.plot$mod = factor(res.plot$mod, 
                      levels=list.model)


jpeg(paste0(path_results,"attribution/iden_model_order2.jpg" ),width = 3000, height = 1800,res = 300)
p <- ggplot(res.plot, aes(fill=mod, y=value, x=series)) + 
  geom_bar(position="dodge", stat="identity")+theme_bw()+ 
  xlab("") + ylab("Percentage of model")+
  theme(axis.text = element_text(size = 14),legend.text=element_text(size=12),
      axis.title = element_text(size=14))
  # theme(
  #   legend.title=element_blank(),
  #   legend.position = c(.5, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )
print(p)
dev.off()


# plot the residual 
dat.bef$gps.era = res.bef
jpeg(paste0(path_results,"attribution/iden_model.jpg" ),width = 2600, height = 1800,res = 300)
p <- ggplot(dat.bef, aes(x=date, y=gps.era)) + 
  geom_line(col="gray")+theme_bw()+ 
  xlab("") + ylab("Residual")+
  theme(axis.text = element_text(size = 14),legend.text=element_text(size=12),
        axis.title = element_text(size=14))
 
print(p)
dev.off()

scr = screen.O(Y = dat.bef , name.var = name.test, method = 'sigma', global.mu = 0, iter = 1, estimator = "Sca", fix.thres = 0, loes = 0, loes.method = 0)
sd0 = unlist(scr$sd.est[[1]])
dat.bef$sd0 = sd0
dat.bef = dat.bef[which(is.na(dat.bef$sd0)==FALSE),]
dat.bef$sd0 =  d
dat.bef$date = Y.with.NA$date

jpeg(paste0(path_results,"attribution/std1.jpg" ),width = 2600, height = 1800,res = 300)
p <- ggplot(dat.bef, aes(x=date, y=sd0)) + 
  geom_line(col="black")+theme_bw()+ 
  xlab("") + ylab("Moving standard deviation of GPS-ERA")+
  theme(axis.text = element_text(size = 16),legend.text=element_text(size=12),
        axis.title = element_text(size=20))

print(p)
dev.off()






# cdf of the coefficients 


ggplot(b, aes(x = value, col = variable))+
  stat_ecdf(lwd = 0.5)+ 
  scale_color_manual(values = c( "#000000", "#FF0000", "#00FF00", "#FFFF00", "#0000FF", "#00FFFF"))+
  labs(y = "CDF", x = "rho(1)")+ 
  theme_bw()






