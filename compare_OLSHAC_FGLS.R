# this prog is used to compare the FPR and TPR of the test of the change in mean between the OLS-HAC and the GLS
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"FGLS.R"))
# source(paste0(path_code_att,"newUsed_functions.R"))

# input: what do you want to test. Ex: TPR of test when data is AR(1) with different rho------------------
one.year = 200
nb.sim = 1000
# n = 
list.param.ar = seq(0, 0.9, 0.15)
list.ma = rep(0, length(list.param.ar))
list.param.sig = seq(0,0.8,0.2)

# specify the condition
mod.sim <- function(heteroscedastic, autocorr, var.inno, list.param.sig, list.param.ar, x.axis, individual, n, T1, noise.name){
  periodic = cos(2*pi*(c(1:n)/T1))
  if(heteroscedastic == 1 & autocorr == 0){
    hetero = 1
    burn.in = 0
    sigma.t = sapply(c(1:length(list.param.sig)), function(x) var.inno - periodic*list.param.sig[x])
    ar = rep(0, length(list.param.sig))
  } else if(heteroscedastic == 0 & autocorr == 1){
    hetero = 0
    burn.in = 1000
    sigma.t = rep(var.inno, length(list.param.ar))
    ar = list.param.ar
  } else if (heteroscedastic == 1 & autocorr == 1){
    hetero = 1
    burn.in = 1000
    if(x.axis == "rho"){
      sigma.0 = var.inno - periodic*list.param.sig[individual]
      sigma.t = matrix( rep(sigma.0,  length(list.param.ar)) , ncol =  length(list.param.ar), byrow = FALSE )
      ar = list.param.ar
    } else if(x.axis == "sig"){
      sigma.t = sapply(c(1:length(list.param.sig)), function(x) var.inno - periodic*list.param.sig[x])
      ar = rep(list.param.ar[individual],length(list.param.sig))
    }
  }
  return(list(hetero = hetero, burn.in = burn.in, sigma.t = sigma.t, ar = ar))
}
simu_performance <- function(off.set, heteroscedast, autocor, x.axis, nb.sim, list.ma, 
                             list.param.ar, list.param.sig, noise.model,noise.name,n,length.wind){
  # generate params
  y.axis = ifelse(off.set !=0, "TPR", "FPR")
  thres =  ifelse(off.set !=0, 0.95, 0.05)
  if(x.axis=="rho"){
    sample = list.param.ar
    individual = 5
  }else{
    sample = list.param.sig
    individual = 5
  }
  gen.test = mod.sim(heteroscedast = heteroscedast, autocorr = autocor, var.inno = 1, list.param.ar = list.param.ar, 
                     list.param.sig = list.param.sig, x.axis = x.axis, individual = individual, n = n, T1 = n/2)
  burn.in = gen.test$burn.in
  hetero = gen.test$hetero
  
  total = list()
  Res.fin = data.frame(matrix(NA, ncol = 4, nrow = length(sample)))
  time.c = c()
  for (l in c(1:length(sample))) {
    tot.res <- list()
    coef.res <- list()
    var.res <- list()
    if(heteroscedast ==0){
      sigma.sim = gen.test$sigma.t[l]
    }else{
      sigma.sim = gen.test$sigma.t[,l]
    }
    ar = gen.test$ar[l]
    ma = list.ma[l]
    for (i in c(1:nb.sim)) {
      set.seed(i)
      if(length(sigma.sim)==1){
        sigma.sim = rep(sigma.sim, n)
      }
      y = simulate.general1(N = n, arma.model = c(ar,ma), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      y[(n/2):n] <- y[(n/2):n] + off.set
      # ggsave(paste0(path_results,"attribution/test.jpg" ), plot = plot(y), width = 5, height = 5, units = "cm", dpi = 1200)
      df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:n])
      Data.mod = construct.design(data.df = df, name.series = "y", break.ind = n/2, one.year)
      Data.mod =  Data.mod[,c(1,10,11)]
      # Test with vcov (HAC)
      list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
      mod.X <-  list.para %>% stringr::str_c(collapse = "+")
      mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
      # ols
      ols.fit = lm(mod.expression, data = Data.mod)
      fit.ols=lmtest::coeftest(ols.fit,df=(n-2))[, ] %>% as.data.frame()
      
      #GLS, HAC with true covariance matrix
      gls.fit = GLS(phi = ar, theta = ma, var.t = sigma.sim, design.matrix = Data.mod)
        
      vcov.para=sandwich::kernHAC(ols.fit, prewhite = TRUE, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
      fit.hac=lmtest::coeftest(ols.fit,df=(n-2),vcov.=vcov.para)[, ] %>% as.data.frame()
      #FGLS
      start_time <- Sys.time()
      fgls.fit = FGLS1(design.m = Data.mod, tol=0.01, day.list = df$date, noise.model = noise.model, length.wind)
      end_time <- Sys.time()
      time.c = c(time.c, (end_time - start_time))
      
      tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fgls = fgls.fit, gls = gls.fit)
      coef.res[[i]] = list( ols = ols.fit$coefficients, fgls = fgls.fit$coefficients,  gls = gls.fit$Coefficients )
      var.res[[i]] = list( ols = vcov(ols.fit), fgls = fgls.fit$varBeta, hac = vcov.para, gls = gls.fit$vcov)
      print(i)
    }
    # significance level
    pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[1,4]))
    pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[1,4]))
    pval.fgls <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fgls$t.table[1,4]))
    pval.gls <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$gls$t.table[1,4]))
    
    p.val = c(length(which(pval.ols>0.05)),
              length(which(pval.fgls>0.05)),
              length(which(pval.gls>0.05)), 
              length(which(pval.hac>0.05)))
    
    Res.fin[l,] = p.val
    total[[l]] = tot.res
  }
  colnames(Res.fin) <- c("OLS", "FGLS", "GLS", "OLS-HAC")
  Res.fin <- Res.fin[c("OLS","OLS-HAC", "FGLS", "GLS")]
  if(y.axis=="FPR"){
    Res.fin$GLS[1] = Res.fin$OLS[1]
  }
  
  Res = list(p = Res.fin, total = total, time= time.c)
  save(Res, file = paste0(path_results,"attribution/",heteroscedast,"auto",autocor, x.axis, y.axis,noise.name,n,"1R.Data"))
  
  res = (nb.sim- Res.fin)/nb.sim
  name.x = x.axis
  if(x.axis == "rho"){
    param.test = list.param.ar
    x.axis1 = expression(Phi)
  }else{
    param.test = list.param.sig*100/0.5
    x.axis1 = "Range of variance(%) "
  }
  res[name.x] = param.test
  dat.plot =reshape2::melt(res, id = name.x)
  face1 = "bold"
  
  
  p2 <- eval(parse(
    text=paste0("ggplot(dat.plot, aes(x =", x.axis, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis) + 
  xlab(x.axis1) + geom_hline(yintercept =", thres,", lwd = 0.2)+
  scale_y_continuous(breaks=seq(0, 1, 0.15), limits =c(0,1))+
  scale_x_continuous(breaks=list.param.ar )+
  theme(axis.text = element_text(size = 5), legend.text=element_text(size=4.5),
        axis.title = element_text(size = 5), legend.title=element_blank())")))
  
  ggsave(paste0(path_results,"attribution/test_sim_h",heteroscedast,"auto",autocor, x.axis, y.axis,noise.name,n,"1.jpg" ), plot = p2, width = 8, height = 5, units = "cm", dpi = 1200)
  
  return(Res)
  
}

simu_performance1 <- function(off.set, heteroscedast, autocor, x.axis, nb.sim, list.ma, 
                             list.param.ar, list.param.sig, noise.model,noise.name,n){
  # generate params
  y.axis = ifelse(off.set !=0, "TPR", "FPR")
  thres =  ifelse(off.set !=0, 0.95, 0.05)
  if(x.axis=="rho"){
    sample = list.param.ar
    individual = 5
  }else{
    sample = list.param.sig
    individual = 5
  }
  gen.test = mod.sim(heteroscedast = heteroscedast, autocorr = autocor, var.inno = 1, list.param.ar = list.param.ar, 
                     list.param.sig = list.param.sig, x.axis = x.axis, individual = individual, n = n, T1 = n/2)
  burn.in = gen.test$burn.in
  hetero = gen.test$hetero
  
  total = list()
  Res.fin = data.frame(matrix(NA, ncol = 4, nrow = length(sample)))
  time.c = c()
  for (l in c(1:length(sample))) {
    tot.res <- list()
    coef.res <- list()
    var.res <- list()
    if(heteroscedast ==0){
      sigma.sim = gen.test$sigma.t[l]
    }else{
      sigma.sim = gen.test$sigma.t[,l]
    }
    ar = gen.test$ar[l]
    ma = list.ma[l]
    for (i in c(1:nb.sim)) {
      set.seed(i)
      if(length(sigma.sim)==1){
        sigma.sim = rep(sigma.sim, n)
      }
      y = simulate.general1(N = n, arma.model = c(ar,ma), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
      y[(n/2):n] <- y[(n/2):n] + off.set
      # ggsave(paste0(path_results,"attribution/test.jpg" ), plot = plot(y), width = 5, height = 5, units = "cm", dpi = 1200)
      df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:n])
      Data.mod = construct.design(data.df = df, name.series = "y", break.ind = n/2, one.year)
      Data.mod =  Data.mod[,c(1,10,11)]
      # Test with vcov (HAC)
      list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
      mod.X <-  list.para %>% stringr::str_c(collapse = "+")
      mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
      # ols
      ols.fit = lm(mod.expression, data = Data.mod)
      fit.ols=lmtest::coeftest(ols.fit,df=(n-2))[, ] %>% as.data.frame()
      
      #GLS, HAC with true covariance matrix
      gls.fit = GLS(phi = ar, theta = ma, var.t = sigma.sim, design.matrix = Data.mod)
      
      vcov.para=sandwich::kernHAC(ols.fit, prewhite = TRUE, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
      fit.hac=lmtest::coeftest(ols.fit,df=(n-2),vcov.=vcov.para)[, ] %>% as.data.frame()
      #FGLS
      start_time <- Sys.time()
      fgls.fit = FGLS1(design.m = Data.mod, tol=0.01, day.list = df$date, noise.model = noise.model)
      end_time <- Sys.time()
      time.c = c(time.c, (end_time - start_time))
      
      tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fgls = fgls.fit, gls = gls.fit)
      coef.res[[i]] = list( ols = ols.fit$coefficients, fgls = fgls.fit$coefficients,  gls = gls.fit$Coefficients )
      var.res[[i]] = list( ols = vcov(ols.fit), fgls = fgls.fit$varBeta, hac = vcov.para, gls = gls.fit$vcov)
      print(i)
    }
    # significance level
    pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[1,4]))
    pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[1,4]))
    pval.fgls <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fgls$t.table[1,4]))
    pval.gls <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$gls$t.table[1,4]))
    
    p.val = c(length(which(pval.ols>0.05)),
              length(which(pval.fgls>0.05)),
              length(which(pval.gls>0.05)), 
              length(which(pval.hac>0.05)))
    
    Res.fin[l,] = p.val
    total[[l]] = tot.res
  }
  colnames(Res.fin) <- c("OLS", "FGLS", "GLS", "OLS-HAC")
  Res.fin <- Res.fin[c("OLS","OLS-HAC", "FGLS", "GLS")]
  if(y.axis=="FPR"){
    Res.fin$GLS[1] = Res.fin$OLS[1]
  }
  
  Res = list(p = Res.fin, total = total, time= time.c)
  save(Res, file = paste0(path_results,"attribution/",heteroscedast,"auto",autocor, x.axis, y.axis,noise.name,n,"1R.Data"))
  
  res = (nb.sim- Res.fin)/nb.sim
  name.x = x.axis
  if(x.axis == "rho"){
    param.test = list.param.ar
    x.axis1 = expression(Phi)
  }else{
    param.test = list.param.sig*100/0.5
    x.axis1 = "Range of variance(%) "
  }
  res[name.x] = param.test
  dat.plot =reshape2::melt(res, id = name.x)
  face1 = "bold"
  
  
  p2 <- eval(parse(
    text=paste0("ggplot(dat.plot, aes(x =", x.axis, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis) + 
  xlab(x.axis1) + geom_hline(yintercept =", thres,", lwd = 0.2)+
  scale_y_continuous(breaks=seq(0, 1, 0.15), limits =c(0,1))+
  scale_x_continuous(breaks=list.param.ar )+
  theme(axis.text = element_text(size = 5), legend.text=element_text(size=4.5),
        axis.title = element_text(size = 5), legend.title=element_blank())")))
  
  ggsave(paste0(path_results,"attribution/test_sim_h",heteroscedast,"auto",autocor, x.axis, y.axis,noise.name,n,"1.jpg" ), plot = p2, width = 8, height = 5, units = "cm", dpi = 1200)
  
  return(Res)
  
}

 # a=simu_performance(off.set, heteroscedast, autocor, x.axis, nb.sim, list.param.ar, list.param.sig, noise.model=c(1,0,0))

# test simulation FGLS-------------------------
case = data.frame(h = rep(c(0,1,1,1),2),
                  a = rep(c(1,0,1,1),2), 
                  offset = rep(c(0, 0.3), each = 4), 
                  x.axis = rep(c("rho", "sig"),4))


for (j in c(3:4)) {
    off.set = case$offset[j]
    heteroscedast = case$h[j]
    autocor = case$a[j]
    x.axis = as.character(case$x.axis[j])
    if(autocor==0){
      noise.list = c(0,0,0)
      model.name = "white"
    }else{
      noise.list = c(1,0,0)
      model.name = "ar1"
    }
    a = simu_performance(off.set, heteroscedast, autocor, x.axis, nb.sim=nb.sim, list.param.ar, list.param.sig,list.ma = list.ma,
                         noise.model=noise.list, noise.name = model.name)
    print(j)
}
noise.list = c(1,0,1)
model.name = "arma1"

a = simu_performance(off.set=0, heteroscedast = 1, autocor = 0, x.axis = "sig", nb.sim=nb.sim, list.param.ar, list.param.sig,list.ma = list.ma,
                     noise.model=c(0,0,0), noise.name = "white", length.wind = 60, n =400)
# test the length of series ------------------------------------

case = data.frame(h = rep(1,4),
                  a = rep(1,4), 
                  n = c(200,400,600,800),
                  x.axis = rep("rho",4))


off.set = 0
heteroscedast = 0
autocor = 1 

for (j in c(3:4)) {
  x.axis = as.character(case$x.axis[j])
  noise.list = c(1,0,0)
  model.name = "ar1"
  n1 = case$n[j]
  one.year = n1/2
  a = simu_performance(off.set, heteroscedast, autocor, x.axis, nb.sim=nb.sim, list.param.ar, list.param.sig,list.ma = list.ma,
                       noise.model=noise.list, noise.name = model.name, n=n1)
  print(j)
}




