# this prog is used to compare the FPR and TPR of the test of the change in mean between the OLS-HAC and the GLS
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"FGLS.R"))
# source(paste0(path_code_att,"newUsed_functions.R"))

# input: what do you want to test. Ex: TPR of test when data is AR(1) with different rho------------------
off.set = 0
heteroscedast = 0
autocor = 1
x.axis = "rho"
one.year = 365

y.axis = ifelse(off.set !=0, "TPR", "FPR")
nb.sim = 1000
n = 500
t = c(1:n)-n/2
list.param.ar = 0.9
list.param.sig = seq(0.1, 0.35, 0.05)

# specify simulation model
trend.sim = 0
# specify the regression model 
trend.reg = 10

# specify the condition
mod.sim <- function(heteroscedastic, autocorr, var.inno, list.param.sig, list.param.ar, x.axis, individual, n, T1){
  a = cos(2*pi*(c(1:n)/T1) -pi)
  if(heteroscedastic == 1 & autocorr == 0){
    hetero = 1
    burn.in = 0
    sigma.t = sapply(c(1:length(list.param.sig)), function(x) var.inno - a*list.param.sig[x])
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
      sigma.0 = var.inno - a*list.param.sig[individual]
      sigma.t = matrix( rep(sigma.0,  length(list.param.ar)) , ncol =  length(list.param.ar), byrow = FALSE )
      ar = list.param.ar
    } else if(x.axis == "sig.v"){
      sigma.t = sapply(c(1:length(list.param.sig)), function(x) var.inno - a*list.param.sig[x])
      ar = rep(list.param.ar[individual],length(list.param.sig))
    }
  }
  return(list(hetero = hetero, burn.in = burn.in, sigma.t = sigma.t, ar = ar))
}
gen.test = mod.sim(heteroscedastic = heteroscedast, autocorr = autocor, var.inno = 0.4, list.param.ar = list.param.ar, 
                   list.param.sig = list.param.sig, x.axis = x.axis, individual = 6, n = n, T1 = n/2)



res.var = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
Res.fin = data.frame(matrix(NA, ncol = 5, nrow = length(gen.test$ar)))

hetero = gen.test$hetero
burn.in = gen.test$burn.in
sigma.sim = gen.test$sigma.t[1]
time.c = c()
time.c1 = c()

# test simulation FGLS
for (l in c(1:length(gen.test$ar))) {
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  ar = gen.test$ar[l]
  
  for (i in c(1:nb.sim)) {
    set.seed(i)
    y = simulate.general(N = n, arma.model = c(ar,0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim),
                         monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + off.set
    df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2015-05-27"), by="days")[1:n])
    Data.mod = construct.design(data.df = df, name.series = "y", break.ind = n/2)
    Data.mod =  Data.mod[,c(1,10,11)]
    # Test with vcov (HAC)
    list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
    mod.X <-  list.para %>% stringr::str_c(collapse = "+")
    mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
    # ols
    ols.fit = lm(mod.expression, data = Data.mod)
    vcov.para=sandwich::kernHAC(ols.fit, prewhite = TRUE, approx = c("AR(1)"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
    
    fit.hac=lmtest::coeftest(ols.fit,df=(n-trend.reg),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.ols=lmtest::coeftest(ols.fit,df=(n-trend.reg))[, ] %>% as.data.frame()
    
    #GLS with true covariance matrix
    if(ar==0){
      gls.fit = ols.fit
    }else{
      gls.fit = GLS(phi = ar, theta = 0, var.t = rep(sigma.sim,n), design.matrix = Data.mod)
    }
    # FGLS
    start_time <- Sys.time()
    fgls.fit = FGLS1(design.m = Data.mod, tol=0.0001, day.list = df$date, noise.model = c(1,0,0))
    end_time <- Sys.time()
    time.c = c(time.c, (end_time - start_time))
    
    start_time <- Sys.time()
    fgls.fit1 = FGLS2(design.m = Data.mod, tol=0.0001, day.list = df$date, noise.model = c(1,0,0))
    end_time <- Sys.time()
    time.c1 = c(time.c1, (end_time - start_time))
    
    tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fgls = fgls.fit, fgls1 = fgls.fit1, gls = gls.fit)
    coef.res[[i]] = list( ols = ols.fit$coefficients, fgls = fgls.fit$coefficients,  fgls1 = fgls.fit1$coefficients, gls = gls.fit$Coefficients )
    var.res[[i]] = list( ols = vcov(ols.fit), fgls = fgls.fit$varBeta, fgls1 = fgls.fit1$varBeta, hac = vcov.para, gls = gls.fit$vcov)
    print(i)
  }
  # significance level
  pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[1,4]))
  pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[1,4]))
  pval.fgls <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fgls$t.table[1,4]))
  pval.fgls1 <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fgls1$t.table[1,4]))
  pval.gls <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$gls$t.table[1,4]))
  
  p.val = c(length(which(pval.ols>0.05)),
            length(which(pval.fgls>0.05)),
            length(which(pval.fgls1>0.05)),
            length(which(pval.gls>0.05)), 
            length(which(pval.hac>0.05)))
  
  Res.fin[l,] = p.val
  
}

save(tot.res, file=paste0(path_result,"attribution/tt.RData"))
a = data.frame(c1 = time.c, c2 = time.c1)
save(tot.res, file=paste0(path_result,"attribution/time.RData"))

colnames(Res.fin) <- c("OLS", "FGLS", "GLS", "OLS-HAC")
Res.fin <- Res.fin[c("OLS","OLS-HAC", "FGLS", "GLS")]

if(off.set == 0){
  res = (nb.sim- Res.fin)/nb.sim
}else{
  res = (nb.sim- Res.fin)/nb.sim
}
name.x = x.axis
if(x.axis == "rho"){
  param.test = list.param.ar
}else{
  param.test = list.param.sig*100/0.4
}
res[name.x] = param.test
dat.plot =reshape2::melt(res, id = name.x)
face1 = "bold"
x.axis1 = "p"

jpeg(paste0(path_results,"attribution/h.", heteroscedast, "a.", autocor, x.axis, y.axis, "trend.reg", trend.reg,"1.jpg" ),
     width = 2600, height = 1800,res = 300)
p2 <- eval(parse(
  text=paste0("ggplot(dat.plot, aes(x =", x.axis, ",y = value, col = variable))+
  geom_point(size=3) + geom_line() +theme_bw() +
  ylab(y.axis) + 
  xlab(x.axis1) +
  scale_y_continuous(breaks=seq(0, 1, 0.15), limits =c(0,1))+
  scale_x_continuous(breaks=list.param.ar )+
  theme(axis.text = element_text(size = 25),legend.text=element_text(size=15),legend.title=element_blank(),
      axis.title = element_text(size=25,face=face1))")))
print(p2)
dev.off()