# this prog is used to compare the FPR and TPR of the test of the change in mean between the OLS-HAC and the GLS
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"newUsed_functions.R"))

# input: what do you want to test. Ex: TPR of test when data is AR(1) with different rho------------------
off.set = 0.3
heteroscedast = 1
autocor = 0
x.axis = "sig.v"

y.axis = ifelse(off.set !=0, "TPR", "FPR")
nb.sim = 1000
n = 200
t = c(1:n)-n/2
list.param.ar = seq(0, 0.9, 0.15)
list.param.sig = seq(0.1, 0.35, 0.05)

# specify simulation model
trend.sim = 0
# specify the regression model 
trend.reg = 0
mod.expression = ifelse(trend.reg == 1, "signal~jump+Xt", "signal~jump")

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
Res.fin = data.frame(matrix(NA, ncol = 4, nrow = length(gen.test$ar)))

hetero = gen.test$hetero
burn.in = gen.test$burn.in

for (l in c(1:length(gen.test$ar))) {
  tot.res <- list()
  coef.res <- list()
  var.res <- list()
  ar = gen.test$ar[l]
  sigma.sim = gen.test$sigma.t[,l]
  for (i in c(1:nb.sim)) {
    set.seed(i)
    y = simulate.general(N = n, arma.model = c(ar,0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim),
                         monthly.var = 0)
    y[(n/2):n] <- y[(n/2):n] + off.set
    Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), weight1 = sigma.sim, t = t, Xt = t)
    
    # Test with vcov (HAC)
    ols.fit = lm(mod.expression, data = Data.mod)
    if(ar ==0 & hetero == 0){
      vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE, kernel = "Quadratic Spectral",adjust = TRUE, sandwich = TRUE)
      gls.fit <- eval(parse(text=paste0("gls(",mod.expression,",data=Data.mod,correlation = NULL, na.action=na.omit,weights=NULL",")")))
    }else if(ar ==0 & hetero != 0){
      vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
      gls.fit <- eval(parse(text=paste0("gls(",mod.expression,",data=Data.mod,correlation = NULL, na.action=na.omit,weights=varFixed(value = ~weight1)",")")))
    }else if(ar != 0 & hetero != 0){
      vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE, approx = c("AR(1)"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
      gls.fit <- eval(parse(text=paste0("gls(",mod.expression,",data=Data.mod,correlation =  corAR1(ar, form = ~ 1), na.action=na.omit,weights=varFixed(value = ~weight1)",")")))
    }else if(ar != 0 & hetero == 0){
      vcov.para=sandwich::kernHAC(ols.fit,prewhite = FALSE, approx = c("AR(1)"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
      gls.fit <- eval(parse(text=paste0("gls(",mod.expression,",data=Data.mod,correlation =  corAR1(ar, form = ~ 1), na.action=na.omit,weights=NULL",")")))
    }
    
    fit.hac=lmtest::coeftest(ols.fit,df=(n-2-trend.reg),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.ols=lmtest::coeftest(ols.fit,df=(n-2-trend.reg))[, ] %>% as.data.frame()
    
    # GLS from nlme package
    fit.gls =lmtest::coeftest(gls.fit,df=(n-2-trend.reg))[, ] %>% as.data.frame()
    
    #GLS with true covariance matrix
    gls.fit.true = gls.true(var.t = Data.mod$weight1, phi = ar, theta = 0, design.matrix = Data.mod, trend = trend.reg)
    
    tot.res[[i]] = list(fit.hac =fit.hac, fit.ols = fit.ols, fit.gls.true = gls.fit.true$fit.gls, fit.gls = fit.gls)
    coef.res[[i]] = list( ols = ols.fit$coefficients, gls = gls.fit$coefficients, gls.t = gls.fit.true$Coefficients)
    var.res[[i]] = list( ols = vcov(ols.fit), gls = gls.fit$varBeta, hac = vcov.para, gls.t = gls.fit.true$vcov)
    
  }
  # significance level
  pval.ols <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.ols[2,4]))
  pval.hac <- unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.hac[2,4]))
  pval.gls.true <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls.true[2,4]))
  pval.gls <-  unlist(sapply(c(1:nb.sim), function(x) tot.res[[x]]$fit.gls[2,4]))
  
  p.val = c(length(which(pval.ols>0.05)),
            length(which(pval.gls.true>0.05)), 
            length(which(pval.gls>0.05)), 
            length(which(pval.hac>0.05)))
  Res.fin[l,] = p.val
}

colnames(Res.fin) <- c("ols", "gls.t","gls", "hac")
if(off.set == 0){
  res = (nb.sim- Res.fin)/nb.sim
}else{
  res = (nb.sim- Res.fin)/nb.sim
}
name.x = x.axis
if(x.axis == "rho"){
  param.test = list.param.ar
}else{
  param.test = list.param.sig
}
res[name.x] = param.test
dat.plot =reshape2::melt(res, id = name.x)
face1 = "bold"
jpeg(paste0(path_results,"attribution/h.", heteroscedast, "a.", autocor, x.axis, y.axis, "trend.reg", trend.reg,"1.jpg" ),
     width = 2600, height = 1800,res = 300)
p2 <- eval(parse(
  text=paste0("ggplot(dat.plot, aes(x =", x.axis, ",y = value, col = variable))+
  geom_point(size=3) + theme_bw() +
  ylab(y.axis) +
  theme(axis.text = element_text(size = 18),legend.text=element_text(size=15),
      axis.title = element_text(size=18,face=face1))")))
print(p2)
dev.off()


a = cos(2*pi*(c(1:n)/T1) )
v = rep(NA,6)
for (l in c(1:6)) {
  s =  0.4 - a*list.param.sig[l]
  v[l] =  1/(sum(1/s))
}

# confusion table for all possible models and regressions -----------------
n = 200
T1 = n/2
nb.sim = 1000
a = cos(2*pi*(c(1:n)/T1))
var.m = 0.4
var.t = var.m - 0.35*a
off.set = 0.3
ar0 = 0.3
ar1 = ar0
trend.reg = 0
t = c(1:n)-n/2
type.hc = "HC3"
kernel1 = "Quadratic Spectral"
approx1 = c("ARMA(1,1)")
mod.expression = ifelse(trend.reg == 1, "signal~jump+Xt", "signal~jump")

name.sim <- function(name.model, var.m, var.t, ar0){
  if(name.model == "IID"){
    ar = 0
    burn.in = 0
    hetero = 0
    sigma = sqrt(var.m)
  } else if (name.model == "AR(1)"){
    ar = ar0 
    burn.in = 1000
    hetero = 0
    sigma = sqrt(var.m)
  } else if (name.model == "hetero"){
    ar = 0
    burn.in = 0
    hetero = 1
    sigma = var.t 
  } else if (name.model == "hetero+AR(1)"){
    ar = ar0 
    burn.in = 1000
    hetero = 1
    sigma = var.t 
  }
  
  return(list(ar = ar, hetero = hetero, burn.in = burn.in, sigma.t = sigma))
}
name.est <- function(name.model, method.lm, mod.expression, Data.mod, ar1, trend.reg){
  fit <- NULL
  fit.test <- NULL
  vcov.para <- NULL
  if(method.lm == "OLS" | method.lm == "OLS-HAC"){
    fit = "lm(mod.expression, data = Data.mod)"
  } else if(method.lm == "GLS-nmle"){
    if(name.model == "IID"){
      fit <- paste0("gls(",mod.expression,",data=Data.mod,correlation = NULL, na.action=na.omit,weights=NULL",")")
    } else if(name.model == "hetero"){
      fit <- paste0("gls(",mod.expression,",data=Data.mod,correlation = NULL, na.action=na.omit,weights=varFixed(value = ~weight1)",")")
    } else if (name.model == "AR(1)"){
      fit <- paste0("gls(",mod.expression,",data=Data.mod,correlation =  corAR1(ar1, form = ~ 1), na.action=na.omit,weights=NULL",")")
    } else if (name.model == "hetero+AR(1)"){
      fit <- paste0("gls(",mod.expression,",data=Data.mod,correlation =  corAR1(ar1, form = ~ 1), na.action=na.omit,weights=varFixed(value = ~weight1)",")")
    }
  } else if (method.lm == "GLS-true"){
    if(name.model == "IID" | name.model == "hetero"){
      ar1 = 0
    } else { ar1 = ar1}
    fit = paste0("gls.true(var.t = Data.mod$weight1, phi = ", ar1, ", theta = 0, design.matrix = Data.mod, trend = trend.reg)")
  }
  if( method.lm == "OLS-HAC"){
    if (name.model == "AR(1)"|name.model == "hetero+AR(1)"){
      vcov.para= paste0("sandwich::kernHAC(fit,prewhite = 1, approx =  approx1, kernel = kernel1 ,adjust = TRUE, sandwich = TRUE)")
    }else if (name.model == "hetero") {
      vcov.para= paste0("sandwich::vcovHC(fit, type = type.hc)")
    }else{
      vcov.para= paste0("vcov(fit)")
    }
    fit.test = "lmtest::coeftest(fit,df=(n-2-trend.reg),vcov.=vcov.para)[, ] %>% as.data.frame()"
  } else if (method.lm == "GLS-true" ){
    fit.test = "fit$fit.gls"
  } else{
    fit.test="lmtest::coeftest(fit,df=(n-2-trend.reg))[, ] %>% as.data.frame()"
  }
  
  return(list(fit.call = fit, vcov.call = vcov.para, fit.test.call = fit.test))
}

#   CHECH FUNCTION IS TRUE BEFORE RUNG ITERATIONS 
# iteration
list.method = c("OLS", "OLS-HAC", "GLS-true","GLS-nmle")
list.model = c( "IID", "AR(1)", "hetero", "hetero+AR(1)")
Total.res = list()
for (l in c(1:4)){
  # specifify the model to simulate 
  sim.mod = name.sim(list.model[l], var.m , var.t, ar0)
  all.res = list()
  for(s in c(1:4)) {
    est.mod = list.model[s]
    tot.res = data.frame(matrix(NA, ncol = 4, nrow = nb.sim))
    for (m in c(2)) {
      # specify model and method to estimate 
      b = name.est(name.model = est.mod, method.lm = list.method[m], mod.expression, Data.mod, ar1 = ar0, trend.reg)
      for (i in c(1:nb.sim)) {
        set.seed(i)
        y = simulate.general(N = n, arma.model = c(sim.mod$ar,0), burn.in = sim.mod$burn.in, hetero = sim.mod$hetero, sigma = sqrt(sim.mod$sigma.t),
                             monthly.var = 0)
        y[(n/2):n] <- y[(n/2):n] + off.set
        Data.mod = data.frame(signal = y, jump = rep(c(0,1), each = n/2), weight1 = sim.mod$sigma.t, t = t, Xt = t)
        fit =  eval(parse(text = b$fit.call))

        if(is.null(b$vcov.call) == FALSE){
          vcov.para = eval(parse(text = b$vcov.call))
        }
        # else{ vcov.para= vcov(fit)}
        test.res = eval(parse(text = b$fit.test.call))
        tot.res[i,m] = test.res[2,4]
      }
    }
    all.res[[s]] = tot.res
  }
  Total.res[[l]] = all.res
}

res.ols = data.frame(matrix(NA, ncol = 4, nrow = 4))
res.olshac = data.frame(matrix(NA, ncol = 4, nrow = 4))
res.glstrue = data.frame(matrix(NA, ncol = 4, nrow = 4))
res.glsnmle = data.frame(matrix(NA, ncol = 4, nrow = 4))
for (j in c(1:4)) {
  for (i in c(1:4)) {
    m = Total.res[[j]][[i]]
    all.met = sapply(c(1:4), function(x) length(which(m[,x] > 0.05)))
    res.ols[j,i] =  all.met[1]/nb.sim
    res.olshac[j,i] =  all.met[2]/nb.sim
    res.glstrue[j,i] =  all.met[3]/nb.sim
    res.glsnmle[j,i] =  all.met[4]/nb.sim
  }
}

write.table( res.ols, 
            file = paste0(path_results, "attribution/sim.OLS", trend.reg, off.set,".txt"),
            sep = "\t", quote = FALSE)
write.table( res.olshac, 
             file = paste0(path_results, "attribution/sim.OLSHAC", trend.reg, off.set,".txt"),
             sep = "\t", quote = FALSE)
write.table( res.glsnmle, 
             file = paste0(path_results, "attribution/sim.GLSnmle", trend.reg, off.set,".txt"),
             sep = "\t", quote = FALSE)
write.table( res.glstrue, 
             file = paste0(path_results, "attribution/sim.GLStrue", trend.reg, off.set,".txt"),
             sep = "\t", quote = FALSE)




