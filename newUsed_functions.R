Test_OLS_vcovhac <- function(Data.mod){
  # OLS estimates
  list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
  mod.X <-  list.para %>% str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X) %>% str_c(collapse = "")
  
  fit.ols <- lm(mod.expression,data=Data.mod)
  names(fit.ols$coefficients)[which(names(fit.ols$coefficients)=="JumpRight")]="Jump"
  
  # covariance matrix HAC
  # L <- round(dim(Data.mod)[1]^(1/4))
  # vcov.para=sandwich::NeweyWest(fit.ols, lag = L)
  vcov.para=sandwich::kernHAC(fit.ols,prewhite = 1,kernel = "Quadratic Spectral",approx = "AR(1)", sandwich = TRUE)
  
  # Test with vcov (HAC)
  fit.hac=lmtest::coeftest(fit.ols,df=Inf,vcov.=vcov.para)[, ] %>% as.data.frame()
  
  
  # Selection parameter by parameter (the intercept will never be removed)
  threshold <- 0.05
  fit.hac.without.NA <- fit.hac[!(rownames(fit.hac) %in% "(Intercept)"),]
  pval.hac.without.intercept=fit.hac.without.NA$`Pr(>|z|)`
  p.max <- max(pval.hac.without.intercept)
  
  while ((p.max>threshold) & (length(list.para)>1)){
    names.max <- rownames(fit.hac.without.NA)[pval.hac.without.intercept==p.max]
    list.para <-list.para[!(list.para %in% names.max)] 
    mod.X <-  list.para %>% str_c(collapse = "+")
    mod.expression <- c("signal","~",mod.X) %>% str_c(collapse = "")
    
    fit.ols <- lm(eval(mod.expression),data=Data.mod)
    vcov.para=sandwich::kernHAC(fit.ols,prewhite = 1,kernel = "Quadratic Spectral",approx = "AR(1)", sandwich = TRUE)
    fit.hac=lmtest::coeftest(fit.ols,df=Inf,vcov.=vcov.para)[, ] %>% as.data.frame()
    # row.names(fit.hac)[which(row.names(fit.hac)=="JumpRight")]="Jump"
    
    fit.hac.without.NA <- fit.hac[!(rownames(fit.hac) %in% "(Intercept)"),]
    pval.hac.without.intercept=fit.hac.without.NA$`Pr(>|z|)`
    p.max <- max(pval.hac.without.intercept)
  }
  return(list(fit.hac = fit.hac, fit.ols = fit.ols, vcov.para = vcov.para, predicted = fit.ols$fitted.values))
}


Test_FGLS <- function(Data.mod){
  list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
  mod.X <-  list.para %>% str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X) %>% str_c(collapse = "")
  
  fit.gls <- eval(parse(text=paste0("gls(",mod.expression,",data=,Data.mod,correlation = corARMA(p=1,q=1), na.action=na.omit,weights=varFunc(~cos1+sin1)",")")))
  vcov.gls=fit.gls$varBeta
  res.gls=lmtest::coeftest(fit.gls,df=Inf,vcov.=vcov.gls)[, ] %>% as.data.frame()
  row.names(res.gls)[which(row.names(res.gls)=="JumpRight")]="Jump"


  # Selection parameter by parameter (the intercept will never be removed)
  threshold <- 0.01
  fit.gls.without.NA <- res.gls[!(rownames(res.gls) %in% "(Intercept)"),]
  pval.gls.without.intercept=fit.gls.without.NA$`Pr(>|z|)`
  p.max <- max(pval.gls.without.intercept)
  
  while ((p.max>threshold) & (length(list.para)>=1)){
    names.max <- rownames(fit.gls.without.NA)[pval.gls.without.intercept==p.max]
    list.para <-list.para[!(list.para %in% names.max)] 
    mod.X <-  list.para %>% str_c(collapse = "+")
    mod.expression <- c("signal","~",mod.X) %>% str_c(collapse = "")
    
    
    fit.gls <- eval(parse(text=paste0("gls(",mod.expression,",data=,Data.mod,correlation = corARMA(p=1,q=1), na.action=na.omit,weights=varFunc(~cos1+sin1)",")")))
    vcov.gls=fit.gls$varBeta
    res.gls=lmtest::coeftest(fit.gls,df=Inf,vcov.=vcov.gls)[, ] %>% as.data.frame()
    row.names(res.gls)[which(row.names(res.gls)=="JumpRight")]="Jump"
    
    fit.gls.without.NA <- res.gls[!(rownames(res.gls) %in% "(Intercept)"),]
    pval.gls.without.intercept=fit.gls.without.NA$`Pr(>|z|)`
    p.max <- max(pval.gls.without.intercept)
    
  }
  return(list(res.gls=res.gls,predicted=fit.gls$fitted))
}

plot.predict <- function(res,Data.mod,one.year,date){
  na.Y <- which(is.na(Data.mod$signal))
  signal.predicted <- rep(NA,dim(Data.mod)[1])
  signal.predicted[!((1:dim(Data.mod)[1]) %in% na.Y)] <- res$predicted

  Left <- 1:one.year
  Right <- (one.year+1):(2*one.year)

  plot(date,Data.mod$signal,type="n",xlab="dates",ylab="gps-era")
  lines(date[Left],Data.mod$signal[Left],col="red", type="l",xlab="time",ylab="signal")
  lines(date[Left],signal.predicted[Left],col="red")
  lines(date[Right],Data.mod$signal[Right],col="blue")
  lines(date[Right],signal.predicted[Right],col="blue")
}


Test_OLS_vcovhac_1step <- function(Data.mod){
  # OLS estimates
  list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
  mod.X <-  list.para %>% str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X) %>% str_c(collapse = "")
  
  fit.ols <- lm(mod.expression,data=Data.mod)
  names(fit.ols$coefficients)[which(names(fit.ols$coefficients)=="JumpRight")]="Jump"
  
  # covariance matrix HAC
  L <- round(dim(Data.mod)[1]^(1/4))
  vcov.para=sandwich::NeweyWest(fit.ols, lag = L)
  
  # Test with vcov (HAC)
  fit.hac=lmtest::coeftest(fit.ols,df=Inf,vcov.=vcov.para)[, ] %>% as.data.frame()
  
  
  # Selection parameter by parameter (the intercept will never be removed)
  # threshold <- 0.01
  # fit.hac.without.NA <- fit.hac[!(rownames(fit.hac) %in% "(Intercept)"),]
  # pval.hac.without.intercept=fit.hac.without.NA$`Pr(>|z|)`
  # p.max <- max(pval.hac.without.intercept)
  
  # while ((p.max>threshold) & (length(list.para)>1)){
  #   names.max <- rownames(fit.hac.without.NA)[pval.hac.without.intercept==p.max]
  #   list.para <-list.para[!(list.para %in% names.max)] 
  #   mod.X <-  list.para %>% str_c(collapse = "+")
  #   mod.expression <- c("signal","~",mod.X) %>% str_c(collapse = "")
  #   
  #   fit.ols <- lm(eval(mod.expression),data=Data.mod)
  #   vcov.para=sandwich::NeweyWest(fit.ols, lag = L)
  #   fit.hac=lmtest::coeftest(fit.ols,df=Inf,vcov.=vcov.para)[, ] %>% as.data.frame()
  #   row.names(fit.hac)[which(row.names(fit.hac)=="JumpRight")]="Jump"
  #   
  #   fit.hac.without.NA <- fit.hac[!(rownames(fit.hac) %in% "(Intercept)"),]
  #   pval.hac.without.intercept=fit.hac.without.NA$`Pr(>|z|)`
  #   p.max <- max(pval.hac.without.intercept)
  # }
  return(list(fit.hac=fit.hac,predicted=fit.ols$fitted.values, fit.ols = fit.ols, vcov.para = vcov.para))
}

plot_HAC <- function(case.name, res.i, data.in, name.var, ver, add.subtitle ){
  Y.with.na = data.in
  Y.without.na = Y.with.na[which(is.na(Y.with.na[name.var])==FALSE),]
  start.day = min(Y.without.na$date)
  end.day = max(Y.without.na$date)
  print(paste(start.day, end.day))
  Y = Y.with.na[which(Y.with.na$date >= start.day & Y.with.na$date <= end.day),]
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  Y$predict = rep(NA, nrow(Y))
  Y$predict[which(is.na(Y[name.var]) == FALSE)] = as.numeric(res.i$predicted)
  data.plot <- data.frame(date = Y$date, name.var = Y[name.var], predicted = Y$predict)
  a = reshape2::melt(data.plot, id = "date")
  subtitle.plot = paste0( paste(rownames(res.i$fit.hac),  
                        round(res.i$fit.hac$Estimate, digits = 4), 
                        round(res.i$fit.hac$`Pr(>|z|)`, digits = 4), 
                        collapse = ", ", sep = "; "), "  VIF: ", add.subtitle)
  
  
  jpeg(paste0(path_results,"attribution/OLS.HAC/", case.name, ver, ".jpeg" ),
       width = 3000, height = 1500,res = 300)
  p <- ggplot(data = a, aes(x = date, y = value, col = variable)) + 
    geom_line()+theme_bw()+
    geom_vline(xintercept = Y$date[which(Y$date == breakpoint)], size = 0.2)+
    ylab(name.var)+
    labs(title = case.name, 
         subtitle = wrapper(subtitle.plot, 130))+ 
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size=15,face="bold"))
  print(p)
  dev.off()
  
}

wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

var_cova <- function(Sig, phi, theta, n){
  M = 9*n
  # compute the variance of et, 
  Sig.k = Sig[2:(M-n)] 
  Sig.k1 = Sig[1:(M-n-1)] 
  k = c(0:(M-n-2))
  
  # to compute the variance 
  c1 = (phi^(2*k))
  c2 = c1*(theta^2)
  c3 = (phi^(2*k+1))*theta
  c4 = (phi^(2*k-1))*theta
  
  term1 = rev(c2) * Sig.k1
  term2 = rev(c1) * Sig.k
  term3 = rev(c3) * Sig.k1
  term4 = rev(c4) * Sig.k
  # variance of point t0-1
  e = term4[length(term4)]
  et = sum(term1+ term2+ term3 + term4) - e
  
  var.t = rep(NA, n)
  for (i in c(1:n)) {
    s = Sig[(M-n+i)]
    s1 = Sig[(M-n+i-1)]
    err = (theta/phi)*s1
    c5 = (s1)*(theta**2) + s + theta*phi*(s1)
    ei = c5 + (phi**2)*(et + s1*theta/phi)
    var.t[i] = ei
    et = ei
  }
  
  # to commpute the covariance 
  var.cov = matrix(NA, ncol = n, nrow = n)
  for (l in c(1:(n-1))) {
    var.et = var.t[l]
    sigma.et = Sig[(M-n+l)]
    for (l1 in c((l+1):(n))) {
      var.cov[l, l1] = (phi^(l1-l))*var.et + (phi^(l1-l-1))*theta*sigma.et
    }
  }
  var.cov[is.na(var.cov)] = 0
  var.cov = var.cov + t(var.cov)
  diag(var.cov) = var.t
  
  return(var.cov)
}


gls.true <- function(var.t, phi, theta, design.matrix, trend){
  if(phi ==0 & theta ==0){
    cov.var = diag(var.t)
  }else{
    Sig = rep(var.t,9)
    cov.var = var_cova(Sig = Sig, phi = phi, theta = theta, n = nrow(design.matrix))
  }
  X = data.frame(intercept = rep(1,n))
  if(trend ==0){
    X = cbind(X, design.matrix[c("jump")])
  }else{
    X = cbind(X,design.matrix[c("jump", "Xt")])
  }
  X = as.matrix(X)
  term1 = t(X) %*% (solve(cov.var)) %*% X
  beta = solve(term1) %*% t(X) %*% (solve(cov.var)) %*% (as.matrix(design.matrix$signal))
  var.beta = solve(term1)
  # form the frame of result as in the ols 
  t.val = beta/sqrt((diag(var.beta)))
  p.val = round(pnorm(-abs(t.val), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 4)
  fit.gls = data.frame(Estimate = beta)
  fit.gls$`Std. Error` = NA
  fit.gls$`t value` = t.val
  fit.gls$`Pr(>|t|)` = p.val
  
  return(list(Coefficients = beta, fit.gls = fit.gls, vcov = var.beta))
}

# selection in the cos+sin 
Test_OLS_vcovhac1 <- function(Data.mod){
  # OLS estimates
  list.para <- colnames(Data.mod)[2:dim(Data.mod)[2]]
  mod.X <-  list.para %>% stringr::str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X) %>% stringr::str_c(collapse = "")
  approx1 = c("AR(1)")
  fit.ols <- lm(mod.expression,data=Data.mod)
  vcov.para<-tryCatch(
    {
      sandwich::kernHAC(fit.ols, prewhite = 1,approx = approx1, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
    }, 
    error = function(e) {
      sandwich::kernHAC(fit.ols , prewhite = 0,approx = approx1, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
    }
  )
  # # options(show.error.messages = TRUE)
  # check.er = try(sandwich::kernHAC(fit.ols, prewhite = 1,approx = c("AR(1)"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE),silent = TRUE)
  # if(class(check.er) == "try-error"){
  #   vcov.para = sandwich::kernHAC(fit.ols, prewhite = 0,approx = c("AR(1)"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
  # }else{
  #   vcov.para=sandwich::kernHAC(fit.ols, prewhite = 1,approx = c("AR(1)"), kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
  # }
  # vcov.para = select1(fit.ols, y1 = 1, y2 = 0)
  # Test with vcov (HAC)
  fit.hac=lmtest::coeftest(fit.ols,df=fit.ols$df.residual,vcov.=vcov.para)[, ] %>% as.data.frame()
  keep.ind1 = which(fit.hac$`Pr(>|t|)` < 0.05)
  keep.ind = keep.ind1[keep.ind1>2]
  list.para.r = c( "Xt", list.para[keep.ind-1])
  mod.X.r <-  list.para.r %>% stringr::str_c(collapse = "+")
  mod.expression.r <- c("signal","~",mod.X.r) %>% stringr::str_c(collapse = "")
  # test with significant variables 
  fit.ols.r <- lm(mod.expression.r,data=Data.mod)
  
  vcov.para.r<-tryCatch(
    {
      sandwich::kernHAC(fit.ols.r, prewhite = 1,approx = approx1, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
    }, 
    error = function(e) {
      sandwich::kernHAC(fit.ols.r , prewhite = 0,approx = approx1, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
    }
  )
  # Test with vcov (HAC)
  fit.hac.r=lmtest::coeftest(fit.ols.r,df=fit.ols.r$df.residual,vcov.=vcov.para.r)[, ] %>% as.data.frame()

  return(list(fit.hac = fit.hac.r, fit.ols = fit.ols.r, vcov.para = vcov.para.r, predicted = fit.ols.r$fitted.values))
}

construct.design <- function(data.df, name.series){
  Data.mod <- data.df %>% dplyr::select(name.series,date) %>%
    rename(signal=name.series) %>% 
    mutate(complete.time=1:nrow(data.df)) %>% 
    dplyr::select(-date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  return(Data.mod)
}

IGLS <- function(design.m, tol, day.list){
  resi0 = rep(NA, nrow(design.m))
  # call expression
  list.para <- colnames(design.m)[2:dim(design.m)[2]]
  mod.X <-  list.para %>% stringr::str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X) %>% stringr::str_c(collapse = "")
  # ols
  ols.fit = lm( mod.expression, data = design.m)
  resi0[which(is.na(design.m$signal)==FALSE)] <- ols.fit$residuals
  old.coef = ols.fit$coefficients
  # estimate initial moving variance 
  Y0 = data.frame(date = day.list, residus = resi0)
  w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
  change1 = 10
  i=0
  while (change1 > tol) {
    design.m$w = w0^2
    gls.fit = eval(parse(text=paste0("gls(",mod.expression,",data=design.m, correlation = NULL, na.action = na.omit, weights=varFixed(value = ~w)",")")))
    change1 = sum((gls.fit$coefficients - old.coef)^2)
    deg =  cbind(rep(1, nrow(design.m)), as.matrix(design.m[,c(2:9)])) 
    fit.val = deg %*% as.matrix(gls.fit$coefficients)
    resi0 = design.m$signal - fit.val
    Y0 = data.frame(date = day.list, residus = resi0)
    w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
    old.coef = gls.fit$coefficients
    i=1+i
  }
  print(i)
  return(list( coefficients = gls.fit$coefficients, var = w0^2, residual = resi0, fit = fit.val))
}

