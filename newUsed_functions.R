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
  vcov.para=sandwich::kernHAC(fit.ols,prewhite = FALSE,kernel = "Quadratic Spectral",approx = "ARMA(1,1)", sandwich = TRUE)
  
  # Test with vcov (HAC)
  fit.hac=lmtest::coeftest(fit.ols,df=Inf,vcov.=vcov.para)[, ] %>% as.data.frame()
  
  
  # Selection parameter by parameter (the intercept will never be removed)
  threshold <- 0.01
  fit.hac.without.NA <- fit.hac[!(rownames(fit.hac) %in% "(Intercept)"),]
  pval.hac.without.intercept=fit.hac.without.NA$`Pr(>|z|)`
  p.max <- max(pval.hac.without.intercept)
  
  while ((p.max>threshold) & (length(list.para)>1)){
    names.max <- rownames(fit.hac.without.NA)[pval.hac.without.intercept==p.max]
    list.para <-list.para[!(list.para %in% names.max)] 
    mod.X <-  list.para %>% str_c(collapse = "+")
    mod.expression <- c("signal","~",mod.X) %>% str_c(collapse = "")
    
    fit.ols <- lm(eval(mod.expression),data=Data.mod)
    vcov.para=sandwich::kernHAC(fit.ols,prewhite = FALSE,kernel = "Quadratic Spectral",approx = "ARMA(1,1)", sandwich = TRUE)
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
