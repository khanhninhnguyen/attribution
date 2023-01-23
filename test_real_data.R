# this function used for the test of the real data
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres, "years_", nearby_ver,"screened.RData")))
source(paste0(path_code_att,"FGLS.R"))
order.arma.l = get(load(file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData")))
# reduced list of cases 
full.list = get(load( file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData")))
reduced.list = full.list[which(full.list$nearby!="pama"),]
full.list$station = paste0(full.list$main,".",as.character(full.list$brp), ".", full.list$nearby)

# run the FGLS 
all.res = list()
for (i in c(1:nrow(full.list))) {
  df = dat[[i]]
  fit.i = list()
  name.i = full.list$station[i]
  if(is.na(full.list$chose[i]) == TRUE){
    list.6 = c(2:6)
  }else{
    list.6 = c(1:6)
  }
  for (j in list.6) {
    name.series = list.test[j]
    df.wtna = remove_na_2sides(df, name.series)
    end.ind = which(df$date==df.wtna$date[nrow(df.wtna)])
    start.ind = 3650-(end.ind-3650)+1
    # print(start.ind)
    # print(end.ind)
    df.test = df[c(start.ind:end.ind),]
    noise.model = unlist(as.data.frame(order.arma.l[[name.series]][[1]])[i,])
    names(noise.model) = NULL
    cor.st = cor.struct.hac(noise.model)
    
    # ols
    Data.mod = construct.design(data.df = df.test, name.series = name.series, break.ind = end.ind-3650)
    # list.para <- colnames(Data.mod)[2:(dim(Data.mod)[2]-1)]
    # mod.X <-  list.para %>% stringr::str_c(collapse = "+")
    # mod.expression <- c("signal","~",mod.X) %>% stringr::str_c(collapse = "")
    # ols.fit = lm( mod.expression, data = Data.mod)
    # 
    # # HAC
    # vcov.para=sandwich::kernHAC(ols.fit, approx = c( cor.st ), prewhite = TRUE, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
    # fit.hac=lmtest::coeftest(ols.fit,df=(ols.fit$df.residual),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.fgls = FGLS1(design.m = Data.mod, tol=0.01, day.list = df.test$date, noise.model = noise.model)
    # fit.fgls1 = FGLS2(design.m = Data.mod, tol=0.0001, day.list = df.test$date, noise.model = noise.model)
    
    fit.i[[name.series]] = fit.fgls
  }
  save(fit.i, file = paste0(path_results,"attribution/FGLS/", name.i, "fgls.RData"))
  all.res[[i]] = fit.i
  
  print(i)
}

save(all.res, file = paste0(path_results,"attribution/all_hac.RData"))
wls_model <- lm(mod.expression , data = Data.mod, weights=var.t)

