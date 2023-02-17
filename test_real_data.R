# this function used for the test of the real data
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution0/data.all_", win.thres, "years_", nearby_ver,"screened.RData")))
source(paste0(path_code_att,"FGLS.R"))
order.arma.l = get(load(file = paste0(path_results,"attribution0/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution0/coef.model.arma", win.thres,".RData")))
# reduced list of cases 
full.list = get(load( file = paste0(path_results, "attribution0/list.segments.selected", win.thres,".RData")))
full.list$station = paste0(full.list$main,".",as.character(full.list$brp), ".", full.list$nearby)
full.list$nbc = sapply(c(1:nrow(full.list)), function(x) min(full.list[x,c(4:5)]))
full.list$nbc.max = sapply(c(1:nrow(full.list)), function(x) max(full.list[x,c(4:5)]))
ind.sel = which(full.list$nearby!="pama" & full.list$min.var>0.002 & full.list$nbc >200 )
ind.sel = ind.sel[ind.sel!=8]
reduced.list = full.list[ind.sel,]
rownames(reduced.list) = NULL

order.arma.l1 = lapply(list.test, function(x) {
  a = as.data.frame(order.arma.l[[x]][[1]])
  return(a[ind.sel,])})
# run the FGLS 
all.res = list()
# a1 = list.files(path = paste0(path_results, "attribution0/FGLS-GE/"))
# list.old = as.character(substr(a1, 1, 20))
# a2 =  as.character(substr(reduced.list$station, 1, 20))
# list.select = which(a2 %in% list.old == FALSE)

list.ind = c(1:494)
ind.sel1 = list.ind[which(list.ind %in% ind.sel == FALSE)]
for (i in c(ind.sel1)) {
  df = dat[[reduced.list$station[i]]]
  fit.i = list()
  name.i = reduced.list$station[i]
  list.6 = c(1:6)
  for (j in list.6) {
    name.series = list.test[j]
    df.test = remove_na_2sides(df, name.series)
    
    if(reduced.list$nbc1[i] > 1000){
      start.day = df$date[3650] - 1000
      df.test = df.test[which(df.test$date>start.day),]
    }
    
    if(reduced.list$nbc2[i] > 1000){
      end.day = df$date[3650] + 1000
      df.test = df.test[which(df.test$date< end.day),]
    }
    ind.brp = which(df.test$date == df$date[3650])
    noise.model = unlist(as.data.frame(order.arma.l1[[j]])[i,])
    names(noise.model) = NULL
    cor.st = cor.struct.hac(noise.model)
    
    # ols
    Data.mod = construct.design(data.df = df.test, name.series = name.series, break.ind = ind.brp, one.year = 365)
    # list.para <- colnames(Data.mod)[2:(dim(Data.mod)[2]-1)]
    # mod.X <-  list.para %>% stringr::str_c(collapse = "+")
    # mod.expression <- c("signal","~",mod.X) %>% stringr::str_c(collapse = "")
    # ols.fit = lm( mod.expression, data = Data.mod)
    # 
    # # HAC
    # vcov.para=sandwich::kernHAC(ols.fit, approx = c( cor.st ), prewhite = TRUE, kernel = "Quadratic Spectral", adjust = TRUE, sandwich = TRUE)
    # fit.hac=lmtest::coeftest(ols.fit,df=(ols.fit$df.residual),vcov.=vcov.para)[, ] %>% as.data.frame()
    fit.fgls = FGLS1(design.m = Data.mod, tol=0.01, day.list = df.test$date, noise.model = noise.model, length.wind0 = 60)
    # fit.fgls1 = FGLS2(design.m = Data.mod, tol=0.0001, day.list = df.test$date, noise.model = noise.model)
    print(noise.model)
    fit.i[[name.series]] = fit.fgls
  }
  save(fit.i, file = paste0(path_results,"attribution0/FGLS/", name.i, "fgls.RData"))
  all.res[[i]] = fit.i
  
  print(i)
}

save(all.res, file = paste0(path_results,"attribution/all_fgls.RData"))
wls_model <- lm(mod.expression , data = Data.mod, weights=var.t)

# run for GPS-ERA
a = aggregate(nbc~main+brp, reduced.list, which.max)
colnames(a)[3] = "ind"
d = left_join(reduced.list, a, by = c("main", "brp"))
d$GE = NA
for (i in c(1:nrow(a))) {
  ind = which(d$main == a$main[i] & d$brp == a$brp[i])
  
  d$GE[ind[a$ind[i]]] = 1
}

reduced.list$GE = d$GE
list.selected = which(is.na(reduced.list$GE)==FALSE)
all.res = list()

list.selected = c(13,15,281,282,283,290,516)
for (i in list.selected) {
  df = dat[[reduced.list$station[i]]]
  fit.i = list()
  name.i = reduced.list$station[i]
  name.series = list.test[1]
  df.test = remove_na_2sides(df, name.series)
  
  if(reduced.list$nbc1[i] > 1000){
    start.day = df$date[3650] - 1000
    df.test = df.test[which(df.test$date>start.day),]
  }
  
  if(reduced.list$nbc2[i] > 1000){
    end.day = df$date[3650] + 1000
    df.test = df.test[which(df.test$date< end.day),]
  }
  ind.brp = which(df.test$date == df$date[3650])
  noise.model = unlist(as.data.frame(order.arma.l1[[1]])[i,])
  names(noise.model) = NULL
  cor.st = cor.struct.hac(noise.model)
  
  # ols
  Data.mod = construct.design(data.df = df.test, name.series = name.series, break.ind = ind.brp)
  fit.fgls = FGLS1(design.m = Data.mod, tol=0.01, day.list = df.test$date, noise.model = noise.model)
  print(noise.model)
  fit.i[[name.series]] = fit.fgls
  save(fit.i, file = paste0(path_results,"attribution/FGLS-GE/", name.i, "fgls.RData"))
  all.res[[i]] = fit.i
  
  print(i)
}
save(all.res, file = paste0(path_results,"attribution/all_fgls.RData"))
