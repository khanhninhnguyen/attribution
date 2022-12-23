# new characterization
# this prog is used to compare the FPR and TPR of the test of the change in mean between the OLS-HAC and the GLS
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"newUsed_functions.R"))
source(paste0(path_code_att,"sliding_variance.R"))
win.thres = 10
one.year=365
L = one.year*win.thres

# choose the longest segment from the screened data ----------------------------


dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
name.series <- "gps.gps"
one.year=365
nb.consecutive <- function(list.day, x){
  a = list.day[which(is.na(x)== FALSE)]
  b = ts(a) 
  y = length(which(diff(b)==1))
  return(y)
}
list.break = data.frame(ref = substr(names(dat), start = 1, stop = 4), 
                        brp = substr(names(dat), start = 6, stop = 15),
                        nb = substr(names(dat), start = 17, stop = 20))
list.break[] <- lapply(list.break, as.character)
list.break$brp = as.Date(list.break$brp , format = "%Y-%m-%d")
list.main = unique((list.break$ref))
length.seg = matrix(NA, ncol = 2, nrow = 0)
list.break[c("len1", "len2")] <- NA
for (i in c(1:length(list.main))) {
  list.s = list.break[which(list.break$ref == list.main[i]),]
  list.nb = split(list.s, list.s$nb)
  for (j in c(1:length(list.nb))) {
    list.ij = paste0(list.nb[[j]]$ref,".",as.character(list.nb[[j]]$brp), ".", list.nb[[j]]$nb)
    data.ij = dat[list.ij]
    length.all <- sapply(c(1:length(data.ij)), function(x){
      seg1 = data.ij[[x]][c(1:L),]
      seg2 = data.ij[[x]][-c(1:L),]
      y = c(nb.consecutive(list.day = seg1$date, x = seg1$gps.gps), nb.consecutive(list.day = seg2$date, x = seg2$gps.gps))
    })
    
    list.break[which(list.break$ref %in% list.nb[[j]]$ref & list.break$nb %in% list.nb[[j]]$nb), c("len1", "len2")] = (t(length.all))
  }
}


res <- data.frame(seg = list.seg, side = list.side)

# add distance
distances <- as.data.frame(get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData"))))
colnames(list.break)[c(1,3)] <- c("main", "nearby")
full.list = left_join(list.break, distances, by = c("main", "nearby"))

r = sapply(c(1:length(list.main)), function(x){
  list.s = full.list[which(full.list$main == list.main[x]),]
  ind.seg = ifelse(c(max(list.s$len1) > max(list.s$len2),  max(list.s$len1) > max(list.s$len2)), c(which.max(list.s$len1),1), c(which.max(list.s$len2),2))
  y = rep(NA, nrow(list.s))
  y[ind.seg[1]] = ind.seg[2]
  print(x)
  return(y)
})
full.list$chose = unlist(r)
save(full.list, file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData"))

# heteroskedasticity ------------------------------------------------------

dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
full.list = get(load( file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData")))
reduced.list = na.omit(full.list)
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

remove_na_2sides <- function(df, name.series){
  a = which(is.na(df[name.series])== FALSE)
  df = df[c(min(a):(max(a))), ]
  return(df)
}

choose_segment <- function(x){
  if(x==1){
    y =c(1:L)
  }else{
    y = c((L+1):(2*L))
  }
}

# run the regression for the whole data
all.coef = list()
all.dat = list()
all.fit = list()
for (i in c(1:nrow(reduced.list))) {
  name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
  dat.i = dat[[name.i]]
  dat.i = dat.i[choose_segment(reduced.list$chose[i]),]
  dat.ij = remove_na_2sides(dat.i, name.series = "gps.gps")
  ind.all = which(is.na(dat.i[["gps.gps"]]) == FALSE)
  print(i)
  for (j in c(1:6)) {
    name.series0 = list.test[j]
    m = construct.design(dat.ij, name.series = name.series0)
    tol0 = 0.000000001
    if(i == 49 & j ==5){ tol0 = 0.0001 }
    fit.igls = IGLS(design.m = m, tol =  tol0, day.list = dat.ij$date)
    dat.i[c(min(ind.all): max(ind.all)), paste0(name.series0, 'var')] <- unlist(fit.igls$var)
    dat.i[c(min(ind.all): max(ind.all)), paste0(name.series0, 'res')] <- unlist(fit.igls$residual)
    dat.i[c(min(ind.all): max(ind.all)), paste0(name.series0, 'fit')] <- unlist(fit.igls$fit)
    all.coef[[name.series0]][[name.i]] = fit.igls$coefficients
    print(j)
  }
  all.dat[[name.i]] = dat.i
}
save(all.coef, file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData"))
save(all.dat, file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData"))

all.coef = get(load( file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData")))
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))

# l = sapply(c(1:length(all.dat)), function(x) length(all.dat[[x]]$gps.era))
length.consecutive = sapply(c(1:nrow(reduced.list)), function(x) max(reduced.list[x, c('len1', 'len2')]))
hist(length.consecutive , breaks = 100, main = "Histogram of the number of consecutive pairs", xlab = "")
range.var <- function(x, day.list, s){
  df = data.frame(date = day.list, x = x)
  df$y = format(df$date, "%Y")
  df[which(is.na(s)==TRUE),] = NA
  if(all(x==1)){
    NA
  }else{
    anu.min = aggregate(x~y, df, function(z) min(z, na.rm=TRUE))[,2]
    anu.max = aggregate(x~y, df, function(z) max(z, na.rm=TRUE))[,2]
    ifelse(length(anu.max)!=1,  mean( (anu.max-anu.min), na.rm = TRUE), NA)
  }
}
diff.range.var <- function(x, day.list,s){
  df = data.frame(date = day.list, x = x)
  df$y = format(df$date, "%Y")
  df[which(is.na(s)==TRUE),] = NA
  if(all(x==1)){
    NA
  }else{
    anu.min = aggregate(x~y, df, function(z) min(z, na.rm=TRUE))[,2]
    anu.max = aggregate(x~y, df, function(z) max(z, na.rm=TRUE))[,2]
    r = (anu.max- anu.min)
    ifelse(length(anu.max)!=1, (max(r, na.rm = TRUE) - min(r, na.rm = TRUE)), NA)
    
  }
}

range.all = list()
range.diff = list()
for (i in c(1:nrow(reduced.list))) {
  name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
  dat.i = all.dat[[name.i]]
  for (j in c(1:6)) {
    name.series = list.test[j]
    var.ij = dat.i[,paste0(name.series, 'var')]
    range.all[[name.series]][[name.i]] = range.var(x = var.ij , day.list = dat.i$date, s = dat.i$gps.gps)
    range.diff[[name.series]][[name.i]] = diff.range.var(x = var.ij , day.list = dat.i$date, s = dat.i$gps.gps)
  }
}

range.all1 = as.data.frame(range.all)
range.diff1 = as.data.frame(range.diff)
range.diff1 = range.diff1/range.all1 
# range.diff1 = range.diff1[which(l>500),]

apply(range.all1, 2, median)
apply(range.diff1, 2, median)
d = reshape2::melt(range.all1)
d = reshape2::melt(range.diff1)

d %>%
  ggplot( aes(x=value, color=variable, fill=variable)) + theme_bw()+
  geom_histogram(alpha=0.6, binwidth = 0.25) +
  ylab("") +
  xlab("Relative difference of the range of moving window variance ") +
  facet_wrap(~variable)

summary(range.all1)
summary(range.diff1)

a= remove_na_2sides(dat.i, name.series = "gps.era")
plot(a$date, a$gps.era1, xlab = "", ylab = "GPS-ERA'", type = 'l', col = "gray")
lines(a$date, a$gps.era1fit, col = "red", xlab = "", ylab = "MW variance of GPS-ERA'")
plot(a$date, a$gps.era1res, col = "red", xlab = "", type = 'l', ylab = "residual of GPS-ERA'")
plot(a$date, a$gps.era1res^2, col = "red", xlab = "", ylab = "MW variance of GPS-ERA'")
lines(a$date, a$gps.era1var, xlab = "", ylab = "MW variance of GPS-ERA'")


# plot an example 
# fit arma model on the residual from the IGLS ------------------------------
# used function  --- ------------
last_signif <- function(signal, pq, alpha, fit.b){  
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
# read residual and fit arima-------------------
full.list = get(load( file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData")))
reduced.list = na.omit(full.list)
all.coef = get(load( file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData")))
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))

fit.arima <- function(signal.test){
  fit.b = forecast::auto.arima(signal.test , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                               max.p = 2, max.q = 2, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
  
  pq <- arimaorder(fit.b)
  # order.init[k, c((testi*3-2): (testi*3))] <- pq
  options(warn = 2)
  
  refit0 = last_signif(signal = signal.test, pq, alpha = significant.level, fit.b = fit.b)
  pq = refit0$pq
  
  if( any(pq > 1)){
    fit.b = forecast::auto.arima( signal.test, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                  max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
    pq = arimaorder(fit.b)
  }
  
  refit1 = last_signif(signal = signal.test, pq, alpha = significant.level, fit.b = fit.b)
  
  pq = refit1$pq
  return(list(pq = pq, coef = refit1$pandcoef$coef, p = refit1$pandcoef$p.value))
}

order.arma.l <- list()
coef.arma.l <- list()
for (testi in c(1:6)) {
  name.test = list.test[testi]
  order.arma = data.frame(matrix(NA, ncol = 3, nrow = length(all.dat)))
  coef.arma = data.frame(matrix(NA, ncol = 4, nrow = length(all.dat)))
  for (i in c(1:nrow(reduced.list))) {
    name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
    dat.i = all.dat[[name.i]]
    arima.fit = fit.arima(dat.i[, paste0(name.test, "res")])
    order.arma[i,] = arima.fit$pq
    coef.arma[i,] = arima.fit$coef
  }
  order.arma.l[[name.test]] <- list(order.arma)
  coef.arma.l[[name.test]] <- list(coef.arma)
}
save(order.arma.l, file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData"))
save(coef.arma.l, file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData"))

# for plot 
list.model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)", "AR(2)", "MA(2)", "ARMA(1,2)", "ARMA(2,1)", "ARMA(2,2)")
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
length.data =nrow(reduced.list)
six.model = data.frame(matrix(NA, ncol = 6, nrow = length.data))
for (i in 1:length(list.test)) {
  name.test = list.test[i]
  six.model[,i] = sapply(c(1:length.data), function(x) model.iden(as.numeric(unlist(order.arma.l[[name.test]][[1]][x,]))))
}
colnames(six.model) <- list.test
six.values = c()
for (i in 1:length(list.test)) {
  value.count = sapply(c(list.model), function(x) length(which(six.model[,i] == x)))
  six.values <- c( six.values, value.count)
}
res.plot = data.frame(series = rep(list.test, each = 9), mod = rep(list.model, 6), value = six.values*100/length.data)
res.plot$series = factor(res.plot$series, 
                         levels=list.test)
res.plot$mod = factor(res.plot$mod, 
                      levels=list.model)

jpeg(paste0(path_results,"attribution/iden_model_longest.jpg" ),width = 3000, height = 1800,res = 300)
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





