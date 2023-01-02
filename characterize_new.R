# new characterization
# this prog is used to compare the FPR and TPR of the test of the change in mean between the OLS-HAC and the GLS
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"newUsed_functions.R"))
source(paste0(path_code_att,"sliding_variance.R"))

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
# length.seg = matrix(NA, ncol = , nrow = 0)
list.break[c("nbc1", "nbc2", "len1", "len2")] <- NA
for (i in c(1:length(list.main))) {
  list.s = list.break[which(list.break$ref == list.main[i]),]
  list.nb = split(list.s, list.s$nb)
  for (j in c(1:length(list.nb))) {
    list.ij = paste0(list.nb[[j]]$ref,".",as.character(list.nb[[j]]$brp), ".", list.nb[[j]]$nb)
    data.ij = dat[list.ij]
    length.all <- sapply(c(1:length(data.ij)), function(x){
      seg1 = data.ij[[x]][c(1:L),]
      seg2 = data.ij[[x]][-c(1:L),]
      y = c(nb.consecutive(list.day = seg1$date, x = seg1$gps.gps), nb.consecutive(list.day = seg2$date, x = seg2$gps.gps),
            length(na.omit(seg1$gps.gps)), length(na.omit(seg2$gps.gps)))
    })
    
    list.break[which(list.break$ref %in% list.nb[[j]]$ref & list.break$nb %in% list.nb[[j]]$nb), c("nbc1", "nbc2", "len1", "len2")] = (t(length.all))
  }
}
list.break$r1 = list.break$nbc1/list.break$len1
list.break$r2 = list.break$nbc2/list.break$len2
list.break <- list.break[which(list.break$nb!= "kaza"),]
# add distance
distances <- as.data.frame(get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData"))))
colnames(list.break)[c(1,3)] <- c("main", "nearby")
full.list = left_join(list.break, distances, by = c("main", "nearby"))

# selection of segment based on its ratio of number of consecutive pairs

r = sapply(c(1:length(list.main)), function(x){
  list.s1 = full.list[which(full.list$main == list.main[x]),]
  ind.1 = which(list.s1$r1>0.9 | list.s1$r2>0.9)
  list.s = list.s1[ind.1,]
  ind.seg = ifelse(c(max(list.s$nbc1) > max(list.s$nbc2),  max(list.s$nbc1) > max(list.s$nbc2)), c(which.max(list.s$nbc1),1), c(which.max(list.s$nbc2),2))
  y = rep(NA, nrow(list.s1))
  y[ind.1[ind.seg[1]]] = ind.seg[2]
  print(x)
  return(y)
})
full.list$chose = unlist(r)
save(full.list, file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData"))

dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres,"years_", nearby_ver,"screened.RData")))
full.list = get(load( file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData")))
reduced.list = na.omit(full.list)

# heteroskedasticity ------------------------------------------------------

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
range.var <- function(x, day.list, s){
  df = data.frame(date = day.list, x = x)
  df$y = format(df$date, "%Y")
  df[which(is.na(s)==TRUE),] = NA
  if(all(x==1)){
    NA
  }else{
    anu.min = min(df$x, na.rm = TRUE)
    anu.max = max(df$x, na.rm = TRUE)
    anu.max - anu.min
  }
} # replaced functionf for 1 year

range.all = list()
range.diff = list()
range.mean = list()
for (i in c(1:nrow(reduced.list))) {
  name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
  dat.i = all.dat[[name.i]]
  for (j in c(1:6)) {
    name.series = list.test[j]
    var.ij = dat.i[,paste0(name.series, 'var')]
    range.all[[name.series]][[name.i]] = range.var(x = var.ij , day.list = dat.i$date, s = dat.i$gps.gps)
    range.diff[[name.series]][[name.i]] = diff.range.var(x = var.ij , day.list = dat.i$date, s = dat.i$gps.gps)
    range.mean[[name.series]][[name.i]] = mean(var.ij, na.rm = TRUE)
  }
}

res.tol = list(range.all, range.mean)
save(res.tol, file = paste0(path_results, "attribution/range_mean_var", win.thres,".RData"))
range.all1 = as.data.frame(range.all)
range.diff1 = as.data.frame(range.diff)
range.diff1 = range.diff1/range.all1 
# range.diff1 = range.diff1[which(l>500),]
# Histogram 
# apply(range.all1, 2, median)
# apply(range.diff1, 2, median)
# d = reshape2::melt(range.all1)
# d = reshape2::melt(range.diff1)
# 
# d %>%
#   ggplot( aes(x=value, color=variable, fill=variable)) + theme_bw()+
#   geom_histogram(alpha=0.6, binwidth = 0.25) +
#   ylab("") +
#   xlab("Range of moving window variance ") +
#   facet_wrap(~variable)

# CDF 
a = rbind(reshape2::melt(range.mean), reshape2::melt(range.all))
a$feature = rep(c("mean", "annual range"), each = nrow(a)/2)
a$series = rep(rep(list.name.test, each = nrow(a)/12),2)

ggplot(a, aes(x = value, col = series ))+ theme_bw()+
  stat_ecdf(lwd = 0.5, aes(linetype=feature))+
  scale_x_continuous(breaks = seq(0, 10, 1), limits = c(0,10))+
  geom_hline(yintercept = 0.5, size = 0.3) +
  scale_color_manual(values = brewer.pal(n = 6, name = 'Dark2'))+
  labs(y = "CDF", x = "Moving window variance", linetype = "")+
  theme(axis.text = element_text(size = 16),legend.text=element_text(size=12),
        axis.title = element_text(size=16))

summary(range.all1)
summary(range.diff1)
# specific case
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
fit.arima.manual <- function(signal.test){
  white = Arima(signal.test, order = c(0, 0 ,0))
  ar1 =  Arima(signal.test, order = c(1, 0 ,0))
  ma1 =  Arima(signal.test, order = c(0, 0 ,1))
  arma1 = Arima(signal.test, order = c(1, 0 ,1))
  
  list.model.name = c("white", "ar1", "ma1", "arma1")
  list.model =  list(white, ar1, ma1, arma1)
  bic.list = c(white$bic, ar1$bic, ma1$bic, arma1$bic)
  ind.chosen = which.min(bic.list)
  model.chosen = list.model[[ind.chosen]]
  pq = arimaorder(model.chosen)
  # choose only significant order and refit
  test.signif = as.data.frame(coeftest(model.chosen)[,])
  if(ncol(test.signif)==1){
    model.s = "white"
  }else{
    test.signif = test.signif[-nrow(test.signif),]
    ind.s = which(test.signif$`Pr(>|z|)`<0.05)
    if(length(ind.s) < 1){
      model.s = "white"
    }else if(length(ind.s) == 1){
      model.s = rownames(test.signif)
    }else{model.s ="arma1"}
  }
  
  model.chosen = list.model[[which(list.model.name %in% model.s)]]
  # refit
  pq = arimaorder(model.chosen)
  coef.list = rep(0, 3)
  coef.list[which(pq>0)] <- model.chosen$coef[-length(model.chosen$coef)]
  test.signif = as.data.frame(coeftest(model.chosen)[,])
  p.list = rep(1, 3)
  p.list[which(pq>0)] <- test.signif$`Pr(>|z|)`[-length(model.chosen$coef)]
  
  return(list(pq = pq, coef = coef.list, p = p.list))
}

order.arma.l <- list()
coef.arma.l <- list()
for (testi in c(1:6)) {
  name.test = list.test[testi]
  order.arma = data.frame(matrix(NA, ncol = 3, nrow = length(all.dat)))
  coef.arma = data.frame(matrix(NA, ncol = 3, nrow = length(all.dat)))
  for (i in c(1:nrow(reduced.list))) {
    name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
    dat.i = all.dat[[name.i]]
    # arima.fit = fit.arima(dat.i[, paste0(name.test, "res")])
    arima.fit = fit.arima.manual(dat.i[, paste0(name.test, "res")])
    order.arma[i,] = arima.fit$pq
    coef.arma[i,] = arima.fit$coef
  }
  order.arma.l[[name.test]] <- list(order.arma)
  coef.arma.l[[name.test]] <- list(coef.arma)
}
save(order.arma.l, file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData"))
save(coef.arma.l, file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData"))

# for plot 
order.arma.l = get(load(file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData")))

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
res.plot = data.frame(series = rep(list.name.test, each = 9), mod = rep(list.model, 6), value = six.values*100/length.data)
res.plot$series = factor(res.plot$series, 
                         levels=list.name.test)
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

x=list( a = six.model$gps.era, b = six.model$gps.gps)
library(ggvenn)

#### Diagramme de Venn (left=before, right=after):

six.model$main = reduced.list$main
six.model.m = reshape2::melt(six.model, id = "main")

name.series = "gps.gps"
x <- list(
  white = six.model.m$main[six.model.m$value=="White"],
  AR1 = six.model.m$main[six.model.m$value=="AR(1)"],
  MA1 = six.model.m$main[six.model.m$value=="MA(1)"],
  ARMA11 = six.model.m$main[six.model.m$value=="ARMA(1,1)"]

)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6
)
# plot the length of data -------
a = na.omit(full.list)
b = sapply(c(1:nrow(a)), function(x) a[x,c(8,9)][a[x,12]])
d = sapply(c(1:nrow(a)), function(x) a[x,c(4,5)][a[x,12]])
df = data.frame(b=unlist(b), d=unlist(d))
scatterplot <- ggplot(df, aes(x = b, y =d)) + theme_bw()+
  geom_point(size = 3, alpha = 0.6) +
  xlab("No. consecutive pairs/Total length")+ylab("No. consecutive pairs")
guides(color = FALSE) +
  theme(plot.margin = margin())

marginal_distribution <- function(x, var) {
  ggplot(x, aes_string(x = var)) +
    geom_histogram(bins = 40, alpha = 0.8) +
    guides(fill = FALSE) +
    theme_void() +
    theme(plot.margin = margin())}

x_hist <- marginal_distribution(df, "b")
y_hist <- marginal_distribution(df, "d") +coord_flip()

library(cowplot)
aligned_x_hist <- align_plots(x_hist, scatterplot, align = "v")[[1]]
aligned_y_hist <- align_plots(y_hist, scatterplot, align = "h")[[1]]


plot_grid(
  aligned_x_hist
  , NULL
  , scatterplot
  , aligned_y_hist
  , ncol = 2
  , nrow = 2
  , rel_heights = c(0.2, 1)
  , rel_widths = c(1, 0.2)
)




# histogram of the mean 
a = sapply(c(1:length(all.dat)), function(x) mean(all.dat[[x]]$gps.era, na.rm = TRUE))




