# last version for the data characterization
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"newUsed_functions.R"))
source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_characterization.R"))

# param
win.thres = 10
one.year=365
L = one.year*win.thres

# data 
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres = 10,"years_", nearby_ver,"screened.RData")))
# list longest segment ---------------------
list.break = data.frame(ref = substr(names(dat), start = 1, stop = 4), 
                        brp = substr(names(dat), start = 6, stop = 15),
                        nb = substr(names(dat), start = 17, stop = 20))
list.break[] <- lapply(list.break, as.character)
list.break$brp = as.Date(list.break$brp , format = "%Y-%m-%d")
list.main = unique((list.break$ref))
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
            nrow(remove_na_2sides(df = seg1, name.series = "gps.gps")), nrow(remove_na_2sides(df = seg2, name.series = "gps.gps")))
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
  ind.1 = which(list.s1$r1>0.5 | list.s1$r2>0.5)
  list.s = list.s1[ind.1,]
  ind.seg = ifelse(c(max(list.s$nbc1) > max(list.s$nbc2),  max(list.s$nbc1) > max(list.s$nbc2)), c(which.max(list.s$nbc1),1), c(which.max(list.s$nbc2),2))
  y = rep(NA, nrow(list.s1))
  y[ind.1[ind.seg[1]]] = ind.seg[2]
  print(x)
  return(y)
})
full.list$chose = unlist(r)
save(full.list, file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData"))
# if limit 1 year, we limit data from 10 year ------------------
full.list = get(load( file = paste0(path_results, "attribution/list.segments.selected", win.thres = 10,".RData")))
reduced.list = na.omit(full.list)
reduced.list$l = sapply(c(1:55), function(x) reduced.list[x, c(4,5)][reduced.list$chose[x]])
reduced.list$r = sapply(c(1:55), function(x) reduced.list[x, c(8,9)][reduced.list$chose[x]])
rownames(reduced.list) = NULL

# compute range and mean of variance from regression IFGLS to see the heteroskedasticity------------------
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

# Match data NO NEED BECAUSE MISSING DATA IN THIS STEP IS DUE TO THE SCREENING -------------------------------------------------------------
# all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))
# for (i in c(1:length(all.dat))) {
#   a = all.dat[[1]]
#   ind.rm = which(is.na(a$gps.era) == TRUE | is.na(a$gps1.era1) == TRUE)
#   a[ind.rm,-c("date")] = NA
# }
# plot the range and mean ---------
all.coef = get(load( file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData")))
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))
range.all = list()
range.mean = list()
for (i in c(1:nrow(reduced.list))) {
  name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
  dat.i = all.dat[[name.i]]
  for (j in c(1:6)) {
    name.series = list.test[j]
    var.ij = dat.i[,paste0(name.series, 'var')]
    range.all[[name.series]][[name.i]] = range.var(x = var.ij , day.list = dat.i$date, s = dat.i$gps.gps)
    range.mean[[name.series]][[name.i]] = mean(var.ij, na.rm = TRUE)
  }
}
a = rbind(reshape2::melt(range.mean), reshape2::melt(range.all))
a$feature = rep(c("mean", "annual range"), each = nrow(a)/2)
a$series = rep(rep(list.name.test, each = nrow(a)/12),2)
a$series = factor(a$series,  levels = reoder.list.name)

jpeg(paste0(path_results,"attribution/heteroskedasticity.jpg" ),width = 2500, height = 1800,res = 300)
library(RColorBrewer)
p <- ggplot(a, aes(x = value, col = series ))+ theme_bw()+
  stat_ecdf(lwd = 0.5, aes(linetype=feature))+
  scale_x_continuous(breaks = seq(0, 10, 1), limits = c(0,10))+
  geom_hline(yintercept = 0.5, size = 0.3) +
  scale_color_manual(values = brewer.pal(n = 6, name = 'Dark2'))+
  labs(y = "CDF", x = "Moving window variance", linetype = "")+
  theme(axis.text = element_text(size = 16),legend.text=element_text(size=12),
        axis.title = element_text(size=16))
print(p)
dev.off()
# Plot specific case to illustrate ------- CONTINUE --------
# construct data 
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))
dat.plot = all.dat$`albh.2015-12-29.sc02`
library(gtable)    
library(grid)
library(gridExtra) 
colors <- c("GPS-ERA" = "gray", "Fitted" = "black")

p1 <- ggplot(data = dat.plot, aes( x = date))+theme_bw()+
  geom_line(aes(y = gps.era, color = "GPS-ERA"))+
  geom_line(aes(y = gps.erafit, color = "Fitted"))+ ylab("GPS-ERA")+
  scale_color_manual(values = colors)+
  theme(axis.text = element_text(size = 14),legend.text=element_text(size=10),
        axis.title = element_text(size=14))+
  theme(
  legend.title=element_blank(),
  legend.background = element_rect(fill=alpha('white', 0.01)),
  legend.position = c(0.1, 0.12)
)
p2 <- ggplot(data = dat.plot, aes( x = date))+theme_bw()+
  geom_line(aes(y = gps.erares), col = "coral")+ ylab("Residual")+ 
  theme(axis.text = element_text(size = 14),legend.text=element_text(size=12),
                                                                         axis.title = element_text(size=14))
p3 <- ggplot(data = dat.plot, aes( x = date))+theme_bw()+
  geom_line(aes(y = gps.eravar), col = "black")+ylab("Moving window variance")+
  theme(axis.text = element_text(size = 14),legend.text=element_text(size=12),
        axis.title = element_text(size=14))
gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
gC <- ggplotGrob(p3)

gA$widths <- gC$widths
gB$widths <- gC$widths

jpeg(paste0(path_results,"figure/", 1,".jpeg"),width = 3000, height = 3000,res = 300) # change
print(grid.arrange(gA, gB, gC, nrow = 3))
dev.off() 
# estimate the arima ------------------------------------------------------
all.coef = get(load( file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData")))
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))

order.arma.l <- list()
coef.arma.l <- list()
for (testi in c(1:6)) {
  print(testi)
  name.test = list.test[testi]
  order.arma = data.frame(matrix(NA, ncol = 3, nrow = length(all.dat)))
  coef.arma = data.frame(matrix(NA, ncol = 4, nrow = length(all.dat)))
  for (i in c(1:nrow(reduced.list))) {
    print(i)
    name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
    dat.i = all.dat[[name.i]]
    y = dat.i[, paste0(name.test, "res")]/sqrt(dat.i[, paste0(name.test, "var")]) 
    arima.fit = fit.arima(dat.i[, paste0(name.test, "res")]/sqrt(dat.i[, paste0(name.test, "var")]) )
    # arima.fit = fit.arima.manual(dat.i[, paste0(name.test, "res")])
    order.arma[i,] = arima.fit$pq
    coef.arma[i,] = arima.fit$coef
  }
  order.arma.l[[name.test]] <- list(order.arma)
  coef.arma.l[[name.test]] <- list(coef.arma)
}
save(order.arma.l, file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData"))
save(coef.arma.l, file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData"))


# plot the histogram of noise model -----------------
order.arma.l = get(load(file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData")))

list.model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)", "AR(2)", "MA(2)", "ARMA(1,2)", "ARMA(2,1)", "ARMA(2,2)")
length.data =nrow(reduced.list)
six.model = data.frame(matrix(NA, ncol = 6, nrow = length.data))
for (i in 1:length(list.test)) {
  name.test = list.test[i]
  six.model[,i] = sapply(c(1:length.data), function(x) model.iden(as.numeric(unlist(order.arma.l[[name.test]][[1]][x,]))))
}
colnames(six.model) <- list.test
save(six.model, file = paste0(path_results,"attribution/six.models", win.thres,".RData"))
six.values = c()
for (i in 1:length(list.test)) {
  value.count = sapply(c(list.model), function(x) length(which(six.model[,i] == x)))
  six.values <- c( six.values, value.count)
}
res.plot = data.frame(series = rep(list.name.test, each = 9), mod = rep(list.model, 6), value = six.values)
res.plot$series = factor(res.plot$series, 
                         levels=reoder.list.name)
res.plot$mod = factor(res.plot$mod, 
                      levels=list.model)
res.plot = res.plot[which(res.plot$value != 0),]
jpeg(paste0(path_results,"attribution/iden_model_longest.jpg" ),width = 3000, height = 1800,res = 300)
p <- ggplot(res.plot, aes(fill=mod, y=value, x=series)) + 
  geom_bar(position="dodge", stat="identity", width = 0.5)+theme_bw()+ 
  xlab("") + ylab("Count")+
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
# Plot coefficients ------------------------
order.arma.l = get(load(file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData")))

n=length(coef.arma.l$gps.era[[1]][,1])

param.list <- c()
model.list <- c()
test.list  <- c()
values <- c()
list.param = c("phi", "theta")
list.name.station = c()
for (testi in c(1:6)) {
  name.test = list.test[testi]
  for (i in c(1:n)) {
    if(is.na(six.model[i, testi]) == FALSE){
      ind.par = unlist(extract_param(six.model[i, testi]))
      param = coef.arma.l[[name.test]][[1]][i,ind.par]
      param.list <- c(param.list, list.param[which(ind.par!=0)])
      values <- c(values, param)
      model.list <- c(model.list, rep(six.model[i, testi], length(param)))
      test.list  <- c(test.list, rep(list.name.test[testi], length(param)))
      list.name.station = c(list.name.station, rep(names(all.dat)[i], length(param)))
    }
  }
}

dat.p = data.frame(name = test.list, param = param.list, model = model.list, value = unlist(values), station = list.name.station)
dat.p$name = factor(dat.p$name,  levels = reoder.list.name)
ggplot(data = dat.p, aes( x = name, y = value, fill = model ,col = param)) + theme_bw()+
  geom_boxplot()+
  xlab("") + ylab(" values of parameters ") +
  theme(axis.text = element_text(size = 16),legend.text=element_text(size=12),
        axis.title = element_text(size=16))+
  scale_color_manual(values=c("red", "green", "purple"))





# estimate using ARMA(1,1) in general -------------------------------------

all.coef = get(load( file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData")))
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))

coef.arma.l <- list()
for (testi in c(1:6)) {
  name.test = list.test[testi]
  coef.arma = data.frame(matrix(NA, ncol = 4, nrow = length(all.dat)))
  for (i in c(1:nrow(reduced.list))) {
    # remove the case of white noise bc it will cause an error
    if(testi ==4 & i ==49){
      coef.arma[i,] = c(0,0,1,1)
    }else{
      name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
      dat.i = all.dat[[name.i]]
      y = dat.i[, paste0(name.test, "res")]/sqrt(dat.i[, paste0(name.test, "var")]) 
      arma11 = Arima( y, order = c(1,0,1), include.mean = FALSE)
      test.signif = as.data.frame(coeftest( arma11 )[,])
      coef.arma[i,] = c(arma11$coef[c(1,2)], test.signif$`Pr(>|z|)`[1:2])
    }
  }
  coef.arma.l[[name.test]] <- list(coef.arma)
}
save(coef.arma.l, file = paste0(path_results,"attribution/coef.model.arma", win.thres,"arma.RData"))
all.coef = get(load( file = paste0(path_results,"attribution/coef.model.arma", win.thres,"arma.RData")))
six.model = get(load(file = paste0(path_results,"attribution/six.models", win.thres,".RData")))

param.list <- c()
model.list <- c()
test.list  <- c()
values <- c()
p.val <- c()
for (i in c(1:6)) {
  name.test = list.test[i]
  model.list = c(model.list, rep(six.model[,i], each = 2))
  coef.df = all.coef[[name.test]][[1]]
  values = c(values, as.vector(rbind(unlist(coef.df$X1), unlist(coef.df$X2))))
  p.val = c(p.val, as.vector(rbind(unlist(coef.df$X3), unlist(coef.df$X4))))
  test.list = c( test.list, rep(name.test, 2*nrow(six.model)))
  param.list = c(param.list, rep(c( "phi", "theta"), nrow(six.model)))
}

dat.p = data.frame(name = test.list, param = param.list, model = model.list, value = unlist(values), p = p.val)
dat.p$name = rep(list.name.test, each = nrow(dat.p)/6)
dat.p$name = factor(dat.p$name,  levels = reoder.list.name)
ggplot(data = dat.p, aes( x = name, y = value, fill = model ,col = param)) + theme_bw()+
  geom_boxplot()+
  xlab("") + ylab(" values of parameters ") +
  theme(axis.text = element_text(size = 16),legend.text=element_text(size=12),
        axis.title = element_text(size=16))+
  scale_color_manual(values=c("red", "green", "purple"))

# impact of gaps on model characterization 
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres= 1,".RData")))
nb.c = sapply(c(1:55), function(x) nb.consecutive(list.day = all.dat[[x]]$date, x = all.dat[[x]]$gps.gpsres))
six.model$nb.c = nb.c
six.model$len = sapply(c(1:55), function(x) nrow(remove_na_2sides(df = all.dat[[x]], name.series = "gps.gpsres")))

d = data.frame(lim1 = unlist(six.model$r), lim10 = unlist(check1$r), l1 = six.model$nb.c, l10 = unlist(check1$L))
d$white1 = sapply(c(1:55), function(x) length(which(six.model[x,c(1:6)] == "White")))
colnames(d) = c("rate1", "rate10", "nb.cons1", "nb.cons10", "nb.white1", "nb.white10")
save(d, file = paste0(path_results, "attribution/length_white_relation.RData"))



# investigate the ARMA(1,1)-----------------------------
# read data case
six.model = get(load(file = paste0(path_results,"attribution/six.models", win.thres,".RData")))
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))
name.series = "era.era"
dat.plot = remove_na_2sides(all.dat$`tidb.2015-09-03.tid1`, name.series = name.series)
y = unlist(dat.plot[paste0(name.series, "res")]/sqrt(dat.plot[paste0(name.series, "var")]))
# plot data and it acf and pacf 
plot(unlist(dat.plot[name.series]), ylab = "raw", type = "l", col = "gray")
lines(unlist(dat.plot[paste0(name.series, "fit")]))
plot(y, ylab = "normalized residual", col = "coral")
plot(unlist(dat.plot[paste0(name.series, "var")]), ylab = "moving variance", col = "blue")
acf(y, na.action = na.exclude)
pacf(y, na.action = na.exclude)
# auto.arima results
fit.b = forecast::auto.arima(y , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                             max.p = 2, max.q = 2, start.p = 0, trace = TRUE, allowdrift = FALSE,  approximation=FALSE)
# fit model
ar1 = forecast::Arima(y , order = c(1,0,0), include.mean = FALSE)
arma1 = forecast::Arima(y , order = c(1,0,1), include.mean = FALSE)
ma1 = forecast::Arima(y , order = c(0,0,1), include.mean = FALSE)
# plot p value of boxtest
p.val = data.frame(matrix(NA, ncol = 2, nrow = 35))
for (l in c(1:35)) {
  p.val[l,] = c(Box.test(ar1$residuals, lag = l)$p.value, Box.test(ma1$residuals, lag = l)$p.value)
}
plot(p.val$X1, ylab = "p.value AR(1)")
plot(p.val$X2, ylab = "p.value ARMA(1,1)")
# ACF, PACF of the residual
acf(ar1$residuals, na.action = na.exclude)
pacf(ar1$residuals, na.action = na.exclude)
acf(arma1$residuals, na.action = na.exclude)
pacf(arma1$residuals, na.action = na.exclude)
# plot BIC
mod_capt <- capture.output(forecast::auto.arima(y , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =TRUE,lambda = NULL,
                                                max.p = 2, max.q = 2, start.p = 0, trace = TRUE, allowdrift = FALSE,  approximation=FALSE))

ind.c = c(2:19)
a = sapply(ind.c, function(x) as.numeric(strsplit(mod_capt[x], split = ":")[[1]][2]))
b = sapply(ind.c, function(x) strsplit(strsplit(mod_capt[x], split = ":")[[1]][1], split = " ")[[1]][2])
d = sapply(ind.c, function(x) strsplit(strsplit(mod_capt[x], split = ":")[[1]][1], split = " ")[[1]][4]) 
data.mod = data.frame(model =b, res =a, mean = d)
ggplot(data = data.mod, aes(x = model, y = res, col = mean))+ theme_bw()+
  geom_point()+ylab("BIC")+ theme(axis.text.y = element_text(size = 10))


# plot theoretical ACF and PACF
a = ARMAacf(ar = 0.24, lag.max = 35, pacf = TRUE)
df <- data.frame(lag = c(1:35), acf = a)
s = c(qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(y))), - qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(y))))
ggplot(data = df, mapping = aes(x = lag, y = acf)) + theme_bw()+
  geom_hline(aes(yintercept = 0)) + ylab("pacf")+
  geom_hline(yintercept = s, linetype = 2)+
  geom_segment(mapping = aes(xend = lag, yend = 0))

# plot periodogram
smooth.spec <- spec.pgram(x)
spec.df <- data.frame(freq = smooth.spec$freq, spec = smooth.spec$spec)
ggplot(data = spec.df) + theme_bw()+
  geom_line(aes(x = freq, y = spec)) + 
  scale_x_log10("Period (years)", 
                breaks = yrs.freqs, labels = yrs.labels) + scale_y_log10()


yrs.period <- rev(c(1/12, 1/6, 1/5, 1/4, 1/3, 0.5, 1, 3, 4))
yrs.labels <- rev(c( "1/12", "1/6", "1/5", "1/4", "1/3", "1/2", "1", "3", "4"))
yrs.freqs <- 1/yrs.period * 1/365



