# run the simulation FGLS with wrong model
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"FGLS.R"))
nb.sim = 1000
length.wind = 60

# Start with true model --------------------------------------------------
## only wrong model: MA(1)
arma.coef = c(0.6, -0.3)
model.est = c(1,0,1)
Tot.res.AR = list()
set.seed(1)
for (i in c(1:nb.sim)) {
  y = simulate.general1(N = N, arma.model = arma.coef, burn.in = 1000, hetero = 0, sigma = 1)
  df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:N])
  Data.mod = construct.design(data.df = df, name.series = "y", break.ind = N/2, one.year = N/2)
  Data.mod =  Data.mod[,c(1,10,11)]
  fgls.fit = FGLS1(design.m = Data.mod, tol=0.01, day.list = df$date, noise.model = model.est, length.wind)
  Tot.res.AR[[i]] =  fgls.fit
}
save(Tot.res.AR, file = paste0(path_results, "attribution0/truemodel.ARMA1.RData"))
## for different length
arma.coef = c(0.3, 0)
model.est = c(1,0,0)
N = 1000
Tot.res.AR = list()
Tot.res = list()
length.list = c(200, 600, 1000, 1400)
for (j in c(1:length(length.list))) {
  N = length.list[j]
  set.seed(1)
  for (i in c(1:nb.sim)) {
    y = simulate.general1(N = N, arma.model = arma.coef, burn.in = 1000, hetero = 0, sigma = 1)
    df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:N])
    Data.mod = construct.design(data.df = df, name.series = "y", break.ind = N/2, one.year = N/2)
    Data.mod =  Data.mod[,c(1,10,11)]
    fgls.fit = FGLS1(design.m = Data.mod, tol=0.01, day.list = df$date, noise.model = model.est, length.wind)
    Tot.res.AR[[i]] =  fgls.fit
  }
  Tot.res[[j]] = Tot.res.AR
}

save(Tot.res, file = paste0(path_results, "attribution0/truemodel.ARMA1.length.RData"))
# Wrong model -------------------------------------------------------------
## only wrong model: a) AR(phi = 0.3) but identified as MA(1)---
arma.coef = c(0.3,0)
model.est = c(1,0,0)
N = 1000
Tot.res.AR = list()
set.seed(1)
for (i in c(1:nb.sim)) {
  y = simulate.general1(N = N, arma.model = arma.coef, burn.in = 1000, hetero = 0, sigma = 1)
  df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:N])
  Data.mod = construct.design(data.df = df, name.series = "y", break.ind = N/2, one.year = N/2)
  Data.mod =  Data.mod[,c(1,10,11)]
  fgls.fit = FGLS1(design.m = Data.mod, tol=0.01, day.list = df$date, noise.model = model.est, length.wind)
  Tot.res.AR[[i]] =  fgls.fit
}
save(Tot.res.AR, file = paste0(path_results, "attribution0/truemodel.AR1.RData"))
## only wrong model: a) MA(phi = 0.2) but identified as AR(1)---
arma.coef = c(0,0.3)
model.est = c(1,0,0)
N = 1000
Tot.res.MA = list()
set.seed(1)
for (i in c(1:nb.sim)) {
  y = simulate.general1(N = N, arma.model = arma.coef, burn.in = 1000, hetero = 0, sigma = 1)
  df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:N])
  Data.mod = construct.design(data.df = df, name.series = "y", break.ind = N/2, one.year = N/2)
  Data.mod =  Data.mod[,c(1,10,11)]
  fgls.fit = FGLS1(design.m = Data.mod, tol=0.01, day.list = df$date, noise.model = model.est, length.wind)
  Tot.res.MA[[i]] =  fgls.fit
}
save(Tot.res.MA, file = paste0(path_results, "attribution0/wrongmodel.MA1a.RData"))

## only wrong model: a) MA(phi = 0.2) but identified as AR(1)---
arma.coef = c(0.6,-0.3)
model.est = c(1,0,0)
N = 1000
Tot.res.MA = list()
set.seed(1)
for (i in c(1:nb.sim)) {
  y = simulate.general1(N = N, arma.model = arma.coef, burn.in = 1000, hetero = 0, sigma = 1)
  df = data.frame(y = y, date = seq(as.Date("2014-01-13"), as.Date("2018-05-27"), by="days")[1:N])
  Data.mod = construct.design(data.df = df, name.series = "y", break.ind = N/2, one.year = N/2)
  Data.mod =  Data.mod[,c(1,10,11)]
  fgls.fit = FGLS1(design.m = Data.mod, tol=0.01, day.list = df$date, noise.model = model.est, length.wind)
  Tot.res.MA[[i]] =  fgls.fit
}
save(Tot.res.MA, file = paste0(path_results, "attribution0/wrongmodel.ARMA1.RData"))


# plot results ------------------------------------------------------------
true.ar = get(load(file = paste0(path_results, "attribution0/truemodel.MA1.RData")))
wrong.ar = get(load(file = paste0(path_results, "attribution0/wrongmodel.MA1.RData")))
true.model = "MA"
t.true.ar = sapply(c(1:nb.sim), function(x) true.ar[[x]]$t.table$`t value`[1] )
t.wrong.ar = sapply(c(1:nb.sim), function(x) wrong.ar[[x]]$t.table$`t value`[1] )

hist(t.true.ar)
hist(t.wrong.ar)

table(abs(t.true.ar)<1.96)
table(abs(t.wrong.ar)<1.96)
d = data.frame(true = t.true.ar, est = t.wrong.ar) %>% reshape2::melt()
ggplot(data = d, aes(x = variable, y = value))+theme_bw()+geom_boxplot()
p1 = ggplot(data = d, aes(x = value, col = variable))+
  theme_bw()+
  stat_ecdf()+
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3,3,1))+
  ylab("CDF")+
  xlab("t-value")
ggsave(filename = paste0(path_results, "attribution0/t_val_",true.model,".jpg"), plot = p1)

phi.true.ar = sapply(c(1:nb.sim), function(x) true.ar[[x]]$coef.arma$theta)
phi.wrong.ar = sapply(c(1:nb.sim), function(x) wrong.ar[[x]]$coef.arma$phi)

d = data.frame(true = phi.true.ar, est = phi.wrong.ar) %>% reshape2::melt()
p2 = ggplot(data = d, aes(x = variable, y = value))+
  theme_bw()+
  geom_boxplot()+
  ylab("Autocorrelation")
ggsave(filename = paste0(path_results, "attribution0/theta_",true.model,".jpg"), plot = p2)

var.true.ar = sapply(c(1:nb.sim), function(x) mean(true.ar[[x]]$t.table$`Std. Error`[1]))
var.wrong.ar = sapply(c(1:nb.sim), function(x) mean(wrong.ar[[x]]$t.table$`Std. Error`[1]))

d = data.frame(true = var.true.ar, est = var.wrong.ar) %>% reshape2::melt()
p3 = ggplot(data = d, aes(x = value, col = variable))+
  theme_bw()+
  stat_ecdf()+
  scale_x_continuous(limits = c(0.07, 0.13), breaks = seq(0.07, 0.13,0.01))+
  ylab("CDF")+
  xlab("Std. Err")
ggsave(filename = paste0(path_results, "attribution0/std_err_",true.model,".jpg"), plot = p3)

est.true.ar = sapply(c(1:nb.sim), function(x) mean(true.ar[[x]]$t.table$Estimate[1]))
est.wrong.ar = sapply(c(1:nb.sim), function(x) mean(wrong.ar[[x]]$t.table$Estimate[1]))

d = data.frame(true = est.true.ar, est = est.wrong.ar) %>% reshape2::melt()
p3 = ggplot(data = d, aes(x = value, col = variable))+
  theme_bw()+
  stat_ecdf()+
  ylab("CDF")+
  xlab("Std. Err")
ggsave(filename = paste0(path_results, "attribution0/beta_",true.model,".jpg"), plot = p3)


