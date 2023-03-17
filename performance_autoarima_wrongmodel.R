# run the simulation FGLS with wrong model
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"FGLS.R"))
nb.sim = 1000
length.wind = 60

# Start with true model --------------------------------------------------
## only wrong model: MA(1)
arma.coef = c(0,0.2)
model.est = c(0,0,1)
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
save(Tot.res.AR, file = paste0(path_results, "attribution0/truemodel.MA1.RData"))
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
arma.coef = c(0,0.2)
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
save(Tot.res.MA, file = paste0(path_results, "attribution0/wrongmodel.MA1.RData"))

