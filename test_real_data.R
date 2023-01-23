# this function used for the test of the real data
win.thres = 10
dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres, "years_", nearby_ver,"screened.RData")))
source(paste0(path_code_att,"FGLS.R"))
order.arma.l = get(load(file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData")))
# reduced list of cases 
full.list = get(load( file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData")))
reduced.list = full.list[which(full.list$nearby!="pama"),]

# run the FGLS 
all.res = list()
for (casei in c(1:length(full.list))) {
  df = dat[[i]]
  noise.model = unlist(as.data.frame(order.arma.l$gps.era[[1]])[i,])
  names(noise.model) = NULL
  Data.mod = construct.design(data.df = df, name.series = "gps.era", break.ind = 3650)
  fit.fgls = FGLS1(design.m = Data.mod, tol=0.01, day.list = df$date, noise.model = noise.model)
  all.res[[i]] = fit.fgls
}

save(all.res, paste0(path_results,"attibution/all.RData"))
wls_model <- lm(mod.expression , data = Data.mod, weights=var.t)

