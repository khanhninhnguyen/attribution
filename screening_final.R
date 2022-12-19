# final procedure for screening 
source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))

window.thres = 10
data.in = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
data.all <- list()
list.outlier <- list()
sd.all <- list()
for (i in c(1:length(data.in))) {
  # read data
  case.name = names(data.in)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  data.i = as.data.frame(data.in[[i]])
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  begin = breakpoint - (window.thres*365-1)
  fin = breakpoint + window.thres*365
  data.i = tidyr::complete(data.i, date = seq(begin, fin, by = "day"))
  
  # screen only series has data 
  before1 = na.omit(data.i[which(data.i$date <= breakpoint),])
  after1 = na.omit(data.i[which(data.i$date > breakpoint),])
  condi = 0
  if(nrow(before1)>100){
    before = data.i[which(data.i$date <= breakpoint),]
    bef.all = list()
    bef.outlier = list()
    bef.sd = list()
    
    for (k in c(1:6)) {
      bef.scr <- screen.O(Y = before, name.var = list.test[k], method = 'sigma', global.mu = 0, iter = 1, estimator = "Sca", fix.thres = 0, loes = 0, loes.method = 0)
      bef.all[[list.test[k]]] <- bef.scr$data
      bef.outlier[[list.test[k]]] <- bef.scr$point.rm
      bef.sd[[list.test[k]]] <- bef.scr$sd.est[[length(bef.scr$sd.est)]]
    }
    bef.all[["date"]] = before$date
    bef.all = as.data.frame(bef.all)
    rownames(bef.all) <- NULL
    condi = condi+1
  }
  if(nrow(after1) > 100){
    after = data.i[which(data.i$date > breakpoint),]
    aft.all = list()
    aft.outlier = list()
    aft.sd = list()
    
    for (k in c(1:6)) {
      aft.scr <- screen.O(Y = after, name.var = list.test[k], method = 'sigma',  global.mu = 0, iter = 1, estimator = "Sca", fix.thres = 0, loes = 0, loes.method = 0)
      aft.all[[list.test[k]]] = aft.scr$data
      aft.outlier[[list.test[k]]] = aft.scr$point.rm
      aft.sd[[list.test[k]]] = aft.scr$sd.est[[length(aft.scr$sd.est)]]
    }
    aft.all[["date"]] = after$date
    aft.all = as.data.frame(aft.all)
    rownames(aft.all) <- NULL
    condi = condi+1
  }
  if(condi == 2){
    data.all[[case.name]] <- rbind(bef.all, aft.all)
    list.outlier[[case.name]] <- list(bef = bef.outlier, aft = aft.outlier)
    sd.all[[case.name]] <- list(bef = bef.sd, aft = aft.sd)
  }
}
save(data.all, file = paste0(path_results,"attribution/data.all_", window.thres,"years_", nearby_ver,"screened.RData"))
save(list.outlier, file = paste0(path_results,"attribution/list.outlier_",  window.thres,"years_", nearby_ver,"screened.RData"))
save(sd.all, file = paste0(path_results,"attribution/sd.all_",  window.thres,"years_", nearby_ver,"screened.RData"))
