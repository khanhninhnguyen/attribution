# distribution of real data after restricted homogeneity
source(paste0(path_code_att,"sliding_variance.R"))

# # Normalize data first:  ------------------------------------------------
window.thres = 2 
data.cr = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
name.var = list.test[1]
dist.mean <- data.frame(matrix(NA, nrow = 0, ncol = (length(seq(-30,30,0.05))-1)))

for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  data.i = data.cr[[i]]
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  before =  data.i[which(data.i$date <= breakpoint),]
  after =  data.i[which(data.i$date > breakpoint),]
  if(nrow(na.omit(before)) > 31){
    bef.norm = two.step.norm(Y = before, name.var)
    if(length(bef.norm)>30){
      hist.b = hist(bef.norm, breaks = seq(-30,30,0.05), plot = FALSE)
      dist.mean <- rbind(dist.mean, hist.b$counts)
    }
  }
  if(nrow(na.omit(after)) > 31){
    aft.norm = two.step.norm(Y = after, name.var)
    if(length(aft.norm) > 30){
      hist.a = hist(aft.norm, breaks = seq(-30,30,0.05), plot = FALSE)
      dist.mean <- rbind(dist.mean, hist.a$counts)
    }
  }
}
save(dist.mean, file = paste0(path_results,"attribution/dist.mean_2years_", nearby_ver,".RData"))




