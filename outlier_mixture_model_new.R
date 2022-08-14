# distribution of real data after restricted homogeneity

# # Normalize data first:  ------------------------------------------------
window.thres = 2 
data.cr = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
name.var = list.test[1]
for (i in c(1:length(data.cr))) {
  case.name = data.cr[[i]]
  std.t <- RobEstiSlidingVariance.S(Y = case.name[c("date", name.var)], name.var, alpha = 0)
  case.name$sd <- std.t
  slide.med <- sliding.median(Y = case.name[c("date", name.var)], name.var)
  # step 1: normalize by std and median
  norm1 <- (case.name[name.var] - slide.med)/std.t
  # step 2: avoid impact of constant in AR(1) by taking the IQR
  # ?should I adding NA on the IQR and median here?
}