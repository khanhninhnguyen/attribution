# THIS PROG USED TO STUDY THE RESULT OF SCREENING ON THE REAL DATA 
window.thres = 1
data.in = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
dat  = get(load(file = paste0(path_results,"attribution/data.all_1year_", nearby_ver,"screened.RData")))
outl =  get(load(file = paste0(path_results,"attribution/list.outlier_1year_", nearby_ver,"screened.RData")))

# check the screening on gps-era 
bef = sapply(c(1:length(outl)), function(x) outl[[x]]$bef$gps.era)
names(bef) <- names(dat)
aft = sapply(c(1:length(outl)), function(x) outl[[x]]$aft$gps.era)
names(aft) <- names(dat)

# extract only gps-era (same for different nearby)
all.cases.name = names(dat)
all.cases.ind = sapply(c(1:length(all.cases.name)), function(x) substr(all.cases.name[x],start = 1, stop = 15))
unique.ind = match(unique(all.cases.ind), all.cases.ind )

gps.era.dat = dat[unique.ind]
bef.outl = bef[unique.ind]
bef.out = sapply(c(1:length(bef.outl)), function(x) length(unlist(bef.outl[[x]])))
aft.outl = aft[unique.ind]
aft.out = sapply(c(1:length(aft.outl)), function(x) length(unlist(aft.outl[[x]])))

which.max(aft.out)
ind.sus =84
ind = unique.ind[ind.sus]
case.name = names(dat)[ind]
y = data.in[[case.name]]
plot(y$gps.era, type = "l")
list.out.bef = unlist(outl[[ind ]]$bef$gps.era)
list.out.aft = unlist(outl[[ind ]]$aft$gps.era) +365
points(list.out.bef, y$gps.era[list.out.bef], col = "red")
points(list.out.aft, y$gps.era[list.out.aft], col = "red")

# plot an individual case to investigate 

plot_screening(case.name =  case.name, data.in = data.in, var.name = "gps.era", side = "aft")
plot_screening(case.name =  case.name, data.in = data.in, var.name = "gps.era", side = "bef")



