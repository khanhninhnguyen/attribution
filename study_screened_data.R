# THIS PROG USED TO STUDY THE RESULT OF SCREENING ON THE REAL DATA 
window.thres = 1
data.in = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
dat  = get(load(file = paste0(path_results,"attribution/data.all_1year_", nearby_ver,"screened.RData")))
outl =  get(load(file = paste0(path_results,"attribution/list.outlier_1year_", nearby_ver,"screened.RData")))

# check the screening on gps-era 
bef = sapply(c(1:length(outl)), function(x) outl[[x]]$bef$gps1.era1)
names(bef) <- names(dat)
aft = sapply(c(1:length(outl)), function(x) outl[[x]]$aft$gps1.era1)
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

which(bef.out > 8)
ind.sus = 94
ind = unique.ind[ind.sus]
case.name = names(dat)[ind]
y = data.in[[case.name]]
plot(y$gps1.era1, type = "l")
list.out.bef = unlist(outl[[ind ]]$bef$gps1.era1)
list.out.aft = unlist(outl[[ind ]]$aft$gps1.era1) +365
points(list.out.bef, y$gps1.era1[list.out.bef], col = "red")
points(list.out.aft, y$gps1.era1[list.out.aft], col = "red")

# plot an individual case to investigate 

plot_screening(case.name =  case.name, data.in = data.in, var.name = "gps1.era1", side = "aft")
plot_screening(case.name =  case.name, data.in = data.in, var.name = "gps1.era1", side = "bef")

# compare between the position of outliers between gps-era and gps'-era'
bef.code = sapply(c(1:length(outl)), function(x) unlist(outl[[x]]$bef$gps.era))
bef.ngl = sapply(c(1:length(outl)), function(x) unlist(outl[[x]]$bef$gps1.era1))
bef.out.ngl = sapply(c(1:length(bef.ngl)), function(x) length(bef.ngl[[x]]))
aft.ngl = sapply(c(1:length(outl)), function(x) unlist(outl[[x]]$aft$gps1.era1))
aft.out.ngl = sapply(c(1:length(aft.ngl)), function(x) length(aft.ngl[[x]]))

coin = sapply(c(1:830), function(x) length(which(bef.ngl[[x]] %in% bef.code[[x]])))
table(coin)
hist(coin, breaks = 10)
ind = which.max(bef.out.ngl)
case.name = names(dat)[ind]





