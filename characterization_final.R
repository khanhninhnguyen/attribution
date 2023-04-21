# last version for the data characterization
source(paste0(path_code_att,"simulate_time_series.R"))
# source(paste0(path_code_att,"newUsed_functions.R"))
source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_characterization.R"))
# UPDATE OR REMOVE PAMA STATION 
# param
win.thres = 10
one.year=365
L = one.year*win.thres

# data 
dat = get(load( file = paste0(path_results,"attribution0/data.all_", win.thres = 10,"years_", nearby_ver,"screened.RData")))
# # list longest segment ---------------------
# list.break = data.frame(ref = substr(names(dat), start = 1, stop = 4),
#                         brp = substr(names(dat), start = 6, stop = 15),
#                         nb = substr(names(dat), start = 17, stop = 20))
# list.break[] <- lapply(list.break, as.character)
# list.break$brp = as.Date(list.break$brp , format = "%Y-%m-%d")
# list.main = unique((list.break$ref))
# list.break[c("nbc1", "nbc2", "len1", "len2", "min.var")] <- NA
# for (i in c(1:length(list.main))) {
#   list.s = list.break[which(list.break$ref == list.main[i]),]
#   list.nb = split(list.s, list.s$nb)
#   for (j in c(1:length(list.nb))) {
#     list.ij = paste0(list.nb[[j]]$ref,".",as.character(list.nb[[j]]$brp), ".", list.nb[[j]]$nb)
#     data.ij = dat[list.ij]
#     length.all <- sapply(c(1:length(data.ij)), function(x){
#       seg1 = data.ij[[x]][c(1:L),]
#       seg2 = data.ij[[x]][-c(1:L),]
#       # compute the minimum of variance too
#       min.var = min(sapply(c(1:6), function(k) sd(unlist(seg1$era.era),na.rm=TRUE)^2))
#       y = c(nb.consecutive(list.day = seg1$date, x = seg1$gps.gps), nb.consecutive(list.day = seg2$date, x = seg2$gps.gps),
#             nrow(remove_na_2sides(df = seg1, name.series = "gps.gps")), nrow(remove_na_2sides(df = seg2, name.series = "gps.gps")), min.var)
#     })
# 
#     list.break[which(list.break$ref %in% list.nb[[j]]$ref & list.break$nb %in% list.nb[[j]]$nb), c("nbc1", "nbc2", "len1", "len2", "min.var")] = (t(length.all))
#   }
# }
# list.break$r1 = list.break$nbc1/list.break$len1
# list.break$r2 = list.break$nbc2/list.break$len2
# list.break <- list.break[which(list.break$nb!= "kaza"),]
# # add distance
# distances <- as.data.frame(get(load(file = paste0(path_results, "attribution/", version_name, nearby_ver, "distances-pairs.RData"))))
# colnames(list.break)[c(1,3)] <- c("main", "nearby")
# full.list = left_join(list.break, distances, by = c("main", "nearby"))
# 
# # selection of segment based on its ratio of number of consecutive pairs
# 
# r = sapply(c(1:length(list.main)), function(x){
#   list.s1 = full.list[which(full.list$main == list.main[x]),]
#   list.s1$ind = c(1:nrow(list.s1))
#   # ind.1 = which(list.s1$r1>0.5 | list.s1$r2>0.5)
#   # ind.2 = which(list.s1$min.var>0.05)
#   ind.1 = which(list.s1$min.var>0.0)
#   list.s = list.s1[ind.1,]
#   list.s$nbc = list.s$nbc1+list.s$nbc2
#   # a = data.frame(nb =  c(list.s$nbc1, list.s$nbc2), rl = c(list.s$r1, list.s$r2), ind = c(list.s$ind, list.s$ind), ind.seg = c(rep(c(1,2), each = nrow(list.s))))
#   # ind2 = a[which(a$rl>0),]
#   # ind3 = ind2[which.max(ind2$nb),]
#   # ind3 = list.s[which.max(list.s$nbc),]
#   # increase nb of nearbys
#   # ind2 = sort(list.s$nbc, index.return=TRUE, decreasing = TRUE)
#   y = rep(NA, nrow(list.s1))
#   if(nrow(list.s)>0){
#     ind0 = aggregate(nbc~nearby, data = list.s, FUN = max) # choose longest segment for each nearby
#     ind0 = inner_join(list.s, ind0)
#     ind2 = sort(ind0$nbc, index.return=TRUE, decreasing = TRUE)
#     if(length(ind2$ix)>0){
#       ind3 = ind0[ind2$ix[1],]
#     }else {
#       ind3 = list.s
#     }
#     # y[ind3$ind] = ind3$ind.seg
#     y[ind3$ind] = ind3$ind
#     print(x)
#   }
#   return(y)
# })
# full.list$chose = unlist(r)
# save(full.list, file = paste0(path_results, "attribution/list.segments.selected", win.thres,".RData"))
# LIST OF STATION   ------------------

full.list = get(load( file = paste0(path_results, "attribution0/list.segments.selected", win.thres = 10,".RData")))
reduced.list = full.list
# reduced.list$chose[51] =1
# reduced.list$l = sapply(c(1:nrow(reduced.list)), function(x) reduced.list[x, c(4,5)][reduced.list$chose[x]])
# reduced.list$r = sapply(c(1:nrow(reduced.list)), function(x) reduced.list[x, c(9,10)][reduced.list$chose[x]])
reduced.list$nbc = reduced.list$nbc1+reduced.list$nbc2
reduced.list$station = paste0(reduced.list$main,".",as.character(reduced.list$brp), ".", reduced.list$nearby)
rownames(reduced.list) = NULL
reduced.list = reduced.list[which(reduced.list$nearby!="pama"),]
# # compute range and mean of variance from regression IFGLS to see the heteroskedasticity------------------
# all.coef = list()
# all.dat = list()
# all.fit = list()
# for (i in c(1:nrow(reduced.list))) {
#   name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
#   dat.i = dat[[name.i]]
#   # six.ser = na.omit(dat.i)
#   # dat.i = tidyr::complete(six.ser, date = seq(dat.i$date[1], dat.i$date[7300], by = "day"))
#   # # dat.i = dat.i[choose_segment(reduced.list$chose[i]),]
#   # dat.ij = remove_na_2sides(dat.i, name.series = "gps.gps")
#   print(i)
#   for (j in c(1:6)) {
#     name.series0 = list.test[j]
#     m = construct.design(dat.i, name.series = name.series0)
#     tol0 = 0.01
#     if(i == 49 & j ==5){ tol0 = 0.0001 }
#     fit.igls = IGLS(design.m = m, tol =  tol0, day.list = dat.i$date)
#     # dat.i[c(min(ind.all): max(ind.all)), paste0(name.series0, 'var')] <- unlist(fit.igls$var)
#     # dat.i[c(min(ind.all): max(ind.all)), paste0(name.series0, 'res')] <- unlist(fit.igls$residual)
#     # dat.i[c(min(ind.all): max(ind.all)), paste0(name.series0, 'fit')] <- unlist(fit.igls$fit)
#     dat.i[, paste0(name.series0, 'var')] <- unlist(fit.igls$var)
#     dat.i[, paste0(name.series0, 'res')] <- unlist(fit.igls$residual)
#     dat.i[, paste0(name.series0, 'fit')] <- unlist(fit.igls$fit)
#     all.coef[[name.series0]][[name.i]] = fit.igls$coefficients
#     print(j)
#   }
#   all.dat[[name.i]] = dat.i
# }
# save(all.coef, file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData"))
# save(all.dat, file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData"))

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
    sd.ij = sqrt(unlist(dat.i[,paste0(name.series, 'var')]))
    range.all[[name.series]][[name.i]] = range.var(x = sd.ij , day.list = dat.i$date, s = dat.i$gps.gps)
    range.mean[[name.series]][[name.i]] = mean(sd.ij, na.rm = TRUE)
  }
}
a = rbind(reshape2::melt(range.mean), reshape2::melt(range.all))
a$feature = rep(c("mean", "annual range"), each = nrow(a)/2)
a$series = rep(rep(list.name.test, each = nrow(a)/12),2)
a$series = factor(a$series,  levels = reoder.list.name)
b= list(mean = range.mean, range = range.all)
save(b, file = paste0(path_results, "moving.var.RData"))

# jpeg(paste0(path_results,"attribution/heteroskedasticity.jpg" ),width = 2500, height = 1800,res = 300)
library(RColorBrewer)
p <- ggplot(a, aes(x = value, col = series ))+ theme_bw()+
  stat_ecdf(lwd = 0.3, aes(linetype=feature))+
  scale_x_continuous(breaks = seq(0, 3, 1), limits = c(0,3))+
  geom_hline(yintercept = 0.5, size = 0.3) +
  scale_color_manual(values = brewer.pal(n = 6, name = 'Dark2'))+
  labs(y = "CDF", x = "Moving window standard deviation", linetype = "")+
  theme(axis.text = element_text(size = 5),legend.text=element_text(size=4.5),
        axis.title = element_text(size=5), legend.key.size = unit(0.15, "cm"),
        legend.box.spacing = unit(0, "pt"), legend.title= element_blank())
# print(p)
# dev.off()
ggsave(paste0(path_results,"attribution/heteroskedasticity.jpg" ), plot = p, width = 8.8, height = 6, units = "cm", dpi = 1200)

# Plot specific case to illustrate ------- CONTINUE --------
# construct data 
all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))
dat.plot = all.dat$`albh.2015-12-29.sc02`
library(gtable)    
library(grid)
library(gridExtra) 
colors <- c("GPS-ERA" = "gray", "Fitted" = "black", "Moving median" = "red")
dat.plot$medi.var = sliding.median(dat.plot, name.var = "gps.erares",length.wind = 60)
p1 <- ggplot(data = dat.plot, aes( x = date))+theme_bw()+
  geom_line(aes(y = gps.era), col = "gray", lwd = 0.3)+
  geom_line(aes(y = gps.erafit, color = "Fitted"), lwd = 0.3)+ ylab("GPS-ERA")+
  scale_color_manual(values = colors)+
  theme(axis.text = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size=5), legend.key.size = unit(0.2, "cm"))+
  theme(
  legend.title=element_blank(),
  legend.background = element_rect(fill=alpha('white', 0.01)),
  legend.position = c(0.1, 0.12)
)
p2 <- ggplot(data = dat.plot, aes( x = date))+theme_bw()+
  geom_line(aes(y = gps.erares), col = "coral", lwd = 0.3)+
  geom_line(aes(y = medi.var, color = "Moving median"), lwd = 0.3) + 
  scale_color_manual(values = colors)+
  ylab("Residual")+ 
  theme(axis.text = element_text(size = 5),legend.text=element_text(size=4),
  axis.title = element_text(size=5), legend.key.size = unit(0.2, "cm"), 
  legend.title=element_blank(),
  legend.background = element_rect(fill=alpha('white', 0.01)),
  legend.position = c(0.1, 0.12))
p3 <- ggplot(data = dat.plot, aes( x = date))+theme_bw()+
  geom_line(aes(y = sqrt(gps.eravar)), col = "black", lwd = 0.3)+ylab("Moving window standard deviation")+
  theme(axis.text = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size=5))
gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
gC <- ggplotGrob(p3)

gA$widths <- gC$widths
gB$widths <- gC$widths

p = grid.arrange(gA, gB, gC, nrow = 3)
ggsave(paste0(path_results,"attribution/specific_case.jpg" ), plot = p, width = 8.8, height = 10, units = "cm", dpi = 1200)

jpeg(paste0(path_results,"figure/", 1,".jpeg"),width = 3000, height = 3000,res = 300) # change
print(grid.arrange(gA, gB, gC, nrow = 3))
dev.off() 
# estimate the arima ------------------------------------------------------
# all.coef = get(load( file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData")))
# all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))
# 
# order.arma.l <- list()
# coef.arma.l <- list()
# for (testi in c(1:6)) {
#   print(testi)
#   name.test = list.test[testi]
#   order.arma = data.frame(matrix(NA, ncol = 3, nrow = length(all.dat)))
#   coef.arma = data.frame(matrix(NA, ncol = 4, nrow = length(all.dat)))
#   for (i in c(1:nrow(reduced.list))) {
#     print(i)
#     name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
#     dat.i = all.dat[[name.i]]
#     y = dat.i[, paste0(name.test, "res")]/sqrt(dat.i[, paste0(name.test, "var")]) 
#     arima.fit = fit.arima(dat.i[, paste0(name.test, "res")]/sqrt(dat.i[, paste0(name.test, "var")]) )
#     # arima.fit = fit.arima.manual(dat.i[, paste0(name.test, "res")])
#     order.arma[i,] = arima.fit$pq
#     coef.arma[i,] = arima.fit$coef
#   }
#   order.arma.l[[name.test]] <- list(order.arma)
#   coef.arma.l[[name.test]] <- list(coef.arma)
# }
# save(order.arma.l, file = paste0(path_results,"attribution/order.model.arma", win.thres,".RData"))
# save(coef.arma.l, file = paste0(path_results,"attribution/coef.model.arma", win.thres,".RData"))


# plot the histogram of noise model -----------------
order.arma.l = get(load(file = paste0(path_results,"attribution0/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution0/coef.model.arma", win.thres,".RData")))

list.model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)")
length.data =nrow(reduced.list)
six.model = data.frame(matrix(NA, ncol = 6, nrow = length.data))
for (i in 1:length(list.test)) {
  name.test = list.test[i]
  six.model[,i] = sapply(c(1:length.data), function(x) model.iden(as.numeric(unlist(order.arma.l[[name.test]][[1]][x,]))))
}
colnames(six.model) <- list.test
# remove dupplicated GPS-ERA
six.model$gps.era[which(is.na(reduced.list$chose)==TRUE)] =NA
# ind1 = which(reduced.list$min.var >0.05 & reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5)
# apply filter on length and variance
ind1 = which(reduced.list$min.var >0.002& reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5 & reduced.list$distances>50)
six.model = six.model[ind1,]
# save(six.model, file = paste0(path_results,"attribution/six.models", win.thres,".RData"))

six.values = c()
for (i in 1:length(list.test)) {
  value.count = sapply(c(list.model), function(x) length(which(six.model[,i] == x)))
  six.values <- c( six.values, value.count)
}
res.plot = data.frame(series = rep(list.name.test, each = 4), mod = rep(list.model, 6), value = six.values, 
                      n = c(rep(length(which(is.na(six.model$gps.era)==FALSE)),4), rep(nrow(six.model),20)))
res.plot$pct = res.plot$value/res.plot$n*100
res.plot$series = factor(res.plot$series, 
                         levels=reoder.list.name)
res.plot$mod = factor(res.plot$mod, 
                      levels=list.model)
# res.plot = res.plot[which(res.plot$value != 0),]
# jpeg(paste0(path_results,"attribution/iden_model_longest.jpg" ),width = 3000, height = 1800,res = 300)
p <- ggplot(res.plot, aes(fill=mod, y=pct, x=series, label = value)) + 
  geom_bar(position="dodge", stat="identity", width = 0.5)+theme_bw()+ 
  xlab("") + ylab("Percentage")+
  geom_text(position = position_dodge(width = .5),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 2)+
  ylim(c(0,100))+
  theme(axis.text = element_text(size = 5),legend.text=element_text(size=4.5),
        axis.title = element_text(size = 5), legend.key.size = unit(0.2, "cm"), 
        legend.title=element_blank())
ggsave(paste0(path_results,"attribution/Datacharacterization_autoarima_test.jpg" ), plot = p, width = 14, height = 5, units = "cm", dpi = 1200)

# Plot coefficients ------------------------
all.dat = get(load(file = paste0(path_results, "attribution0/all.dat.longest", win.thres,".RData")))
order.arma.l = get(load(file = paste0(path_results,"attribution0/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution0/coef.model.arma", win.thres,".RData")))
moving.var = get(load(file = paste0(path_results, "moving.var.RData")))

n=length(coef.arma.l$gps.era[[1]][,1])

param.list <- c()
model.list <- c()
test.list  <- c()
values <- c()
list.param = c("", "theta")
list.name.station = c()
distance = c()
mean.var = c()
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
      distance = c(distance, rep(reduced.list$distances[i],length(param)))
      mean.var = c(mean.var, rep(moving.var$mean[[name.test]][i],length(param)))
    }
  }
}

dat.p = data.frame(name = test.list, param = param.list, model = model.list, value = unlist(values), 
                   station = list.name.station, distance = distance, var = mean.var)
dat.p$name = factor(dat.p$name,  levels = reoder.list.name)
# sort data
ind1 = which(reduced.list$min.var>0.002&reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5& reduced.list$distances<50)

name.selected = reduced.list$station[ind1]
dat.p = dat.p[which(dat.p$station %in% name.selected),]

p <- ggplot(data = dat.p, aes( x = name, y = value, fill = model ,col = param)) + theme_bw()+
  geom_boxplot(lwd=0.2, outlier.size=0.2)+
  xlab("") + ylab(" values of parameters ") +
  theme(axis.text = element_text(size = 5),legend.text=element_text(size=3),
        axis.title = element_text(size = 5), legend.key.size = unit(0.2, "cm"), 
        legend.title=element_blank())+
  scale_color_manual(values = c("green", "deepskyblue4"))+
  scale_fill_manual(values = c("#D6604D", "#FDDBC7", "#92C5DE"))
ggsave(paste0(path_results,"attribution/Datacharacterization_test1.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 1200)


# # estimate using ARMA(1,1) in general -------------------------------------
# 
# all.coef = get(load( file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData")))
# all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))
# 
# coef.arma.l <- list()
# for (testi in c(1:6)) {
#   name.test = list.test[testi]
#   coef.arma = data.frame(matrix(NA, ncol = 4, nrow = length(all.dat)))
#   for (i in c(1:nrow(reduced.list))) {
#     # remove the case of white noise bc it will cause an error
#     if(testi ==4 & i ==49){
#       coef.arma[i,] = c(0,0,1,1)
#     }else{
#       name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
#       dat.i = all.dat[[name.i]]
#       y = dat.i[, paste0(name.test, "res")]/sqrt(dat.i[, paste0(name.test, "var")]) 
#       arma11 = Arima( y, order = c(1,0,1), include.mean = FALSE)
#       test.signif = as.data.frame(coeftest( arma11 )[,])
#       coef.arma[i,] = c(arma11$coef[c(1,2)], test.signif$`Pr(>|z|)`[1:2])
#     }
#   }
#   coef.arma.l[[name.test]] <- list(coef.arma)
# }
# save(coef.arma.l, file = paste0(path_results,"attribution/coef.model.arma", win.thres,"arma.RData"))
# all.coef = get(load( file = paste0(path_results,"attribution/coef.model.arma", win.thres,"arma.RData")))
# six.model = get(load(file = paste0(path_results,"attribution/six.models", win.thres,".RData")))
# 
# param.list <- c()
# model.list <- c()
# test.list  <- c()
# values <- c()
# p.val <- c()
# for (i in c(1:6)) {
#   name.test = list.test[i]
#   model.list = c(model.list, rep(six.model[,i], each = 2))
#   coef.df = all.coef[[name.test]][[1]]
#   values = c(values, as.vector(rbind(unlist(coef.df$X1), unlist(coef.df$X2))))
#   p.val = c(p.val, as.vector(rbind(unlist(coef.df$X3), unlist(coef.df$X4))))
#   test.list = c( test.list, rep(name.test, 2*nrow(six.model)))
#   param.list = c(param.list, rep(c( "phi", "theta"), nrow(six.model)))
# }
# 
# dat.p = data.frame(name = test.list, param = param.list, model = model.list, value = unlist(values), p = p.val)
# dat.p$name = rep(list.name.test, each = nrow(dat.p)/6)
# dat.p$name = factor(dat.p$name,  levels = reoder.list.name)
# ggplot(data = dat.p, aes( x = name, y = value, fill = model ,col = param)) + theme_bw()+
#   geom_boxplot()+
#   xlab("") + ylab(" values of parameters ") +
#   theme(axis.text = element_text(size = 16),legend.text=element_text(size=12),
#         axis.title = element_text(size=16))+
#   scale_color_manual(values=c("red", "green", "purple"))
# 
# # impact of gaps on model characterization 
# all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres= 1,".RData")))
# nb.c = sapply(c(1:52), function(x) nb.consecutive(list.day = all.dat[[x]]$date, x = all.dat[[x]]$gps.gpsres))
# six.model$nb.c = nb.c
# six.model$len = sapply(c(1:52), function(x) nrow(remove_na_2sides(df = all.dat[[x]], name.series = "gps.gpsres")))
# 
# d = data.frame(lim1 = unlist(six.model$r), lim10 = unlist(check1$r), l1 = six.model$nb.c, l10 = unlist(check1$L))
# d$white1 = sapply(c(1:52), function(x) length(which(six.model[x,c(1:6)] == "White")))
# colnames(d) = c("rate1", "rate10", "nb.cons1", "nb.cons10", "nb.white1", "nb.white10")
# save(d, file = paste0(path_results, "attribution/length_white_relation.RData"))
# 
# 

# # investigate the ARMA(1,1)-----------------------------
# # read data case
# six.model = get(load(file = paste0(path_results,"attribution/six.models", win.thres,".RData")))
# all.dat = get(load(file = paste0(path_results, "attribution/all.dat.longest", win.thres,".RData")))
# name.series = "gps.gps"
# dat.plot = remove_na_2sides(all.dat$`albh.2015-12-29.sc02`, name.series = name.series)
# y = unlist(dat.plot[paste0(name.series, "res")]/sqrt(dat.plot[paste0(name.series, "var")]))
# # plot data and it acf and pacf 
# plot(unlist(dat.plot[name.series]), ylab = "raw", type = "l", col = "gray")
# lines(unlist(dat.plot[paste0(name.series, "fit")]))
# plot(y, ylab = "normalized residual", col = "coral")
# plot(unlist(dat.plot[paste0(name.series, "var")]), ylab = "moving variance", col = "blue")
# acf(y, na.action = na.exclude)
# pacf(y, na.action = na.exclude)
# # auto.arima results
# fit.b = forecast::auto.arima(y , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
#                              max.p = 1, max.q = 1, start.p = 0, trace = TRUE, allowdrift = FALSE,  approximation=FALSE)
# # fit model
# ar1 = forecast::Arima(y , order = c(1,0,0), include.mean = FALSE)
# arma1 = forecast::Arima(y , order = c(1,0,1), include.mean = FALSE)
# ma1 = forecast::Arima(y , order = c(0,0,1), include.mean = FALSE)
# # plot p value of boxtest
# p.val = data.frame(matrix(NA, ncol = 2, nrow = 35))
# for (l in c(1:35)) {
#   p.val[l,] = c(Box.test(ar1$residuals, lag = l)$p.value, Box.test(arma1$residuals, lag = l)$p.value)
# }
# plot(p.val$X1, ylab = "p.value AR(1)")
# abline(h = 0.05)
# plot(p.val$X2, ylab = "p.value MA(1,1)")
# abline(h = 0.05)
# 
# # ACF, PACF of the residual
# acf(ar1$residuals, na.action = na.exclude)
# pacf(ar1$residuals, na.action = na.exclude)
# acf(arma1$residuals, na.action = na.exclude)
# pacf(arma1$residuals, na.action = na.exclude)
# # plot BIC
# mod_capt <- capture.output(forecast::auto.arima(y , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean = FALSE,lambda = NULL,
#                                                 max.p = 2, max.q = 2, start.p = 0, trace = TRUE, allowdrift = FALSE,  approximation=FALSE))
# 
# ind.c = c(2:10)
# a = sapply(ind.c, function(x) as.numeric(strsplit(mod_capt[x], split = ":")[[1]][2]))
# b = sapply(ind.c, function(x) strsplit(strsplit(mod_capt[x], split = ":")[[1]][1], split = " ")[[1]][2])
# d = sapply(ind.c, function(x) strsplit(strsplit(mod_capt[x], split = ":")[[1]][1], split = " ")[[1]][4]) 
# data.mod = data.frame(model =b, res =a, mean = d)
# ggplot(data = data.mod, aes(x = model, y = res, col = mean))+ theme_bw()+
#   geom_point()+ylab("BIC")+ theme(axis.text = element_text(size = 10))
# 
# 
# # plot theoretical ACF and PACF
# a = ARMAacf(ar = 0.2, ma =0.28, lag.max = 35, pacf = FALSE)
# df <- data.frame(lag = c(0:35), acf = a)
# s = c(qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(y))), - qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(y))))
# ggplot(data = df, mapping = aes(x = lag, y = acf)) + theme_bw()+
#   geom_hline(aes(yintercept = 0)) + ylab("acf")+
#   geom_hline(yintercept = s, linetype = 2)+
#   geom_segment(mapping = aes(xend = lag, yend = 0))+
#   theme(axis.title=element_text(size=15,face="bold"),
#     axis.text = element_text(size = 15))
# 
# # plot periodogram
# smooth.spec <- spec.pgram(x)
# spec.df <- data.frame(freq = smooth.spec$freq, spec = smooth.spec$spec)
# ggplot(data = spec.df) + theme_bw()+
#   geom_line(aes(x = freq, y = spec)) + 
#   scale_x_log10("Period (years)", 
#                 breaks = yrs.freqs, labels = yrs.labels) + scale_y_log10()
# 
# 
# yrs.period <- rev(c(1/12, 1/6, 1/5, 1/4, 1/3, 0.5, 1, 3, 4))
# yrs.labels <- rev(c( "1/12", "1/6", "1/5", "1/4", "1/3", "1/2", "1", "3", "4"))
# yrs.freqs <- 1/yrs.period * 1/365
# 
# 
# # fit ARMA for 4 series  --------------------------------------------------
# data.test = dat$`guam.2017-09-26.guug`
# be = data.test[c(1:3650),]
# af = data.test[c(3651:7300),]
# s = "gps.gps"
# 
# m = construct.design(be, name.series = s)
# fit.igls = IGLS(design.m = m, tol =  tol0, day.list = be$date)
# y1 = fit.igls$residual/sqrt(fit.igls$var)
# fit.arima(y1)
# m = construct.design(af, name.series = s)
# fit.igls = IGLS(design.m = m, tol =  tol0, day.list = af$date)
# y = fit.igls$residual/sqrt(fit.igls$var)
# fit.arima(y)
# 
# data.test$res = c(y1, y)
# m = construct.design(data.test, name.series = "res")
# fit.igls = IGLS(design.m = m, tol =  tol0, day.list = data.test$date)
# y = fit.igls$residual/sqrt(fit.igls$var)
# fit.arima(y)
# 
# 
# res = data.frame(matrix(NA, ncol = 3, nrow = 1000))
# res.coef = data.frame(matrix(NA, ncol = 4, nrow = 1000))
# list.ser = c("gps","gps1","era","era1")
# for (i in c(1:nrow(reduced.list))) {
#   name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
#   dat.i = dat[[name.i]]
#   dat.i = dat.i[choose_segment(reduced.list$chose[i]),]
#   dat.ij = remove_na_2sides(dat.i, name.series = "gps.gps")
#   data.test = dat.ij
#   a = read.series(path_series_main, station = reduced.list$main[i], na.rm = 0, add.full = 1)
#   b = left_join(data.test, a, by = "date")
#   data.test$gps = b$GPS
#   data.test$era = b$ERAI
#   data.test$gps1 = data.test$gps - data.test$gps.gps
#   data.test$era1 = data.test$era - data.test$era.era
#   for (j in c(1:4)) {
#     name.series0 = list.ser[j]
#     m = construct.design(data.test, name.series = name.series0)
#     tol0 = 0.00001
#     fit.igls = IGLS(design.m = m, tol =  tol0, day.list = data.test$date)
#     y = fit.igls$residual/sqrt(fit.igls$var)
#     fit.arma = fit.arima(y)
#     res[(52*(j-1)+i),] = fit.arma$pq
#     res.coef[(52*(j-1)+i),] = fit.arma$coef
#     print(j)
#   }
# }
# 
# save(res, file=paste0(path_results,"attribution/model.4series.RData"))
# 
# save(res.coef, file=paste0(path_results,"attribution/coef.4series.RData"))
# 
# res = get(load(file=paste0(path_results,"attribution/model.4series.RData")))
# res.coef = get(load(file=paste0(path_results,"attribution/coef.4series.RData")))
# 
# a = sapply(c(1:208), function(x) model.iden(as.numeric(unlist(b[x,]))))
# 
# b = na.omit(res)
# d = na.omit(res.coef)
# 
# e = d[which(d$model=="AR(1)"),]
# f = d[which(d$model=="ARMA(1,1)"),]
# 
# d$model = a
# d$series = rep(list.ser, each = 52)
# 
# 
# 
# 
# # study seasonal bias -----------------------------------------------------
# all.coef = get(load( file = paste0(path_results, "attribution/all.coef.longest", win.thres,".RData")))
# tot.s = list()
# for (i in c(1:6)) {
#   name.series = list.test[i]
#   dat.coef = as.data.frame(t(as.data.frame(all.coef[[name.series]])))
#   dat.coef$oneyear = sqrt(dat.coef$cos1^2 + dat.coef$sin1^2)
#   dat.coef$halfyear = sqrt(dat.coef$cos2^2 + dat.coef$sin2^2)
#   dat.coef$thirdyear = sqrt(dat.coef$cos3^2 + dat.coef$sin3^2)
#   dat.coef$quartyear = sqrt(dat.coef$cos4^2 + dat.coef$sin4^2)
#   tot.s[[name.series]] = dat.coef
# }
# a = data.frame(value = sapply(c(1:6), function(x) tot.s[[x]]$haftyear))
# colnames(a)[1:6] = list.name.test
# b = reshape2::melt(a)
# b$series = factor(b$variable,  levels = reoder.list.name)
# 
# library(RColorBrewer)
# p <- ggplot(b, aes(x = value, col = series ))+ theme_bw()+
#   stat_ecdf(lwd = 0.3)+
#   scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))+
#   geom_hline(yintercept = 0.5, size = 0.3) +
#   scale_color_manual(values = brewer.pal(n = 6, name = 'Dark2'))+
#   labs(y = "CDF", x = "Moving window standard deviation", linetype = "")+
#   theme(axis.text = element_text(size = 5),legend.text=element_text(size=4.5),
#         axis.title = element_text(size=5), legend.key.size = unit(0.15, "cm"),
#         legend.box.spacing = unit(0, "pt"), legend.title= element_blank())
# print(p)
# # dev.off()
# ggsave(paste0(path_results,"attribution/coef_haftyear.jpg" ), plot = p, width = 8.8, height = 6, units = "cm", dpi = 1200)
# # study the variance 
# a = get(load( file = paste0(path_results, "moving.var.RData")))
# 
# 
# 
# 
# 
# 
# 
# # WHAT HAPPEN IF WE CONCATENATE ALL DATA BEFORE AND AFTER  ----------------
# 
# all.coef = list()
# all.dat = list()
# all.fit = list()
# # for (i in c(1:nrow(reduced.list))) {
# #   name.i =  paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
# for (i in c(1:length(dat))) {
#   name.i = names(dat)[i]
#   dat.i = dat[[name.i]]
#   # dat.ij = remove_na_2sides(dat.i, name.series = "gps.gps")
#   six.ser = na.omit(dat.i)
#   dat.i = tidyr::complete(six.ser, date = seq(dat.i$date[1], dat.i$date[7300], by = "day"))
#   print(i)
#   for (j in c(1:6)) {
#     name.series0 = list.test[j]
#     dat.ib = dat.i[c(1:3650),]
#     dat.ia = dat.i[c(3651:7300),]
#     
#     mb = construct.design(dat.ib, name.series = name.series0)
#     ma = construct.design(dat.ia, name.series = name.series0)
#     
#     tol0 = 0.000001
#     if(i %in% c(30,521)){
#       tol0 = 0.001 
#     }
#     
#     fit.iglsa = IGLS(design.m = ma, tol =  tol0, day.list = dat.ia$date)
#     fit.iglsb = IGLS(design.m = mb, tol =  tol0, day.list = dat.ib$date)
#     
#     dat.i[, paste0(name.series0, 'var')] <- c( unlist(fit.iglsb$var), unlist(fit.iglsa$var))
#     dat.i[, paste0(name.series0, 'res')] <- c( unlist(fit.iglsb$residual), unlist(fit.iglsa$residual))
#     dat.i[, paste0(name.series0, 'fit')] <- c( unlist(fit.iglsb$fit), unlist(fit.iglsa$fit))
#     all.coef[[name.series0]][[name.i]] = list(bef = fit.iglsb$coefficients, aft = fit.iglsa$coefficients)
#     print(j)
#   }
#   all.dat[[name.i]] = dat.i
# }
# save(all.coef, file = paste0(path_results, "attribution/all.coef.longest", win.thres,"full.RData"))
# save(all.dat, file = paste0(path_results, "attribution/all.dat.longest", win.thres,"full.RData"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


# PLOT FOR PAPER ----------------------------------------------------------
# put 4 plots in the same figures 
library(ggplot2)   
library(gtable)    
library(grid)
library(gridExtra) 
unicode_minus = function(x) sub('^-', '\U2212', format(x))
# unicode_minus(format(round(SE,  digits = 3), nsmall = 3))

all.dat = get(load(file = paste0(path_results, "attribution0/all.dat.longest", win.thres,".RData")))
order.arma.l = get(load(file = paste0(path_results,"attribution0/order.model.arma", win.thres,".RData")))
coef.arma.l = get(load(file = paste0(path_results,"attribution0/coef.model.arma", win.thres,".RData")))

list.model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)")

text1 = "Distance < 50 km"
text2 = "Distance > 50 km"

# plot1-----
length.data =nrow(reduced.list)
six.model = data.frame(matrix(NA, ncol = 6, nrow = length.data))
for (i in 1:length(list.test)) {
  name.test = list.test[i]
  six.model[,i] = sapply(c(1:length.data), function(x) model.iden(as.numeric(unlist(order.arma.l[[name.test]][[1]][x,]))))
}
colnames(six.model) <- list.test
six.model$gps.era[which(is.na(reduced.list$chose)==TRUE)] =NA
ind1 = which(reduced.list$min.var >0.002& reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5 & reduced.list$distances<50)
six.model = six.model[ind1,]

six.values = c()
for (i in 1:length(list.test)) {
  value.count = sapply(c(list.model), function(x) length(which(six.model[,i] == x)))
  six.values <- c( six.values, value.count)
}
res.plot = data.frame(series = rep(list.name.test, each = 4), mod = rep(list.model, 6), value = six.values, 
                      n = c(rep(length(which(is.na(six.model$gps.era)==FALSE)),4), rep(nrow(six.model),20)))
res.plot$pct = res.plot$value/res.plot$n*100
res.plot$series = factor(res.plot$series, 
                         levels=reoder.list.name)
res.plot$mod = factor(res.plot$mod, 
                      levels=list.model)

p1 <- ggplot(res.plot, aes(fill=mod, y=pct, x=series, label = value)) + 
  geom_bar(position="dodge", stat="identity", width = 0.5)+theme_bw()+ 
  xlab("") + ylab("Percentage")+
  labs(tag = "(a)", subtitle = text1) +
  geom_text(position = position_dodge(width = .5),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 1)+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.2, "cm"), 
        plot.tag = element_text(size = 6) , plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
# plot2-----

length.data =nrow(reduced.list)
six.model = data.frame(matrix(NA, ncol = 6, nrow = length.data))
for (i in 1:length(list.test)) {
  name.test = list.test[i]
  six.model[,i] = sapply(c(1:length.data), function(x) model.iden(as.numeric(unlist(order.arma.l[[name.test]][[1]][x,]))))
}
colnames(six.model) <- list.test
# remove dupplicated GPS-ERA
six.model$gps.era[which(is.na(reduced.list$chose)==TRUE)] =NA
ind1 = which(reduced.list$min.var >0.002& reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5 & reduced.list$distances>50)
six.model = six.model[ind1,]

six.values = c()
for (i in 1:length(list.test)) {
  value.count = sapply(c(list.model), function(x) length(which(six.model[,i] == x)))
  six.values <- c( six.values, value.count)
}
res.plot = data.frame(series = rep(list.name.test, each = 4), mod = rep(list.model, 6), value = six.values, 
                      n = c(rep(length(which(is.na(six.model$gps.era)==FALSE)),4), rep(nrow(six.model),20)))
res.plot$pct = res.plot$value/res.plot$n*100
res.plot$series = factor(res.plot$series, 
                         levels=reoder.list.name)
res.plot$mod = factor(res.plot$mod, 
                      levels=list.model)

p2 <- ggplot(res.plot, aes(fill=mod, y=pct, x=series, label = value)) + 
  geom_bar(position="dodge", stat="identity", width = 0.5)+theme_bw()+ 
  xlab("") + ylab("Percentage")+
  labs(tag = "(b)", subtitle = text2) + 
  geom_text(position = position_dodge(width = .5),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 1)+
  ylim(c(0,100))+
  theme(axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.2, "cm"), plot.tag = element_text(size = 6),
        legend.title=element_blank(), plot.subtitle = element_text(size = 6), legend.box.spacing = unit(0, "pt"), 
        plot.margin = rep(unit(0,"null"),4))

# plot3 -----
length.data =nrow(reduced.list)
six.model = data.frame(matrix(NA, ncol = 6, nrow = length.data))
for (i in 1:length(list.test)) {
  name.test = list.test[i]
  six.model[,i] = sapply(c(1:length.data), function(x) model.iden(as.numeric(unlist(order.arma.l[[name.test]][[1]][x,]))))
}
colnames(six.model) <- list.test
six.model$gps.era[which(is.na(reduced.list$chose)==TRUE)] =NA
n=length(coef.arma.l$gps.era[[1]][,1])
param.list <- c()
model.list <- c()
test.list  <- c()
values <- c()
list.param = c("Phi", "theta")

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
# sort data
ind1 = which(reduced.list$min.var>0.002&reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5& reduced.list$distances<50)

name.selected = reduced.list$station[ind1]
dat.p = dat.p[which(dat.p$station %in% name.selected),]

p3 <- ggplot(data = dat.p, aes( x = name, y = value, fill = model ,col = param)) + theme_bw()+
  geom_boxplot(lwd=0.2, outlier.size=0.2)+
  xlab("") + ylab(" values of coefficients ") + ylim(-1,1)+ labs(tag = "(c)") + 
  theme(axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.2, "cm"), plot.tag = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))+
  scale_color_manual(values = c("green", "deepskyblue4"), labels = expression(Phi, theta))+
  scale_fill_manual(values = c("#D6604D", "#FDDBC7", "#92C5DE"))

# plot4 ---------
dat.p = data.frame(name = test.list, param = param.list, model = model.list, value = unlist(values), station = list.name.station)
dat.p$name = factor(dat.p$name,  levels = reoder.list.name)
# sort data
ind1 = which(reduced.list$min.var>0.002&reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5& reduced.list$distances>50)

name.selected = reduced.list$station[ind1]
dat.p = dat.p[which(dat.p$station %in% name.selected),]

p4 <- ggplot(data = dat.p, aes( x = name, y = value, fill = model ,col = param)) + theme_bw()+
  geom_boxplot(lwd=0.2, outlier.size=0.2)+
  xlab("") + ylab(" values of coefficients ") + ylim(-1,1)+ labs(tag = "(d)") + 
  theme(axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.2, "cm"), plot.tag = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))+
  scale_color_manual(values = c("green", "deepskyblue4"), labels = expression(Phi, theta))+
  scale_fill_manual(values = c("#D6604D", "#FDDBC7", "#92C5DE"))

# grid of 4 ------
gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
gC <- ggplotGrob(p3)
gD <- ggplotGrob(p4)

gB$widths <- gA$widths
gC$widths <- gA$widths
gD$widths <- gA$widths
# Arrange the two charts.

grid.newpage()
p = (grid.arrange(gA, gB, gC, gD,nrow = 2))

ggsave(paste0(path_results,"attribution/Autoarima.jpg" ), plot = p, width = 14.4, height = 9, units = "cm", dpi = 1200)




# table of variance for paper ---------------------------------------------

all.coef = get(load( file = paste0(path_results, "attribution0/all.coef.longest", win.thres,".RData")))
all.dat = get(load(file = paste0(path_results, "attribution0/all.dat.longest", win.thres,".RData")))
mean.res = data.frame(matrix(NA, ncol = 0, nrow = nrow(reduced.list)))
range.res = data.frame(matrix(NA, ncol = 0, nrow = nrow(reduced.list)))

for (i in c(1:nrow(reduced.list))) {
  name.i = paste0(reduced.list$main[i],".",as.character(reduced.list$brp[i]), ".", reduced.list$nearby[i])
  dat.i = all.dat[[name.i]]
  for (j in c(1:6)) {
    name.series = list.test[j]
    sd.ij = sqrt(unlist(dat.i[,paste0(name.series, 'var')]))
    mean.res[i,name.series] = mean(sd.ij, na.rm = TRUE)
    range.res[i,name.series]  = range.var(x = sd.ij , day.list = dat.i$date, s = dat.i[name.series])
  }
}
mean.res$gps.era[which(is.na(reduced.list$chose)==TRUE)] =NA
range.res$gps.era[which(is.na(reduced.list$chose)==TRUE)] =NA
tot.res = cbind(mean.res, range.res)

ind1 = which(reduced.list$min.var>0.002&reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5)
ind2 = which(reduced.list$min.var>0.002&reduced.list$nbc>1000 & reduced.list$r1>0.5 & reduced.list$r2>0.5 & reduced.list$distances >50)
ind = ind1
res = data.frame(mean = colMeans(tot.res[ind,c(1:6)], na.rm = TRUE),
                 sd.mean = sapply(c(1:6), function(x) sd(tot.res[ind,x], na.rm = TRUE)), 
                 range = colMeans(tot.res[ind,c(7:12)]/ tot.res[ind,c(1:6)], na.rm = TRUE),
                 sd.range = sapply(c(1:6), function(x) sd(tot.res[ind,(x+6)], na.rm = TRUE)),
                 name = list.name.test)

a = res[match(reoder.list.name, res$name),]
a[,c(1:4)]= round(a[,c(1:4)], digits = 2)
b = a

ind = ind2
res = data.frame(mean = colMeans(tot.res[ind,c(1:6)], na.rm = TRUE),
                 sd.mean = sapply(c(1:6), function(x) sd(tot.res[ind,x], na.rm = TRUE)), 
                 range = colMeans(tot.res[ind,c(7:12)], na.rm = TRUE),
                 sd.range = sapply(c(1:6), function(x) sd(tot.res[ind,(x+6)], na.rm = TRUE)),
                 name = list.name.test)
a = res[match(reoder.list.name, res$name),]
a[,c(1:4)]= round(a[,c(1:4)], digits = 2)

out = data.frame(cbind(b[,c(1:2)], a[,c(1:2)]), b[,c(3:4)], a[,c(3:4)], b[,5])
write.table(a, file = paste0(path_results,"attribution/sd_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# percentage pf range 

for (i in c(1:6)) {
  print(sd(range.res[,i]/mean.res[,i], na.rm=TRUE ))
}




