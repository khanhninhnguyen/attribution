# study the performance of auto.arima 
source(paste0(path_code_att,"FGLS.R"))
source(paste0(path_code_att,"simulate_time_series.R"))
source(paste0(path_code_att,"support_characterization.R"))

# As a function of length -------------------------------------------------
length.list = seq(200, 2000, 200)
coef.list = seq(0, 0.8, 0.1)

nb.sim = 10000
tot.res = data.frame(n = length.list, 
                     TPR.white = rep(NA, length(length.list)),
                     TPR.ar = rep(NA, length(length.list)),
                     TPR.ma = rep(NA, length(length.list)),
                     TPR.arma = rep(NA, length(length.list)))
burn.in = 1000
hetero = 0
sigma.sim = 1
## As a function of length -------------------------------------------------
set.seed(1)
tot.fit = list()
for (i in c(1:length(length.list))) {
  n = length.list[i]
  TPR = data.frame(matrix(NA, ncol = 4, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.white = rnorm(n, mean = 0, sd = 1)
    y.ar = simulate.general1(N = n, arma.model = c(ar=0.3,ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.ma = simulate.general1(N = n, arma.model = c(ar=0,ma=0.2), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.arma = simulate.general1(N = n, arma.model = c(ar=0.6,ma=-0.3), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit
    fit.ar = fit.arima(y.ar)
    fit.ma = fit.arima(y.ma)
    fit.arma = fit.arima(y.arma)
    fit.white = fit.arima(y.white)
    TPR[j,] = c( model.iden(fit.ar$pq), model.iden(fit.ma$pq), model.iden(fit.arma$pq), model.iden(fit.white$pq))
    fit.i[[j]] = list( ar = fit.ar, ma = fit.ma, arma = fit.arma, white = fit.white)
  }
  colnames(TPR) = c("ar", "ma", "arma", "white")
  tot.res[i,c(2:5)] = c(length(which(TPR$white == "White")), 
                  length(which(TPR$ar == "AR(1)")), 
                  length(which(TPR$ma == "MA(1)")),
                  length(which(TPR$arma == "ARMA(1,1)")))
  tot.fit[[i]] = fit.i
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_length.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_length_all.RData"))


## As a function of coefficients -------------------------------------------------

tot.res = data.frame(coef = coef.list, 
                     TPR.ar = rep(NA, length(coef.list)),
                     TPR.ma = rep(NA, length(coef.list)),
                     TPR.arma = rep(NA, length(coef.list),
                     TPR.arma1 = rep(NA, length(coef.list)),
                     ))
set.seed(1)
tot.fit = list()
n=1000
for (i in c(1:length(coef.list))) {
  phi = coef.list[i]
  theta = 0.3 - phi
  TPR = data.frame(matrix(NA, ncol = 4, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.ar = simulate.general1(N = n, arma.model = c(ar=phi, ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.ma = simulate.general1(N = n, arma.model = c(ar=0, ma=phi), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.arma = simulate.general1(N = n, arma.model = c(ar=phi, ma=theta), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    y.arma1 = simulate.general1(N = n, arma.model = c(ar=phi, ma=(phi -0.1)), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit
    fit.ar = fit.arima(y.ar)
    fit.ma = fit.arima(y.ma)
    fit.arma = fit.arima(y.arma)
    fit.arma1 = fit.arima(y.arma1)
    
    TPR[j,] = c(model.iden(fit.ar$pq), model.iden(fit.ma$pq), model.iden(fit.arma$pq), model.iden(fit.arma1$pq))
    fit.i[[j]] = list( ar = fit.ar, ma = fit.ma, arma = fit.arma, arma1 = fit.arma1)
  }
  colnames(TPR) = c("ar", "ma", "arma", "arma1")
  tot.res[i,c(2:5)] = c(length(which(TPR$ar == "AR(1)")), 
                        length(which(TPR$ma == "MA(1)")),
                        length(which(TPR$arma == "ARMA(1,1)")),
                        length(which(TPR$arma1 == "ARMA(1,1)")))
  tot.fit[[i]] = fit.i
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_coef_all.RData"))


# a = get(load( file = paste0(path_results, "attribution0/performance_autoarima_length.RData")))
# colnames(a) = c("N", "AR(1)", "MAR(1)", "ARMA(1,1)")
# dat.p = reshape2::melt(a, id="N")
# p = ggplot(data = dat.p, aes(x = N, y = value/nb.sim, col = variable))+
#   theme_bw()+
#   geom_point(size = 0.3)+
#   geom_hline(yintercept = 0.95, lwd = 0.3)+
#   ylab("TPR")+
#   theme(axis.text.x = element_text(size = 5), 
#       axis.text.y = element_text(size = 5),
#       legend.text=element_text(size=4),
#       axis.title = element_text(size = 5), 
#       legend.key.size = unit(0.3, "cm"), 
#       plot.tag = element_text(size = 5), 
#       plot.subtitle = element_text(size = 5),
#       legend.title=element_blank())
# 
# ggsave(paste0(path_results,"attribution0/TPR.auto.arima.identification_length.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)


# coefficient dependence AR(1) --------------------------------------------------
coef.list = seq(0, 0.8, 0.1)
set.seed(1)
TPR.ar <- rep(NA, nb.sim)
tot.res <- list()
tot.fit = list()

for (i in c(1:length(coef.list))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.ar = simulate.general1(N = n, arma.model = c(ar=coef.list[i],ma=0), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.ar = fit.arima(y.ar)
    TPR.ar[j] = model.iden(fit.ar$pq)
    fit.i[[j]] = fit.ar
  }
  tot.res[[i]] = TPR.ar
  tot.fit[[i]] = fit.i
}


save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_ar.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_coef_ar_all.RData"))

# b = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_ar.RData")))
# coef.list = seq(0.1, 0.8, 0.1)
# 
# tpr = sapply(c(1:length(b)), function(x) length(which(b[[x]]=="AR(1)")))[-1]/nb.sim
# dat.p = data.frame(ar = coef.list, tpr = tpr)
# p = ggplot(data = dat.p, aes(x = ar, y = tpr))+
#   theme_bw()+
#   geom_point(size = 0.3)+
#   geom_hline(yintercept = 0.95, lwd = 0.3)+
#   ylab("TPR")+
#   theme(axis.text.x = element_text(size = 5), 
#         axis.text.y = element_text(size = 5),
#         legend.text=element_text(size=4),
#         axis.title = element_text(size = 5), 
#         legend.key.size = unit(0.3, "cm"), 
#         plot.tag = element_text(size = 5), 
#         plot.subtitle = element_text(size = 5),
#         legend.title=element_blank())
# ggsave(paste0(path_results,"attribution0/TPR.auto.arima.identification_ar.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

# coefficient dependence MA(1)  --------------------------------------------------
coef.list = seq(0, 0.8, 0.1)
set.seed(1)
TPR.ma <- rep(NA, nb.sim)
tot.res <- list()
tot.fit = list()
for (i in c(1:length(coef.list))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.ma = simulate.general1(N = n, arma.model = c(ar=0,ma=coef.list[i]), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.ma = fit.arima(y.ma)
    TPR.ma[j] = model.iden(fit.ma$pq)
    fit.i[[j]] = fit.ma
  }
  tot.res[[i]] = TPR.ma
  tot.fit[[i]] = fit.i
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_ma.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_coef_ma_all.RData"))

# coefficient dependence ARMA(1,1) --------------------------------------------------
coef.list.ar = c(0.1, 0.2, 0.3)
coef.list.ma = c(0.1, 0.2, 0.3)
tot.fit = list()
set.seed(1)
TPR.ar <- rep(NA, nb.sim)
tot.res <- list()
for (i in c(1:length(coef.list.ar))) {
  n = 1000
  TPR = data.frame(matrix(NA, ncol = 3, nrow = nb.sim)) 
  fit.i = list()
  for (j in c(1:nb.sim)) {
    y.arma = simulate.general1(N = n, arma.model = c(ar=coef.list.ar[i],ma=coef.list.ma[i]), burn.in = burn.in, hetero = hetero, sigma = sqrt(sigma.sim))
    # fit 
    fit.arma = fit.arima(y.arma)
    TPR.ar[j] = model.iden(fit.arma$pq)
    fit.i[[j]] = fit.arma
  }
  tot.res[[i]] = TPR.ar
  tot.fit[[i]] = fit.i
}

save(tot.res, file = paste0(path_results, "attribution0/performance_autoarima_coef_arma1.RData"))
save(tot.fit, file = paste0(path_results, "attribution0/performance_autoarima_coef_arma_all1.RData"))



# plot --------------------------------------------------------------------
all.length = get(load(file = paste0(path_results, "attribution0/performance_autoarima_length_all.RData")))
all.ar1 = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_ar_all.RData")))
all.ma1 = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_ma_all.RData")))
all.arma11 = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_arma_all.RData")))
all.coef = get(load( file = paste0(path_results, "attribution0/performance_autoarima_coef_all.RData")))

extract_arima <- function(df.res, model, true.param, phi, length.c){
  # choose model
  if(model == "AR(1)"){name = "ar"}
  if(model == "MA(1)"){name = "ma"}
  if(model == "ARMA(1,1)"){name = "arma"}
  
  # choose param 
  if(phi==1){
    p = 1
  }else{p=3}
  if(model == "White"){
    name = "white"
    p = 1
  }
  # read all info
  n.iter = length(df.res)
  if(length(true.param)==1){ true.param = rep(true.param, n.iter)}
  
  if(length.c==1){
    model.iden = sapply(c(1:n.iter), function(y) sapply(c(1:length(df.res[[y]])), function(x) model.iden(df.res[[y]][[x]][[name]]$pq))) 
    phi = sapply(c(1:n.iter), function(y) sapply(c(1:length(df.res[[y]])), function(x) df.res[[y]][[x]][[name]]$coef[p])) 
  }else{
    model.iden = sapply(c(1:n.iter), function(y) sapply(c(1:length(df.res[[y]])), function(x) model.iden(df.res[[y]][[x]]$pq))) 
    phi = sapply(c(1:n.iter), function(y) sapply(c(1:length(df.res[[y]])), function(x) df.res[[y]][[x]]$coef[p])) 
  }
  
  res = data.frame(matrix(NA, ncol = 2, nrow = n.iter))
  ests = list()
  for (i in c(1:n.iter)) {
    ind.ar = which(model.iden[,i] == model)
    if(length(ind.ar)==0){
      rms = NA
    }else{
      rms = sqrt(sum((phi[ind.ar,i]-true.param[i])^2)/length(ind.ar))
    }
    res[i,] = c(length(ind.ar), rms)
    ests[[i]] = phi[ind.ar,i]
  }
  colnames(res) = c("tpr", "rmsd")
  return(list(tpr = res, est = ests))
}

## length---------------
#NEED TO BE CHECKED ONLY WHEN THEY IDENTIRY THE TRUE MODEL

white = extract_arima(all.length, model = "White", true.param = 0, phi = 1, length.c = 1)
ar1 = extract_arima(all.length, model = "AR(1)", true.param = 0.3, phi = 1, length.c = 1)
ma1 = extract_arima(all.length, model = "MA(1)", true.param = 0.2, phi = 0, length.c = 1)
arma.phi = extract_arima(all.length, model = "ARMA(1,1)", true.param = 0.6, phi = 1, length.c = 1)
arma.theta = extract_arima(all.length, model = "ARMA(1,1)", true.param = -0.3, phi = 0, length.c = 1)

### plot TPR -----------

TPR = data.frame(n = length.list, 
                 white = white$tpr$tpr,
                 ar = ar1$tpr$tpr, 
                 ma1 = ma1$tpr$tpr, 
                 arma1 = arma.phi$tpr$tpr)

colnames(TPR) = c("N", "White", "AR:0.3", "MA:0.2", "ARMA:0.6,-0.3")
dat.p = reshape2::melt(TPR, id="N")
p = ggplot(data = dat.p, aes(x = N, y = value/nb.sim, col = variable))+
  theme_bw()+
  geom_point(size = 0.3)+
  geom_hline(yintercept = 0.95, lwd = 0.3)+
  ylab("TPR")+
  scale_x_continuous(breaks = length.list, 
                     limits = c(200, 2000))+
  theme(axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      legend.text=element_text(size=4),
      axis.title = element_text(size = 5),
      legend.key.size = unit(0.3, "cm"),
      plot.tag = element_text(size = 5),
      plot.subtitle = element_text(size = 5),
      legend.title=element_blank())
# 
ggsave(paste0(path_results,"attribution0/auto.arima.identification_length.tpr.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

### plot coefficients ----------------
dat.ar = as.data.frame(plyr::ldply(ar1$est, rbind))  %>% 
  mutate(length = as.factor(length.list)) %>% 
  reshape2::melt(id = "length")
dat.ma = as.data.frame(plyr::ldply(ma1$est, rbind)) %>% 
  mutate(length = as.factor(length.list)) %>% 
  reshape2::melt(id = "length")
dat.arma.phi = as.data.frame(plyr::ldply(arma.phi$est, rbind)) %>% 
  mutate(length = as.factor(length.list)) %>% 
  reshape2::melt(id = "length")
dat.arma.theta = as.data.frame(plyr::ldply(arma.theta$est, rbind)) %>% 
  mutate(length = as.factor(length.list)) %>% 
  reshape2::melt(id = "length")
dat.arma = rbind(dat.arma.phi, dat.arma.theta) %>% 
  mutate(param = as.factor(rep(c("phi", "theta"), each = nrow(dat.arma.phi))))
  

p1 = ggplot(data = dat.ar, aes(x = length, y = value))+ theme_bw()+
  geom_boxplot(lwd = 0.3, outlier.size = 0.3)
p2 = ggplot(data = dat.ma, aes(x = length, y = value))+ theme_bw()+
  geom_boxplot(lwd = 0.3, outlier.size = 0.3)
p3 = ggplot(data = dat.arma, aes(x = length, y = value, fill = param))+ theme_bw()+
  geom_boxplot(lwd = 0.15, outlier.size = 0.15)

p = p3 + ylab("Coefficients")+
  geom_hline(yintercept = c(0.7,-0.4), lwd = 0.3, col = "gray")+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text=element_text(size=4),
        axis.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        plot.tag = element_text(size = 5),
        plot.subtitle = element_text(size = 5),
        legend.title=element_blank())

ggsave(paste0(path_results,"attribution0/estimates_length_arma1.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)


## coefficicents -----------------
coeff = seq(0, 0.8, 0.1)
ar1 = extract_arima(all.ar1, model = "AR(1)", true.param = coeff, phi = 1, length.c = 0)
ma1 = extract_arima(all.ma1, model = "MA(1)", true.param = coeff, phi = 0, length.c = 0)
arma.phi = extract_arima(all.arma11, model = "ARMA(1,1)", true.param = coeff, phi = 1, length.c = 0)
arma.theta = extract_arima(all.arma11, model = "ARMA(1,1)", true.param = coeff, phi = 0, length.c = 0)

ar1 = extract_arima(all.coef, model = "AR(1)", true.param = coeff, phi = 1, length.c = 1)
ma1 = extract_arima(all.coef, model = "MA(1)", true.param = coeff, phi = 0, length.c = 1)
arma.phi = extract_arima(all.coef, model = "ARMA(1,1)", true.param = coeff, phi = 1, length.c = 1)
arma.theta = extract_arima(all.coef, model = "ARMA(1,1)", true.param = coeff, phi = 0, length.c = 1)

### TPR------
TPR = data.frame(coef = coeff, 
                 ar = ar1$tpr$tpr, 
                 ma1 = ma1$tpr$tpr, 
                 arma1 = arma.phi$tpr$tpr)
TPR = TPR[-1,]
colnames(TPR) = c("coef", "AR(1)", "MA(1)", "ARMA(1,1)")
dat.p = reshape2::melt(TPR, id="coef")
p = ggplot(data = dat.p, aes(x = coef, y = value/nb.sim, col = variable))+
  theme_bw()+
  geom_point(size = 0.3)+
  geom_hline(yintercept = 0.95, lwd = 0.3)+
  ylab("TPR")+
  xlab("Coefficient")+
  # xlab(expression(phi))+
  scale_x_continuous(breaks = coeff[-1], 
                     limits = c(0.1, 0.8))+
                     # sec.axis = sec_axis(~ (0.3-.x), breaks = seq(0.2, -0.5, -0.1)), name = expression(theta))+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text=element_text(size=4),
        axis.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        plot.tag = element_text(size = 5),
        plot.subtitle = element_text(size = 5),
        legend.title=element_blank())
# 
ggsave(paste0(path_results,"attribution0/auto.arima.identification_coef.tpr.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)

### estimates impact ----------

dat.ar = as.data.frame(plyr::ldply(ar1$est, rbind))  %>% 
  mutate(coef = as.factor(coeff)) %>% 
  reshape2::melt(id = "coef")
dat.ma = as.data.frame(plyr::ldply(ma1$est, rbind)) %>% 
  mutate(coef = as.factor(coeff)) %>% 
  reshape2::melt(id = "coef")
dat.arma.phi = as.data.frame(plyr::ldply(arma.phi$est[4:9], rbind)) %>% 
  mutate(coef = as.factor(coeff[4:9])) %>% 
  reshape2::melt(id = "coef")
dat.arma.theta = as.data.frame(plyr::ldply(arma.theta$est[4:9], rbind)) %>% 
  mutate(coef = as.factor(coeff[4:9])) %>% 
  reshape2::melt(id = "coef")
dat.arma = rbind(dat.arma.phi, dat.arma.theta) %>% 
  mutate(param = as.factor(rep(c("phi", "theta"), each = nrow(dat.arma.phi)))) %>% 
  mutate(true = as.numeric(as.character(coef)))
dat.arma$true[which(dat.arma$param!= "phi")] = 0.3-dat.arma$true[which(dat.arma$param!= "phi")]
p1 = ggplot(data = dat.ar, aes(x = coef, y = value))+ theme_bw()+
  geom_boxplot(lwd = 0.15, outlier.size = 0.15, width=0.5 )
p2 = ggplot(data = dat.ma, aes(x = coef, y = value))+ theme_bw()+
  geom_boxplot(lwd = 0.15, outlier.size = 0.15, width=0.5 )
p3 = ggplot(data = dat.arma, aes(x = coef, fill = param))+ theme_bw()+
  geom_boxplot(aes( y = value), lwd = 0.15, outlier.size = 0.15, width=0.5 )+
  geom_point(aes( y = true, col = param), size = 0.15)

p = p3 + ylab("Estimates")+
  scale_y_continuous(breaks = seq(-0.8, 0.8, 0.2), limits = c(-0.8,0.8))+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text=element_text(size=4),
        axis.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        plot.tag = element_text(size = 5),
        plot.subtitle = element_text(size = 5),
        legend.title=element_blank())

ggsave(paste0(path_results,"attribution0/estimates_coef_arma1.jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)




# plot a specific case 
chosen = 4
print(paste0("Coefficient is chosen:", coeff[chosen]))
list.model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)")
ar1.iden = sapply(c(1:length(all.ar1[[chosen]])), function(x) model.iden(all.ar1[[chosen]][[x]]$pq))
ma1.iden = sapply(c(1:length(all.ma1[[chosen]])), function(x) model.iden(all.ma1[[chosen]][[x]]$pq))
arma.iden = sapply(c(1:length(all.arma11[[chosen]])), function(x) model.iden(all.arma11[[chosen]][[x]]$pq))
dat.all = data.frame(ar = ar1.iden, ma = ma1.iden, arma = arma.iden)
count = sapply(c(1:3), function(y) sapply(list.model, function(x) length(which(dat.all[,y] == x))))
dat.p = reshape2::melt(count)
dat.p$Var2 = rep(list.model[-1], each = 4)
colnames(dat.p) = c("predict", "truth", "value")
dat.p$predict = as.factor(dat.p$predict)

p = ggplot(data = dat.p, aes(x = truth, y = value/nb.sim, fill = predict))+
  theme_bw()+
  geom_bar(position="stack", stat="identity", width = 0.5)+
  ylab("Percentage")+
  labs(subtitle = paste0("coeff = ",  coeff[chosen], ", N = 1000"))+
  scale_y_continuous(labels = scales::percent) +
  # geom_text(aes(label = paste0(value*100/nb.sim,"%")), 
  #           position = position_stack(vjust = 0.5), size = 1)+
  ggrepel::geom_text_repel(aes(label = paste0(value*100/nb.sim,"%")), 
            position = position_stack(vjust = 0.5), size = 1, direction = "y", 
            box.padding = unit(0.01, "lines"))+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text=element_text(size=4),
        axis.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        plot.tag = element_text(size = 5),
        plot.subtitle = element_text(size = 5),
        legend.title=element_blank())
  
ggsave(paste0(path_results,"attribution0/TPR.auto.arima.identification_coef_specific", coeff[chosen], ".jpg" ), plot = p, width = 8.8, height = 5, units = "cm", dpi = 600)




# RUN THE ESTIMATION OF COEFFICIENTS REPARATELY  --------------------------
sim_est <- function(model.order, param, x.axis, sigma.sim, phi, theta, nb.sim){
  res <- list()
  for (i in c(1:length(param))) {
    set.seed(1)
    res.j = data.frame(matrix(NA, ncol = 2, nrow = nb.sim))
    for (j in c(1:nb.sim)) {
      if(x.axis == "n"){
        n = param[i]
        y = simulate.general1(N = n, arma.model = c(ar=phi, ma=theta), burn.in = 1000, hetero = 0, sigma = sigma.sim)
      }else if(x.axis == "phi"){
        n = 1000
        y = simulate.general1(N = n, arma.model = c(ar=phi[i], ma=theta), burn.in = 1000, hetero = 0, sigma = sigma.sim)
      }else if(x.axis == "theta"){
        n = 1000
        y = simulate.general1(N = n, arma.model = c(ar=phi, ma=theta[i]), burn.in = 1000, hetero = 0, sigma = sigma.sim)
      }else{
        n = 1000
        y = simulate.general1(N = n, arma.model = c(ar=phi[i], ma=theta[i]), burn.in = 1000, hetero = 0, sigma = sigma.sim)
      }
      fit.arma = arima(y, order = model.order, method="ML")
      
      if(model.order[1] !=0){
        res.j[j,1] = fit.arma$coef[which(names(fit.arma$coef) == "ar1")]
      }
      if(model.order[3] !=0){
        res.j[j,2] = fit.arma$coef[which(names(fit.arma$coef) == "ma1")]
      }
    }
    res[[i]] = res.j
  }
  return(res)
}
## with length
sim.ar = sim_est(model.order = c(1,0,0), param = length.list, x.axis = "n", sigma.sim = 1, phi = 0.3, theta = 0, nb.sim = 10000)
save(sim.ar, file = paste0(path_results, "attribution0/performance_autoarima_coef_AR_length.RData"))
sim.ma = sim_est(model.order = c(0,0,1), param = length.list, x.axis = "n", sigma.sim = 1, phi = 0, theta = 0.2, nb.sim = 10000)
save(sim.ma, file = paste0(path_results, "attribution0/performance_autoarima_coef_MA_length.RData"))
sim.arma = sim_est(model.order = c(1,0,1), param = length.list, x.axis = "n", sigma.sim = 1, phi = 0.6, theta = -0.3, nb.sim = 10000)
save(sim.arma, file = paste0(path_results, "attribution0/performance_autoarima_coef_ARMA_length.RData"))

# with coefficient 
sim.ar = sim_est(model.order = c(1,0,0), param = coef.list, x.axis = "phi", sigma.sim = 1, phi = coef.list, theta = 0, nb.sim = 10000)
save(sim.ar, file = paste0(path_results, "attribution0/performance_autoarima_coef_AR.RData"))
sim.ma = sim_est(model.order = c(0,0,1), param = coef.list, x.axis = "theta", sigma.sim = 1, phi = 0, theta = coef.list, nb.sim = 10000)
save(sim.ma, file = paste0(path_results, "attribution0/performance_autoarima_coef_MA.RData"))
sim.arma = sim_est(model.order = c(1,0,1), param = coef.list, x.axis = "", sigma.sim = 1, phi = coef.list, theta = (0.3-coef.list), nb.sim = 10000)
save(sim.arma, file = paste0(path_results, "attribution0/performance_autoarima_coef_ARMA.RData"))

dat.ar = as.data.frame(plyr::ldply(sim.ar, rbind))  %>% 
  select("X1") %>% 
  dplyr::rename( "phi" = "X1") %>% 
  mutate(length =rep(as.factor(coef.list), nb.sim)) %>% 
  reshape2::melt(id = "length")

